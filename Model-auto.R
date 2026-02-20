rm(list = ls())
library(tidyverse)
library(rms)
library(forestplot)
library(ggDCA)

dat_expr_train <- readRDS('input/GEO_combined_dataset.Rds')
dat_group_train <- readRDS('input/GEO_combined_group.Rds')
dat_expr_vali <- readRDS('input/dat_GSE106724_expr.Rds')
dat_group_vali <- readRDS('input/dat_GSE106724_group.Rds')
key_genes <- read.csv('output/hub genes.CSV')$x

# key_genes <- readRDS('input/key genes.Rds')$X
# saveRDS(key_genes, file = 'input/key genes.Rds')
# write.table(hub_genes,"input/key genes.csv",row.names = F)

# 保存验证集表达数据
# write.table(dat_expr_vali,"input/dat_expr_vali.csv",row.names = T)

# 整理数据
dat_expr_train_key <- dat_expr_train[key_genes, ] %>%
  na.omit() %>% t() %>% as.data.frame()
# dat_expr_train_key$group <- factor(dat_group_train$group, levels = c('Control', 'PCOS'))
dat_expr_train_key$group <- ifelse(dat_group_train$group == 'Control', 0, 1)

# 将数据打包
ddist <- datadist(dat_expr_train_key)
options(datadist = 'ddist')

# 4A 诊断列线图 ---------------------------------------------------------------

# 用lrm方法再构建一个Logistic回归模型
fit_log_lrm = lrm(as.formula(paste0(
  # 公式为所有筛选好的基因与分组之间的关系
  'group ~ ', paste(key_genes, collapse = ' + ')
)), data = dat_expr_train_key, x = T, y = T)

# 诊断列线图
nomo <- nomogram(
  # 这个函数绘制诊断列线图要用lrm构建的回归模型
  fit_log_lrm,
  # 进行logit转换
  fun = plogis,
  # 概率坐标轴刻度
  fun.at = c(0.01, 0.1, 0.5, 0.9, 0.99),
  # 显示预测值
  lp = T,
  funlabel = 'Risk'
)

pdf(file = 'output/4A_diag_nomogram.pdf', width = 12, height = 10)
# 网络颜色
plot(nomo, col.grid = c('grey50', 'lightgrey'))
dev.off()


# 4B 诊断校准曲线 ---------------------------------------------------------------

# 校准分析
set.seed(2024)
dat_cal <- calibrate(
  # 这个函数绘制诊断校准曲线要用lrm构建的回归模型
  fit_log_lrm,
  # 抽样方法
  method = 'boot',
  # 抽样次数
  B = 500)

# 提取结果
dat_cal <- dat_cal[, 3:1] %>% as.data.frame()
colnames(dat_cal) <- c('bias_actual', 'apparent_actual', 'pre')

# 整理结果
dat_cal <- data.frame(
  actual = c(dat_cal$bias_actual, dat_cal$apparent_actual),
  pre = rep(dat_cal$pre, 2),
  group = c(rep('Bias Corrected', nrow(dat_cal)), rep('Apparent', nrow(dat_cal)))
)

# 校准曲线
p <- ggplot(dat_cal, aes(pre, actual, color = group)) +
  # 添加对角线
  geom_abline(
    # 斜率
    slope = 1,
    # 截距
    intercept = 0,
    color = 'grey',
    lty = 2
  ) +
  geom_line() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(x = 'Predicted Probability', y = 'Actual Probability') +
  # 坐标轴范围
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  # 固定xy轴比例
  coord_fixed() +
  scale_color_manual(values = c('#4DBBD5', '#E64B35'))

ggsave(file = 'output/4B_diag_calibration.pdf', p, height = 8, width = 8)


# 4C 诊断DCA ---------------------------------------------------------------
# 1.16版本的R包data.table与可能与R包ggDCA不兼容，如果发生这种情况需要卸载并重新安装R包data.table
# remove.packages('data.table')
# devtools::install_version('data.table',version = '1.14.10')

# 表达矩阵的分组改为字符串
# dat_expr_train_key$group <- as.character(dat_expr_train_key$group)

# DCA分析
p <- dca(fit_log_lrm) %>%
  # DCA曲线
  ggplot(aes(thresholds, NB, color = model)) +
  theme_bw() +
  scale_color_manual(
    values = c('#4DBBD5', '#E64B35', '#3C5488'),
    labels = c('Logistic Model', 'All', 'None')
  ) +
  scale_linetype_manual(values = rep(1, 3)) +
  # 不展示线型图例
  guides(lty = 'none') +
  labs(color = element_blank())

ggsave(file = 'output/4C_diag_DCA.pdf', p, height = 6, width = 7)

# 4D 临床影响曲线 ---------------------------------------------------------------
# 需要用rmda包中的函数绘制
# install.packages('rmda')
library(rmda)

# 首先用rmda包中的dca函数构建dca模型用于临床影响曲线的绘制
dca <- decision_curve(as.formula(paste0(
  # 公式为所有筛选好的基因与分组之间的关系
  'group ~ ', paste(key_genes, collapse = ' + ')
)), dat_expr_train_key)


pdf(file = 'output/4D_CIC.pdf', height = 6, width = 7)
plot_clinical_impact(dca, col = c('#E64B35', '#4DBBD5'))
dev.off()

# 4E、F ROC ---------------------------------------------------------------

# 需要用nomogramFormula包计算患者的nomogram得分
# install.packages('nomogramFormula')
library(pROC)
library(nomogramFormula)

# 计算nomogram得分
nomo_formula <- formula_rd(nomo)
nomo_point_train <- points_cal(nomo_formula$formula, as.matrix(dat_expr_train_key))
nomo_point_vali <- points_cal(nomo_formula$formula, t(dat_expr_vali))

# Fig 4E
# 整理ROC数据
dat_roc_nomo_train <- data.frame(
  row.names = rownames(dat_group_train), 
  value = nomo_point_train, 
  group = factor(dat_group_train$group, levels = c('Control', 'PCOS'))
)

# ROC分析
roc_res_nomo_train <- roc(group ~ value, 
                          data = dat_roc_nomo_train,
                          # 计算AUC和置信区间
                          auc = T,
                          ci = T)

# 整理敏感性和特异性数据
dat_roc_plot_nomo_train <- data.frame(specificity = roc_res_nomo_train$specificities, 
                                      sensitivity = roc_res_nomo_train$sensitivities) %>%
  dplyr::arrange(desc(specificity), sensitivity)

# ROC曲线
p <- ggplot(dat_roc_plot_nomo_train, aes(1 - specificity, sensitivity)) + 
  geom_line(color = '#4DBBD5') +
  # 添加对角线
  geom_abline(
    # 斜率
    slope = 1,
    # 截距
    intercept = 0,
    color = 'grey',
    lty = 'dashed'
  ) +
  # 展示曲线下面积
  geom_area(fill = '#4DBBD5', alpha = 0.2) +
  # 展示AUC值和置信区间
  annotate(
    'text',
    label = paste0(
      'AUC = ',
      round(roc_res_nomo_train$auc, 3),
      ' (',
      round(roc_res_nomo_train$ci[1], 3),
      '-',
      round(roc_res_nomo_train$ci[3], 3),
      ')'
    ),
    x = 0.8,
    y = 0.1,
    color = '#4DBBD5'
  ) +
  theme_bw() +
  # 固定xy轴比例
  coord_fixed() +
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)')

ggsave(file = 'output/4E_ROC_nomo.pdf', p, height = 6, width = 6)

# Fig 4F
# 整理ROC数据
dat_roc_nomo_vali <- data.frame(
  row.names = rownames(dat_group_vali), 
  value = nomo_point_vali, 
  group = factor(dat_group_vali$group, levels = c('Control', 'PCOS'))
)

# ROC分析
roc_res_nomo_vali <- roc(group ~ value, 
                         data = dat_roc_nomo_vali,
                         # 计算AUC和置信区间
                         auc = T,
                         ci = T)

# 整理敏感性和特异性数据
dat_roc_plot_nomo_vali <- data.frame(specificity = roc_res_nomo_vali$specificities, 
                                     sensitivity = roc_res_nomo_vali$sensitivities) %>%
  dplyr::arrange(desc(specificity), sensitivity)

# ROC曲线
p <- ggplot(dat_roc_plot_nomo_vali, aes(1 - specificity, sensitivity)) + 
  geom_line(color = '#4DBBD5') + 
  # 添加对角线
  geom_abline(
    # 斜率
    slope = 1,
    # 截距
    intercept = 0,
    color = 'grey',
    lty = 'dashed'
  ) +
  # 展示曲线下面积
  geom_area(fill = '#4DBBD5', alpha = 0.2) + 
  # 展示AUC值和置信区间
  annotate(
    'text',
    label = paste0(
      'AUC = ',
      round(roc_res_nomo_vali$auc, 3),
      ' (',
      round(roc_res_nomo_vali$ci[1], 3),
      '-',
      round(roc_res_nomo_vali$ci[3], 3),
      ')'
    ),
    x = 0.8,
    y = 0.1,
    color = '#4DBBD5'
  ) +
  theme_bw() + 
  # 固定xy轴比例
  coord_fixed() + 
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)')

ggsave(file = 'output/4F_ROC_nomo.pdf', p, height = 6, width = 6)
