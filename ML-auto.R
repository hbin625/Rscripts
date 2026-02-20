rm(list = ls())
library(tidyverse)
library(caret)
library(DALEX)
library(pROC)
library(ggvenn)
dat_expr <- readRDS('input/GEO_combined_dataset.Rds')
dat_group <- readRDS('input/GEO_combined_group.Rds')

# 去除第一列,将第二列列名X改为id
ch_genes <- read.csv('output/geoDEGs.CSV') 
ch_genes <- ch_genes[, -1]
write.csv(ch_genes, file='output/geoDEGs-M.csv', row.names = FALSE)
ch_genes <- read.csv('output/geoDEGs-M.CSV')
colnames(ch_genes)[1] <- "id"


dat_expr<-dat_expr[,dat_group$sample]

# 构建训练集数据和测试集数据 ---------------------------------------------------------------

# 整理数据
dat_expr_ch <- dat_expr[ch_genes$id, ] %>% t() %>% as.data.frame()
dat_expr_ch$group <- ifelse(dat_group$group == 'Control', '0', '1')

# 创建训练集和验证集
set.seed(2024)
idx <- createDataPartition(dat_expr_ch$group, p = 0.7, list = F) %>% as.vector()
dat_train <- dat_expr_ch[idx, ]
dat_test <- dat_expr_ch[-idx, ]


# 构建模型 ---------------------------------------------------------------

# 创建模型训练配置对象
tr_control <- trainControl(
  # 使用重复交叉验证的方法
  method = 'repeatedcv',
  # 五折交叉验证
  repeats = 5)

# 使用不同的机器学习算法构建模型
# SVM
model_svm <- train(
  group ~ .,
  data = dat_train,
  # 训练方法
  method = 'svmRadial',
  # 需要保存建模时的概率模型
  prob.model = T,
  # 只对参数进行一次调优
  tuneLength = 1,
  # 模型训练配置对象
  trControl = tr_control
)

# RF
model_rf <- train(
  group ~ .,
  data = dat_train,
  # 训练方法
  method = 'rf',
  # 设置森林中树的数量
  ntree = 100,
  # 需要保存建模时的概率模型
  prob.model = T,
  # 只对参数进行一次调优
  tuneLength = 1,
  # 模型训练配置对象
  trControl = tr_control
)
# xgb
model_xgb <- train(
  group ~ .,
  data = dat_train,
  # 训练方法
  method = 'xgbTree',
  # 只对参数进行一次调优
  tuneLength = 1,
  # 模型训练配置对象
  trControl = tr_control
)
# glm
model_glm <- train(
  group ~ .,
  data = dat_train,
  # 训练方法
  method = 'glm',
  # 结局为二分类变量
  family = 'binomial',
  # 只对参数进行一次调优
  tuneLength = 1,
  # 模型训练配置对象
  trControl = tr_control
)


# 对算法模型进行验证分析 ---------------------------------------------------------------

# 对测试集进行预测，并返回预测值
predicted_fun <- function(object,newdata){
  predict(object,newdata,type = 'prob')[,2]
}

# 解释分析
# SVM
explainer_svm <- DALEX::explain(
  model_svm,
  label = 'SVM',
  data = dat_test,
  # 观测值
  y = as.numeric(dat_test$group),
  # 预测值
  predict_function = predicted_fun
)
# RF
explainer_rf <- DALEX::explain(
  model_rf,
  label = 'RF',
  data = dat_test,
  # 观测值
  y = as.numeric(dat_test$group),
  # 预测值
  predict_function = predicted_fun
)
# XGB
explainer_xgb <- DALEX::explain(
  model_xgb,
  label = 'XGB',
  data = dat_test,
  # 观测值
  y = as.numeric(dat_test$group),
  # 预测值
  predict_function = predicted_fun
)
# GLM
explainer_glm <- DALEX::explain(
  model_glm,
  label = 'GLM',
  data = dat_test,
  # 观测值
  y = as.numeric(dat_test$group),
  # 预测值
  predict_function = predicted_fun
)


# 残差 ---------------------------------------------------------------

# 计算残差
# SVM
mp_svm = model_performance(explainer_svm)
# RF
mp_rf = model_performance(explainer_rf)
# XGB
mp_xgb = model_performance(explainer_xgb)
# glm
mp_glm = model_performance(explainer_glm)

# 残差图
pdf(file = 'Output/8b-residual.pdf', width = 6,height = 5)
plot(mp_svm, mp_rf, mp_xgb, mp_glm)
dev.off()
pdf(file = 'Output/8c-residual_boxplot.pdf', width = 6,height = 5)
plot(mp_svm, mp_rf, mp_xgb, mp_glm, geom = 'boxplot')
dev.off()


# ROC ---------------------------------------------------------------

# ROC数据
dat_roc <- data.frame(
  group = dat_test$group,
  # SVM模型的预测值
  SVM = predict(model_svm, dat_test, type = 'prob')$`1`,
  # RF模型的预测值
  RF = predict(model_rf, dat_test, type = 'prob')$`1`,
  # XGB模型的预测值
  XGB = predict(model_xgb, dat_test, type = 'prob')$`1`,
  # GLM模型的预测值
  GLM = predict(model_glm, dat_test, type = 'prob')$`1`
)

# ROC分析
roc_res_model <- roc(group ~ SVM + RF + XGB + GLM,
                     data = dat_roc,
                     # 计算AUC和置信区间
                     auc = T,
                     ci = T)

# 整理敏感性和特异性数据
dat_roc_plot_model <- rbind(
  data.frame(
    specificity = roc_res_model$SVM$specificities,
    sensitivity = roc_res_model$SVM$sensitivities,
    model = 'SVM'
  ) %>%
    dplyr::arrange(desc(specificity), sensitivity),
  data.frame(
    specificity = roc_res_model$RF$specificities,
    sensitivity = roc_res_model$RF$sensitivities,
    model = 'RF'
  ) %>%
    dplyr::arrange(desc(specificity), sensitivity),
  data.frame(
    specificity = roc_res_model$XGB$specificities,
    sensitivity = roc_res_model$XGB$sensitivities,
    model = 'XGB'
  ) %>%
    dplyr::arrange(desc(specificity), sensitivity),
  data.frame(
    specificity = roc_res_model$GLM$specificities,
    sensitivity = roc_res_model$GLM$sensitivities,
    model = 'GLM'
  ) %>%
    dplyr::arrange(desc(specificity), sensitivity)
)
dat_roc_plot_model$model <- factor(dat_roc_plot_model$model, levels = c('SVM', 'RF', 'XGB', 'GLM'))

# ROC曲线
p <- ggplot(dat_roc_plot_model,
            aes(1 - specificity, sensitivity, color = model)) +
  geom_line() +
  # 添加对角线
  geom_abline(
    # 斜率
    slope = 1,
    # 截距
    intercept = 0,
    color = 'grey',
    lty = 'dashed'
  ) +
  # 展示AUC值和置信区间
  annotate(
    'text',
    label = paste0(
      'AUC(SVM) = ',
      round(roc_res_model$SVM$auc, 3),
      ' (',
      round(roc_res_model$SVM$ci[1], 3),
      '-',
      round(roc_res_model$SVM$ci[3], 3),
      ')'
    ),
    x = 0.75,
    y = 0.2,
    color = '#4DBBD5'
  ) +
  annotate(
    'text',
    label = paste0(
      'AUC(RF) = ',
      round(roc_res_model$RF$auc, 3),
      ' (',
      round(roc_res_model$RF$ci[1], 3),
      '-',
      round(roc_res_model$RF$ci[3], 3),
      ')'
    ),
    x = 0.75,
    y = 0.15,
    color = '#E64B35'
  ) +
  annotate(
    'text',
    label = paste0(
      'AUC(XGB) = ',
      round(roc_res_model$XGB$auc, 3),
      ' (',
      round(roc_res_model$XGB$ci[1], 3),
      '-',
      round(roc_res_model$XGB$ci[3], 3),
      ')'
    ),
    x = 0.75,
    y = 0.1,
    color = '#3C5488'
  ) +
  annotate(
    'text',
    label = paste0(
      'AUC(GLM) = ',
      round(roc_res_model$GLM$auc, 3),
      ' (',
      round(roc_res_model$GLM$ci[1], 3),
      '-',
      round(roc_res_model$GLM$ci[3], 3),
      ')'
    ),
    x = 0.75,
    y = 0.05,
    color = '#F39B7F'
  ) +
  theme_bw() +
  # 固定xy轴比例
  coord_fixed() +
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)') +
  scale_color_manual(values = c('#4DBBD5', '#E64B35', '#3C5488', '#F39B7F'))

ggsave(file = 'Output/8d-ROC_model.pdf', p, height = 6, width = 7)


# 变量重要性分析 ---------------------------------------------------------------

# 计算模型各基因的重要性
# SVM
vi_svm <- variable_importance(explainer_svm,
                              # 计算均方根
                              loss_function = loss_root_mean_square)
# RF
vi_rf <- variable_importance(explainer_rf,
                             # 计算均方根
                             loss_function = loss_root_mean_square)
# XGB
vi_xgb <- variable_importance(explainer_xgb,
                              # 计算均方根
                              loss_function = loss_root_mean_square)
# GLM
vi_glm <- variable_importance(explainer_glm,
                              # 计算均方根
                              loss_function = loss_root_mean_square)

# 特征变量重要性柱状图
pdf(file = 'Output/8a-variable_importance.pdf', width = 6,height = 8)
plot(vi_svm, vi_rf, vi_xgb, vi_glm,
     # 展示前10个基因
     max_vars = 10)
dev.off()

# 选择基因
# SVM
top_gene_svm <- vi_svm %>%
  dplyr::filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(m = median(dropout_loss)) %>%
  arrange(desc(m))
top_gene_svm <- top_gene_svm$variable[1:20]
# RF
top_gene_rf <- vi_rf %>%
  dplyr::filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(m = median(dropout_loss)) %>%
  arrange(desc(m))
top_gene_rf <- top_gene_rf$variable[1:20]
# XGB
top_gene_xgb <- vi_xgb %>%
  dplyr::filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(m = median(dropout_loss)) %>%
  arrange(desc(m))
top_gene_xgb <- top_gene_xgb$variable[1:20]
# GLM
top_gene_glm <- vi_glm %>%
  dplyr::filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(m = median(dropout_loss)) %>%
  arrange(desc(m))
top_gene_glm <- top_gene_glm$variable[1:20]

# 交集韦恩图
p <- ggvenn(
  list(
    SVM = top_gene_svm,
    RF = top_gene_rf,
    XGB = top_gene_xgb,
    GLM = top_gene_glm
  ),
  # 下面要与列表中的命名一致
  c('SVM', 'RF', 'XGB', 'GLM'),
  # 不展示比例
  show_percentage = F,
  fill_alpha = 0.5,
  stroke_color = NA,
  fill_color = c('#4DBBD5', '#E64B35', '#3C5488', '#F39B7F')
)

ggsave(file = 'Output/8e-model_gene_venn.pdf', p, width = 5, height = 5)

# 交集情况
gene_model <- Reduce(intersect,
                     list(top_gene_svm, top_gene_rf, top_gene_xgb, top_gene_glm))
# 保存交集基因
write.table(gene_model,"output/hub genes.csv",row.names = F)
saveRDS(gene_model, file = 'output/hub genes.Rds')

# 诊断ROC 以ATG2B为例 ---------------------------------------------------------------
# 整理分组数据
dat_test1 <- dat_test[,c('ATG2B','group')]

# ROC分析
roc_res_ATG2B <- roc(group ~ ATG2B,
                     data = dat_test1,
                     # 计算AUC和置信区间
                     auc = T,
                     ci = T)

# 找最佳cutoff值
# dat_roc_cutoff_ATG2B <- cutoff::roc(dat_test1$ATG2B, dat_test1$group)

# 整理敏感性和特异性数据
dat_roc_plot_ATG2B <- data.frame(specificity = roc_res_ATG2B$specificities,
                                 sensitivity = roc_res_ATG2B$sensitivities) %>%
  dplyr::arrange(desc(specificity), sensitivity)

# ROC曲线
p <- ggplot(dat_roc_plot_ATG2B, aes(1 - specificity, sensitivity)) +
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
  # geom_area(fill = '#4DBBD5', alpha = 0.2) +
  # 展示最佳cutoff值点
  # geom_point(data = dat_roc_cutoff_ATG2B,
  #            aes(x = 1 - specificity, y = sensitivity),
  #            color = '#4DBBD5') +
  # 展示最佳cutoff值点标签
  # ggrepel::geom_text_repel(
  #   data = dat_roc_cutoff_ATG2B,
  #   aes(
  #     x = 1 - specificity,
  #     y = sensitivity,
  #     label = paste0(
  #       round(cutoff, 3),
  #       ' (',
  #       round(1 - specificity, 3),
  #       ',',
  #       round(sensitivity, 3),
  #       ')'
  #     )
  #   ),
  #   color = '#4DBBD5'
  # ) +
  # 展示AUC值和置信区间
  annotate(
    'text',
    label = paste0(
      'AUC = ',
      round(roc_res_ATG2B$auc, 3),
      ' (',
      round(roc_res_ATG2B$ci[1], 3),
      '-',
      round(roc_res_ATG2B$ci[3], 3),
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

ggsave(file = 'output/8f-ROC_ATG2B.pdf', p, height = 6, width = 6)

