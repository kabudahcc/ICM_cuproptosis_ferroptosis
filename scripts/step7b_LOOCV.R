# 留一数据集交叉验证（LOOCV）：用两个数据集训练，第三个数据集测试
# 评估 LASSO-LR / RandomForest / SVM 在跨数据集场景下的泛化性能

rm(list = ls())
library(glmnet)
library(randomForest)
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)

datasets <- c("GSE16499", "GSE5406", "GSE57338")

# 读取特征池（铜凋亡 + 铁凋亡 DEGs 并集）
crg_genes <- unique(unlist(lapply(paste0("CRG_genes_", datasets, ".txt"), function(f) {
  if (!file.exists(f)) return(character(0))
  read.table(f, header = F, sep = "\t", stringsAsFactors = F)[, 1]
})))
frg_genes <- unique(unlist(lapply(paste0("FRG_genes_", datasets, ".txt"), function(f) {
  if (!file.exists(f)) return(character(0))
  read.table(f, header = F, sep = "\t", stringsAsFactors = F)[, 1]
})))
feature_pool <- unique(c(crg_genes, frg_genes))
cat("Feature pool:", length(feature_pool), "genes\n")

# 稳定特征基因
stable_genes <- read.table("ML_stable_signature_genes.txt", header = T, sep = "\t", stringsAsFactors = F)$gene

core_genes <- read.table("cell_death_intersect_three.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]

# 数据读取函数
prepare_data <- function(ds, feature_genes) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = T, row.names = 1, check.names = F)
  group_df <- read.csv(paste0(ds, "_group.csv"), header = T, check.names = F)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  
  common_samples <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common_samples]
  group <- group_df[common_samples, "group"]
  
  available_features <- intersect(feature_genes, rownames(exp))
  X <- t(exp[available_features, , drop = FALSE])
  rownames(X) <- common_samples
  colnames(X) <- available_features
  
  y <- ifelse(group == "Disease", 1, 0)
  names(y) <- common_samples
  
  return(list(X = X, y = y, group = group, features = available_features))
}

# 模型评估
evaluate_model <- function(prob, y_true) {
  roc_obj <- roc(y_true, prob, quiet = TRUE)
  auc <- as.numeric(auc(roc_obj))
  pred <- ifelse(prob > 0.5, 1, 0)
  cm <- confusionMatrix(factor(pred, levels = c(0, 1)), factor(y_true, levels = c(0, 1)))
  list(
    AUC = auc,
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    roc = roc_obj
  )
}

# 标准化：用训练集参数标准化测试集
scale_with_params <- function(X_train, X_test) {
  center <- attr(X_train, "scaled:center")
  scale <- attr(X_train, "scaled:scale")
  X_test_scaled <- scale(X_test, center = center, scale = scale)
  # 若测试集某些特征方差为 0，scale 会产生 NaN，替换为 0
  X_test_scaled[is.nan(X_test_scaled)] <- 0
  X_test_scaled
}

# 留一数据集 CV
loocv_results <- list()
roc_data_list <- list()
selected_genes_list <- list()

for (test_ds in datasets) {
  train_datasets <- setdiff(datasets, test_ds)
  cat("\n========== Leave-one-out:", test_ds, "as test ==========\n")
  
  # 读取训练数据
  train_list <- lapply(train_datasets, prepare_data, feature_genes = feature_pool)
  # 取训练集共同特征
  common_train_features <- Reduce(intersect, lapply(train_list, function(x) colnames(x$X)))
  cat("Common features in training datasets:", length(common_train_features), "\n")
  
  X_train <- do.call(rbind, lapply(train_list, function(x) x$X[, common_train_features, drop = FALSE]))
  y_train <- do.call(c, lapply(train_list, function(x) x$y))
  
  # 测试数据
  test_dat <- prepare_data(test_ds, feature_genes = common_train_features)
  X_test <- test_dat$X
  y_test <- test_dat$y
  
  # 只保留训练集和测试集都有的特征
  final_features <- intersect(colnames(X_train), colnames(X_test))
  if (length(final_features) < 2) {
    cat("Too few overlapping features, skipping\n")
    next
  }
  X_train <- X_train[, final_features, drop = FALSE]
  X_test <- X_test[, final_features, drop = FALSE]
  
  # 标准化训练集，并记录参数
  X_train_scaled <- scale(X_train)
  X_test_scaled <- scale_with_params(X_train_scaled, X_test)
  
  # 移除训练集方差为 0 的特征
  var_train <- apply(X_train_scaled, 2, var, na.rm = TRUE)
  nonzero_var <- var_train > 0
  X_train_scaled <- X_train_scaled[, nonzero_var, drop = FALSE]
  X_test_scaled <- X_test_scaled[, nonzero_var, drop = FALSE]
  
  # LASSO 特征选择
  set.seed(123)
  cv_lasso <- cv.glmnet(X_train_scaled, y_train, family = "binomial", alpha = 1,
                        nfolds = 5, type.measure = "auc")
  coef_lasso <- coef(cv_lasso, s = "lambda.min")
  selected <- setdiff(rownames(coef_lasso)[coef_lasso[, 1] != 0], "(Intercept)")
  if (length(selected) == 0) selected <- colnames(X_train_scaled)
  selected_genes_list[[test_ds]] <- selected
  cat("Selected genes:", length(selected), "\n")
  
  X_train_sel <- X_train_scaled[, selected, drop = FALSE]
  X_test_sel <- X_test_scaled[, selected, drop = FALSE]
  
  # LASSO-LR
  model_lr <- glmnet(X_train_scaled, y_train, family = "binomial", alpha = 1,
                     lambda = cv_lasso$lambda.min)
  prob_lr <- as.numeric(predict(model_lr, newx = X_test_scaled,
                                s = cv_lasso$lambda.min, type = "response"))
  
  # Random Forest
  rf_data <- as.data.frame(X_train_sel)
  rf_data$y <- as.factor(y_train)
  model_rf <- randomForest(y ~ ., data = rf_data, ntree = 500)
  prob_rf <- as.numeric(predict(model_rf, newdata = as.data.frame(X_test_sel), type = "prob")[, 2])
  
  # SVM
  model_svm <- svm(y ~ ., data = rf_data, probability = TRUE, kernel = "radial")
  pred_svm <- predict(model_svm, newdata = as.data.frame(X_test_sel), probability = TRUE)
  prob_svm <- as.numeric(attr(pred_svm, "probabilities")[, 2])
  
  # 评估
  eval_lr <- evaluate_model(prob_lr, y_test)
  eval_rf <- evaluate_model(prob_rf, y_test)
  eval_svm <- evaluate_model(prob_svm, y_test)
  
  cat("LASSO-LR AUC:", round(eval_lr$AUC, 3), "\n")
  cat("RF AUC:", round(eval_rf$AUC, 3), "\n")
  cat("SVM AUC:", round(eval_svm$AUC, 3), "\n")
  
  loocv_results[[test_ds]] <- data.frame(
    TestDataset = test_ds,
    TrainDatasets = paste(train_datasets, collapse = "+"),
    Model = c("LASSO-LR", "RandomForest", "SVM"),
    AUC = c(eval_lr$AUC, eval_rf$AUC, eval_svm$AUC),
    Sensitivity = c(eval_lr$Sensitivity, eval_rf$Sensitivity, eval_svm$Sensitivity),
    Specificity = c(eval_lr$Specificity, eval_rf$Specificity, eval_svm$Specificity),
    Accuracy = c(eval_lr$Accuracy, eval_rf$Accuracy, eval_svm$Accuracy),
    stringsAsFactors = F
  )
  
  roc_data_list[[test_ds]] <- rbind(
    data.frame(TestDataset = test_ds, Model = "LASSO-LR",
               Specificity = 1 - eval_lr$roc$specificities,
               Sensitivity = eval_lr$roc$sensitivities),
    data.frame(TestDataset = test_ds, Model = "RandomForest",
               Specificity = 1 - eval_rf$roc$specificities,
               Sensitivity = eval_rf$roc$sensitivities),
    data.frame(TestDataset = test_ds, Model = "SVM",
               Specificity = 1 - eval_svm$roc$specificities,
               Sensitivity = eval_svm$roc$sensitivities)
  )
}

# 汇总结果
perf_df <- bind_rows(loocv_results)
write.table(perf_df, "LOOCV_model_performance.txt", sep = "\t", quote = F, row.names = F)
cat("\nSaved: LOOCV_model_performance.txt\n")
print(perf_df)

# 保存每个测试集选中的基因
selected_df <- bind_rows(lapply(names(selected_genes_list), function(ds) {
  data.frame(TestDataset = ds, Gene = selected_genes_list[[ds]], stringsAsFactors = F)
}))
write.table(selected_df, "LOOCV_selected_genes.txt", sep = "\t", quote = F, row.names = F)

# 绘制 ROC
roc_all <- bind_rows(roc_data_list)
roc_all$Model <- factor(roc_all$Model, levels = c("LASSO-LR", "RandomForest", "SVM"))

p <- ggplot(roc_all, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ TestDataset, ncol = 2) +
  scale_color_manual(values = c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000", "SVM" = "#CD534C")) +
  theme_bw() +
  labs(title = "Leave-one-dataset cross-validation ROC",
       x = "1 - Specificity", y = "Sensitivity", color = "Model")
ggsave("LOOCV_ROC_curve.pdf", p, width = 10, height = 8)
cat("Saved: LOOCV_ROC_curve.pdf\n")
