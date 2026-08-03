# 机器学习诊断模型：基于铜凋亡/铁凋亡基因区分 ICM vs Normal
# 输入：GSE*_symbol.csv, GSE*_group.csv, cell_death_intersect_three.txt
# 输出：ML_signature_genes.txt, model_performance.txt, ROC_curve_*.pdf, feature_importance.pdf

rm(list = ls())
library(glmnet)
library(randomForest)
library(e1071)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# 1. 参数设置 -----------------------------------------------------------
datasets <- c("GSE16499", "GSE5406", "GSE57338")
use_core_only <- FALSE   # FALSE = 用全部 54 个 CRG+FRG DEGs；TRUE = 只用 5 个核心基因

# 读取核心基因和扩展基因
 core_genes <- read.table("cell_death_intersect_three.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]

crg_genes <- unique(unlist(lapply(paste0("CRG_genes_", datasets, ".txt"), function(f) {
  if (!file.exists(f)) return(character(0))
  read.table(f, header = F, sep = "\t", stringsAsFactors = F)[, 1]
})))

frg_genes <- unique(unlist(lapply(paste0("FRG_genes_", datasets, ".txt"), function(f) {
  if (!file.exists(f)) return(character(0))
  read.table(f, header = F, sep = "\t", stringsAsFactors = F)[, 1]
})))

feature_pool <- if (use_core_only) core_genes else unique(c(crg_genes, frg_genes))
cat("Feature pool:", length(feature_pool), "genes\n")
print(feature_pool)

# 2. 数据读取与预处理函数 -----------------------------------------------
prepare_data <- function(ds, feature_genes) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = T, row.names = 1, check.names = F)
  group_df <- read.csv(paste0(ds, "_group.csv"), header = T, check.names = F)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  
  common_samples <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common_samples]
  group <- group_df[common_samples, "group"]
  
  # 提取特征基因
  available_features <- intersect(feature_genes, rownames(exp))
  X <- t(exp[available_features, , drop = FALSE])
  X <- scale(X)  # 标准化
  rownames(X) <- common_samples
  colnames(X) <- available_features
  
  y <- ifelse(group == "Disease", 1, 0)
  names(y) <- common_samples
  
  return(list(X = X, y = y, group = group, features = available_features))
}

# 3. 模型评估函数 -------------------------------------------------------
evaluate_model <- function(prob, y_true) {
  roc_obj <- roc(y_true, prob, quiet = TRUE)
  auc <- auc(roc_obj)
  
  # 使用 0.5 阈值计算敏感性/特异性/准确率（CV 场景下更稳定）
  pred <- ifelse(prob > 0.5, 1, 0)
  cm <- confusionMatrix(factor(pred, levels = c(0, 1)), factor(y_true, levels = c(0, 1)))
  
  return(list(
    AUC = as.numeric(auc),
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Accuracy = cm$overall["Accuracy"],
    Threshold = 0.5,
    roc = roc_obj
  ))
}

# 4. 交叉验证建模 -------------------------------------------------------
all_results <- list()
signature_genes_by_dataset <- list()
roc_data_list <- list()

for (ds in datasets) {
  cat("\n==========", ds, "==========\n")
  dat <- prepare_data(ds, feature_pool)
  X <- dat$X
  y <- dat$y
  
  cat("Samples:", length(y), "| Disease:", sum(y == 1), "| Normal:", sum(y == 0), "\n")
  cat("Available features:", ncol(X), "\n")
  
  if (ncol(X) < 2) {
    cat("Too few features, skipping.\n")
    next
  }
  
  # 设置 10 折交叉验证
  set.seed(123)
  folds <- createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
  
  # 存储每折的预测概率
  prob_lr <- prob_rf <- prob_svm <- numeric(length(y))
  names(prob_lr) <- names(prob_rf) <- names(prob_svm) <- names(y)
  selected_genes_list <- list()
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test <- X[test_idx, , drop = FALSE]
    y_test <- y[test_idx]
    
    # 移除训练集中方差为 0 的特征
    var_train <- apply(X_train, 2, var, na.rm = TRUE)
    X_train <- X_train[, var_train > 0, drop = FALSE]
    X_test <- X_test[, colnames(X_train), drop = FALSE]
    
    # LASSO 特征选择
    if (ncol(X_train) >= 2) {
      cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc")
      coef_lasso <- coef(cv_lasso, s = "lambda.min")
      selected <- rownames(coef_lasso)[coef_lasso[, 1] != 0]
      selected <- setdiff(selected, "(Intercept)")
      selected_genes_list[[i]] <- selected
    } else {
      selected <- colnames(X_train)
      selected_genes_list[[i]] <- selected
    }
    
    # 如果没有选中任何基因，用全部基因
    if (length(selected) == 0) {
      selected <- colnames(X_train)
      selected_genes_list[[i]] <- selected
    }
    
    X_train_sel <- X_train[, selected, drop = FALSE]
    X_test_sel <- X_test[, selected, drop = FALSE]
    
    # Logistic Regression (with LASSO coefficients as probability)
    model_lr <- glmnet(X_train, y_train, family = "binomial", alpha = 1, lambda = cv_lasso$lambda.min)
    prob_lr[test_idx] <- as.numeric(predict(model_lr, newx = X_test, s = cv_lasso$lambda.min, type = "response"))
    
    # Random Forest
    rf_data <- as.data.frame(X_train_sel)
    rf_data$y <- as.factor(y_train)
    model_rf <- randomForest(y ~ ., data = rf_data, ntree = 500, importance = TRUE)
    prob_rf[test_idx] <- as.numeric(predict(model_rf, newdata = as.data.frame(X_test_sel), type = "prob")[, 2])
    
    # SVM
    model_svm <- svm(y ~ ., data = rf_data, probability = TRUE, kernel = "radial")
    pred_svm <- predict(model_svm, newdata = as.data.frame(X_test_sel), probability = TRUE)
    prob_svm[test_idx] <- as.numeric(attr(pred_svm, "probabilities")[, 2])
  }
  
  # 汇总每折选中的基因
  selected_freq <- table(unlist(selected_genes_list))
  signature_genes_by_dataset[[ds]] <- data.frame(
    gene = names(selected_freq),
    freq = as.numeric(selected_freq),
    selected_in_folds = paste0(as.numeric(selected_freq), "/10"),
    stringsAsFactors = F
  ) %>% arrange(desc(freq))
  
  # 评估模型
  eval_lr <- evaluate_model(prob_lr, y)
  eval_rf <- evaluate_model(prob_rf, y)
  eval_svm <- evaluate_model(prob_svm, y)
  
  cat("LASSO-LR AUC:", round(eval_lr$AUC, 3), "\n")
  cat("RF AUC:", round(eval_rf$AUC, 3), "\n")
  cat("SVM AUC:", round(eval_svm$AUC, 3), "\n")
  
  all_results[[ds]] <- data.frame(
    Dataset = ds,
    Model = c("LASSO-LR", "RandomForest", "SVM"),
    AUC = c(eval_lr$AUC, eval_rf$AUC, eval_svm$AUC),
    Sensitivity = c(eval_lr$Sensitivity, eval_rf$Sensitivity, eval_svm$Sensitivity),
    Specificity = c(eval_lr$Specificity, eval_rf$Specificity, eval_svm$Specificity),
    Accuracy = c(eval_lr$Accuracy, eval_rf$Accuracy, eval_svm$Accuracy),
    stringsAsFactors = F
  )
  
  # ROC 数据
  roc_data_list[[ds]] <- rbind(
    data.frame(Dataset = ds, Model = "LASSO-LR", 
               Specificity = 1 - eval_lr$roc$specificities, 
               Sensitivity = eval_lr$roc$sensitivities),
    data.frame(Dataset = ds, Model = "RandomForest", 
               Specificity = 1 - eval_rf$roc$specificities, 
               Sensitivity = eval_rf$roc$sensitivities),
    data.frame(Dataset = ds, Model = "SVM", 
               Specificity = 1 - eval_svm$roc$specificities, 
               Sensitivity = eval_svm$roc$sensitivities)
  )
}

# 5. 保存性能结果 -------------------------------------------------------
performance_df <- bind_rows(all_results)
write.table(performance_df, "model_performance.txt", sep = "\t", quote = F, row.names = F)
cat("\nModel performance:\n")
print(performance_df)

# 6. 保存特征基因 -------------------------------------------------------
sig_df <- bind_rows(lapply(names(signature_genes_by_dataset), function(ds) {
  df <- signature_genes_by_dataset[[ds]]
  df$Dataset <- ds
  return(df)
}))
write.table(sig_df, "ML_signature_genes.txt", sep = "\t", quote = F, row.names = F)

# 找出跨数据集稳定被选中的基因
stable_genes <- sig_df %>%
  filter(freq >= 5) %>%
  group_by(gene) %>%
  summarise(datasets = paste(Dataset, collapse = ", "), total_freq = sum(freq), .groups = "drop") %>%
  arrange(desc(total_freq))
write.table(stable_genes, "ML_stable_signature_genes.txt", sep = "\t", quote = F, row.names = F)
cat("\nStable signature genes (selected in >= 5 folds in at least one dataset):\n")
print(stable_genes)

# 7. 绘制 ROC 曲线 ------------------------------------------------------
roc_all <- bind_rows(roc_data_list)
roc_all$Model <- factor(roc_all$Model, levels = c("LASSO-LR", "RandomForest", "SVM"))

# 计算 AUC 标签
auc_labels <- performance_df %>%
  mutate(label = paste0(Model, " AUC=", round(AUC, 3)))

p_roc <- ggplot(roc_all, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Dataset, ncol = 3) +
  scale_color_manual(values = c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000", "SVM" = "#CD534C")) +
  theme_bw() +
  labs(title = "ROC curves for ICM diagnosis models",
       x = "1 - Specificity", y = "Sensitivity", color = "Model") +
  theme(legend.position = "bottom")

ggsave(p_roc, file = "ROC_curve_all_datasets.pdf", width = 10, height = 4)
cat("\nSaved: ROC_curve_all_datasets.pdf\n")

# 8. 特征选择频率热图 ---------------------------------------------------
sig_wide <- sig_df %>%
  select(gene, Dataset, freq) %>%
  pivot_wider(names_from = Dataset, values_from = freq, values_fill = 0)
sig_wide_mat <- as.matrix(sig_wide[, -1])
rownames(sig_wide_mat) <- sig_wide$gene

if (nrow(sig_wide_mat) > 1) {
  pdf("feature_selection_frequency.pdf", width = 8, height = max(4, nrow(sig_wide_mat) * 0.3))
  pheatmap(sig_wide_mat,
           scale = "none",
           color = colorRampPalette(c("white", "#EFC000", "#CD534C"))(10),
           cluster_cols = F,
           cluster_rows = T,
           display_numbers = TRUE,
           number_color = "black",
           fontsize_number = 7,
           main = "Feature selection frequency across 10-fold CV")
  dev.off()
  cat("Saved: feature_selection_frequency.pdf\n")
}

cat("\nDone!\n")
cat("Output files:\n")
cat("- model_performance.txt\n")
cat("- ML_signature_genes.txt\n")
cat("- ML_stable_signature_genes.txt\n")
cat("- ROC_curve_all_datasets.pdf\n")
cat("- feature_selection_frequency.pdf\n")
