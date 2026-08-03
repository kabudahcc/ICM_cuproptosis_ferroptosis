# 外部验证数据集 GSE116250 分析
# 输入：GSE116250_rpkm.txt.gz, GSE116250_series_matrix.txt.gz
# 输出：GSE116250_diff.txt, external_validation_core_genes.pdf, external_validation_ROC.pdf

rm(list = ls())
library(dplyr)
library(ggplot2)
library(limma)
library(glmnet)
library(randomForest)
library(pROC)
library(tidyr)
library(tibble)

# 1. 读取 GSE116250 RPKM 数据 ---------------------------------------------
exp <- read.table("GSE116250_rpkm.txt.gz", header = T, sep = "\t", check.names = F, quote = "")
exp <- exp %>% filter(Common_name != "")
exp <- exp[!duplicated(exp$Common_name), ]
rownames(exp) <- exp$Common_name
exp <- exp[, -(1:2)]
exp <- as.matrix(exp)

cat("GSE116250 expression matrix:", nrow(exp), "genes x", ncol(exp), "samples\n")

# 2. 读取分组信息 ---------------------------------------------------------
library(GEOquery)
gse <- getGEO("GSE116250", destdir = ".", getGPL = FALSE)
pd <- pData(phenoData(gse[[1]]))
disease_col <- grep("disease", names(pd), value = TRUE)
pd$disease <- pd[[disease_col]]
pd <- pd %>% filter(disease %in% c("non-failing", "ischemic cardiomyopathy"))
pd$group <- ifelse(pd$disease == "non-failing", "Normal", "Disease")
pd$geo_accession <- rownames(pd)

cat("External validation samples:", nrow(pd), "\n")
print(table(pd$group))

# 3. 提取 ICM vs NF 的表达矩阵 -------------------------------------------
# RPKM 文件列名是 title（如 NF10, ICM47），需要映射到 geo_accession
title_to_geo <- data.frame(
  title = pd$title,
  geo_accession = rownames(pd),
  group = pd$group,
  stringsAsFactors = F
)
# 去除 RPKM 列名中可能的空格或特殊字符
colnames(exp) <- trimws(colnames(exp))
# 找到匹配的样本
matched <- title_to_geo %>% filter(title %in% colnames(exp))
common_titles <- matched$title
exp_val <- exp[, common_titles, drop = FALSE]
group_val <- factor(matched$group, levels = c("Normal", "Disease"))
# 将列名改为 geo_accession 以便后续统一
colnames(exp_val) <- matched$geo_accession

cat("Common samples:", ncol(exp_val), "\n")

# 4. 差异表达分析（limma）-------------------------------------------------
# RPKM 数据已经是 log2 转换后的值（看数值范围），但稳妥起见做 log2(x+1)
exp_log <- log2(exp_val + 1)

design <- model.matrix(~ 0 + group_val)
colnames(design) <- c("Normal", "Disease")
contrast.matrix <- makeContrasts(Disease - Normal, levels = design)

fit <- lmFit(exp_log, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff_val <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
diff_val$symbol <- rownames(diff_val)
diff_val <- diff_val %>% select(symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

write.table(diff_val, "GSE116250_diff.txt", sep = "\t", quote = F, row.names = F)
cat("GSE116250 DEGs (|logFC|>0.5 & P<0.05):", sum(abs(diff_val$logFC) > 0.5 & diff_val$P.Value < 0.05), "\n")

# 5. 验证 5 个核心基因 ----------------------------------------------------
core_genes <- read.table("cell_death_intersect_three.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]
core_in_val <- intersect(core_genes, rownames(exp_log))
cat("Core genes in validation set:", core_in_val, "\n")

if (length(core_in_val) > 0) {
  core_exp <- exp_log[core_in_val, , drop = FALSE]
  core_df <- as.data.frame(t(core_exp))
  core_df$sample <- rownames(core_df)
  core_df$group <- as.character(group_val)
  
  core_long <- core_df %>%
    pivot_longer(cols = all_of(core_in_val), names_to = "Gene", values_to = "Expression")
  core_long$Gene <- factor(core_long$Gene, levels = core_in_val)
  
  # Wilcoxon 检验
  test_results <- core_long %>%
    group_by(Gene) %>%
    summarise(pvalue = wilcox.test(Expression ~ group)$p.value,
              logFC = mean(Expression[group == "Disease"]) - mean(Expression[group == "Normal"]),
              .groups = "drop") %>%
    mutate(signif = ifelse(pvalue < 0.001, "***",
                           ifelse(pvalue < 0.01, "**",
                                  ifelse(pvalue < 0.05, "*", "ns"))))
  write.table(test_results, "external_validation_core_genes_test.txt", sep = "\t", quote = F, row.names = F)
  
  # 箱线图
  p1 <- ggplot(core_long, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~ Gene, scales = "free_y", ncol = length(core_in_val)) +
    scale_fill_manual(values = c("Normal" = "#0073C2", "Disease" = "#CD534C")) +
    labs(title = "Core gene expression in external validation (GSE116250)",
         x = "", y = "log2(RPKM+1)", fill = "Group") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(p1, file = "external_validation_core_genes.pdf", width = 12, height = 4)
  cat("Saved: external_validation_core_genes.pdf\n")
}

# 6. 用稳定特征基因在验证集内训练模型 ------------------------------------
stable_genes <- read.table("ML_stable_signature_genes.txt", header = T, sep = "\t", stringsAsFactors = F)$gene
stable_in_val <- intersect(stable_genes, rownames(exp_log))
cat("Stable ML genes in validation set:", length(stable_in_val), "\n")

if (length(stable_in_val) >= 2) {
  X_val <- t(exp_log[stable_in_val, , drop = FALSE])
  X_val <- scale(X_val)
  y_val <- ifelse(group_val == "Disease", 1, 0)
  
  # 10-fold CV
  set.seed(123)
  folds <- caret::createFolds(y_val, k = 5, list = TRUE)  # 样本少，用 5-fold
  
  prob_lr <- prob_rf <- numeric(length(y_val))
  names(prob_lr) <- names(prob_rf) <- colnames(exp_log)
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y_val), test_idx)
    
    X_train <- X_val[train_idx, , drop = FALSE]
    y_train <- y_val[train_idx]
    X_test <- X_val[test_idx, , drop = FALSE]
    
    # 移除训练集方差为 0 的特征
    var_train <- apply(X_train, 2, var, na.rm = TRUE)
    X_train <- X_train[, var_train > 0, drop = FALSE]
    X_test <- X_test[, colnames(X_train), drop = FALSE]
    
    # LASSO-LR
    cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 3, type.measure = "auc")
    model_lr <- glmnet(X_train, y_train, family = "binomial", alpha = 1, lambda = cv_lasso$lambda.min)
    prob_lr[test_idx] <- as.numeric(predict(model_lr, newx = X_test, s = cv_lasso$lambda.min, type = "response"))
    
    # Random Forest
    rf_data <- as.data.frame(X_train)
    rf_data$y <- as.factor(y_train)
    model_rf <- randomForest(y ~ ., data = rf_data, ntree = 500)
    prob_rf[test_idx] <- as.numeric(predict(model_rf, newdata = as.data.frame(X_test), type = "prob")[, 2])
  }
  
  # 评估
  roc_lr <- roc(y_val, prob_lr, quiet = TRUE)
  roc_rf <- roc(y_val, prob_rf, quiet = TRUE)
  
  auc_lr <- as.numeric(auc(roc_lr))
  auc_rf <- as.numeric(auc(roc_rf))
  
  cat("External validation AUC - LASSO-LR:", round(auc_lr, 3), "\n")
  cat("External validation AUC - RandomForest:", round(auc_rf, 3), "\n")
  
  # ROC 曲线
  roc_df <- rbind(
    data.frame(Model = paste0("LASSO-LR (AUC=", round(auc_lr, 3), ")"),
               Specificity = 1 - roc_lr$specificities,
               Sensitivity = roc_lr$sensitivities),
    data.frame(Model = paste0("RandomForest (AUC=", round(auc_rf, 3), ")"),
               Specificity = 1 - roc_rf$specificities,
               Sensitivity = roc_rf$sensitivities)
  )
  
  p2 <- ggplot(roc_df, aes(x = Specificity, y = Sensitivity, color = Model)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("#0073C2", "#EFC000")) +
    theme_bw() +
    labs(title = "External validation ROC (GSE116250)",
         x = "1 - Specificity", y = "Sensitivity", color = "Model") +
    theme(legend.position = "bottom")
  
  ggsave(p2, file = "external_validation_ROC.pdf", width = 6, height = 6)
  cat("Saved: external_validation_ROC.pdf\n")
  
  # 保存性能
  perf_df <- data.frame(
    Dataset = "GSE116250",
    Model = c("LASSO-LR", "RandomForest"),
    AUC = c(auc_lr, auc_rf),
    stringsAsFactors = F
  )
  write.table(perf_df, "external_validation_performance.txt", sep = "\t", quote = F, row.names = F)
}

cat("\nDone!\n")
