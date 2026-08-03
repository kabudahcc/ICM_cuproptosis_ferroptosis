# 外部验证：GSE55296（RNA-seq）和 GSE42955（Affymetrix 微阵列）
# 使用跨发现队列稳定特征基因评估 ICM vs Normal 分类性能

rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(glmnet)
library(randomForest)
library(pROC)
library(caret)
library(limma)
library(edgeR)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)

# 1. 读取签名基因 ---------------------------------------------------------
stable_genes <- read.table("ML_stable_signature_genes.txt", header = T, sep = "\t", stringsAsFactors = F)$gene
core_genes <- read.table("cell_death_intersect_three.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]

cat("Stable signature genes:", length(stable_genes), "\n")
cat("Core genes:", length(core_genes), "\n")

# 评估函数
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

# 2. GSE55296 RNA-seq -----------------------------------------------------
cat("\n========== GSE55296 ==========\n")

# 读取 count 数据
counts <- read.delim(gzfile("GSE55296_count_data.txt.gz"), check.names = FALSE, stringsAsFactors = FALSE)
# 第一列 Ensembl，第二列 Symbol，之后为样本列
colnames(counts)[1:2] <- c("Ensembl", "Symbol")
colnames(counts) <- make.names(colnames(counts), unique = TRUE)
# 去除文件末尾多余的空列
counts <- counts[, colSums(!is.na(counts) & counts != "") > 0]
counts <- counts %>% filter(Symbol != "") %>% filter(!duplicated(Symbol))
rownames(counts) <- counts$Symbol
exp_counts <- as.matrix(counts[, -(1:2)])
# 文件中个别计数位置混入了样本 ID 字符串，强制转数值（错误值变 NA）
storage.mode(exp_counts) <- "numeric"

# 读取样本映射与分组
readme <- read.delim("GSE55296_processed_data_readme.txt", check.names = FALSE, stringsAsFactors = FALSE)
colnames(readme) <- c("geo_accession", "title", "col_header")
readme$group <- ifelse(grepl("Control", readme$title, ignore.case = TRUE), "Normal",
                       ifelse(grepl("Ischemic", readme$title, ignore.case = TRUE), "Disease", NA))
readme <- readme %>% filter(!is.na(group))
readme <- readme[readme$col_header %in% colnames(exp_counts), ]

cat("GSE55296 samples used:", nrow(readme), "\n")
print(table(readme$group))

exp55296 <- exp_counts[, readme$col_header, drop = FALSE]
colnames(exp55296) <- readme$geo_accession

# 数据已是标准化后的表达量（类似 RPKM/FPKM），直接 log2(x+1) 转换
exp_log <- log2(exp55296 + 1)

# 选择稳定签名基因中可映射的基因
features <- intersect(stable_genes, rownames(exp_log))
cat("Stable features in GSE55296:", length(features), "\n")

perf_55296 <- list()
roc_55296 <- NULL
if (length(features) >= 2) {
  X <- t(exp_log[features, , drop = FALSE])
  y <- ifelse(readme$group == "Disease", 1, 0)
  names(y) <- readme$geo_accession
  
  # 5-fold CV（样本少）
  set.seed(123)
  folds <- createFolds(factor(y), k = 5, list = TRUE)
  prob_lr <- prob_rf <- numeric(length(y))
  names(prob_lr) <- names(prob_rf) <- rownames(X)
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    X_train <- scale(X[train_idx, , drop = FALSE])
    X_test <- scale(X[test_idx, , drop = FALSE],
                    center = attr(X_train, "scaled:center"),
                    scale = attr(X_train, "scaled:scale"))
    X_test[is.nan(X_test)] <- 0
    y_train <- y[train_idx]
    
    var_train <- apply(X_train, 2, var, na.rm = TRUE)
    keep <- var_train > 0
    X_train <- X_train[, keep, drop = FALSE]
    X_test <- X_test[, keep, drop = FALSE]
    
    cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1,
                          nfolds = length(y_train), type.measure = "deviance")
    model_lr <- glmnet(X_train, y_train, family = "binomial", alpha = 1,
                       lambda = cv_lasso$lambda.min)
    prob_lr[test_idx] <- as.numeric(predict(model_lr, newx = X_test,
                                            s = cv_lasso$lambda.min, type = "response"))
    
    rf_data <- as.data.frame(X_train)
    rf_data$y <- as.factor(y_train)
    model_rf <- randomForest(y ~ ., data = rf_data, ntree = 500)
    prob_rf[test_idx] <- as.numeric(predict(model_rf, newdata = as.data.frame(X_test), type = "prob")[, 2])
  }
  
  eval_lr <- evaluate_model(prob_lr, y)
  eval_rf <- evaluate_model(prob_rf, y)
  
  cat("GSE55296 LASSO-LR AUC:", round(eval_lr$AUC, 3), "\n")
  cat("GSE55296 RF AUC:", round(eval_rf$AUC, 3), "\n")
  
  perf_55296 <- data.frame(
    Dataset = "GSE55296",
    Model = c("LASSO-LR", "RandomForest"),
    AUC = c(eval_lr$AUC, eval_rf$AUC),
    Sensitivity = c(eval_lr$Sensitivity, eval_rf$Sensitivity),
    Specificity = c(eval_lr$Specificity, eval_rf$Specificity),
    Accuracy = c(eval_lr$Accuracy, eval_rf$Accuracy),
    stringsAsFactors = F
  )
  
  roc_55296 <- rbind(
    data.frame(Dataset = "GSE55296", Model = "LASSO-LR",
               Specificity = 1 - eval_lr$roc$specificities,
               Sensitivity = eval_lr$roc$sensitivities),
    data.frame(Dataset = "GSE55296", Model = "RandomForest",
               Specificity = 1 - eval_rf$roc$specificities,
               Sensitivity = eval_rf$roc$sensitivities)
  )
}

# 3. GSE42955 Affymetrix 微阵列 -------------------------------------------
cat("\n========== GSE42955 ==========\n")

# 解析 series matrix：定位表格起始行
lines <- readLines(gzfile("GSE42955_series_matrix.txt.gz"))
table_start <- which(grepl("!series_matrix_table_begin", lines))
header_line <- lines[table_start + 1]
# 用 read.delim 读取，跳过 metadata 行
expr <- read.delim(gzfile("GSE42955_series_matrix.txt.gz"),
                   skip = table_start, comment.char = "!",
                   check.names = FALSE, stringsAsFactors = FALSE)
# 第一列 ID_REF
rownames(expr) <- expr[, 1]
expr <- expr[, -1, drop = FALSE]
expr <- as.matrix(expr)

# 探针注释到 Symbol（GPL6244: Affymetrix Human Gene 1.0 ST）
probe_symbols <- mapIds(hugene10sttranscriptcluster.db, keys = rownames(expr), column = "SYMBOL",
                        keytype = "PROBEID", multiVals = "first")
probe_df <- data.frame(probe = rownames(expr), symbol = probe_symbols, stringsAsFactors = F)
probe_df <- probe_df %>% filter(!is.na(symbol) & symbol != "")
expr_sym <- expr[probe_df$probe, , drop = FALSE]
rownames(expr_sym) <- probe_df$symbol

# 去重：同一 symbol 取最大均值探针
mean_expr <- rowMeans(expr_sym, na.rm = TRUE)
expr_sym <- expr_sym[order(mean_expr, decreasing = TRUE), ]
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), ]

cat("GSE42955 genes after annotation:", nrow(expr_sym), "\n")

# 解析表型
meta_lines <- lines[1:table_start]
source_line <- meta_lines[grep("^!Sample_source_name_ch1", meta_lines)]
geo_line <- meta_lines[grep("^!Sample_geo_accession", meta_lines)]

# 提取字段（按 tab 分割，去掉引号）
parse_meta <- function(line) {
  vals <- strsplit(line, "\t")[[1]]
  gsub('"', "", vals[-1])
}
source_vals <- parse_meta(source_line)
geo_vals <- parse_meta(geo_line)

pd42955 <- data.frame(geo_accession = geo_vals,
                      source = source_vals,
                      stringsAsFactors = F)
pd42955$group <- ifelse(grepl("Normal heart", pd42955$source, ignore.case = TRUE), "Normal",
                        ifelse(grepl("Ischemic cardiomyopathy", pd42955$source, ignore.case = TRUE), "Disease", NA))
pd42955 <- pd42955 %>% filter(!is.na(group))
pd42955 <- pd42955[pd42955$geo_accession %in% colnames(expr_sym), ]

cat("GSE42955 samples used:", nrow(pd42955), "\n")
print(table(pd42955$group))

exp42955 <- expr_sym[, pd42955$geo_accession, drop = FALSE]

features42955 <- intersect(stable_genes, rownames(exp42955))
cat("Stable features in GSE42955:", length(features42955), "\n")

perf_42955 <- list()
roc_42955 <- NULL
if (length(features42955) >= 2) {
  X <- t(exp42955[features42955, , drop = FALSE])
  y <- ifelse(pd42955$group == "Disease", 1, 0)
  names(y) <- pd42955$geo_accession
  
  set.seed(123)
  folds <- createFolds(factor(y), k = 5, list = TRUE)
  prob_lr <- prob_rf <- numeric(length(y))
  names(prob_lr) <- names(prob_rf) <- rownames(X)
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    X_train <- scale(X[train_idx, , drop = FALSE])
    X_test <- scale(X[test_idx, , drop = FALSE],
                    center = attr(X_train, "scaled:center"),
                    scale = attr(X_train, "scaled:scale"))
    X_test[is.nan(X_test)] <- 0
    y_train <- y[train_idx]
    
    var_train <- apply(X_train, 2, var, na.rm = TRUE)
    keep <- var_train > 0
    X_train <- X_train[, keep, drop = FALSE]
    X_test <- X_test[, keep, drop = FALSE]
    
    cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1,
                          nfolds = length(y_train), type.measure = "deviance")
    model_lr <- glmnet(X_train, y_train, family = "binomial", alpha = 1,
                       lambda = cv_lasso$lambda.min)
    prob_lr[test_idx] <- as.numeric(predict(model_lr, newx = X_test,
                                            s = cv_lasso$lambda.min, type = "response"))
    
    rf_data <- as.data.frame(X_train)
    rf_data$y <- as.factor(y_train)
    model_rf <- randomForest(y ~ ., data = rf_data, ntree = 500)
    prob_rf[test_idx] <- as.numeric(predict(model_rf, newdata = as.data.frame(X_test), type = "prob")[, 2])
  }
  
  eval_lr <- evaluate_model(prob_lr, y)
  eval_rf <- evaluate_model(prob_rf, y)
  
  cat("GSE42955 LASSO-LR AUC:", round(eval_lr$AUC, 3), "\n")
  cat("GSE42955 RF AUC:", round(eval_rf$AUC, 3), "\n")
  
  perf_42955 <- data.frame(
    Dataset = "GSE42955",
    Model = c("LASSO-LR", "RandomForest"),
    AUC = c(eval_lr$AUC, eval_rf$AUC),
    Sensitivity = c(eval_lr$Sensitivity, eval_rf$Sensitivity),
    Specificity = c(eval_lr$Specificity, eval_rf$Specificity),
    Accuracy = c(eval_lr$Accuracy, eval_rf$Accuracy),
    stringsAsFactors = F
  )
  
  roc_42955 <- rbind(
    data.frame(Dataset = "GSE42955", Model = "LASSO-LR",
               Specificity = 1 - eval_lr$roc$specificities,
               Sensitivity = eval_lr$roc$sensitivities),
    data.frame(Dataset = "GSE42955", Model = "RandomForest",
               Specificity = 1 - eval_rf$roc$specificities,
               Sensitivity = eval_rf$roc$sensitivities)
  )
}

# 4. 汇总结果 -------------------------------------------------------------
perf_all <- bind_rows(perf_55296, perf_42955)
write.table(perf_all, "external_validation_GSE55296_GSE42955_performance.txt",
            sep = "\t", quote = F, row.names = F)
cat("\nSaved: external_validation_GSE55296_GSE42955_performance.txt\n")
print(perf_all)

# 5. ROC 曲线 -------------------------------------------------------------
roc_all <- bind_rows(roc_55296, roc_42955)
if (!is.null(roc_all) && nrow(roc_all) > 0) {
  p <- ggplot(roc_all, aes(x = Specificity, y = Sensitivity, color = Model)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ Dataset, ncol = 2) +
    scale_color_manual(values = c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000")) +
    theme_bw() +
    labs(title = "External validation ROC (GSE55296 & GSE42955)",
         x = "1 - Specificity", y = "Sensitivity", color = "Model")
  ggsave("external_validation_GSE55296_GSE42955_ROC.pdf", p, width = 10, height = 5)
  cat("Saved: external_validation_GSE55296_GSE42955_ROC.pdf\n")
}

# 6. 合并三个外部验证数据集的 ROC 性能摘要（若 GSE116250 结果存在）----
if (file.exists("external_validation_performance.txt")) {
  perf_116250 <- read.table("external_validation_performance.txt", header = T, sep = "\t", stringsAsFactors = F)
  perf_combined <- bind_rows(perf_116250, perf_all)
  write.table(perf_combined, "external_validation_combined_performance.txt",
              sep = "\t", quote = F, row.names = F)
  cat("Saved: external_validation_combined_performance.txt\n")
}

cat("\nDone!\n")
