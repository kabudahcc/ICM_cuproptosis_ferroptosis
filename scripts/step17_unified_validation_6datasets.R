# step17: 统一外部验证（六数据集完整版）
# 复用 step9c 流程：合并发现队列 -> ComBat -> z-score -> LASSO/RF -> 六个外部数据集预测
# 新增 GSE26887 / GSE46224 / GSE52601，补齐 Results 3.6 的六队列数字来源

rm(list = ls())

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
proj_dir <- dirname(script_dir)
data_dir <- file.path(proj_dir, "Data")
tab_dir  <- file.path(proj_dir, "Tables")
fig_dir  <- file.path(proj_dir, "Figures")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr); library(tibble)
  library(glmnet); library(randomForest); library(pROC); library(caret)
  library(PRROC); library(sva); library(GEOquery)
  library(hugene10sttranscriptcluster.db); library(AnnotationDbi)
})

datasets <- c("GSE16499", "GSE5406", "GSE57338")

# ---- 特征池 ----
read_genes <- function(f) if (file.exists(f)) read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1] else character(0)
crg_genes <- unique(unlist(lapply(file.path(data_dir, paste0("CRG_genes_", datasets, ".txt")), read_genes)))
frg_genes <- unique(unlist(lapply(file.path(data_dir, paste0("FRG_genes_", datasets, ".txt")), read_genes)))
feature_pool <- unique(c(crg_genes, frg_genes))
cat("Feature pool:", length(feature_pool), "genes\n")

# ---- 发现队列 ----
prepare_data <- function(ds, feature_genes) {
  exp <- read.csv(file.path(data_dir, paste0(ds, "_symbol.csv")), header = TRUE, row.names = 1, check.names = FALSE)
  group_df <- read.csv(file.path(data_dir, paste0(ds, "_group.csv")), header = TRUE, check.names = FALSE)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  common_samples <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common_samples, drop = FALSE]
  group <- group_df[common_samples, "group"]
  available_features <- intersect(feature_genes, rownames(exp))
  X <- t(exp[available_features, , drop = FALSE])
  rownames(X) <- common_samples
  y <- ifelse(group == "Disease", 1, 0); names(y) <- common_samples
  list(X = X, y = y)
}

# ---- 评估函数（与 step9c 一致） ----
evaluate_model <- function(prob, y_true) {
  roc_obj <- roc(y_true, prob, quiet = TRUE)
  auc <- as.numeric(auc(roc_obj))
  pr_obj <- pr.curve(scores.class0 = prob[y_true == 1], scores.class1 = prob[y_true == 0], curve = TRUE)
  auprc <- pr_obj$auc.integral
  best_coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"),
                        best.method = "youden", transpose = FALSE)
  opt_threshold <- as.numeric(best_coords$threshold[1])
  pred_05 <- ifelse(prob > 0.5, 1, 0)
  cm_05 <- confusionMatrix(factor(pred_05, levels = c(0, 1)), factor(y_true, levels = c(0, 1)), positive = "1")
  pred_opt <- ifelse(prob > opt_threshold, 1, 0)
  cm_opt <- confusionMatrix(factor(pred_opt, levels = c(0, 1)), factor(y_true, levels = c(0, 1)), positive = "1")
  cal_df <- data.frame(prob = prob, y = y_true)
  if (length(unique(prob)) <= 1) {
    cal_summary <- data.frame(bin = "all", n = length(prob), mean_pred = prob[1], obs_rate = mean(y_true))
  } else {
    breaks <- unique(quantile(prob, probs = seq(0, 1, length.out = 11)))
    if (length(breaks) < 3) breaks <- seq(min(prob), max(prob), length.out = 11)
    cal_df$bin <- cut(cal_df$prob, breaks = breaks, include.lowest = TRUE, dig.lab = 4)
    cal_summary <- cal_df %>% group_by(bin) %>%
      summarise(n = n(), mean_pred = mean(prob), obs_rate = mean(y), .groups = "drop") %>% filter(!is.na(bin))
  }
  list(AUC = auc, AUPRC = auprc, OptThreshold = opt_threshold,
       Sensitivity_0.5 = as.numeric(cm_05$byClass["Sensitivity"]),
       Specificity_0.5 = as.numeric(cm_05$byClass["Specificity"]),
       Accuracy_0.5 = as.numeric(cm_05$overall["Accuracy"]),
       Sensitivity_opt = as.numeric(cm_opt$byClass["Sensitivity"]),
       Specificity_opt = as.numeric(cm_opt$byClass["Specificity"]),
       Accuracy_opt = as.numeric(cm_opt$overall["Accuracy"]),
       roc = roc_obj, pr_curve = pr_obj, calibration = cal_summary)
}

# ---- 训练 ----
train_list <- lapply(datasets, prepare_data, feature_genes = feature_pool)
names(train_list) <- datasets
common_train_features <- Reduce(intersect, lapply(train_list, function(x) colnames(x$X)))
cat("Common features across discovery datasets:", length(common_train_features), "\n")

X_train <- do.call(rbind, lapply(train_list, function(x) x$X[, common_train_features, drop = FALSE]))
y_train <- do.call(c, lapply(train_list, function(x) x$y))
batch_train <- do.call(c, lapply(names(train_list), function(ds) rep(ds, nrow(train_list[[ds]]$X))))

combat_out <- ComBat(dat = t(X_train), batch = factor(batch_train),
                     mod = model.matrix(~ y_train), par.prior = TRUE, prior.plots = FALSE)
X_train_scaled <- scale(t(combat_out))
train_center <- attr(X_train_scaled, "scaled:center")
train_scale <- attr(X_train_scaled, "scaled:scale")
keep_features <- apply(X_train_scaled, 2, var, na.rm = TRUE) > 0
X_train_scaled <- X_train_scaled[, keep_features, drop = FALSE]
train_center <- train_center[keep_features]
train_scale <- train_scale[keep_features]

set.seed(123)
cv_lasso <- cv.glmnet(X_train_scaled, y_train, family = "binomial", alpha = 1, nfolds = 10, type.measure = "auc")
coef_lasso <- coef(cv_lasso, s = "lambda.min")
selected <- setdiff(rownames(coef_lasso)[coef_lasso[, 1] != 0], "(Intercept)")
if (length(selected) == 0) selected <- colnames(X_train_scaled)
cat("Final selected genes:", length(selected), "\n")

model_lr <- glmnet(X_train_scaled, y_train, family = "binomial", alpha = 1, lambda = cv_lasso$lambda.min)
rf_data <- as.data.frame(X_train_scaled[, selected, drop = FALSE]); rf_data$y <- as.factor(y_train)
model_rf <- randomForest(y ~ ., data = rf_data, ntree = 500)

# ---- 外部验证函数 ----
validate_external <- function(name, X_ext, y_ext) {
  cat("\n==========", name, "==========\n")
  cat("Samples:", length(y_ext), "Disease:", sum(y_ext == 1), "Normal:", sum(y_ext == 0), "\n")
  train_features <- colnames(X_train_scaled)
  common_features <- intersect(colnames(X_ext), train_features)
  cat("Overlapping features with training:", length(common_features), "/", length(train_features), "\n")
  X_full <- matrix(0, nrow = nrow(X_ext), ncol = length(train_features),
                   dimnames = list(rownames(X_ext), train_features))
  X_full[, common_features] <- X_ext[, common_features, drop = FALSE]
  X_ext_scaled <- scale(X_full, center = train_center, scale = train_scale)
  X_ext_scaled[is.nan(X_ext_scaled)] <- 0
  prob_lr <- as.numeric(predict(model_lr, newx = X_ext_scaled, s = cv_lasso$lambda.min, type = "response"))
  sel_present <- intersect(selected, common_features)
  X_ext_rf_full <- matrix(0, nrow = nrow(X_ext_scaled), ncol = length(selected),
                          dimnames = list(NULL, selected))
  if (length(sel_present) > 0) X_ext_rf_full[, sel_present] <- X_ext_scaled[, sel_present, drop = FALSE]
  prob_rf <- as.numeric(predict(model_rf, newdata = as.data.frame(X_ext_rf_full), type = "prob")[, 2])
  eval_lr <- evaluate_model(prob_lr, y_ext)
  eval_rf <- evaluate_model(prob_rf, y_ext)
  cat("LASSO-LR AUC:", round(eval_lr$AUC, 3), "AUPRC:", round(eval_lr$AUPRC, 3), "\n")
  cat("RF AUC:", round(eval_rf$AUC, 3), "AUPRC:", round(eval_rf$AUPRC, 3), "\n")
  list(Dataset = name, prob_lr = prob_lr, prob_rf = prob_rf, y = y_ext, eval_lr = eval_lr, eval_rf = eval_rf)
}

# ---- 通用：symbol+group csv（已 log2 转换的矩阵） ----
validate_symbol_csv <- function(name) {
  exp <- read.csv(file.path(data_dir, paste0(name, "_symbol.csv")), header = TRUE, row.names = 1, check.names = FALSE)
  group_df <- read.csv(file.path(data_dir, paste0(name, "_group.csv")), header = TRUE, check.names = FALSE)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  common <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common, drop = FALSE]
  if (max(exp, na.rm = TRUE) > 50) exp <- log2(exp + 1)   # 保险：非 log 矩阵才转换
  feats <- intersect(names(train_center), rownames(exp))
  X <- t(exp[feats, , drop = FALSE])
  y <- ifelse(group_df[common, "group"] == "Disease", 1, 0)
  validate_external(name, X, y)
}

# ---- GSE116250 (RNA-seq RPKM) ----
exp116250 <- read.table(gzfile(file.path(data_dir, "GSE116250_rpkm.txt.gz")), header = TRUE, sep = "\t", check.names = FALSE, quote = "")
exp116250 <- exp116250 %>% filter(Common_name != "")
exp116250 <- exp116250[!duplicated(exp116250$Common_name), ]
rownames(exp116250) <- exp116250$Common_name
exp116250 <- as.matrix(exp116250[, -(1:2)])
gse <- getGEO(filename = file.path(data_dir, "GSE116250_series_matrix.txt.gz"), getGPL = FALSE)
pd <- pData(phenoData(gse))
disease_col <- grep("disease", names(pd), value = TRUE)[1]
pd$disease <- pd[[disease_col]]
pd <- pd %>% filter(disease %in% c("non-failing", "ischemic cardiomyopathy"))
pd$group <- ifelse(pd$disease == "non-failing", "Normal", "Disease")
title_to_geo <- data.frame(title = pd$title, geo_accession = rownames(pd), group = pd$group, stringsAsFactors = FALSE)
colnames(exp116250) <- trimws(colnames(exp116250))
matched <- title_to_geo %>% filter(title %in% colnames(exp116250))
exp_val <- exp116250[, matched$title, drop = FALSE]
exp_log <- log2(exp_val + 1)
X116250 <- t(exp_log[intersect(names(train_center), rownames(exp_log)), , drop = FALSE])
y116250 <- ifelse(matched$group == "Disease", 1, 0)
res116250 <- validate_external("GSE116250", X116250, y116250)

# ---- GSE55296 (RNA-seq counts) ----
counts <- read.delim(gzfile(file.path(data_dir, "GSE55296_count_data.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts)[1:2] <- c("Ensembl", "Symbol")
colnames(counts) <- make.names(colnames(counts), unique = TRUE)
counts <- counts[, colSums(!is.na(counts) & counts != "") > 0]
counts <- counts %>% filter(Symbol != "") %>% filter(!duplicated(Symbol))
rownames(counts) <- counts$Symbol
exp_counts <- as.matrix(counts[, -(1:2)]); storage.mode(exp_counts) <- "numeric"
readme <- read.delim(file.path(data_dir, "GSE55296_processed_data_readme.txt"), check.names = FALSE, stringsAsFactors = FALSE)
colnames(readme) <- c("geo_accession", "title", "col_header")
readme$group <- ifelse(grepl("Control", readme$title, ignore.case = TRUE), "Normal",
                       ifelse(grepl("Ischemic", readme$title, ignore.case = TRUE), "Disease", NA))
readme <- readme %>% filter(!is.na(group)) %>% filter(col_header %in% colnames(exp_counts))
exp55296_log <- log2(exp_counts[, readme$col_header, drop = FALSE] + 1)
X55296 <- t(exp55296_log[intersect(names(train_center), rownames(exp55296_log)), , drop = FALSE])
y55296 <- ifelse(readme$group == "Disease", 1, 0)
res55296 <- validate_external("GSE55296", X55296, y55296)

# ---- GSE42955 (microarray, 探针注释) ----
lines <- readLines(gzfile(file.path(data_dir, "GSE42955_series_matrix.txt.gz")))
table_start <- which(grepl("!series_matrix_table_begin", lines))
expr <- read.delim(gzfile(file.path(data_dir, "GSE42955_series_matrix.txt.gz")),
                   skip = table_start, comment.char = "!", check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr) <- expr[, 1]; expr <- as.matrix(expr[, -1, drop = FALSE])
probe_symbols <- mapIds(hugene10sttranscriptcluster.db, keys = rownames(expr), column = "SYMBOL",
                        keytype = "PROBEID", multiVals = "first")
probe_df <- data.frame(probe = rownames(expr), symbol = probe_symbols, stringsAsFactors = FALSE) %>%
  filter(!is.na(symbol) & symbol != "")
expr_sym <- expr[probe_df$probe, , drop = FALSE]; rownames(expr_sym) <- probe_df$symbol
expr_sym <- expr_sym[order(rowMeans(expr_sym, na.rm = TRUE), decreasing = TRUE), ]
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), ]
parse_meta <- function(line) gsub('"', "", strsplit(line, "\t")[[1]][-1])
source_vals <- parse_meta(grep("^!Sample_source_name_ch1", lines[1:table_start], value = TRUE))
geo_vals <- parse_meta(grep("^!Sample_geo_accession", lines[1:table_start], value = TRUE))
pd42955 <- data.frame(geo_accession = geo_vals, source = source_vals, stringsAsFactors = FALSE)
pd42955$group <- ifelse(grepl("Normal heart", pd42955$source, ignore.case = TRUE), "Normal",
                        ifelse(grepl("Ischemic cardiomyopathy", pd42955$source, ignore.case = TRUE), "Disease", NA))
pd42955 <- pd42955 %>% filter(!is.na(group)) %>% filter(geo_accession %in% colnames(expr_sym))
exp42955 <- expr_sym[, pd42955$geo_accession, drop = FALSE]
X42955 <- t(exp42955[intersect(names(train_center), rownames(exp42955)), , drop = FALSE])
y42955 <- ifelse(pd42955$group == "Disease", 1, 0)
res42955 <- validate_external("GSE42955", X42955, y42955)

# ---- GSE26887 / GSE46224 / GSE52601 (symbol+group csv) ----
res26887 <- validate_symbol_csv("GSE26887")
res46224 <- validate_symbol_csv("GSE46224")
res52601 <- validate_symbol_csv("GSE52601")

# ---- 汇总 ----
results <- list(res116250, res55296, res42955, res26887, res46224, res52601)
results <- results[!sapply(results, is.null)]
saveRDS(results, file.path(data_dir, "unified6_validation_results.rds"))

metric <- function(r, f) c(r$eval_lr[[f]], r$eval_rf[[f]])
perf_df <- bind_rows(lapply(results, function(r) data.frame(
  Dataset = r$Dataset, Model = c("LASSO-LR", "RandomForest"),
  AUC = metric(r, "AUC"), AUPRC = metric(r, "AUPRC"),
  Sensitivity_0.5 = metric(r, "Sensitivity_0.5"), Specificity_0.5 = metric(r, "Specificity_0.5"),
  Accuracy_0.5 = metric(r, "Accuracy_0.5"),
  Sensitivity_opt = metric(r, "Sensitivity_opt"), Specificity_opt = metric(r, "Specificity_opt"),
  Accuracy_opt = metric(r, "Accuracy_opt"),
  Threshold_opt = metric(r, "OptThreshold"), stringsAsFactors = FALSE)))
write.table(perf_df, file.path(tab_dir, "external_validation_unified6_performance.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n"); print(perf_df)

# ---- 六分面图 ----
model_colors <- c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000")
ds_levels <- c("GSE116250", "GSE55296", "GSE42955", "GSE26887", "GSE46224", "GSE52601")

roc_df <- bind_rows(lapply(results, function(r) rbind(
  data.frame(Dataset = r$Dataset, Model = "LASSO-LR", Specificity = 1 - r$eval_lr$roc$specificities, Sensitivity = r$eval_lr$roc$sensitivities),
  data.frame(Dataset = r$Dataset, Model = "RandomForest", Specificity = 1 - r$eval_rf$roc$specificities, Sensitivity = r$eval_rf$roc$sensitivities))))
roc_df$Model <- factor(roc_df$Model, levels = c("LASSO-LR", "RandomForest"))
roc_df$Dataset <- factor(roc_df$Dataset, levels = ds_levels)
p_roc <- ggplot(roc_df, aes(Specificity, Sensitivity, color = Model)) +
  geom_line(linewidth = 0.9) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Dataset, ncol = 3) + scale_color_manual(values = model_colors) + theme_bw() +
  labs(title = "Unified external validation ROC", x = "1 - Specificity", y = "Sensitivity", color = "Model")
ggsave(file.path(fig_dir, "external_validation_unified6_ROC.pdf"), p_roc, width = 11, height = 7)

pr_df <- bind_rows(lapply(results, function(r) rbind(
  data.frame(Dataset = r$Dataset, Model = "LASSO-LR", Recall = r$eval_lr$pr_curve$curve[, 1], Precision = r$eval_lr$pr_curve$curve[, 2]),
  data.frame(Dataset = r$Dataset, Model = "RandomForest", Recall = r$eval_rf$pr_curve$curve[, 1], Precision = r$eval_rf$pr_curve$curve[, 2]))))
pr_df$Model <- factor(pr_df$Model, levels = c("LASSO-LR", "RandomForest"))
pr_df$Dataset <- factor(pr_df$Dataset, levels = ds_levels)
p_pr <- ggplot(pr_df, aes(Recall, Precision, color = Model)) +
  geom_line(linewidth = 0.9) + facet_wrap(~ Dataset, ncol = 3) +
  scale_color_manual(values = model_colors) + theme_bw() +
  labs(title = "Unified external validation PR curve", x = "Recall", y = "Precision", color = "Model")
ggsave(file.path(fig_dir, "external_validation_unified6_PR_curve.pdf"), p_pr, width = 11, height = 7)

cal_df <- bind_rows(lapply(results, function(r) rbind(
  data.frame(Dataset = r$Dataset, Model = "LASSO-LR", r$eval_lr$calibration),
  data.frame(Dataset = r$Dataset, Model = "RandomForest", r$eval_rf$calibration))))
cal_df$Model <- factor(cal_df$Model, levels = c("LASSO-LR", "RandomForest"))
cal_df$Dataset <- factor(cal_df$Dataset, levels = ds_levels)
p_cal <- ggplot(cal_df, aes(mean_pred, obs_rate, color = Model)) +
  geom_point(aes(size = n), alpha = 0.7) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Dataset, ncol = 3) + scale_color_manual(values = model_colors) + theme_bw() +
  labs(title = "Unified external validation calibration plot",
       x = "Mean predicted probability", y = "Observed event rate", color = "Model", size = "n") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
ggsave(file.path(fig_dir, "external_validation_unified6_calibration.pdf"), p_cal, width = 11, height = 7)

cat("\nUnified 6-dataset external validation done!\n")
