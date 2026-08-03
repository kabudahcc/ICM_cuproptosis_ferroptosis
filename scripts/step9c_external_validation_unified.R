# 统一外部验证：用全部三个发现队列训练模型，独立验证 GSE116250 / GSE55296 / GSE42955
# 流程：合并发现队列 -> ComBat 批次校正 -> z-score -> LASSO/RF 训练 -> 外部数据集预测
# 新增评估：AUC、AUPRC、混淆矩阵（0.5 与 Youden 阈值）、校准曲线

rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(glmnet)
library(randomForest)
library(pROC)
library(caret)
library(PRROC)
library(sva)
library(GEOquery)
library(limma)
library(edgeR)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)

datasets <- c("GSE16499", "GSE5406", "GSE57338")

# 读取铜凋亡 + 铁凋亡特征池
CRG_files <- paste0("CRG_genes_", datasets, ".txt")
FRG_files <- paste0("FRG_genes_", datasets, ".txt")
crg_genes <- unique(unlist(lapply(CRG_files, function(f) {
  if (!file.exists(f)) return(character(0))
  read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1]
})))
frg_genes <- unique(unlist(lapply(FRG_files, function(f) {
  if (!file.exists(f)) return(character(0))
  read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1]
})))
feature_pool <- unique(c(crg_genes, frg_genes))
cat("Feature pool:", length(feature_pool), "genes\n")

# 发现队列数据读取
prepare_data <- function(ds, feature_genes) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  group_df <- read.csv(paste0(ds, "_group.csv"), header = TRUE, check.names = FALSE)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  
  common_samples <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common_samples, drop = FALSE]
  group <- group_df[common_samples, "group"]
  
  available_features <- intersect(feature_genes, rownames(exp))
  X <- t(exp[available_features, , drop = FALSE])
  rownames(X) <- common_samples
  colnames(X) <- available_features
  
  y <- ifelse(group == "Disease", 1, 0)
  names(y) <- common_samples
  
  list(X = X, y = y, group = group, features = available_features)
}

# 评估函数
evaluate_model <- function(prob, y_true) {
  roc_obj <- roc(y_true, prob, quiet = TRUE)
  auc <- as.numeric(auc(roc_obj))
  
  pr_obj <- pr.curve(scores.class0 = prob[y_true == 1],
                     scores.class1 = prob[y_true == 0],
                     curve = TRUE)
  auprc <- pr_obj$auc.integral
  
  best_coords <- coords(roc_obj, "best",
                        ret = c("threshold", "specificity", "sensitivity"),
                        best.method = "youden", transpose = FALSE)
  opt_threshold <- as.numeric(best_coords$threshold[1])
  
  pred_05 <- ifelse(prob > 0.5, 1, 0)
  cm_05 <- confusionMatrix(factor(pred_05, levels = c(0, 1)),
                           factor(y_true, levels = c(0, 1)), positive = "1")
  pred_opt <- ifelse(prob > opt_threshold, 1, 0)
  cm_opt <- confusionMatrix(factor(pred_opt, levels = c(0, 1)),
                            factor(y_true, levels = c(0, 1)), positive = "1")
  
  cal_df <- data.frame(prob = prob, y = y_true)
  if (length(unique(prob)) <= 1) {
    cal_summary <- data.frame(bin = "all", n = length(prob), mean_pred = prob[1], obs_rate = mean(y_true), stringsAsFactors = FALSE)
  } else {
    breaks <- unique(quantile(prob, probs = seq(0, 1, length.out = 11)))
    if (length(breaks) < 3) breaks <- seq(min(prob), max(prob), length.out = 11)
    cal_df$bin <- cut(cal_df$prob, breaks = breaks, include.lowest = TRUE, dig.lab = 4)
    cal_summary <- cal_df %>%
      group_by(bin) %>%
      summarise(n = n(), mean_pred = mean(prob), obs_rate = mean(y), .groups = "drop") %>%
      filter(!is.na(bin))
  }
  
  list(
    AUC = auc, AUPRC = auprc, OptThreshold = opt_threshold,
    Sensitivity_0.5 = as.numeric(cm_05$byClass["Sensitivity"]),
    Specificity_0.5 = as.numeric(cm_05$byClass["Specificity"]),
    Accuracy_0.5 = as.numeric(cm_05$overall["Accuracy"]),
    Kappa_0.5 = as.numeric(cm_05$overall["Kappa"]),
    PPV_0.5 = as.numeric(cm_05$byClass["Pos Pred Value"]),
    NPV_0.5 = as.numeric(cm_05$byClass["Neg Pred Value"]),
    Sensitivity_opt = as.numeric(cm_opt$byClass["Sensitivity"]),
    Specificity_opt = as.numeric(cm_opt$byClass["Specificity"]),
    Accuracy_opt = as.numeric(cm_opt$overall["Accuracy"]),
    Kappa_opt = as.numeric(cm_opt$overall["Kappa"]),
    PPV_opt = as.numeric(cm_opt$byClass["Pos Pred Value"]),
    NPV_opt = as.numeric(cm_opt$byClass["Neg Pred Value"]),
    roc = roc_obj, pr_curve = pr_obj, calibration = cal_summary
  )
}

# 1. 构建发现队列训练集 -----------------------------------------------------
train_list <- lapply(datasets, prepare_data, feature_genes = feature_pool)
names(train_list) <- datasets
common_train_features <- Reduce(intersect, lapply(train_list, function(x) colnames(x$X)))
cat("Common features across discovery datasets:", length(common_train_features), "\n")

X_train <- do.call(rbind, lapply(train_list, function(x) x$X[, common_train_features, drop = FALSE]))
y_train <- do.call(c, lapply(train_list, function(x) x$y))
batch_train <- do.call(c, lapply(names(train_list), function(ds) rep(ds, nrow(train_list[[ds]]$X))))

# ComBat 批次校正
combat_input <- t(X_train)
mod <- model.matrix(~ y_train)
combat_out <- ComBat(dat = combat_input, batch = factor(batch_train), mod = mod, par.prior = TRUE, prior.plots = FALSE)
X_train_corrected <- t(combat_out)

# z-score（保存参数供外部验证使用）
X_train_scaled <- scale(X_train_corrected)
train_center <- attr(X_train_scaled, "scaled:center")
train_scale <- attr(X_train_scaled, "scaled:scale")

# 移除方差为 0 的特征
var_train <- apply(X_train_scaled, 2, var, na.rm = TRUE)
keep_features <- var_train > 0
X_train_scaled <- X_train_scaled[, keep_features, drop = FALSE]
train_center <- train_center[keep_features]
train_scale <- train_scale[keep_features]

# 2. 训练最终模型 ------------------------------------------------------------
set.seed(123)
cv_lasso <- cv.glmnet(X_train_scaled, y_train, family = "binomial", alpha = 1,
                      nfolds = 10, type.measure = "auc")
coef_lasso <- coef(cv_lasso, s = "lambda.min")
selected <- setdiff(rownames(coef_lasso)[coef_lasso[, 1] != 0], "(Intercept)")
if (length(selected) == 0) selected <- colnames(X_train_scaled)
cat("Final selected genes:", length(selected), "\n")

model_lr <- glmnet(X_train_scaled, y_train, family = "binomial", alpha = 1,
                   lambda = cv_lasso$lambda.min)

rf_data <- as.data.frame(X_train_scaled[, selected, drop = FALSE])
rf_data$y <- as.factor(y_train)
model_rf <- randomForest(y ~ ., data = rf_data, ntree = 500)

# 3. 外部验证数据集处理与预测函数 ------------------------------------------
validate_external <- function(name, X_ext, y_ext) {
  cat("\n==========", name, "==========\n")
  cat("Samples:", length(y_ext), "Disease:", sum(y_ext == 1), "Normal:", sum(y_ext == 0), "\n")
  
  train_features <- colnames(X_train_scaled)
  common_features <- intersect(colnames(X_ext), train_features)
  cat("Overlapping features with training:", length(common_features), "/", length(train_features), "\n")
  if (length(common_features) < 2) {
    cat("Too few overlapping features, skipping\n")
    return(NULL)
  }
  
  # 构建与训练集同维度的矩阵，缺失特征用 0（训练集均值）填充
  X_full <- matrix(0, nrow = nrow(X_ext), ncol = length(train_features))
  colnames(X_full) <- train_features
  rownames(X_full) <- rownames(X_ext)
  X_full[, common_features] <- X_ext[, common_features, drop = FALSE]
  
  X_ext_scaled <- scale(X_full, center = train_center, scale = train_scale)
  X_ext_scaled[is.nan(X_ext_scaled)] <- 0
  
  # LR 使用全部训练特征
  prob_lr <- as.numeric(predict(model_lr, newx = X_ext_scaled,
                                s = cv_lasso$lambda.min, type = "response"))
  # RF 需要与训练时完全一致的列；缺失列用 0 填充
  sel_present <- intersect(selected, common_features)
  X_ext_rf_full <- matrix(0, nrow = nrow(X_ext_scaled), ncol = length(selected))
  colnames(X_ext_rf_full) <- selected
  if (length(sel_present) > 0) {
    X_ext_rf_full[, sel_present] <- X_ext_scaled[, sel_present, drop = FALSE]
  }
  X_ext_sel <- as.data.frame(X_ext_rf_full)
  prob_rf <- as.numeric(predict(model_rf, newdata = X_ext_sel, type = "prob")[, 2])
  
  eval_lr <- evaluate_model(prob_lr, y_ext)
  eval_rf <- evaluate_model(prob_rf, y_ext)
  
  cat("LASSO-LR AUC:", round(eval_lr$AUC, 3), "AUPRC:", round(eval_lr$AUPRC, 3), "\n")
  cat("RF AUC:", round(eval_rf$AUC, 3), "AUPRC:", round(eval_rf$AUPRC, 3), "\n")
  
  list(
    Dataset = name,
    prob_lr = prob_lr, prob_rf = prob_rf, y = y_ext,
    eval_lr = eval_lr, eval_rf = eval_rf
  )
}

# 4. GSE116250 ---------------------------------------------------------------
exp116250 <- read.table("GSE116250_rpkm.txt.gz", header = TRUE, sep = "\t", check.names = FALSE, quote = "")
exp116250 <- exp116250 %>% filter(Common_name != "")
exp116250 <- exp116250[!duplicated(exp116250$Common_name), ]
rownames(exp116250) <- exp116250$Common_name
exp116250 <- exp116250[, -(1:2)]
exp116250 <- as.matrix(exp116250)

gse <- getGEO("GSE116250", destdir = ".", getGPL = FALSE)
pd <- pData(phenoData(gse[[1]]))
disease_col <- grep("disease", names(pd), value = TRUE)
pd$disease <- pd[[disease_col]]
pd <- pd %>% filter(disease %in% c("non-failing", "ischemic cardiomyopathy"))
pd$group <- ifelse(pd$disease == "non-failing", "Normal", "Disease")
pd$geo_accession <- rownames(pd)

title_to_geo <- data.frame(title = pd$title, geo_accession = rownames(pd), group = pd$group, stringsAsFactors = FALSE)
colnames(exp116250) <- trimws(colnames(exp116250))
matched <- title_to_geo %>% filter(title %in% colnames(exp116250))
exp_val <- exp116250[, matched$title, drop = FALSE]
colnames(exp_val) <- matched$geo_accession
group_val <- factor(matched$group, levels = c("Normal", "Disease"))

exp_log <- log2(exp_val + 1)
features116250 <- intersect(names(train_center), rownames(exp_log))
X116250 <- t(exp_log[features116250, , drop = FALSE])
y116250 <- ifelse(group_val == "Disease", 1, 0)
res116250 <- validate_external("GSE116250", X116250, y116250)

# 5. GSE55296 ----------------------------------------------------------------
counts <- read.delim(gzfile("GSE55296_count_data.txt.gz"), check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts)[1:2] <- c("Ensembl", "Symbol")
colnames(counts) <- make.names(colnames(counts), unique = TRUE)
counts <- counts[, colSums(!is.na(counts) & counts != "") > 0]
counts <- counts %>% filter(Symbol != "") %>% filter(!duplicated(Symbol))
rownames(counts) <- counts$Symbol
exp_counts <- as.matrix(counts[, -(1:2)])
storage.mode(exp_counts) <- "numeric"

readme <- read.delim("GSE55296_processed_data_readme.txt", check.names = FALSE, stringsAsFactors = FALSE)
colnames(readme) <- c("geo_accession", "title", "col_header")
readme$group <- ifelse(grepl("Control", readme$title, ignore.case = TRUE), "Normal",
                       ifelse(grepl("Ischemic", readme$title, ignore.case = TRUE), "Disease", NA))
readme <- readme %>% filter(!is.na(group))
readme <- readme[readme$col_header %in% colnames(exp_counts), ]

exp55296 <- exp_counts[, readme$col_header, drop = FALSE]
colnames(exp55296) <- readme$geo_accession
exp55296_log <- log2(exp55296 + 1)

features55296 <- intersect(names(train_center), rownames(exp55296_log))
X55296 <- t(exp55296_log[features55296, , drop = FALSE])
y55296 <- ifelse(readme$group == "Disease", 1, 0)
res55296 <- validate_external("GSE55296", X55296, y55296)

# 6. GSE42955 ----------------------------------------------------------------
lines <- readLines(gzfile("GSE42955_series_matrix.txt.gz"))
table_start <- which(grepl("!series_matrix_table_begin", lines))
expr <- read.delim(gzfile("GSE42955_series_matrix.txt.gz"),
                   skip = table_start, comment.char = "!",
                   check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr) <- expr[, 1]
expr <- expr[, -1, drop = FALSE]
expr <- as.matrix(expr)

probe_symbols <- mapIds(hugene10sttranscriptcluster.db, keys = rownames(expr), column = "SYMBOL",
                        keytype = "PROBEID", multiVals = "first")
probe_df <- data.frame(probe = rownames(expr), symbol = probe_symbols, stringsAsFactors = FALSE)
probe_df <- probe_df %>% filter(!is.na(symbol) & symbol != "")
expr_sym <- expr[probe_df$probe, , drop = FALSE]
rownames(expr_sym) <- probe_df$symbol
mean_expr <- rowMeans(expr_sym, na.rm = TRUE)
expr_sym <- expr_sym[order(mean_expr, decreasing = TRUE), ]
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), ]

meta_lines <- lines[1:table_start]
source_line <- meta_lines[grep("^!Sample_source_name_ch1", meta_lines)]
geo_line <- meta_lines[grep("^!Sample_geo_accession", meta_lines)]
parse_meta <- function(line) {
  vals <- strsplit(line, "\t")[[1]]
  gsub('"', "", vals[-1])
}
source_vals <- parse_meta(source_line)
geo_vals <- parse_meta(geo_line)
pd42955 <- data.frame(geo_accession = geo_vals, source = source_vals, stringsAsFactors = FALSE)
pd42955$group <- ifelse(grepl("Normal heart", pd42955$source, ignore.case = TRUE), "Normal",
                        ifelse(grepl("Ischemic cardiomyopathy", pd42955$source, ignore.case = TRUE), "Disease", NA))
pd42955 <- pd42955 %>% filter(!is.na(group))
pd42955 <- pd42955[pd42955$geo_accession %in% colnames(expr_sym), ]

exp42955 <- expr_sym[, pd42955$geo_accession, drop = FALSE]
features42955 <- intersect(names(train_center), rownames(exp42955))
X42955 <- t(exp42955[features42955, , drop = FALSE])
y42955 <- ifelse(pd42955$group == "Disease", 1, 0)
res42955 <- validate_external("GSE42955", X42955, y42955)

# 7. 汇总结果 ----------------------------------------------------------------
results <- list(res116250, res55296, res42955)
results <- results[!sapply(results, is.null)]

perf_rows <- lapply(results, function(r) {
  data.frame(
    Dataset = r$Dataset,
    Model = c("LASSO-LR", "RandomForest"),
    AUC = c(r$eval_lr$AUC, r$eval_rf$AUC),
    AUPRC = c(r$eval_lr$AUPRC, r$eval_rf$AUPRC),
    Sensitivity_0.5 = c(r$eval_lr$Sensitivity_0.5, r$eval_rf$Sensitivity_0.5),
    Specificity_0.5 = c(r$eval_lr$Specificity_0.5, r$eval_rf$Specificity_0.5),
    Accuracy_0.5 = c(r$eval_lr$Accuracy_0.5, r$eval_rf$Accuracy_0.5),
    Kappa_0.5 = c(r$eval_lr$Kappa_0.5, r$eval_rf$Kappa_0.5),
    PPV_0.5 = c(r$eval_lr$PPV_0.5, r$eval_rf$PPV_0.5),
    NPV_0.5 = c(r$eval_lr$NPV_0.5, r$eval_rf$NPV_0.5),
    Sensitivity_opt = c(r$eval_lr$Sensitivity_opt, r$eval_rf$Sensitivity_opt),
    Specificity_opt = c(r$eval_lr$Specificity_opt, r$eval_rf$Specificity_opt),
    Accuracy_opt = c(r$eval_lr$Accuracy_opt, r$eval_rf$Accuracy_opt),
    Kappa_opt = c(r$eval_lr$Kappa_opt, r$eval_rf$Kappa_opt),
    PPV_opt = c(r$eval_lr$PPV_opt, r$eval_rf$PPV_opt),
    NPV_opt = c(r$eval_lr$NPV_opt, r$eval_rf$NPV_opt),
    Threshold_opt = c(r$eval_lr$OptThreshold, r$eval_rf$OptThreshold),
    stringsAsFactors = FALSE
  )
})
perf_df <- bind_rows(perf_rows)
write.table(perf_df, "external_validation_unified_performance.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSaved: external_validation_unified_performance.txt\n")
print(perf_df)

# 8. 绘图 --------------------------------------------------------------------
# ROC
roc_df <- bind_rows(lapply(results, function(r) {
  rbind(
    data.frame(Dataset = r$Dataset, Model = "LASSO-LR",
               Specificity = 1 - r$eval_lr$roc$specificities,
               Sensitivity = r$eval_lr$roc$sensitivities),
    data.frame(Dataset = r$Dataset, Model = "RandomForest",
               Specificity = 1 - r$eval_rf$roc$specificities,
               Sensitivity = r$eval_rf$roc$sensitivities)
  )
}))
roc_df$Model <- factor(roc_df$Model, levels = c("LASSO-LR", "RandomForest"))
p_roc <- ggplot(roc_df, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Dataset, ncol = 3) +
  scale_color_manual(values = c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000")) +
  theme_bw() +
  labs(title = "Unified external validation ROC",
       x = "1 - Specificity", y = "Sensitivity", color = "Model")
ggsave("external_validation_unified_ROC.pdf", p_roc, width = 12, height = 4)

# PR
pr_df <- bind_rows(lapply(results, function(r) {
  rbind(
    data.frame(Dataset = r$Dataset, Model = "LASSO-LR",
               Recall = r$eval_lr$pr_curve$curve[, 1],
               Precision = r$eval_lr$pr_curve$curve[, 2]),
    data.frame(Dataset = r$Dataset, Model = "RandomForest",
               Recall = r$eval_rf$pr_curve$curve[, 1],
               Precision = r$eval_rf$pr_curve$curve[, 2])
  )
}))
pr_df$Model <- factor(pr_df$Model, levels = c("LASSO-LR", "RandomForest"))
p_pr <- ggplot(pr_df, aes(x = Recall, y = Precision, color = Model)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Dataset, ncol = 3) +
  scale_color_manual(values = c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000")) +
  theme_bw() +
  labs(title = "Unified external validation PR curve",
       x = "Recall", y = "Precision", color = "Model")
ggsave("external_validation_unified_PR_curve.pdf", p_pr, width = 12, height = 4)

# Calibration
cal_df <- bind_rows(lapply(results, function(r) {
  rbind(
    data.frame(Dataset = r$Dataset, Model = "LASSO-LR", r$eval_lr$calibration),
    data.frame(Dataset = r$Dataset, Model = "RandomForest", r$eval_rf$calibration)
  )
}))
cal_df$Model <- factor(cal_df$Model, levels = c("LASSO-LR", "RandomForest"))
p_cal <- ggplot(cal_df, aes(x = mean_pred, y = obs_rate, color = Model)) +
  geom_point(aes(size = n), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Dataset, ncol = 3) +
  scale_color_manual(values = c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000")) +
  theme_bw() +
  labs(title = "Unified external validation calibration plot",
       x = "Mean predicted probability", y = "Observed event rate",
       color = "Model", size = "n") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
ggsave("external_validation_unified_calibration.pdf", p_cal, width = 12, height = 4)

cat("Saved: external_validation_unified_ROC.pdf\n")
cat("Saved: external_validation_unified_PR_curve.pdf\n")
cat("Saved: external_validation_unified_calibration.pdf\n")
cat("\nUnified external validation done!\n")
