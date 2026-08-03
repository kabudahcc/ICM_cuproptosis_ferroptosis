# 跨数据集 Meta-analysis：整合多个 ICM vs Normal 数据集发现稳健 DEG
# 方法 1：ComBat 批次校正 + limma
# 方法 2：RankProd（基于每个数据集内排序的跨数据集一致性分析）

rm(list = ls())
library(dplyr)
library(ggplot2)
library(ggrepel)
library(limma)
library(sva)
library(RankProd)
library(pheatmap)
library(GEOquery)
library(edgeR)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)

# ---------- 1. 数据集加载函数 ----------
load_csv <- function(ds) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  group_df <- read.csv(paste0(ds, "_group.csv"), header = TRUE, check.names = FALSE)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  common <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common, drop = FALSE]
  group <- group_df[common, "group"]
  list(expr = as.matrix(exp), group = group)
}

load_GSE116250 <- function() {
  exp <- read.table("GSE116250_rpkm.txt.gz", header = TRUE, sep = "\t", check.names = FALSE, quote = "")
  exp <- exp %>% filter(Common_name != "")
  exp <- exp[!duplicated(exp$Common_name), ]
  rownames(exp) <- exp$Common_name
  exp <- exp[, -(1:2)]
  exp <- as.matrix(exp)
  
  gse <- getGEO("GSE116250", destdir = ".", getGPL = FALSE)
  pd <- pData(phenoData(gse[[1]]))
  disease_col <- grep("disease", names(pd), value = TRUE)
  pd$disease <- pd[[disease_col]]
  pd <- pd %>% filter(disease %in% c("non-failing", "ischemic cardiomyopathy"))
  pd$group <- ifelse(pd$disease == "non-failing", "Normal", "Disease")
  
  title_to_geo <- data.frame(title = pd$title, geo_accession = rownames(pd), group = pd$group, stringsAsFactors = FALSE)
  colnames(exp) <- trimws(colnames(exp))
  matched <- title_to_geo %>% filter(title %in% colnames(exp))
  expr <- exp[, matched$title, drop = FALSE]
  colnames(expr) <- matched$geo_accession
  list(expr = log2(expr + 1), group = matched$group)
}

load_GSE55296 <- function() {
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
  
  expr <- exp_counts[, readme$col_header, drop = FALSE]
  colnames(expr) <- readme$geo_accession
  list(expr = log2(expr + 1), group = readme$group)
}

load_GSE42955 <- function() {
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
  pd <- data.frame(geo_accession = geo_vals, source = source_vals, stringsAsFactors = FALSE)
  pd$group <- ifelse(grepl("Normal heart", pd$source, ignore.case = TRUE), "Normal",
                     ifelse(grepl("Ischemic cardiomyopathy", pd$source, ignore.case = TRUE), "Disease", NA))
  pd <- pd %>% filter(!is.na(group))
  pd <- pd[pd$geo_accession %in% colnames(expr_sym), ]
  
  expr <- expr_sym[, pd$geo_accession, drop = FALSE]
  list(expr = expr, group = pd$group)
}

# ---------- 2. 读取所有数据集 ----------
datasets <- list(
  GSE16499 = function() load_csv("GSE16499"),
  GSE5406 = function() load_csv("GSE5406"),
  GSE57338 = function() load_csv("GSE57338"),
  GSE116250 = load_GSE116250,
  GSE55296 = load_GSE55296,
  GSE42955 = load_GSE42955,
  GSE26887 = function() load_csv("GSE26887"),
  GSE46224 = function() load_csv("GSE46224")
)

dat_list <- list()
for (ds in names(datasets)) {
  cat("\nLoading", ds, "...\n")
  dat <- datasets[[ds]]()
  cat("  samples:", length(dat$group), "Disease:", sum(dat$group == "Disease"), "Normal:", sum(dat$group == "Normal"), "\n")
  cat("  genes:", nrow(dat$expr), "\n")
  # 去除 NA/Inf
  dat$expr[is.na(dat$expr)] <- 0
  dat$expr[is.infinite(dat$expr)] <- 0
  dat_list[[ds]] <- dat
}

# 取共同基因
common_genes <- Reduce(intersect, lapply(dat_list, function(x) rownames(x$expr)))
cat("\nCommon genes across all datasets:", length(common_genes), "\n")

# 构建合并表达矩阵
expr_all_list <- lapply(dat_list, function(x) x$expr[common_genes, , drop = FALSE])
expr_all <- do.call(cbind, expr_all_list)

batch_all <- do.call(c, lapply(names(dat_list), function(ds) rep(ds, ncol(dat_list[[ds]]$expr))))
group_all <- do.call(c, lapply(dat_list, function(x) x$group))
cl_all <- ifelse(group_all == "Disease", 1, 0)

# 过滤表达量过低的基因（至少在 20% 样本中 > 0）
keep_expr <- rowSums(expr_all > 0) > 0.2 * ncol(expr_all)
expr_all <- expr_all[keep_expr, ]
common_genes <- rownames(expr_all)
cat("Genes after filtering:", length(common_genes), "\n")

# ---------- 3. ComBat + limma Meta-analysis ----------
cat("\nRunning ComBat + limma meta-analysis...\n")
combat_out <- ComBat(dat = expr_all, batch = factor(batch_all),
                     mod = model.matrix(~ cl_all), par.prior = TRUE, prior.plots = FALSE)

# limma：校正批次后只保留疾病效应
design <- model.matrix(~ cl_all)
fit <- lmFit(combat_out, design)
fit <- eBayes(fit)
limma_res <- topTable(fit, coef = "cl_all", number = Inf, adjust.method = "BH", sort.by = "P")
limma_res$symbol <- rownames(limma_res)
limma_res <- limma_res %>% dplyr::select(symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
write.table(limma_res, "meta_analysis_DEG_limma.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: meta_analysis_DEG_limma.txt\n")

# ---------- 4. RankProd Meta-analysis ----------
cat("\nRunning RankProd meta-analysis...\n")

# 辅助：基于置换的跨数据集 rank product（当 RPadvance 因内存/大数失败时作为稳健回退）
run_manual_rankprod <- function(expr, group, batch, nperm = 200, seed = 123) {
  set.seed(seed)
  datasets <- unique(batch)
  genes <- rownames(expr)
  ngenes <- length(genes)
  nds <- length(datasets)

  calc_logfc <- function(e, g) {
    rowMeans(e[, g == "Disease", drop = FALSE]) - rowMeans(e[, g == "Normal", drop = FALSE])
  }

  obs_fc <- matrix(NA_real_, nrow = ngenes, ncol = nds,
                   dimnames = list(genes, datasets))
  for (d in datasets) {
    idx <- which(batch == d)
    obs_fc[, d] <- calc_logfc(expr[, idx, drop = FALSE], group[idx])
  }

  rank_up <- apply(-obs_fc, 2, rank, ties.method = "random")
  rank_down <- apply(obs_fc, 2, rank, ties.method = "random")
  obs_logrp_up <- rowSums(log(rank_up))
  obs_logrp_down <- rowSums(log(rank_down))

  perm_logrp_up <- matrix(NA_real_, nrow = ngenes, ncol = nperm)
  perm_logrp_down <- matrix(NA_real_, nrow = ngenes, ncol = nperm)
  for (b in seq_len(nperm)) {
    perm_fc <- matrix(NA_real_, nrow = ngenes, ncol = nds,
                      dimnames = list(genes, datasets))
    for (d in datasets) {
      idx <- which(batch == d)
      perm_g <- sample(group[idx])
      perm_fc[, d] <- calc_logfc(expr[, idx, drop = FALSE], perm_g)
    }
    pr_up <- apply(-perm_fc, 2, rank, ties.method = "random")
    pr_down <- apply(perm_fc, 2, rank, ties.method = "random")
    perm_logrp_up[, b] <- rowSums(log(pr_up))
    perm_logrp_down[, b] <- rowSums(log(pr_down))
  }

  p_up <- (rowSums(perm_logrp_up <= obs_logrp_up) + 1) / (nperm + 1)
  p_down <- (rowSums(perm_logrp_down <= obs_logrp_down) + 1) / (nperm + 1)

  data.frame(
    gene.names = genes,
    ave.fold.change = rowMeans(obs_fc),
    RP_up = exp(obs_logrp_up),
    pval_up = p_up,
    FDR_up = p.adjust(p_up, method = "BH"),
    RP_down = exp(obs_logrp_down),
    pval_down = p_down,
    FDR_down = p.adjust(p_down, method = "BH"),
    stringsAsFactors = FALSE
  )
}

origin <- as.numeric(factor(batch_all))
cl <- cl_all
rp_res <- tryCatch({
  RPadvance(data = expr_all, cl = cl, origin = origin, logged = TRUE,
            gene.names = rownames(expr_all), num.perm = 100, rand = 100, huge = TRUE)
}, error = function(e) {
  cat("RPadvance failed:", conditionMessage(e), "\n")
  cat("Falling back to manual permutation-based rank product...\n")
  run_manual_rankprod(expr_all, group_all, batch_all, nperm = 200)
})

if ("RP_up" %in% names(rp_res)) {
  # 来自手动实现
  rp_up <- data.frame(
    symbol = rp_res$gene.names,
    logFC = rp_res$ave.fold.change,
    RP_score = rp_res$RP_up,
    pvalue = rp_res$pval_up,
    FDR = rp_res$FDR_up,
    stringsAsFactors = FALSE
  ) %>% arrange(pvalue)
  rp_down <- data.frame(
    symbol = rp_res$gene.names,
    logFC = rp_res$ave.fold.change,
    RP_score = rp_res$RP_down,
    pvalue = rp_res$pval_down,
    FDR = rp_res$FDR_down,
    stringsAsFactors = FALSE
  ) %>% arrange(pvalue)
} else {
  rp_up <- data.frame(
    symbol = rp_res$AveFC$gene.names,
    logFC = rp_res$AveFC$ave.fold.change,
    RP_score = rp_res$RPs[, 1],
    pvalue = rp_res$pvals[, 1],
    FDR = rp_res$pfp[, 1]
  ) %>% arrange(pvalue)

  rp_down <- data.frame(
    symbol = rp_res$AveFC$gene.names,
    logFC = rp_res$AveFC$ave.fold.change,
    RP_score = rp_res$RPs[, 2],
    pvalue = rp_res$pvals[, 2],
    FDR = rp_res$pfp[, 2]
  ) %>% arrange(pvalue)
}

rp_up$direction <- "up"
rp_down$direction <- "down"
rp_all <- bind_rows(rp_up, rp_down) %>% arrange(pvalue)
write.table(rp_all, "meta_analysis_DEG_rankprod.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: meta_analysis_DEG_rankprod.txt\n")

# ---------- 5. 与当前 5 个核心基因比较 ----------
core_genes <- read.table("cell_death_intersect_three.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1]

if ("RP_up" %in% names(rp_res)) {
  rp_fc <- rp_res$ave.fold.change
  rp_names <- rp_res$gene.names
} else {
  rp_fc <- rp_res$AveFC$ave.fold.change
  rp_names <- rp_res$AveFC$gene.names
}

comparison <- data.frame(
  symbol = core_genes,
  limma_logFC = limma_res$logFC[match(core_genes, limma_res$symbol)],
  limma_pvalue = limma_res$P.Value[match(core_genes, limma_res$symbol)],
  rankprod_logFC = rp_fc[match(core_genes, rp_names)],
  rankprod_pvalue_up = rp_up$pvalue[match(core_genes, rp_up$symbol)],
  rankprod_pvalue_down = rp_down$pvalue[match(core_genes, rp_down$symbol)]
)
write.table(comparison, "meta_analysis_core_genes_comparison.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: meta_analysis_core_genes_comparison.txt\n")
print(comparison)

# ---------- 6. Top 稳健 DEG 摘要 ----------
# RankProd FDR < 0.05 且方向一致
robust_deg <- rp_all %>% filter(FDR < 0.05) %>% arrange(direction, pvalue)
cat("\nRankProd robust DEGs (FDR<0.05):", nrow(robust_deg), "\n")
write.table(robust_deg, "meta_analysis_robust_DEG_rankprod.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Limma adj.P.Val < 0.05
robust_limma <- limma_res %>% filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val)
cat("Limma robust DEGs (adj.P.Val<0.05):", nrow(robust_limma), "\n")
write.table(robust_limma, "meta_analysis_robust_DEG_limma.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- 7. 可视化 ----------
# Volcano plot
df_volc <- limma_res %>%
  mutate(
    sig = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 0.5 ~ "sig",
      TRUE ~ "ns"
    ),
    label = ifelse(symbol %in% core_genes, symbol, "")
  )

p_volc <- ggplot(df_volc, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(aes(label = label), color = "black", max.overlaps = 30) +
  scale_color_manual(values = c("ns" = "grey70", "sig" = "#CD534C")) +
  theme_bw() +
  labs(title = "Meta-analysis volcano plot (ComBat + limma)",
       x = "log2 fold change (Disease vs Normal)", y = "-log10(p-value)")
ggsave("meta_analysis_volcano.pdf", p_volc, width = 8, height = 6)
cat("Saved: meta_analysis_volcano.pdf\n")

# Top 50 RankProd 基因热图
top_genes <- head(rp_all$symbol[rp_all$FDR < 0.05], 50)
if (length(top_genes) > 10) {
  heat_mat <- expr_all[top_genes, ]
  heat_mat <- heat_mat[apply(heat_mat, 1, var) > 0, ]
  anno_col <- data.frame(Dataset = batch_all, Group = group_all, row.names = colnames(expr_all))
  pdf("meta_analysis_top50_heatmap.pdf", width = 12, height = 10)
  pheatmap(heat_mat, annotation_col = anno_col, cluster_cols = TRUE, cluster_rows = TRUE,
           scale = "row", show_colnames = FALSE, color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "Top 50 robust DEGs across datasets")
  dev.off()
  cat("Saved: meta_analysis_top50_heatmap.pdf\n")
}

cat("\nMeta-analysis done!\n")
