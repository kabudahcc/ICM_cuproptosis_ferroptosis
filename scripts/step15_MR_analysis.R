# step15_MR_analysis.R
# 基于 meta-analysis  top DEGs，使用 eQTLGen 血液 cis-eQTL 作为暴露、HERMES 心力衰竭 GWAS 作为结局
# 进行两样本孟德尔随机化（TwoSample MR）。暴露数据来自 eQTLGen 显著 cis-eQTL 文件，结局数据来自
# GWAS Catalog GCST009541（ShahS_31919418_HeartFailure.gz）。

.libPaths(c('Rlib', .libPaths()))

library(data.table)
library(ggplot2)
library(dplyr)

# ---------------------------- 参数 ----------------------------
n_top       <- 100          # 从 meta DEG 中取 top 基因
eqtl_file   <- "MR_data/eqtlgen_cis_fdr0.05.txt.gz"
gwas_file   <- "MR_data/ShahS_31919418_HeartFailure.gz"
deg_file    <- "meta_analysis_DEG_limma.txt"
clump_kb    <- 500          # 基因内 distance clumping 窗口 (kb)
min_F       <- 10           # 弱工具变量过滤
out_dir     <- "MR_results"
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------- 读取 meta DEGs ----------------------------
deg <- fread(deg_file)
deg <- deg[order(adj.P.Val)]
top_genes <- deg$symbol[1:n_top]
writeLines(top_genes, file.path(out_dir, "MR_top_genes.txt"))
cat("Top genes selected:", length(top_genes), "\n")

# ---------------------------- 从 eQTLGen 提取目标基因 ----------------------------
eqtl_filtered <- "MR_data/eqtlgen_filtered.txt.gz"
gene_list_file <- "MR_data/top_genes.txt"
writeLines(top_genes, gene_list_file)
if (!file.exists(eqtl_filtered)) {
  cat("Filtering eQTLGen file for target genes ...\n")
  cmd <- paste0("awk -F'\\t' 'NR==FNR{a[$1];next} FNR==1 || ($9 in a)' ",
                gene_list_file, " <(zcat ", eqtl_file, ") | gzip > ", eqtl_filtered)
  system(paste("bash -c", shQuote(cmd)))
}

eqtl <- fread(eqtl_filtered)
setnames(eqtl, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                 "Zscore","Gene","GeneSymbol","GeneChr","GenePos",
                 "NrCohorts","NrSamples","FDR","BonferroniP"))
# 仅保留目标基因且 P 显著的行（FDR<0.05 文件本身已过滤）
eqtl <- eqtl[GeneSymbol %in% top_genes]
eqtl[, AssessedAllele := toupper(AssessedAllele)]
eqtl[, OtherAllele := toupper(OtherAllele)]
# 由 Zscore 推导 beta/se（标准化表达表型）
eqtl[, beta_eqtl := Zscore / sqrt(NrSamples)]
eqtl[, se_eqtl   := 1 / sqrt(NrSamples)]
eqtl[, F_stat    := beta_eqtl^2 / se_eqtl^2]
cat("eQTL rows after gene filter:", nrow(eqtl), "\n")

# ---------------------------- 读取结局 GWAS ----------------------------
gwas <- fread(gwas_file)
setnames(gwas, c("SNP","CHR","BP","A1","A2","freq","b","se","p","N"))
gwas[, A1 := toupper(A1)]
gwas[, A2 := toupper(A2)]

# 只保留暴露中存在的 SNP（减少内存）
gwas <- gwas[SNP %in% eqtl$SNP]
cat("HF GWAS SNPs matching exposure:", nrow(gwas), "\n")

# ---------------------------- 合并与等位基因协调 ----------------------------
m <- merge(eqtl, gwas, by = "SNP", suffixes = c("_eqtl", "_gwas"))

is_palindromic <- function(a1, a2) {
  paste0(a1, a2) %in% c("AT","TA","CG","GC")
}

harmonize <- function(df) {
  df[, beta_eqtl_harm := as.numeric(NA)]
  df[, keep := FALSE]
  # 方向一致
  df[AssessedAllele == A1 & OtherAllele == A2, `:=`(beta_eqtl_harm = beta_eqtl, keep = TRUE)]
  # 方向相反
  df[AssessedAllele == A2 & OtherAllele == A1, `:=`(beta_eqtl_harm = -beta_eqtl, keep = TRUE)]
  # 剔除回文 SNP（无 MAF 信息，保守处理）
  df[keep == TRUE & is_palindromic(AssessedAllele, OtherAllele), keep := FALSE]
  df[keep == TRUE]
}

mh <- harmonize(m)
mh <- mh[F_stat >= min_F]
cat("SNPs after harmonization and F>=", min_F, ":", nrow(mh), "\n")

if (nrow(mh) == 0) stop("No valid instruments after harmonization.")

# ---------------------------- 按基因 distance clumping ----------------------------
setorder(mh, GeneSymbol, Pvalue)
clump_gene <- function(df, kb = clump_kb) {
  df <- as.data.frame(df)
  kept <- rep(FALSE, nrow(df))
  for (i in seq_len(nrow(df))) {
    if (!kept[i]) {
      kept_pos <- df$SNPPos[kept]
      if (length(kept_pos) == 0 || all(abs(df$SNPPos[i] - kept_pos) > kb * 1000)) {
        kept[i] <- TRUE
      }
    }
  }
  df[kept, ]
}

clumped <- mh[, clump_gene(.SD, kb = clump_kb), by = GeneSymbol]
clumped <- as.data.table(clumped)
cat("Instruments after clumping:", nrow(clumped), "across", length(unique(clumped$GeneSymbol)), "genes\n")

# ---------------------------- MR 分析 ----------------------------
run_mr <- function(df) {
  df <- as.data.frame(df)
  k <- nrow(df)
  if (k == 0) return(data.frame(nsnp = 0L))
  bx <- df$beta_eqtl_harm
  bxse <- df$se_eqtl
  by <- df$b
  byse <- df$se
  F_min <- min(df$F_stat)

  if (k == 1) {
    b_mr <- by / bx
    se_mr <- byse / abs(bx)
    p_mr <- 2 * pnorm(-abs(b_mr / se_mr))
    return(data.frame(
      nsnp = 1L, method = "Wald_ratio",
      b = b_mr, se = se_mr, p = p_mr,
      Q = NA_real_, p_het = NA_real_,
      b_egger = NA_real_, se_egger = NA_real_, p_egger = NA_real_,
      egger_intercept = NA_real_, p_pleio = NA_real_,
      F_min = F_min, stringsAsFactors = FALSE
    ))
  }

  # IVW (fixed-effect)
  w <- bx^2 / byse^2
  b_ivw <- sum(bx * by / byse^2) / sum(w)
  se_ivw <- sqrt(1 / sum(w))
  p_ivw <- 2 * pnorm(-abs(b_ivw / se_ivw))
  Q <- sum(((by - b_ivw * bx) / byse)^2)
  p_het <- pchisq(Q, df = k - 1, lower.tail = FALSE)

  if (k >= 3) {
    fit <- summary(lm(by ~ bx, weights = 1 / byse^2))
    b_eg <- coef(fit)[2, 1]
    se_eg <- coef(fit)[2, 2]
    p_eg <- coef(fit)[2, 4]
    int_eg <- coef(fit)[1, 1]
    p_int <- coef(fit)[1, 4]
  } else {
    b_eg <- se_eg <- p_eg <- int_eg <- p_int <- NA_real_
  }

  data.frame(
    nsnp = k, method = "IVW",
    b = b_ivw, se = se_ivw, p = p_ivw,
    Q = Q, p_het = p_het,
    b_egger = b_eg, se_egger = se_eg, p_egger = p_eg,
    egger_intercept = int_eg, p_pleio = p_int,
    F_min = F_min, stringsAsFactors = FALSE
  )
}

mr_res <- clumped[, run_mr(.SD), by = GeneSymbol]
mr_res <- as.data.table(mr_res)
mr_res <- mr_res[nsnp > 0]

# 合并表达方向、计算 OR 与校正
mr_res <- merge(mr_res, deg[, .(symbol, logFC)], by.x = "GeneSymbol", by.y = "symbol", all.x = TRUE)
mr_res[, OR := exp(b)]
mr_res[, OR_lower := exp(b - 1.96 * se)]
mr_res[, OR_upper := exp(b + 1.96 * se)]
mr_res[, p_bonf := pmin(p * .N, 1)]
mr_res[, p_fdr := p.adjust(p, method = "BH")]
setorder(mr_res, p)

fwrite(mr_res, file.path(out_dir, "MR_results.txt"), sep = "\t")
fwrite(clumped, file.path(out_dir, "MR_instruments.txt"), sep = "\t")
cat("Saved:", file.path(out_dir, "MR_results.txt"), "\n")
cat("Significant genes (p < 0.05):", sum(mr_res$p < 0.05, na.rm = TRUE), "\n")
cat("Bonferroni-significant genes:", sum(mr_res$p_bonf < 0.05, na.rm = TRUE), "\n")

# ---------------------------- 可视化 ----------------------------
sig <- mr_res[p < 0.05]
if (nrow(sig) > 0) {
  p <- ggplot(sig, aes(x = reorder(GeneSymbol, b), y = b, ymin = b - 1.96 * se, ymax = b + 1.96 * se)) +
    geom_pointrange(aes(color = p < 0.05)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    coord_flip() +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2")) +
    theme_bw() +
    labs(title = "MR estimates for top DEGs on heart failure",
         subtitle = paste0("IVW / Wald ratio; nsnp>=", min_F, " F-stat"),
         x = "Gene", y = "MR beta (log OR per SD expression)", color = "p < 0.05") +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "MR_forest_plot.pdf"), p, width = 6, height = max(4, nrow(sig) * 0.3))
  cat("Saved:", file.path(out_dir, "MR_forest_plot.pdf"), "\n")
}

cat("\nMR analysis done!\n")
