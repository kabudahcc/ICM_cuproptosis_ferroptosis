# 铜凋亡 / 铁凋亡通路层面的 GSEA 与 ssGSEA 评分分析
# 1) GSEA：在每个数据集的 limma t 统计量排序基因列表上检测通路富集
# 2) ssGSEA：计算每个样本的通路活性评分，比较 Disease vs Normal

rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(limma)
library(GEOquery)
library(edgeR)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(pheatmap)

# ---------- 1. 读取基因集 ----------
cup <- read.delim("cuproptosis_genes.txt", stringsAsFactors = FALSE)
fer <- read.delim("ferroptosis_genes.txt", stringsAsFactors = FALSE)
fer <- fer %>% filter(confidence == "Validated", organism == "Human")

gene_sets <- list(
  Cuproptosis_core = cup$symbol[cup$source == "Tsvetkov2022_core"],
  Cuproptosis_all = unique(cup$symbol),
  Ferroptosis_driver = unique(fer$symbol[fer$category == "driver"]),
  Ferroptosis_suppressor = unique(fer$symbol[fer$category == "suppressor"]),
  Ferroptosis_marker = unique(fer$symbol[fer$category == "marker"]),
  Ferroptosis_all = unique(fer$symbol)
)
# 过滤空集合
gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
cat("Gene set sizes:\n")
print(sapply(gene_sets, length))

term2gene <- bind_rows(lapply(names(gene_sets), function(nm) {
  data.frame(term = nm, gene = gene_sets[[nm]], stringsAsFactors = FALSE)
}))

# ---------- 2. 数据集加载函数 ----------
load_discovery <- function(ds) {
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = TRUE, row.names = 1, check.names = FALSE)
  group_df <- read.csv(paste0(ds, "_group.csv"), header = TRUE, check.names = FALSE)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  common <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common, drop = FALSE]
  group <- group_df[common, "group"]
  list(expr = as.matrix(exp), group = group,
       disease = ifelse(group == "Disease", 1, 0))
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
  
  title_to_geo <- data.frame(title = pd$title, geo_accession = rownames(pd),
                             group = pd$group, stringsAsFactors = FALSE)
  colnames(exp) <- trimws(colnames(exp))
  matched <- title_to_geo %>% filter(title %in% colnames(exp))
  expr <- exp[, matched$title, drop = FALSE]
  colnames(expr) <- matched$geo_accession
  list(expr = log2(expr + 1), group = matched$group,
       disease = ifelse(matched$group == "Disease", 1, 0))
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
  list(expr = log2(expr + 1), group = readme$group,
       disease = ifelse(readme$group == "Disease", 1, 0))
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
  list(expr = expr, group = pd$group,
       disease = ifelse(pd$group == "Disease", 1, 0))
}

# ---------- 3. 对每个数据集做 GSEA 和 ssGSEA ----------
ds_loaders <- list(
  GSE16499 = load_discovery,
  GSE5406 = load_discovery,
  GSE57338 = load_discovery,
  GSE116250 = load_GSE116250,
  GSE55296 = load_GSE55296,
  GSE42955 = load_GSE42955
)

gsea_results <- list()
gsva_scores_list <- list()
gsva_comparison <- list()

for (ds in names(ds_loaders)) {
  cat("\n==========", ds, "==========\n")
  if (ds %in% c("GSE16499", "GSE5406", "GSE57338")) {
    dat <- ds_loaders[[ds]](ds)
  } else {
    dat <- ds_loaders[[ds]]()
  }
  expr <- dat$expr
  expr[is.na(expr)] <- 0
  group <- dat$group
  
  # --- GSEA：基于 limma t 统计量 ---
  # 只有 discovery 和 GSE116250 已有差异分析文件；其余直接跑 limma
  diff_file <- paste0(ds, "_diff.txt")
  if (file.exists(diff_file)) {
    diff_df <- read.delim(diff_file, stringsAsFactors = FALSE)
  } else {
    design <- model.matrix(~ 0 + factor(group, levels = c("Normal", "Disease")))
    colnames(design) <- c("Normal", "Disease")
    cont <- makeContrasts(Disease - Normal, levels = design)
    fit <- lmFit(expr, design)
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- eBayes(fit2)
    diff_df <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
    diff_df$symbol <- rownames(diff_df)
  }
  
  geneList <- setNames(diff_df$t, diff_df$symbol)
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- geneList[!is.na(geneList)]
  
  ego <- tryCatch(
    GSEA(geneList = geneList, TERM2GENE = term2gene,
         pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE),
    error = function(e) {
      cat("GSEA error in", ds, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego)
    ego_df$Dataset <- ds
    gsea_results[[ds]] <- ego_df
    cat("GSEA significant pathways (p<0.05):", sum(ego_df$pvalue < 0.05), "\n")
  } else {
    cat("No GSEA result for", ds, "\n")
  }
  
  # --- ssGSEA 评分 ---
  # 只保留基因集中在表达矩阵里有的基因
  gsets_in_expr <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr)))
  gsets_in_expr <- gsets_in_expr[sapply(gsets_in_expr, length) >= 3]
  if (length(gsets_in_expr) == 0) next
  
  param <- gsvaParam(expr, gsets_in_expr, assay = NA_character_, minSize = 1)
  gsva_res <- gsva(param, verbose = FALSE)
  
  score_df <- as.data.frame(t(gsva_res))
  score_df$sample <- rownames(score_df)
  score_df$group <- group
  score_df$dataset <- ds
  gsva_scores_list[[ds]] <- score_df
  
  # 组间比较
  comp <- lapply(names(gsets_in_expr), function(path) {
    x <- score_df[[path]][score_df$group == "Disease"]
    y <- score_df[[path]][score_df$group == "Normal"]
    if (length(x) < 2 || length(y) < 2) return(NULL)
    w <- wilcox.test(x, y)
    data.frame(
      Dataset = ds,
      Pathway = path,
      Disease_mean = mean(x),
      Normal_mean = mean(y),
      logFC = mean(x) - mean(y),
      pvalue = w$p.value,
      stringsAsFactors = FALSE
    )
  })
  gsva_comparison[[ds]] <- bind_rows(comp)
}

# ---------- 4. 保存结果 ----------
# GSEA 汇总
gsea_all <- bind_rows(gsea_results)
if (nrow(gsea_all) > 0) {
  write.table(gsea_all, "GSEA_cuproptosis_ferroptosis_summary.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nSaved: GSEA_cuproptosis_ferroptosis_summary.txt\n")
}

# ssGSEA 评分
gsva_scores_all <- bind_rows(gsva_scores_list)
write.table(gsva_scores_all, "GSVA_pathway_scores_all_datasets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: GSVA_pathway_scores_all_datasets.txt\n")

# ssGSEA 组间比较
gsva_comp_all <- bind_rows(gsva_comparison)
gsva_comp_all$FDR <- p.adjust(gsva_comp_all$pvalue, method = "BH")
gsva_comp_all <- gsva_comp_all %>% arrange(pvalue)
write.table(gsva_comp_all, "GSVA_pathway_comparison.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: GSVA_pathway_comparison.txt\n")
print(head(gsva_comp_all, 20))

# ---------- 5. 可视化 ----------
if (nrow(gsva_scores_all) > 0) {
  # 箱线图
  plot_df <- gsva_scores_all %>%
    pivot_longer(cols = all_of(names(gene_sets)), names_to = "Pathway", values_to = "Score")
  plot_df$group <- factor(plot_df$group, levels = c("Normal", "Disease"))
  
  p_box <- ggplot(plot_df, aes(x = group, y = Score, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) +
    facet_grid(Pathway ~ dataset, scales = "free_y") +
    scale_fill_manual(values = c("Normal" = "#0073C2", "Disease" = "#CD534C")) +
    theme_bw() +
    labs(title = "ssGSEA pathway scores across datasets",
         x = "", y = "ssGSEA score", fill = "Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("GSVA_pathway_boxplot.pdf", p_box, width = 16, height = 12)
  cat("Saved: GSVA_pathway_boxplot.pdf\n")
}

# GSEA NES 热图
if (nrow(gsea_all) > 0) {
  nes_mat <- gsea_all %>%
    dplyr::select(Dataset, ID, NES) %>%
    pivot_wider(names_from = Dataset, values_from = NES, values_fill = NA)
  nes_mat <- tibble::column_to_rownames(nes_mat, "ID")
  nes_mat <- as.matrix(nes_mat)
  
  pdf("GSEA_pathway_NES_heatmap.pdf", width = 8, height = 5)
  pheatmap(nes_mat, color = colorRampPalette(c("blue", "white", "red"))(50),
           display_numbers = TRUE, main = "GSEA NES across datasets")
  dev.off()
  cat("Saved: GSEA_pathway_NES_heatmap.pdf\n")
}

cat("\nPathway scoring done!\n")
