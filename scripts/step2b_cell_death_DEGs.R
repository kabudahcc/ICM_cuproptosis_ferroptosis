# 铜凋亡 / 铁凋亡相关差异基因筛选与可视化
# 输入：GSE16499_diff.txt, GSE5406_diff.txt, GSE57338_diff.txt
#        cuproptosis_genes.txt, ferroptosis_genes.txt
# 输出：
#   - CRG_DEG_GSE*.txt / FRG_DEG_GSE*.txt（各数据集铜凋亡/铁凋亡差异基因）
#   - CRG_intersect_three.txt / FRG_intersect_three.txt（三数据集共有）
#   - CRG_FRG_crosstalk_three.txt（同时与铜凋亡和铁凋亡相关的共有差异基因）
#   - venn_CRG_three.pdf / venn_FRG_three.pdf
#   - cell_death_DEG_summary.txt（汇总表）

rm(list = ls())
library(grid)
library(VennDiagram)
library(dplyr)
library(data.table)

# 1. 设置阈值 -----------------------------------------------------------
logFC_filter <- 0.5    # log2 倍数变化阈值
pval_filter  <- 0.05   # P 值阈值（也可用 adj.P.Val < 0.05）

datasets <- c("GSE16499", "GSE5406", "GSE57338")
diff_files <- paste0(datasets, "_diff.txt")

# 2. 读取铜凋亡、铁凋亡基因列表 -----------------------------------------
cup_genes <- fread("cuproptosis_genes.txt", header = T, sep = "\t") %>%
  pull(symbol) %>% unique() %>% as.character() %>% trimws()

fer_genes <- fread("ferroptosis_genes.txt", header = T, sep = "\t") %>%
  pull(symbol) %>% unique() %>% as.character() %>% trimws()

cat("Cuproptosis genes loaded:", length(cup_genes), "\n")
cat("Ferroptosis genes loaded:", length(fer_genes), "\n")

# 3. 读取并筛选差异基因 -------------------------------------------------
get_deg <- function(file) {
  df <- read.table(file, header = T, sep = "\t", check.names = F)
  df$symbol <- trimws(as.character(df$symbol))
  deg <- df %>% dplyr::filter(abs(logFC) > logFC_filter, P.Value < pval_filter)
  return(deg)
}

deg_list <- lapply(diff_files, get_deg)
names(deg_list) <- datasets

for (ds in datasets) {
  cat(ds, "total DEGs:", nrow(deg_list[[ds]]), "\n")
}

# 4. 筛选铜凋亡 / 铁凋亡差异基因 ----------------------------------------
get_celldeath_deg <- function(deg_df, gene_set, set_name) {
  res <- deg_df %>% dplyr::filter(symbol %in% gene_set)
  return(res)
}

cup_deg_list <- lapply(deg_list, get_celldeath_deg, gene_set = cup_genes, set_name = "cuproptosis")
fer_deg_list <- lapply(deg_list, get_celldeath_deg, gene_set = fer_genes, set_name = "ferroptosis")

# 输出各数据集结果
for (ds in datasets) {
  fwrite(cup_deg_list[[ds]], paste0("CRG_DEG_", ds, ".txt"), sep = "\t", quote = F)
  fwrite(fer_deg_list[[ds]], paste0("FRG_DEG_", ds, ".txt"), sep = "\t", quote = F)
  cat(ds, "cuproptosis DEGs:", nrow(cup_deg_list[[ds]]),
      "| ferroptosis DEGs:", nrow(fer_deg_list[[ds]]), "\n")
}

# 提取基因名向量
cup_gene_list <- lapply(cup_deg_list, function(x) unique(x$symbol))
fer_gene_list <- lapply(fer_deg_list, function(x) unique(x$symbol))

# 5. 三数据集共有铜凋亡 / 铁凋亡差异基因 --------------------------------
cup_inter_three <- Reduce(intersect, cup_gene_list)
fer_inter_three <- Reduce(intersect, fer_gene_list)

cat("\nThree-way cuproptosis DEG intersection:", length(cup_inter_three), "\n")
if (length(cup_inter_three) > 0) cat(cup_inter_three, "\n")

cat("Three-way ferroptosis DEG intersection:", length(fer_inter_three), "\n")
if (length(fer_inter_three) > 0) cat(fer_inter_three, "\n")

write.table(cup_inter_three, "CRG_intersect_three.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(fer_inter_three, "FRG_intersect_three.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# 6. 铜凋亡-铁凋亡 crosstalk（同时属于两个基因集） ----------------------
cup_fer_overlap <- intersect(cup_genes, fer_genes)   # 本身既是铜凋亡又是铁凋亡的基因
cat("\nGenes annotated as both cuproptosis and ferroptosis:", length(cup_fer_overlap), "\n")
if (length(cup_fer_overlap) > 0) cat(cup_fer_overlap, "\n")

# 在三数据集共有差异基因中，同时与铜凋亡和铁凋亡相关的基因
cup_fer_deg_inter <- intersect(cup_inter_three, fer_inter_three)
cat("Three-way DEGs with both cuproptosis & ferroptosis annotation:", length(cup_fer_deg_inter), "\n")
if (length(cup_fer_deg_inter) > 0) cat(cup_fer_deg_inter, "\n")
write.table(cup_fer_deg_inter, "CRG_FRG_crosstalk_three.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# 7. 绘制韦恩图 ---------------------------------------------------------
plot_venn <- function(gene_list, title, outpdf, colors = c("#0073C2FF", "#EFC000FF", "#CD534CFF")) {
  venn.plot <- venn.diagram(
    x = gene_list,
    filename = NULL,
    fill = colors,
    scaled = FALSE,
    cat.pos = c(-20, 20, 180),
    cat.col = colors,
    cat.cex = 1.2,
    category = names(gene_list),
    main = title,
    main.cex = 1.5
  )
  pdf(file = outpdf, width = 6, height = 6)
  grid.draw(venn.plot)
  dev.off()
}

if (sum(sapply(cup_gene_list, length)) > 0) {
  plot_venn(cup_gene_list, "Cuproptosis-related DEGs", "venn_CRG_three.pdf")
}
if (sum(sapply(fer_gene_list, length)) > 0) {
  plot_venn(fer_gene_list, "Ferroptosis-related DEGs", "venn_FRG_three.pdf")
}

# 8. 生成汇总表 ---------------------------------------------------------
summary_df <- data.frame(
  Dataset = datasets,
  Total_DEGs = sapply(deg_list, nrow),
  Cuproptosis_DEGs = sapply(cup_deg_list, nrow),
  Ferroptosis_DEGs = sapply(fer_deg_list, nrow),
  Cuproptosis_unique_genes = sapply(cup_gene_list, length),
  Ferroptosis_unique_genes = sapply(fer_gene_list, length)
)

sink("cell_death_DEG_summary.txt")
cat("=== 铜凋亡 / 铁凋亡差异基因分析汇总 ===\n")
cat("筛选阈值: |logFC| >", logFC_filter, ", P.Value <", pval_filter, "\n\n")
print(summary_df, row.names = F)
cat("\n三数据集共有铜凋亡差异基因 (", length(cup_inter_three), "):\n")
cat(ifelse(length(cup_inter_three) > 0, paste(cup_inter_three, collapse = ", "), "无"), "\n")
cat("\n三数据集共有铁凋亡差异基因 (", length(fer_inter_three), "):\n")
cat(ifelse(length(fer_inter_three) > 0, paste(fer_inter_three, collapse = ", "), "无"), "\n")
cat("\n三数据集共有且同时标注为铜凋亡+铁凋亡的差异基因 (", length(cup_fer_deg_inter), "):\n")
cat(ifelse(length(cup_fer_deg_inter) > 0, paste(cup_fer_deg_inter, collapse = ", "), "无"), "\n")
sink()

# 9. 保存每个数据集的铜凋亡/铁凋亡差异基因名 ---------------------------
for (ds in datasets) {
  write.table(cup_gene_list[[ds]], paste0("CRG_genes_", ds, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(fer_gene_list[[ds]], paste0("FRG_genes_", ds, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
}

cat("\nDone! Output files:\n")
cat("- CRG_DEG_GSE*.txt / FRG_DEG_GSE*.txt\n")
cat("- CRG_intersect_three.txt / FRG_intersect_three.txt\n")
cat("- CRG_FRG_crosstalk_three.txt\n")
cat("- venn_CRG_three.pdf / venn_FRG_three.pdf\n")
cat("- cell_death_DEG_summary.txt\n")
