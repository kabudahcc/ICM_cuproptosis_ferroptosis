# 提取三数据集共有细胞死亡相关差异基因的表达模式并绘制热图
# 输入：CRG_intersect_three.txt, FRG_intersect_three.txt, GSE*_diff.txt
# 输出：cell_death_logFC_matrix.txt, cell_death_heatmap.pdf

rm(list = ls())
library(dplyr)
library(data.table)
library(pheatmap)

# 读取细胞死亡相关交集基因
crg <- read.table("CRG_intersect_three.txt", header = F, sep = "\t")$V1
frg <- read.table("FRG_intersect_three.txt", header = F, sep = "\t")$V1
cell_death_genes <- sort(unique(c(crg, frg)))

cat("Cell death intersect genes:", length(cell_death_genes), "\n")
print(cell_death_genes)

# 读取三个差异基因文件并提取 logFC
datasets <- c("GSE16499", "GSE5406", "GSE57338")
mat <- matrix(NA, nrow = length(cell_death_genes), ncol = length(datasets))
rownames(mat) <- cell_death_genes
colnames(mat) <- datasets

for (i in seq_along(datasets)) {
  ds <- datasets[i]
  df <- read.table(paste0(ds, "_diff.txt"), header = T, sep = "\t", check.names = F)
  df$symbol <- trimws(as.character(df$symbol))
  df <- df[!duplicated(df$symbol), ]
  row_idx <- match(cell_death_genes, df$symbol)
  mat[, i] <- df$logFC[row_idx]
}

# 保存 logFC 矩阵
write.table(mat, "cell_death_logFC_matrix.txt", sep = "\t", quote = F, row.names = T, col.names = NA)

# 绘制热图
if (nrow(mat) > 1) {
  bk <- c(seq(min(mat, na.rm = T), -0.01, length = 50),
          seq(0.01, max(mat, na.rm = T), length = 50))
  col <- colorRampPalette(c("#0073C2", "white", "#CD534C"))(99)
  
  pdf("cell_death_heatmap.pdf", width = 5, height = 4)
  pheatmap(mat,
           scale = "none",
           color = col,
           breaks = bk,
           cluster_cols = F,
           cluster_rows = T,
           display_numbers = TRUE,
           number_color = "black",
           fontsize_number = 8,
           main = "Cell death-related DEGs across GEO datasets",
           filename = NA)
  dev.off()
} else if (nrow(mat) == 1) {
  # 只有一行时用条形图
  pdf("cell_death_logFC_barplot.pdf", width = 5, height = 3)
  barplot(mat[1, ], col = ifelse(mat[1, ] > 0, "#CD534C", "#0073C2"),
          main = paste("logFC of", rownames(mat)),
          ylab = "logFC", ylim = range(c(0, mat[1, ]), na.rm = T) * 1.2)
  dev.off()
}

# 输出每个基因在三个数据集中的表达方向
sink("cell_death_logFC_direction.txt")
cat("Gene\tGSE16499\tGSE5406\tGSE57338\tConsistency\n")
for (g in rownames(mat)) {
  dirs <- ifelse(mat[g, ] > 0, "UP", ifelse(mat[g, ] < 0, "DOWN", "NA"))
  consist <- ifelse(length(unique(dirs[!is.na(dirs)])) == 1, "Consistent", "Inconsistent")
  cat(g, "\t", paste(dirs, collapse = "\t"), "\t", consist, "\n")
}
sink()

cat("\nDone!\n")
cat("- cell_death_logFC_matrix.txt\n")
cat("- cell_death_heatmap.pdf\n")
cat("- cell_death_logFC_direction.txt\n")
