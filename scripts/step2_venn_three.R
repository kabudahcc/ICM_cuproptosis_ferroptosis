# 三个 GEO 数据集差异基因韦恩图
# 输入：GSE16499_diff.txt, GSE5406_diff.txt, GSE57338_diff.txt
# 输出：venn_three.pdf, interGenes_three.txt

rm(list = ls())
library(grid)
library(futile.logger)
library(VennDiagram)
library(tidyverse)

# 1. 设置阈值 -----------------------------------------------------------
logFC_filter <- 0.5    # log2 倍数变化阈值
pval_filter <- 0.05    # P 值阈值

# 2. 读取三个差异基因文件 -----------------------------------------------
file1 <- "GSE16499_diff.txt"
file2 <- "GSE5406_diff.txt"
file3 <- "GSE57338_diff.txt"

deg1 <- read.table(file1, header = T, sep = "\t", check.names = F)
deg2 <- read.table(file2, header = T, sep = "\t", check.names = F)
deg3 <- read.table(file3, header = T, sep = "\t", check.names = F)

# 3. 筛选差异基因 -------------------------------------------------------
# 根据 logFC 和 P.Value 筛选
genes1 <- deg1 %>%
  dplyr::filter(abs(logFC) > logFC_filter, P.Value < pval_filter) %>%
  dplyr::pull(symbol) %>%
  as.character() %>%
  gsub("^ | $", "", .) %>%
  unique()

genes2 <- deg2 %>%
  dplyr::filter(abs(logFC) > logFC_filter, P.Value < pval_filter) %>%
  dplyr::pull(symbol) %>%
  as.character() %>%
  gsub("^ | $", "", .) %>%
  unique()

genes3 <- deg3 %>%
  dplyr::filter(abs(logFC) > logFC_filter, P.Value < pval_filter) %>%
  dplyr::pull(symbol) %>%
  as.character() %>%
  gsub("^ | $", "", .) %>%
  unique()

# 查看每个数据集的差异基因数量
cat("GSE16499 DEGs:", length(genes1), "\n")
cat("GSE5406 DEGs:", length(genes2), "\n")
cat("GSE57338 DEGs:", length(genes3), "\n")

# 4. 准备韦恩图数据 -----------------------------------------------------
venn_list <- list(
  GSE16499 = genes1,
  GSE5406 = genes2,
  GSE57338 = genes3
)

# 5. 绘制韦恩图 ----------------------------------------------------------
color <- c("#0073C2FF", "#EFC000FF", "#CD534CFF")

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = color,
  scaled = FALSE,
  cat.pos = c(-20, 20, 180),      # 三个类别标签的位置
  cat.col = color,
  cat.cex = 1.2,
  category = c("GSE16499", "GSE5406", "GSE57338"),
  main = "DEGs Venn Diagram",
  main.cex = 1.5
)

pdf(file = "venn_three.pdf", width = 6, height = 6)
grid.draw(venn.plot)
dev.off()

# 6. 输出各种交集结果 ---------------------------------------------------
# 三个数据集共同交集
inter_three <- Reduce(intersect, venn_list)
cat("Three-way intersection:", length(inter_three), "\n")
write.table(inter_three, file = "interGenes_three.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# 两两交集
inter_12 <- intersect(genes1, genes2)
inter_13 <- intersect(genes1, genes3)
inter_23 <- intersect(genes2, genes3)

cat("GSE16499 & GSE5406:", length(inter_12), "\n")
cat("GSE16499 & GSE57338:", length(inter_13), "\n")
cat("GSE5406 & GSE57338:", length(inter_23), "\n")

# 保存所有两两交集和三界交集
write.table(inter_12, file = "interGenes_GSE16499_GSE5406.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(inter_13, file = "interGenes_GSE16499_GSE57338.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(inter_23, file = "interGenes_GSE5406_GSE57338.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# 7. 也可以保存每个数据集的差异基因列表 --------------------------------
write.table(genes1, file = "DEG_GSE16499.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(genes2, file = "DEG_GSE5406.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(genes3, file = "DEG_GSE57338.txt", sep = "\t", quote = F, col.names = F, row.names = F)

cat("\nDone! Output files:\n")
cat("- venn_three.pdf\n")
cat("- interGenes_three.txt (three-way intersection)\n")
cat("- interGenes_*_*.txt (pairwise intersections)\n")
cat("- DEG_GSE*.txt (individual DEG lists)\n")
