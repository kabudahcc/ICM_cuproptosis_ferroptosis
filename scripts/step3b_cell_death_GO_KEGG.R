# 铜凋亡 / 铁凋亡共有差异基因的 GO 与 KEGG 富集分析
# 输入：cell_death_intersect_three.txt（或自行替换为 CRG_intersect_three.txt / FRG_intersect_three.txt）
# 输出：GO.txt, KEGG.txt, GO_barplot.pdf, GO_dotplot.pdf, KEGG_barplot.pdf, KEGG_dotplot.pdf

rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 1. 设置参数 -----------------------------------------------------------
inputFile <- "cell_death_intersect_three.txt"   # 可替换：CRG_intersect_three.txt / FRG_intersect_three.txt
pvalueFilter <- 0.05
qvalueFilter <- 1

# 2. 读取基因 -----------------------------------------------------------
genes <- read.table(inputFile, header = F, sep = "\t", check.names = F)[, 1]
genes <- unique(as.character(genes))
cat("Input genes:", length(genes), "\n")
print(genes)

# 3. 转换为 ENTREZ ID ---------------------------------------------------
gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_ids <- gene_df$ENTREZID
cat("Mapped ENTREZ IDs:", length(entrez_ids), "\n")

# 4. GO 富集分析 --------------------------------------------------------
if (length(entrez_ids) >= 1) {
  kk <- enrichGO(
    gene = entrez_ids,
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    ont = "all",
    readable = TRUE,
    minGSSize = 1
  )
  GO <- as.data.frame(kk)
  GO_sig <- GO[GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter, ]
  
  write.table(GO, file = "GO_cell_death.txt", sep = "\t", quote = F, row.names = F)
  cat("GO terms (p<0.05):", nrow(GO_sig), "\n")
  
  if (nrow(GO_sig) > 0) {
    # 柱状图
    pdf("GO_cell_death_barplot.pdf", width = 9, height = 7)
    p <- barplot(kk, drop = TRUE, showCategory = 10, label_format = 130, split = "ONTOLOGY", color = "pvalue") +
      facet_grid(ONTOLOGY ~ ., scales = "free")
    print(p)
    dev.off()
    
    # 气泡图
    pdf("GO_cell_death_dotplot.pdf", width = 9, height = 7)
    p <- dotplot(kk, showCategory = 10, orderBy = "GeneRatio", label_format = 130, split = "ONTOLOGY", color = "pvalue") +
      facet_grid(ONTOLOGY ~ ., scales = "free")
    print(p)
    dev.off()
  }
}

# 5. KEGG 富集分析 ------------------------------------------------------
if (length(entrez_ids) >= 1) {
  kk_kegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 1
  )
  KEGG <- as.data.frame(kk_kegg)
  KEGG_sig <- KEGG[KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter, ]
  
  write.table(KEGG, file = "KEGG_cell_death.txt", sep = "\t", quote = F, row.names = F)
  cat("KEGG pathways (p<0.05):", nrow(KEGG_sig), "\n")
  
  if (nrow(KEGG_sig) > 0) {
    pdf("KEGG_cell_death_barplot.pdf", width = 8, height = 10)
    print(barplot(kk_kegg, showCategory = 20, color = "pvalue"))
    dev.off()
    
    pdf("KEGG_cell_death_dotplot.pdf", width = 8, height = 10)
    print(dotplot(kk_kegg, showCategory = 20))
    dev.off()
  }
}

cat("\nDone!\n")
