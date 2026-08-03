# compareCluster KEGG：比较铜凋亡 vs 铁凋亡相关差异基因的通路富集
# 使用本地 KEGG 注释文件（避免在线 API 不稳定）
# 输入：CRG_DEG_GSE*.txt / FRG_DEG_GSE*.txt, kegg_pathway_names.txt, kegg_pathway_genes.txt
# 输出：kegg_compareCluster_dotplot.pdf, kegg_compareCluster_alluvial.pdf, KEGG_compareCluster.txt

rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(tidyr)

# 1. 读取本地 KEGG 注释 --------------------------------------------------
# 通路 ID -> 通路名称
term2name <- read.table("kegg_pathway_names.txt", header = F, sep = "\t", stringsAsFactors = F, quote = "")
colnames(term2name) <- c("Term", "Name")
term2name$Term <- gsub("^path:", "", term2name$Term)

# 通路 ID -> hsa:ENTREZID
term2gene <- read.table("kegg_pathway_genes.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(term2gene) <- c("Term", "Gene")
term2gene$Term <- gsub("^path:", "", term2gene$Term)
term2gene$Gene <- gsub("^hsa:", "", term2gene$Gene)

cat("KEGG pathways:", nrow(term2name), "\n")
cat("KEGG pathway-gene pairs:", nrow(term2gene), "\n")

# 2. 读取各数据集铜凋亡 / 铁凋亡差异基因 ---------------------------------
datasets <- c("GSE16499", "GSE5406", "GSE57338")

read_deg_symbols <- function(file) {
  df <- read.table(file, header = T, sep = "\t", check.names = F)
  df$symbol <- trimws(as.character(df$symbol))
  return(unique(df$symbol))
}

# 合并三个数据集的铜凋亡 / 铁凋亡 DEGs
crg_union <- unique(unlist(lapply(paste0("CRG_DEG_", datasets, ".txt"), read_deg_symbols)))
frg_union <- unique(unlist(lapply(paste0("FRG_DEG_", datasets, ".txt"), read_deg_symbols)))

cat("CRG union DEGs:", length(crg_union), "\n")
cat("FRG union DEGs:", length(frg_union), "\n")

# 3. 转换为 ENTREZ ID ----------------------------------------------------
to_entrez <- function(symbols) {
  df <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  return(df$ENTREZID)
}

geneList <- list(
  CRG_DEGs = to_entrez(crg_union),
  FRG_DEGs = to_entrez(frg_union)
)

print(sapply(geneList, length))

# 4. compareCluster 使用本地 enricher ------------------------------------
kegg <- compareCluster(
  geneList,
  fun = "enricher",
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 1,
  maxGSSize = 500
)

kegg_readable <- setReadable(kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
kegg_result <- as.data.frame(kegg_readable)

write.table(kegg_result, "KEGG_compareCluster.txt", sep = "\t", quote = F, row.names = F)
cat("KEGG compareCluster terms:", nrow(kegg_result), "\n")

# 5. 点图 ----------------------------------------------------------------
if (nrow(kegg_result) > 0) {
  # 每个 Cluster 取 top 10 通路
  kegg_top <- kegg_result %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # 统一通路顺序：按整体显著性
  kegg_top$Description <- factor(kegg_top$Description,
                                  levels = unique(kegg_top$Description[order(kegg_top$p.adjust, decreasing = T)]))
  
  p1 <- ggplot(kegg_top, aes(x = Cluster, y = Description, size = Count, color = p.adjust)) +
    geom_point(stroke = 0.5) +
    scale_color_gradient(low = "red", high = "blue", trans = "log10") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 10, angle = 0),
          panel.grid = element_blank()) +
    labs(x = "", y = "KEGG Pathway", size = "Gene Count", color = "p.adjust",
         title = "KEGG pathway enrichment: CRG vs FRG DEGs")
  
  ggsave(p1, file = "kegg_compareCluster_dotplot.pdf", width = 8.5, height = 7)
  cat("Saved: kegg_compareCluster_dotplot.pdf\n")
}

# 6. Alluvial / Sankey 图：基因-通路关系 ---------------------------------
# 取每个 Cluster top 5 通路
alluvial_top <- kegg_result %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  separate_rows(geneID, sep = "/") %>%
  dplyr::select(Cluster, Description, geneID) %>%
  mutate(geneID = trimws(geneID))

if (nrow(alluvial_top) > 0) {
  # 统计每个 Cluster-Description 的基因数
  alluvial_count <- alluvial_top %>%
    group_by(Cluster, Description) %>%
    summarise(n = n(), .groups = "drop")
  
  # 设定通路顺序与点图一致（按 Cluster 内 p.adjust）
  pathway_order <- kegg_result %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  alluvial_count$Description <- factor(alluvial_count$Description, levels = pathway_order)
  
  p2 <- ggplot(alluvial_count,
               aes(axis1 = Cluster, axis2 = Description, y = n)) +
    geom_alluvium(aes(fill = Cluster), width = 1/12, alpha = 0.6) +
    geom_stratum(width = 1/12, fill = "grey90", color = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
    scale_x_discrete(limits = c("Gene Set", "KEGG Pathway"), expand = c(0.05, 0.05)) +
    theme_void() +
    theme(legend.position = "right") +
    labs(title = "Gene-to-KEGG pathway alluvial diagram",
         y = "Gene count")
  
  ggsave(p2, file = "kegg_compareCluster_alluvial.pdf", width = 10, height = 7)
  cat("Saved: kegg_compareCluster_alluvial.pdf\n")
}

cat("\nDone!\n")
