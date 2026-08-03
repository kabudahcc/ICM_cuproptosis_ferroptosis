# 整理汇总表（Excel）并绘制组合图（Figure Panel）
# 输入：前面所有分析结果文件
# 输出：ICM_cell_death_summary.xlsx, Figure_panel.pdf

rm(list = ls())
library(writexl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(pheatmap)
library(igraph)
library(ggraph)
library(stringr)

# 1. 创建 Excel 汇总表 ---------------------------------------------------

# Sheet 1: 分析概览
datasets <- c("GSE16499", "GSE5406", "GSE57338")
overview_list <- list()
for (ds in datasets) {
  group_df <- read.csv(paste0(ds, "_group.csv"), header = T, check.names = F)
  diff_df <- read.table(paste0(ds, "_diff.txt"), header = T, sep = "\t", check.names = F)
  crg_df <- read.table(paste0("CRG_DEG_", ds, ".txt"), header = T, sep = "\t", check.names = F)
  frg_df <- read.table(paste0("FRG_DEG_", ds, ".txt"), header = T, sep = "\t", check.names = F)
  
  overview_list[[ds]] <- data.frame(
    Dataset = ds,
    Platform = c("GPL5175", "GPL96", "GPL11532")[match(ds, datasets)],
    Grouping = c("Normal vs Ischemic HF", "nonfailing vs failing", "Non-failing vs Ischemic/Idiopathic CMP")[match(ds, datasets)],
    Total_samples = nrow(group_df),
    Normal_samples = sum(group_df$group == "Normal"),
    Disease_samples = sum(group_df$group == "Disease"),
    Total_DEGs = sum(abs(diff_df$logFC) > 0.5 & diff_df$P.Value < 0.05),
    CRG_DEGs = nrow(crg_df),
    FRG_DEGs = nrow(frg_df),
    stringsAsFactors = F
  )
}
overview_df <- bind_rows(overview_list)

# Sheet 2: 核心 5 基因 logFC
core_logfc <- read.table("cell_death_logFC_matrix.txt", header = T, sep = "\t", check.names = F, row.names = 1)
core_logfc$Gene <- rownames(core_logfc)
core_logfc <- core_logfc %>% select(Gene, everything())
core_logfc$Direction <- ifelse(rowMeans(core_logfc[, datasets], na.rm = T) > 0, "UP", "DOWN")
core_logfc$Category <- c("Cuproptosis", "Ferroptosis", "Cuproptosis", "Ferroptosis", "Ferroptosis")

# Sheet 3: PPI Hub 基因
hub_genes <- read.table("hub_genes.txt", header = T, sep = "\t", check.names = F)
hub_genes$Rank <- seq_len(nrow(hub_genes))
hub_genes <- hub_genes %>% select(Rank, everything())

# Sheet 4: ML 稳定特征基因
ml_stable <- read.table("ML_stable_signature_genes.txt", header = T, sep = "\t", check.names = F)
annotations <- c(
  PRDX6 = "Peroxiredoxin 6, antioxidant/lipid peroxidation",
  BCL6 = "B-cell lymphoma 6, transcription repressor",
  ATF4 = "Activating transcription factor 4, stress response",
  BRD4 = "Bromodomain-containing protein 4, inflammation/fibrosis",
  CBS = "Cystathionine beta-synthase, H2S production",
  STEAP3 = "STEAP3 metalloreductase, iron metabolism",
  TGM2 = "Transglutaminase 2, fibrosis/apoptosis",
  LTBP2 = "Latent TGF-beta binding protein 2, fibrosis",
  FBLN1 = "Fibulin-1, extracellular matrix/fibrosis",
  HSPD1 = "Heat shock protein D1, mitochondrial/cuproptosis",
  BEX1 = "Brain-expressed X-linked protein 1",
  CD38 = "CD38, NAD+ metabolism/immune",
  PLTP = "Phospholipid transfer protein",
  PTPRC = "Protein tyrosine phosphatase receptor type C, immune",
  SLC2A1 = "Glucose transporter 1",
  SLC40A1 = "Ferroportin, iron export",
  STC1 = "Stanniocalcin-1, calcium/metabolism",
  TIGAR = "TP53-inducible glycolysis regulator",
  ENPP2 = "Ectonucleotide pyrophosphatase 2",
  LPCAT3 = "Lysophosphatidylcholine acyltransferase 3, ferroptosis",
  NQO1 = "NADPH quinone oxidoreductase 1, antioxidant/ferroptosis",
  CDKN1A = "Cyclin-dependent kinase inhibitor 1A, p21",
  NOX4 = "NADPH oxidase 4, ROS production",
  PDK4 = "Pyruvate dehydrogenase kinase 4",
  POR = "Cytochrome P450 oxidoreductase",
  CYBB = "Cytochrome b-245 beta chain, NOX2",
  STAT3 = "Signal transducer and activator of transcription 3"
)
ml_stable$Annotation <- annotations[ml_stable$gene]

# Sheet 5: 免疫细胞差异（合并三个数据集）
immune_sig_list <- list()
for (ds in datasets) {
  df <- read.table(paste0("immune_wilcox_", ds, ".txt"), header = T, sep = "\t", check.names = F)
  df$Dataset <- ds
  immune_sig_list[[ds]] <- df
}
immune_sig <- bind_rows(immune_sig_list)
immune_sig$Direction <- ifelse(immune_sig$signif != "ns", "Up in Disease", "Not significant")

# Sheet 6: 基因-免疫平均相关
cor_avg <- read.table("immune_correlation_average.txt", header = T, sep = "\t", check.names = F, row.names = 1)
cor_avg$Gene <- rownames(cor_avg)
cor_avg_long <- cor_avg %>% pivot_longer(cols = -Gene, names_to = "CellType", values_to = "Spearman_r")

# Sheet 7: KEGG compareCluster top
kegg_comp <- read.table("KEGG_compareCluster.txt", header = T, sep = "\t", check.names = F, quote = "")
kegg_top <- kegg_comp %>%
  group_by(Cluster) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  select(Cluster, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID)

# Sheet 8: GO top
go_df <- read.table("GO_cell_death.txt", header = T, sep = "\t", check.names = F, quote = "")
go_top <- go_df %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>%
  head(20) %>%
  select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID)

# Sheet 9: 文件清单
file_list <- data.frame(
  File = c(
    "interGenes_three.txt", "CRG_intersect_three.txt", "FRG_intersect_three.txt",
    "cell_death_intersect_three.txt", "cell_death_logFC_matrix.txt",
    "CRG_DEG_GSE*.txt", "FRG_DEG_GSE*.txt",
    "GO_cell_death.txt", "KEGG_cell_death.txt", "KEGG_compareCluster.txt",
    "PPI_edges.txt", "PPI_nodes.txt", "hub_genes.txt",
    "immune_scores_all_datasets.txt", "immune_correlation_average.txt",
    "ML_signature_genes.txt", "ML_stable_signature_genes.txt", "model_performance.txt"
  ),
  Description = c(
    "三数据集总 DEGs 交集（54 个）",
    "三数据集共有铜凋亡 DEGs（MT1M）",
    "三数据集共有铁凋亡 DEGs（BCL6/LTBP2/STAT3/MYC）",
    "合并的细胞死亡相关交集基因（5 个）",
    "5 个核心基因在三个数据集中的 logFC",
    "各数据集铜凋亡差异基因",
    "各数据集铁凋亡差异基因",
    "5 核心基因 GO 富集结果",
    "5 核心基因 KEGG 富集结果",
    "CRG vs FRG compareCluster KEGG 结果",
    "PPI 网络边",
    "PPI 网络节点属性",
    "Top 10 Hub 基因",
    "三个数据集免疫细胞浸润评分",
    "三数据集平均免疫相关性",
    "ML 每折 CV 特征选择频率",
    "ML 跨数据集稳定特征基因",
    "三种模型性能指标"
  ),
  stringsAsFactors = F
)

# 写入 Excel
excel_list <- list(
  "分析概览" = overview_df,
  "核心5基因" = core_logfc,
  "PPI_Hub基因" = hub_genes,
  "ML稳定特征基因" = ml_stable,
  "免疫细胞差异" = immune_sig,
  "基因免疫相关" = cor_avg_long,
  "KEGG_compareCluster" = kegg_top,
  "GO_top" = go_top,
  "文件清单" = file_list
)
write_xlsx(excel_list, "ICM_cell_death_summary.xlsx")
cat("Saved: ICM_cell_death_summary.xlsx\n")

# 2. 绘制组合图 ----------------------------------------------------------

# Panel A: 核心 5 基因 logFC 热图
core_mat <- as.matrix(core_logfc[, datasets])
rownames(core_mat) <- core_logfc$Gene
core_mat <- core_mat[order(rowMeans(core_mat)), , drop = FALSE]

# 将热图转换为 ggplot 可用的数据
heatmap_df <- core_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Dataset", values_to = "logFC")
heatmap_df$Gene <- factor(heatmap_df$Gene, levels = rownames(core_mat))
heatmap_df$Dataset <- factor(heatmap_df$Dataset, levels = datasets)

pA <- ggplot(heatmap_df, aes(x = Dataset, y = Gene, fill = logFC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", logFC)), color = "black", size = 3) +
  scale_fill_gradient2(low = "#0073C2", mid = "white", high = "#CD534C", midpoint = 0,
                       limits = c(-1.5, 1.5), oob = scales::squish) +
  theme_minimal() +
  labs(title = "A. Core genes logFC", x = "", y = "") +
  theme(legend.position = "right")

# Panel B: 免疫细胞浸润显著性（GSE57338 最显著，选它展示）
immune_plot <- read.table("immune_wilcox_GSE57338.txt", header = T, sep = "\t", check.names = F)
immune_plot$logP <- -log10(immune_plot$pvalue)
immune_plot$CellType <- factor(immune_plot$CellType, levels = immune_plot$CellType[order(immune_plot$pvalue, decreasing = TRUE)])

pB <- ggplot(immune_plot, aes(x = CellType, y = logP, fill = signif)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("***" = "#CD534C", "**" = "#EFC000", "*" = "#0073C2", "ns" = "grey70")) +
  coord_flip() +
  theme_bw() +
  labs(title = "B. Immune infiltration (GSE57338)", x = "", y = "-log10(pvalue)", fill = "Signif") +
  theme(legend.position = "right")

# Panel C: 基因-免疫平均相关热图（简化为 top 相关条形图）
cor_top <- cor_avg_long %>%
  mutate(abs_r = abs(Spearman_r)) %>%
  arrange(desc(abs_r)) %>%
  head(15)
cor_top$Pair <- paste0(cor_top$Gene, " - ", cor_top$CellType)
cor_top$Pair <- factor(cor_top$Pair, levels = cor_top$Pair[order(cor_top$Spearman_r)])

pC <- ggplot(cor_top, aes(x = Pair, y = Spearman_r, fill = Spearman_r > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#CD534C", "FALSE" = "#0073C2"),
                    labels = c("Negative", "Positive")) +
  coord_flip() +
  theme_bw() +
  labs(title = "C. Top gene-immune correlations", x = "", y = "Spearman r", fill = "Direction") +
  theme(legend.position = "right")

# Panel D: ROC 曲线（用之前保存的性能数据重新绘图）
roc_data <- read.table("model_performance.txt", header = T, sep = "\t", check.names = F)
# 为了 ROC 曲线需要重新构建，这里用 AUC 画柱状图代替
d_perf <- roc_data %>%
  select(Dataset, Model, AUC) %>%
  mutate(Model = factor(Model, levels = c("LASSO-LR", "RandomForest", "SVM")))

pD <- ggplot(d_perf, aes(x = Dataset, y = AUC, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", AUC)), position = position_dodge(width = 0.8), vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c("LASSO-LR" = "#0073C2", "RandomForest" = "#EFC000", "SVM" = "#CD534C")) +
  ylim(0, 1.1) +
  theme_bw() +
  labs(title = "D. Diagnostic model performance (AUC)", x = "", y = "AUC", fill = "Model") +
  theme(legend.position = "right")

# Panel E: ML 稳定特征基因频率
ml_plot <- ml_stable %>%
  head(15) %>%
  mutate(gene = factor(gene, levels = gene[order(total_freq)]))

pE <- ggplot(ml_plot, aes(x = gene, y = total_freq, fill = total_freq)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#EFC000", high = "#CD534C") +
  coord_flip() +
  theme_bw() +
  labs(title = "E. Stable ML signature genes", x = "", y = "Selection frequency") +
  theme(legend.position = "right")

# Panel F: KEGG compareCluster top 通路（气泡图）
kegg_fig <- kegg_comp %>%
  group_by(Cluster) %>%
  arrange(pvalue) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(Description = str_remove(Description, " - Homo sapiens \\(human\\)"),
         Description = factor(Description, levels = unique(Description[order(pvalue, decreasing = TRUE)])))

pF <- ggplot(kegg_fig, aes(x = Cluster, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", trans = "log10") +
  theme_bw() +
  labs(title = "F. KEGG pathways (CRG vs FRG)", x = "", y = "", size = "Gene count", color = "p.adjust") +
  theme(legend.position = "right")

# 组合图
combined_figure <- (pA + pB) / (pC + pD) / (pE + pF) +
  plot_annotation(title = "Figure: Cuproptosis/ferroptosis-related genes in ischemic cardiomyopathy",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))

ggsave(combined_figure, file = "Figure_panel.pdf", width = 14, height = 16)
cat("Saved: Figure_panel.pdf\n")

cat("\nDone!\n")
