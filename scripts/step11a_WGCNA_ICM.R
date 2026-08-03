# WGCNA 分析：GSE57338（ICM + DCM + Normal）
# 目标：识别与心肌病相关的共表达模块
# 输出：WGCNA 结果、模块-性状相关性、组合图

rm(list = ls())
options(stringsAsFactors = FALSE)

# 第一部分：只加载 WGCNA 和基础包，避免与 Bioconductor 包的 cor 函数冲突
library(WGCNA)
library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)

# 1. 读取数据 -------------------------------------------------------------
expr_file <- "GSE57338_symbol.csv"
group_file <- "GSE57338_group.csv"

expr <- read.csv(expr_file, check.names = FALSE, row.names = 1)
expr <- as.matrix(expr)
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", ]

group <- read.csv(group_file, check.names = FALSE, row.names = 1)
common_samples <- intersect(colnames(expr), rownames(group))
expr <- expr[, common_samples]
group <- group[common_samples, ]

cat("Samples:", ncol(expr), "Genes:", nrow(expr), "\n")
cat("Group distribution:\n")
print(table(group$group))

# 2. 提取 Ischemic/Idiopathic 信息 ---------------------------------------
group$ischemic <- ifelse(grepl("[Ii]schem", group$title), "Ischemic", 
                         ifelse(grepl("[Ii]diopathic", group$title), "Idiopathic", "Normal"))
cat("Ischemic status:\n")
print(table(group$ischemic))

# 3. 过滤基因 -------------------------------------------------------------
expr <- expr[apply(expr, 1, function(x) sum(x > 0) > 0.5 * ncol(expr)), ]
expr <- avereps(expr, ID = rownames(expr))
mads <- apply(expr, 1, mad)
expr_top <- expr[order(mads, decreasing = TRUE)[1:min(5000, nrow(expr))], ]

cat("After filtering - Genes:", nrow(expr_top), "\n")

datExpr <- t(expr_top)

# 4. 构建性状矩阵 ---------------------------------------------------------
trait_disease <- ifelse(group$group == "Disease", 1, 0)
trait_ischemic <- ifelse(group$ischemic == "Ischemic", 1, 0)
trait_idiopathic <- ifelse(group$ischemic == "Idiopathic", 1, 0)

datTraits <- data.frame(
  Disease = trait_disease,
  Ischemic = trait_ischemic,
  Idiopathic = trait_idiopathic,
  row.names = rownames(group)
)

# 5. 样本聚类 -------------------------------------------------------------
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  datTraits <- datTraits[gsg$goodSamples, ]
}

sampleTree <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = TRUE)

pdf("WGCNA_A_sample_dendrogram_traits.pdf", width = 10, height = 6)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    dendroLabels = FALSE)
dev.off()
cat("Saved: WGCNA_A_sample_dendrogram_traits.pdf\n")

# 6. 选择软阈值 -----------------------------------------------------------
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

pdf("WGCNA_soft_threshold.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
softPower <- sft$fitIndices[which(fit >= 0.9)[1], 1]
if (is.na(softPower) || softPower < 1) softPower <- 6
cat("Selected soft power:", softPower, "\n")

# 7. 构建网络 -------------------------------------------------------------
net <- blockwiseModules(datExpr,
                        power = softPower,
                        corType = "pearson",
                        maxPOutliers = 0,
                        TOMType = "signed",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

cat("Modules found:", length(unique(moduleColors)), "\n")
print(table(moduleColors))

# 保存模块基因
module_df <- data.frame(gene = colnames(datExpr), module = moduleColors, stringsAsFactors = FALSE)
write.table(module_df, "WGCNA_module_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 8. 模块-性状相关性 ------------------------------------------------------
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

write.table(round(moduleTraitCor, 3), "WGCNA_module_trait_correlation.txt", sep = "\t", quote = FALSE)
write.table(signif(moduleTraitPvalue, 3), "WGCNA_module_trait_pvalue.txt", sep = "\t", quote = FALSE)

# 9. 识别 top 模块 --------------------------------------------------------
disease_cor <- moduleTraitCor[, "Disease"]
top_modules <- names(disease_cor)[order(abs(disease_cor), decreasing = TRUE)][1:6]
cat("Top modules associated with Disease:\n")
print(disease_cor[top_modules])

# 10. 保存中间结果 --------------------------------------------------------
save(MEs, moduleColors, module_df, moduleTraitCor, moduleTraitPvalue, datExpr, datTraits,
     expr_top, group, top_modules, net, file = "WGCNA_results_part1.RData")
cat("Saved: WGCNA_results_part1.RData\n")

# 11. 生成可视化（不含 GO）------------------------------------------------
# 图 E：模块-性状相关性热图
cor_df <- as.data.frame(moduleTraitCor)
cor_df$module <- rownames(cor_df)
cor_long <- gather(cor_df, trait, correlation, -module)

pval_df <- as.data.frame(moduleTraitPvalue)
pval_df$module <- rownames(pval_df)
pval_long <- gather(pval_df, trait, pvalue, -module)

plot_df <- left_join(cor_long, pval_long, by = c("module", "trait"))
plot_df$star <- ifelse(plot_df$pvalue < 0.001, "***",
                       ifelse(plot_df$pvalue < 0.01, "**",
                              ifelse(plot_df$pvalue < 0.05, "*", "")))
plot_df$label <- paste0(round(plot_df$correlation, 2), "\n", plot_df$star)
plot_df$label[plot_df$star == ""] <- ""

pE <- ggplot(plot_df, aes(x = trait, y = module, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limit = c(-1, 1)) +
  theme_minimal() +
  labs(title = "Module-trait relationships", x = "", y = "Module") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("WGCNA_E_module_trait_correlation.pdf", pE, width = 6, height = 8)

# 图 C：Top 模块 eigengene 箱线图
ME_plot <- MEs[, top_modules, drop = FALSE]
ME_plot$Group <- factor(group$group, levels = c("Normal", "Disease"))
ME_long <- gather(ME_plot, module, eigengene, -Group)

pC <- ggplot(ME_long, aes(x = Group, y = eigengene, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5) +
  facet_wrap(~ module, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#0073C2", "Disease" = "#CD534C")) +
  theme_bw() +
  labs(title = "Module eigengenes by group", y = "Module eigengene") +
  theme(legend.position = "none")
ggsave("WGCNA_C_module_eigengene_boxplot.pdf", pC, width = 10, height = 6)

# 图 D：Top 模块基因表达热图
# 将 ME 数字标签映射为颜色
ME_to_color <- sapply(top_modules, function(me) {
  mod_num <- as.numeric(gsub("ME", "", me))
  unique(moduleColors[net$colors == mod_num])[1]
})

selected_genes <- c()
for (i in seq_along(top_modules)) {
  mod_color <- ME_to_color[i]
  g <- module_df$gene[module_df$module == mod_color]
  selected_genes <- c(selected_genes, head(g, 50))
}
selected_genes <- selected_genes[selected_genes %in% rownames(expr_top)]
heat_mat <- expr_top[selected_genes, , drop = FALSE]
heat_mat <- heat_mat[apply(heat_mat, 1, function(x) all(is.finite(x))), ]

# 去除方差为 0 的行，避免 scale = "row" 产生 Inf
heat_mat <- heat_mat[apply(heat_mat, 1, sd) > 0, ]

if (nrow(heat_mat) > 10) {
  heat_data_sub <- t(heat_mat)
  heat_anno <- data.frame(Group = group$group, row.names = rownames(heat_data_sub))
  
  pdf("WGCNA_D_module_heatmap.pdf", width = 10, height = 8)
  pheatmap(t(heat_data_sub),
           annotation_col = heat_anno,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           show_rownames = FALSE,
           show_colnames = FALSE,
           scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "Module gene expression heatmap")
  dev.off()
  cat("Saved: WGCNA_D_module_heatmap.pdf\n")
} else {
  cat("Warning: Too few valid genes for module heatmap\n")
}

# 图 F：模块网络 PCA（基因维度）------------------------------------------
# datExpr: 行是样本，列是基因；转置后 PCA 看基因间关系
pca <- prcomp(t(datExpr), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  module = moduleColors,
  gene = colnames(datExpr)
)
MM <- as.data.frame(cor(datExpr, MEs, use = "p"))
hub_genes_idx <- c()
for (mod in top_modules) {
  mm_col <- paste0("ME", mod)
  if (mm_col %in% colnames(MM)) {
    idx <- which(moduleColors == mod)
    scores <- MM[idx, mm_col]
    top <- idx[order(scores, decreasing = TRUE)[1:5]]
    hub_genes_idx <- c(hub_genes_idx, top)
  }
}
pca_df$is_hub <- pca_df$gene %in% colnames(datExpr)[hub_genes_idx]

pF <- ggplot(pca_df, aes(x = PC1, y = PC2, color = module)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_text_repel(data = pca_df[pca_df$is_hub, ],
                  aes(label = gene), size = 3, max.overlaps = 30) +
  theme_minimal() +
  labs(title = "Module network (PCA)", x = "PC1", y = "PC2")
ggsave("WGCNA_F_module_network.pdf", pF, width = 8, height = 6)

# 12. Hub genes -----------------------------------------------------------
MM_full <- signedKME(datExpr, MEs)
top_hub_genes <- data.frame()
for (me in colnames(MEs)) {
  mod_num <- as.numeric(gsub("ME", "", me))
  mod_color <- unique(moduleColors[net$colors == mod_num])[1]
  col_name <- paste0("kME", mod_num)
  if (col_name %in% colnames(MM_full)) {
    idx <- order(MM_full[, col_name], decreasing = TRUE)[1:10]
    top_hub_genes <- rbind(top_hub_genes,
                           data.frame(module = paste0(me, "_", mod_color),
                                      gene = rownames(MM_full)[idx],
                                      kME = MM_full[idx, col_name],
                                      stringsAsFactors = FALSE))
  }
}
write.table(top_hub_genes, "WGCNA_top10_hub_genes_per_module.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nWGCNA Part 1 done!\n")
