# WGCNA 6 panel 组合图（纯矢量版，直接绘制后拼接）
rm(list = ls())
options(stringsAsFactors = FALSE)

library(WGCNA)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)
library(ggplotify)

# 加载 WGCNA 中间结果（包含 datExpr、MEs、moduleColors 等）
load("WGCNA_results_part1.RData")

cat("Top modules:", paste(top_modules, collapse = ", "), "\n")

# 数字标签 -> 颜色映射
ME_to_color <- setNames(WGCNA::labels2colors(as.numeric(gsub("ME", "", top_modules))),
                        top_modules)

# ---------- Panel A：样本聚类 + 性状热图（base graphics） ----------
pA <- as.ggplot(function() {
  sampleTree <- hclust(dist(datExpr), method = "average")
  traitColors <- numbers2colors(datTraits, signed = TRUE)
  plotDendroAndColors(sampleTree, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap",
                      dendroLabels = FALSE,
                      marAll = c(1, 5, 3, 1))
})

# ---------- Panel B：Top 模块 GO-BP 富集横向条形图 ----------
go_plot_list <- list()
for (me in top_modules) {
  f <- paste0("WGCNA_GO_module_", me, ".txt")
  if (!file.exists(f)) next
  df <- read.delim(f, check.names = FALSE)
  if (nrow(df) == 0) next
  df <- head(df[order(df$pvalue), ], 5)
  df$logp <- -log10(df$pvalue)
  mod_color <- ME_to_color[me]
  p <- ggplot(df, aes(x = logp, y = reorder(Description, logp))) +
    geom_bar(stat = "identity", fill = mod_color, alpha = 0.7) +
    theme_bw() +
    labs(title = paste(me, "(", mod_color, ")"),
         x = expression(-log[10]~(p-value)), y = "") +
    theme(axis.text.y = element_text(size = 7),
          title = element_text(size = 9))
  go_plot_list[[me]] <- p
}
pB <- wrap_plots(go_plot_list, ncol = 2)

# ---------- Panel C：Top 模块 eigengene 箱线图 ----------
ME_plot <- MEs[, top_modules, drop = FALSE]
ME_plot$Group <- factor(group$group, levels = c("Normal", "Disease"))
ME_long <- gather(ME_plot, module, eigengene, -Group)

pC <- ggplot(ME_long, aes(x = Group, y = eigengene, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) +
  facet_wrap(~ module, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#0073C2", "Disease" = "#CD534C")) +
  theme_bw() +
  labs(title = "Module eigengenes by group", y = "Module eigengene") +
  theme(legend.position = "none",
        axis.title.x = element_blank())

# ---------- Panel D：Top 模块基因表达热图 ----------
selected_genes <- c()
for (i in seq_along(top_modules)) {
  mod_color <- ME_to_color[i]
  g <- module_df$gene[module_df$module == mod_color]
  selected_genes <- c(selected_genes, head(g, 50))
}
selected_genes <- selected_genes[selected_genes %in% rownames(expr_top)]
heat_mat <- expr_top[selected_genes, , drop = FALSE]
heat_mat <- heat_mat[apply(heat_mat, 1, function(x) all(is.finite(x))), ]
heat_mat <- heat_mat[apply(heat_mat, 1, sd) > 0, ]

heat_data_sub <- t(heat_mat)
heat_anno <- data.frame(Group = group$group, row.names = rownames(heat_data_sub))

pD <- as.ggplot(~ pheatmap(t(heat_data_sub),
                           annotation_col = heat_anno,
                           cluster_cols = TRUE,
                           cluster_rows = TRUE,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           scale = "row",
                           color = colorRampPalette(c("blue", "white", "red"))(50),
                           main = "Module gene expression heatmap"))

# ---------- Panel E：模块-性状相关性热图 ----------
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
                       limit = c(-1, 1), name = "r") +
  theme_minimal() +
  labs(title = "Module-trait relationships", x = "", y = "Module") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------- Panel F：模块网络 PCA ----------
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
  labs(title = "Module network (PCA)", x = "PC1", y = "PC2") +
  theme(legend.position = "right")

# ---------- 组合 6 panel ----------
# 把 Panel B（patchwork）封装为一个整体 grob，避免内部子图也被打标签
pB_element <- wrap_elements(full = pB)

# 使用 3 行 x 2 列，添加 A-F 标签
p_all <- wrap_plots(pA, pB_element, pC, pD, pE, pF,
                    ncol = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# 输出纯矢量 PDF（无 raster 转换）
ggsave("WGCNA_6panel_combined_vector.pdf", p_all,
       width = 16, height = 22, device = cairo_pdf)
# 同时输出高分辨率 PNG 供预览
ggsave("WGCNA_6panel_combined_vector.png", p_all,
       width = 16, height = 22, dpi = 300)

cat("Saved: WGCNA_6panel_combined_vector.pdf\n")
cat("Saved: WGCNA_6panel_combined_vector.png\n")
