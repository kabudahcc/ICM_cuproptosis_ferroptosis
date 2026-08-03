# 核心基因与免疫细胞浸润的相关性分析
# 输入：step6 生成的 immune_scores_*.txt 和 *_symbol.csv
# 输出：相关性热图、气泡图、统计表格

rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 1. 设置 -----------------------------------------------------------------
datasets <- c("GSE16499", "GSE5406", "GSE57338")
core_genes <- c("BCL6", "LTBP2", "MT1M", "MYC", "STAT3")

# 2. 读取免疫评分 ---------------------------------------------------------
immune_all <- data.frame()
for (ds in datasets) {
  imm <- read.table(paste0("immune_scores_", ds, ".txt"), header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  immune_all <- rbind(immune_all, imm)
}

# 3. 读取表达矩阵并提取核心基因 ------------------------------------------
expr_all <- data.frame()
for (ds in datasets) {
  expr <- read.csv(paste0(ds, "_symbol.csv"), check.names = FALSE, row.names = 1)
  expr <- as.matrix(expr)
  
  present <- intersect(core_genes, rownames(expr))
  if (length(present) == 0) next
  
  expr_sub <- expr[present, , drop = FALSE]
  df <- as.data.frame(t(expr_sub))
  df$sample <- rownames(df)
  df$dataset <- ds
  expr_all <- rbind(expr_all, df)
}

# 4. 合并数据 -------------------------------------------------------------
merged <- merge(immune_all, expr_all, by = c("sample", "dataset"))

# 5. 计算相关性 -----------------------------------------------------------
immune_cols <- setdiff(colnames(immune_all), c("sample", "group", "dataset"))
cor_results <- data.frame()

for (ds in datasets) {
  sub <- merged[merged$dataset == ds, ]
  if (nrow(sub) == 0) next
  
  for (gene in core_genes) {
    if (!(gene %in% colnames(sub))) next
    for (cell in immune_cols) {
      x <- as.numeric(sub[, cell])
      y <- as.numeric(sub[, gene])
      if (sd(x) == 0 || sd(y) == 0) next
      
      ct <- cor.test(x, y, method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        dataset = ds,
        gene = gene,
        cell = cell,
        cor = ct$estimate,
        pvalue = ct$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cor_results$padj <- p.adjust(cor_results$pvalue, method = "BH")
write.table(cor_results, "core_gene_immune_correlation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 6. 计算跨数据集平均相关性 ----------------------------------------------
avg_cor <- cor_results %>%
  group_by(gene, cell) %>%
  summarise(mean_cor = mean(cor),
            mean_p = mean(pvalue),
            n = n(),
            .groups = "drop")

# 标记显著性
avg_cor$significance <- ifelse(avg_cor$mean_p < 0.001, "***",
                               ifelse(avg_cor$mean_p < 0.01, "**",
                                      ifelse(avg_cor$mean_p < 0.05, "*", "ns")))
write.table(avg_cor, "core_gene_immune_correlation_average.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 7. 可视化 ---------------------------------------------------------------
# 图 1：每个数据集的相关性热图
plot_list <- list()
for (ds in datasets) {
  sub <- cor_results[cor_results$dataset == ds, ]
  if (nrow(sub) == 0) next
  
  p <- ggplot(sub, aes(x = gene, y = cell, fill = cor)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limit = c(-1, 1), name = "Spearman r") +
    geom_text(aes(label = ifelse(pvalue < 0.05, "*", "")), size = 4, color = "black") +
    theme_minimal() +
    labs(title = ds, x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_list[[ds]] <- p
}

combined_heatmap <- wrap_plots(plot_list, ncol = 1)
ggsave("core_gene_immune_correlation_heatmap.pdf", combined_heatmap, width = 8, height = 14)

# 图 2：平均相关性热图
p_avg <- ggplot(avg_cor, aes(x = gene, y = cell, fill = mean_cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limit = c(-1, 1), name = "Mean\nSpearman r") +
  geom_text(aes(label = significance), size = 4, color = "black") +
  theme_minimal() +
  labs(title = "Average correlation across three datasets", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("core_gene_immune_correlation_average_heatmap.pdf", p_avg, width = 7, height = 6)

# 图 3：气泡图（跨数据集平均）
p_bubble <- ggplot(avg_cor, aes(x = gene, y = cell, size = abs(mean_cor), color = mean_cor)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        limit = c(-1, 1), name = "Mean r") +
  scale_size_continuous(range = c(2, 10), name = "|r|") +
  geom_text(aes(label = significance), color = "black", size = 3, vjust = -0.8) +
  theme_bw() +
  labs(title = "Core genes vs immune infiltration", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("core_gene_immune_correlation_bubble.pdf", p_bubble, width = 8, height = 6)

cat("Done!\n")
cat("Saved: core_gene_immune_correlation.txt\n")
cat("Saved: core_gene_immune_correlation_average.txt\n")
cat("Saved: core_gene_immune_correlation_heatmap.pdf\n")
cat("Saved: core_gene_immune_correlation_average_heatmap.pdf\n")
cat("Saved: core_gene_immune_correlation_bubble.pdf\n")
