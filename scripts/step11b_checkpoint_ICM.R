# 免疫检查点基因在 ICM vs Normal 中的差异表达
# 输入：三个数据集的 symbol 表达矩阵和 group 文件
# 输出：箱线图 + 统计结果

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

# 1. 设置数据 -------------------------------------------------------------
datasets <- c("GSE16499", "GSE5406", "GSE57338")

# 免疫检查点基因（可扩展）
checkpoint_genes <- c("CD274", "CTLA4", "HAVCR2", "TIGIT", "LAG3",
                      "PDCD1", "PDCD1LG2", "SIGLEC15", "ICOS", "ICOSLG",
                      "CD80", "CD86", "CD276", "VTCN1")

# 2. 读取并合并数据 -------------------------------------------------------
all_data <- data.frame()

for (ds in datasets) {
  expr <- read.csv(paste0(ds, "_symbol.csv"), check.names = FALSE, row.names = 1)
  expr <- as.matrix(expr)
  
  group <- read.csv(paste0(ds, "_group.csv"), check.names = FALSE, row.names = 1)
  common <- intersect(colnames(expr), rownames(group))
  expr <- expr[, common]
  group <- group[common, ]
  
  # 只保留 Normal 和 Disease
  keep <- group$group %in% c("Normal", "Disease")
  expr <- expr[, keep]
  group <- group[keep, ]
  
  # 提取检查点基因
  present <- intersect(checkpoint_genes, rownames(expr))
  if (length(present) == 0) next
  
  expr_sub <- expr[present, , drop = FALSE]
  
  for (gene in present) {
    df <- data.frame(
      dataset = ds,
      sample = colnames(expr_sub),
      gene = gene,
      expression = as.numeric(expr_sub[gene, ]),
      group = group$group,
      stringsAsFactors = FALSE
    )
    all_data <- rbind(all_data, df)
  }
}

all_data$group <- factor(all_data$group, levels = c("Normal", "Disease"))
all_data$dataset <- factor(all_data$dataset, levels = datasets)

# 3. 统计检验 -------------------------------------------------------------
stat_results <- data.frame()
for (ds in datasets) {
  for (gene in checkpoint_genes) {
    sub <- all_data[all_data$dataset == ds & all_data$gene == gene, ]
    if (nrow(sub) == 0) next
    normal <- sub$expression[sub$group == "Normal"]
    disease <- sub$expression[sub$group == "Disease"]
    if (length(normal) < 2 || length(disease) < 2) next
    
    test <- wilcox.test(disease, normal)
    logFC <- log2(mean(disease) / mean(normal))
    stat_results <- rbind(stat_results, data.frame(
      dataset = ds,
      gene = gene,
      normal_mean = mean(normal),
      disease_mean = mean(disease),
      logFC = logFC,
      pvalue = test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}
stat_results$padj <- p.adjust(stat_results$pvalue, method = "BH")
write.table(stat_results, "immune_checkpoint_stat_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 4. 绘制各数据集图 -------------------------------------------------------
plot_list <- list()
for (ds in datasets) {
  sub <- all_data[all_data$dataset == ds, ]
  if (nrow(sub) == 0) next
  
  p <- ggboxplot(sub, x = "gene", y = "expression", fill = "group",
                 palette = c("Normal" = "#0073C2", "Disease" = "#CD534C"),
                 xlab = "", ylab = "Expression") +
    rotate_x_text(45) +
    stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif") +
    facet_wrap(~ dataset) +
    theme(legend.position = "right")
  plot_list[[ds]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 1)
ggsave("immune_checkpoint_boxplot_by_dataset.pdf", combined_plot, width = 12, height = 14)

# 5. 汇总图：每个基因在三数据集中的平均表达 -----------------------------
summary_data <- all_data %>%
  group_by(dataset, gene, group) %>%
  summarise(mean_expr = mean(expression), .groups = "drop")

p_summary <- ggplot(summary_data, aes(x = gene, y = mean_expr, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#0073C2", "Disease" = "#CD534C")) +
  theme_bw() +
  rotate_x_text(45) +
  labs(title = "Immune checkpoint gene expression across datasets",
       y = "Mean expression", x = "")

ggsave("immune_checkpoint_summary_barplot.pdf", p_summary, width = 12, height = 10)

cat("Done!\n")
cat("Saved: immune_checkpoint_stat_results.txt\n")
cat("Saved: immune_checkpoint_boxplot_by_dataset.pdf\n")
cat("Saved: immune_checkpoint_summary_barplot.pdf\n")
