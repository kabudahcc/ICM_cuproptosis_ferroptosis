# WGCNA 模块 GO 富集（第二部分，单独运行避免包冲突）
rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 读取 WGCNA 结果
load("WGCNA_results_part1.RData")

cat("Top modules:", paste(top_modules, collapse = ", "), "\n")

# top_modules 是 ME0, ME1 等 WGCNA 数字标签，使用 labels2colors 映射为颜色
ME_to_color <- setNames(WGCNA::labels2colors(as.numeric(gsub("ME", "", top_modules))),
                        top_modules)
cat("Module color mapping:\n")
print(ME_to_color)

# GO 富集
Module_GO <- list()
for (i in seq_along(top_modules)) {
  me <- top_modules[i]
  mod_color <- ME_to_color[me]
  genes <- module_df$gene[module_df$module == mod_color]
  cat("Module", me, "(", mod_color, "):", length(genes), "genes\n")
  
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (is.null(entrez) || nrow(entrez) < 5) {
    cat("  Skipping - too few mappable genes\n")
    next
  }
  
  ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1,
                  readable = TRUE)
  
  ego_df <- ego@result
  ego_df <- ego_df[ego_df$pvalue < 0.05 & !is.na(ego_df$pvalue), ]
  if (!is.null(ego) && nrow(ego_df) > 0) {
    Module_GO[[me]] <- ego_df
    write.table(ego_df, paste0("WGCNA_GO_module_", me, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Significant GO terms:", nrow(ego_df), "\n")
  } else {
    cat("  No significant GO terms\n")
  }
}

# 合并所有模块 GO 结果（Top 10），方便论文附表
GO_summary <- dplyr::bind_rows(lapply(names(Module_GO), function(me) {
  df <- Module_GO[[me]]
  df <- df[order(df$pvalue), ]
  df$Module <- me
  df$ModuleColor <- ME_to_color[me]
  head(df, 10)
}), .id = NULL)

if (nrow(GO_summary) > 0) {
  write.table(GO_summary, "WGCNA_GO_modules_top10_summary.txt",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved: WGCNA_GO_modules_top10_summary.txt\n")
}

# 绘制 GO 富集图
go_plot_list <- list()
for (me in names(Module_GO)) {
  df <- head(Module_GO[[me]], 5)
  if (nrow(df) == 0) next
  df$logp <- -log10(df$pvalue)
  mod_color <- ME_to_color[me]
  p <- ggplot(df, aes(x = logp, y = reorder(Description, logp))) +
    geom_bar(stat = "identity", fill = mod_color, alpha = 0.7) +
    theme_bw() +
    labs(title = paste(me, "(", mod_color, ")"), x = "-log10(p-value)", y = "") +
    theme(axis.text.y = element_text(size = 8))
  go_plot_list[[me]] <- p
}

if (length(go_plot_list) > 0) {
  pB <- wrap_plots(go_plot_list, ncol = 2)
  ggsave("WGCNA_B_GO_enrichment.pdf", pB, width = 12, height = 10)
  cat("Saved: WGCNA_B_GO_enrichment.pdf\n")
} else {
  cat("No GO terms to plot\n")
}

cat("\nWGCNA GO enrichment done!\n")
