# 免疫浸润分析（GSVA + MCPcounter signatures）
# 输入：GSE*_symbol.csv, GSE*_group.csv, cell_death_intersect_three.txt
# 输出：immune_scores_*.txt, immune_boxplot_*.pdf, immune_correlation_*.pdf

rm(list = ls())
library(GSVA)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# 1. MCPcounter 免疫细胞标志物（已发表标志物，直接内置于脚本） -----------
mcpsig <- list(
  T_cells = c("CD3D", "CD3E", "CD2", "CD5", "CD6"),
  CD8_T_cells = c("CD8A", "CD8B"),
  Cytotoxic_lymphocytes = c("EOMES", "GZMA", "GZMB", "PRF1", "GNLY"),
  NK_cells = c("NCR1", "NCR2", "NCR3", "KIR2DL3", "KIR3DL1", "KIR3DL2"),
  B_lineage = c("CD19", "CD79A", "CD79B", "MS4A1"),
  Monocytic_lineage = c("CD14", "CD163", "CD68", "CSF1R", "FCGR1A", "ITGAM"),
  Myeloid_dendritic_cells = c("CD1C", "CLEC9A", "CD209"),
  Neutrophils = c("CEACAM3", "FCGR3B", "S100A12", "CSF3R"),
  Endothelial_cells = c("VWF", "CDH5", "PECAM1", "ENG"),
  Fibroblasts = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "LUM", "FN1")
)

cat("Immune signatures:\n")
print(sapply(mcpsig, length))

# 2. 读取表达矩阵和分组 ---------------------------------------------------
datasets <- c("GSE16499", "GSE5406", "GSE57338")
key_genes <- read.table("cell_death_intersect_three.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]

all_immune_results <- list()
all_cor_results <- list()

for (ds in datasets) {
  cat("\n==========", ds, "==========\n")
  
  # 读取表达矩阵
  exp <- read.csv(paste0(ds, "_symbol.csv"), header = T, row.names = 1, check.names = F)
  exp <- as.matrix(exp)
  
  # 读取分组
  group_df <- read.csv(paste0(ds, "_group.csv"), header = T, check.names = F)
  group_df$group <- trimws(as.character(group_df$group))
  rownames(group_df) <- group_df$geo_accession
  
  # 确保分组和表达矩阵样本一致
  common_samples <- intersect(colnames(exp), group_df$geo_accession)
  exp <- exp[, common_samples]
  group_df <- group_df[common_samples, ]
  
  cat("Samples:", length(common_samples), "\n")
  cat("Group table:\n")
  print(table(group_df$group))
  
  # 3. GSVA 免疫浸润评分 -------------------------------------------------
  # GSVA >= 1.50 新版 API
  gsva_param <- gsvaParam(exprData = exp, geneSets = mcpsig, kcdf = "Gaussian")
  gsva_result <- gsva(gsva_param, verbose = FALSE)
  
  immune_scores <- as.data.frame(t(gsva_result))
  immune_scores$sample <- rownames(immune_scores)
  immune_scores$group <- group_df[rownames(immune_scores), "group"]
  immune_scores$dataset <- ds
  
  write.table(immune_scores, paste0("immune_scores_", ds, ".txt"), sep = "\t", quote = F, row.names = F)
  all_immune_results[[ds]] <- immune_scores
  
  # 4. 比较 Normal vs Disease --------------------------------------------
  cell_types <- setdiff(colnames(immune_scores), c("sample", "group", "dataset"))
  
  # 转换为长格式
  immune_long <- immune_scores %>%
    select(-dataset) %>%
    pivot_longer(cols = all_of(cell_types), names_to = "CellType", values_to = "Score")
  
  # Wilcoxon 检验
  test_results <- data.frame(CellType = cell_types, pvalue = NA, stringsAsFactors = F)
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    x <- immune_scores[immune_scores$group == "Normal", ct]
    y <- immune_scores[immune_scores$group == "Disease", ct]
    if (length(x) > 1 && length(y) > 1) {
      test_results$pvalue[i] <- wilcox.test(x, y)$p.value
    }
  }
  test_results$signif <- ifelse(test_results$pvalue < 0.001, "***",
                                 ifelse(test_results$pvalue < 0.01, "**",
                                        ifelse(test_results$pvalue < 0.05, "*", "ns")))
  write.table(test_results, paste0("immune_wilcox_", ds, ".txt"), sep = "\t", quote = F, row.names = F)
  
  # 箱线图
  p1 <- ggplot(immune_long, aes(x = CellType, y = Score, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(position = position_dodge(width = 0.75), size = 0.8, alpha = 0.5) +
    scale_fill_manual(values = c("Normal" = "#0073C2", "Disease" = "#CD534C")) +
    labs(title = paste(ds, "immune infiltration"), x = "", y = "GSVA score", fill = "Group") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = test_results, aes(x = CellType, y = max(immune_long$Score) * 1.1,
                                        label = signif, fill = NULL),
              inherit.aes = FALSE, size = 3)
  
  ggsave(p1, file = paste0("immune_boxplot_", ds, ".pdf"), width = 10, height = 6)
  cat("Saved: immune_boxplot_", ds, ".pdf\n", sep = "")
  
  # 5. 关键基因与免疫细胞相关性 ------------------------------------------
  # 提取关键基因表达
  key_exp <- exp[key_genes, , drop = FALSE]
  key_exp <- key_exp[rownames(key_exp) %in% rownames(exp), ]
  
  if (nrow(key_exp) > 0) {
    # 合并免疫评分和关键基因表达
    combined <- cbind(t(key_exp), immune_scores[, cell_types])
    colnames(combined)[1:nrow(key_exp)] <- rownames(key_exp)
    
    # 计算相关性
    cor_mat <- cor(combined, method = "spearman", use = "pairwise.complete.obs")
    cor_key_immune <- cor_mat[rownames(key_exp), cell_types, drop = FALSE]
    
    # p 值矩阵
    p_mat <- matrix(NA, nrow = nrow(cor_key_immune), ncol = ncol(cor_key_immune))
    rownames(p_mat) <- rownames(cor_key_immune)
    colnames(p_mat) <- colnames(cor_key_immune)
    for (g in rownames(cor_key_immune)) {
      for (c in colnames(cor_key_immune)) {
        test <- cor.test(combined[, g], combined[, c], method = "spearman")
        p_mat[g, c] <- test$p.value
      }
    }
    
    write.table(cor_key_immune, paste0("immune_correlation_", ds, ".txt"), sep = "\t", quote = F)
    write.table(p_mat, paste0("immune_correlation_pvalue_", ds, ".txt"), sep = "\t", quote = F)
    all_cor_results[[ds]] <- list(cor = cor_key_immune, p = p_mat)
    
    # 热图
    if (nrow(cor_key_immune) > 0 && ncol(cor_key_immune) > 0) {
      bk <- c(seq(-1, -0.01, length = 50), seq(0.01, 1, length = 50))
      col <- colorRampPalette(c("#0073C2", "white", "#CD534C"))(99)
      
      pdf(paste0("immune_correlation_heatmap_", ds, ".pdf"), width = 8, height = 4)
      pheatmap(cor_key_immune,
               scale = "none",
               color = col,
               breaks = bk,
               cluster_cols = T,
               cluster_rows = T,
               display_numbers = TRUE,
               number_color = "black",
               fontsize_number = 7,
               main = paste(ds, "correlation: key genes vs immune cells"))
      dev.off()
      cat("Saved: immune_correlation_heatmap_", ds, ".pdf\n", sep = "")
    }
  }
}

# 6. 合并三个数据集的免疫评分做总览 --------------------------------------
all_immune_df <- bind_rows(all_immune_results)
write.table(all_immune_df, "immune_scores_all_datasets.txt", sep = "\t", quote = F, row.names = F)

# 合并三个数据集的相关性（取平均）
if (length(all_cor_results) > 0) {
  cor_mats <- lapply(all_cor_results, function(x) x$cor)
  common_genes <- Reduce(intersect, lapply(cor_mats, rownames))
  common_cells <- Reduce(intersect, lapply(cor_mats, colnames))
  
  if (length(common_genes) > 0 && length(common_cells) > 0) {
    avg_cor <- matrix(0, nrow = length(common_genes), ncol = length(common_cells))
    rownames(avg_cor) <- common_genes
    colnames(avg_cor) <- common_cells
    
    for (ds in names(cor_mats)) {
      avg_cor <- avg_cor + cor_mats[[ds]][common_genes, common_cells]
    }
    avg_cor <- avg_cor / length(cor_mats)
    
    write.table(avg_cor, "immune_correlation_average.txt", sep = "\t", quote = F)
    
    pdf("immune_correlation_heatmap_average.pdf", width = 8, height = 4)
    pheatmap(avg_cor,
             scale = "none",
             color = colorRampPalette(c("#0073C2", "white", "#CD534C"))(99),
             breaks = c(seq(-1, -0.01, length = 50), seq(0.01, 1, length = 50)),
             cluster_cols = T,
             cluster_rows = T,
             display_numbers = TRUE,
             number_color = "black",
             fontsize_number = 7,
             main = "Average correlation across 3 datasets")
    dev.off()
    cat("Saved: immune_correlation_heatmap_average.pdf\n")
  }
}

cat("\nDone!\n")
cat("Output files:\n")
cat("- immune_scores_*.txt\n")
cat("- immune_wilcox_*.txt\n")
cat("- immune_boxplot_*.pdf\n")
cat("- immune_correlation_*.txt\n")
cat("- immune_correlation_heatmap_*.pdf\n")
cat("- immune_scores_all_datasets.txt\n")
cat("- immune_correlation_average.txt\n")
cat("- immune_correlation_heatmap_average.pdf\n")
