# WGCNA Top 10 Hub 基因（基于已保存的 WGCNA 结果）
rm(list = ls())

library(WGCNA)
library(dplyr)

load("WGCNA_results_part1.RData")

# 计算基因模块成员度（kME）
MM_full <- signedKME(datExpr, MEs)

cat("kME columns:", paste(colnames(MM_full), collapse = ", "), "\n")

# 将 ME 名称（ME0, ME1...）映射到颜色
ME_to_color <- setNames(labels2colors(as.numeric(gsub("ME", "", colnames(MEs)))),
                        colnames(MEs))

top_hub_genes <- data.frame()
for (me in colnames(MEs)) {
  mod_num <- as.numeric(gsub("ME", "", me))
  mod_color <- ME_to_color[me]
  kme_col <- paste0("kME", mod_num)
  if (!kme_col %in% colnames(MM_full)) {
    cat("  Skip", me, "- kME column not found\n")
    next
  }
  idx <- order(MM_full[, kme_col], decreasing = TRUE)[1:10]
  top_hub_genes <- rbind(
    top_hub_genes,
    data.frame(module = me,
               module_color = mod_color,
               gene = rownames(MM_full)[idx],
               kME = round(unname(MM_full[idx, kme_col]), 4),
               stringsAsFactors = FALSE)
  )
}

write.table(top_hub_genes, "WGCNA_top10_hub_genes_per_module.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: WGCNA_top10_hub_genes_per_module.txt\n")
print(top_hub_genes)
