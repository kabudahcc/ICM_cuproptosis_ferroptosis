# 单细胞分析：GSE145154 人心脏单细胞数据
# 输入：GSE145154_RAW/ 目录下 10x 文件
# 输出：single_cell_UMAP.pdf, single_cell_dotplot.pdf, single_cell_violin.pdf,
#       single_cell_marker_genes.txt, single_cell_metadata.txt

rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. 设置数据目录和样本 ---------------------------------------------------
data_dir <- "GSE145154_RAW"
core_genes <- read.table("cell_death_intersect_three.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]

# 选择左心室（LV）样本：1 个 Normal + 2 个 DCM + 2 个 ICM
selected_samples <- c(
  "GSM4307515_N-1-LVP", "GSM4307516_N-1-LVN",
  "GSM4307520_DCM-2-LVP", "GSM4307521_DCM-2-LVN",
  "GSM4307525_DCM-3-LVP", "GSM4307526_DCM-3-LVN",
  "GSM4307535_ICM-2-LVP", "GSM4307536_ICM-2-LVN",
  "GSM4307540_ICM-3-LVP", "GSM4307541_ICM-3-LVN"
)

sample_meta <- data.frame(
  sample = selected_samples,
  patient = c("N-1", "N-1", "DCM-2", "DCM-2", "DCM-3", "DCM-3", "ICM-2", "ICM-2", "ICM-3", "ICM-3"),
  condition = c("Normal", "Normal", "DCM", "DCM", "DCM", "DCM", "ICM", "ICM", "ICM", "ICM"),
  region = c("LVP", "LVN", "LVP", "LVN", "LVP", "LVN", "LVP", "LVN", "LVP", "LVN"),
  stringsAsFactors = F
)

cat("Selected samples:", length(selected_samples), "\n")

# 2. 读取每个样本的 10x 数据 ----------------------------------------------
library(Matrix)

read_10x_sample <- function(prefix, data_dir) {
  barcodes <- read.table(gzfile(file.path(data_dir, paste0(prefix, "_barcodes.tsv.gz"))), header = F, stringsAsFactors = F)[, 1]
  features <- read.table(gzfile(file.path(data_dir, paste0(prefix, "_features.tsv.gz"))), header = F, stringsAsFactors = F, sep = "\t")
  gene_names <- features[, 2]  # 第二列是 gene symbol
  
  # 读取 matrix.mtx
  mat <- readMM(gzfile(file.path(data_dir, paste0(prefix, "_matrix.mtx.gz"))))
  rownames(mat) <- make.unique(gene_names)
  colnames(mat) <- barcodes
  return(mat)
}

seurat_list <- list()
for (s in selected_samples) {
  cat("Reading sample:", s, "\n")
  counts <- read_10x_sample(s, data_dir)
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = counts, project = s, min.cells = 3, min.features = 200)
  meta <- sample_meta[sample_meta$sample == s, ]
  seurat_obj$patient <- meta$patient
  seurat_obj$condition <- meta$condition
  seurat_obj$region <- meta$region
  seurat_obj$sample <- s
  
  seurat_list[[s]] <- seurat_obj
}

# 3. 合并样本 --------------------------------------------------------------
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = selected_samples, project = "GSE145154")
combined$condition <- factor(combined$condition, levels = c("Normal", "DCM", "ICM"))

cat("Combined cells before QC:", ncol(combined), "\n")

# 4. 质控 ------------------------------------------------------------------
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA > 500 & percent.mt < 15)

cat("Cells after QC:", ncol(combined), "\n")

# Seurat v5: 合并后 counts 按样本分层，需 JoinLayers 再下游分析
combined <- JoinLayers(combined)
cat("Joined layers. Features:", length(rownames(combined)), "\n")

# 保存中间结果，避免失败后从头开始
saveRDS(combined, "single_cell_after_QC.rds")

# 5. 标准化、找高变基因、降维、聚类 ----------------------------------------
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.8)
combined <- RunUMAP(combined, dims = 1:20)

# 6. 细胞类型注释（基于已知标志物）-----------------------------------------
# 心脏主要细胞类型标志物
cell_markers <- list(
  Cardiomyocytes = c("MYH6", "MYH7", "TNNT2", "TTN"),
  Fibroblasts = c("DCN", "COL1A1", "COL1A2", "PDGFRA"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "ENG"),
  Pericytes = c("RGS5", "PDGFRB", "CSPG4"),
  Macrophages = c("CD14", "CD68", "CSF1R", "CD163"),
  T_cells = c("CD3D", "CD3E", "TRAC"),
  B_cells = c("CD79A", "CD79B", "MS4A1"),
  NK_cells = c("NKG7", "GNLY", "KLRD1"),
  Neutrophils = c("S100A8", "S100A9", "CSF3R"),
  Mast_cells = c("TPSAB1", "KIT"),
  Smooth_muscle = c("ACTA2", "MYH11"),
  Adipocytes = c("PLIN1", "ADIPOQ")
)

# 为每个细胞计算各细胞类型的评分
cat("\nMarker presence check:\n")
for (ct in names(cell_markers)) {
  markers <- cell_markers[[ct]]
  markers_present <- intersect(markers, rownames(combined))
  cat("  ", ct, ":", length(markers_present), "/", length(markers), "present\n")
}

score_added <- c()
for (ct in names(cell_markers)) {
  markers <- cell_markers[[ct]]
  markers_present <- intersect(markers, rownames(combined))
  if (length(markers_present) >= 2) {
    combined <- tryCatch(
      AddModuleScore(combined, features = list(markers_present), name = ct, ctrl = 5),
      error = function(e) {
        cat("  AddModuleScore failed for", ct, ":", conditionMessage(e), "\n")
        return(combined)
      }
    )
    score_added <- c(score_added, ct)
  } else {
    cat("  Skipping", ct, ": insufficient markers present\n")
  }
}

# 基于最高评分分配细胞类型
if (length(score_added) > 0) {
  score_cols <- paste0(score_added, "1")
  score_cols <- score_cols[score_cols %in% colnames(combined@meta.data)]
  scores <- combined@meta.data[, score_cols, drop = FALSE]
  colnames(scores) <- gsub("1$", "", colnames(scores))
  combined$cell_type <- apply(scores, 1, function(x) {
    if (all(is.na(x))) return("Unknown")
    return(colnames(scores)[which.max(x)])
  })
} else {
  combined$cell_type <- paste0("Cluster_", combined$seurat_clusters)
  cat("Warning: No module scores added; using cluster IDs as cell type labels.\n")
}

# 7. 可视化 ----------------------------------------------------------------
# UMAP 按细胞类型和 condition 着色
p1 <- DimPlot(combined, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  labs(title = "Cell types") +
  theme(legend.position = "right")

p2 <- DimPlot(combined, reduction = "umap", group.by = "condition", label = FALSE) +
  labs(title = "Condition") +
  theme(legend.position = "right")

p_combined <- p1 + p2

ggsave(p_combined, file = "single_cell_UMAP.pdf", width = 14, height = 6)
cat("Saved: single_cell_UMAP.pdf\n")

# 8. 核心基因表达点图 ------------------------------------------------------
core_present <- intersect(core_genes, rownames(combined))
cat("Core genes present in scRNA-seq:", core_present, "\n")

if (length(core_present) > 0) {
  p3 <- DotPlot(combined, features = core_present, group.by = "cell_type") +
    RotatedAxis() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Core gene expression across cardiac cell types") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(p3, file = "single_cell_dotplot.pdf", width = 10, height = 6)
  cat("Saved: single_cell_dotplot.pdf\n")
  
  # 小提琴图
  p4 <- VlnPlot(combined, features = core_present, group.by = "cell_type", pt.size = 0, ncol = length(core_present)) +
    labs(title = "Core gene expression distribution")
  
  ggsave(p4, file = "single_cell_violin.pdf", width = 14, height = 8)
  cat("Saved: single_cell_violin.pdf\n")
}

# 9. 按 condition 比较核心基因在主要细胞类型中的表达 -------------------------
major_cell_types <- c("Cardiomyocytes", "Fibroblasts", "Endothelial", "Macrophages", "T_cells")
for (gene in core_present) {
  for (ct in major_cell_types) {
    cells <- subset(combined, subset = cell_type == ct)
    if (ncol(cells) < 10) next
    
    expr_df <- data.frame(
      expression = as.vector(GetAssayData(cells, layer = "data")[gene, ]),
      condition = cells$condition,
      cell_type = ct,
      gene = gene,
      stringsAsFactors = F
    )
    
    # 保存到文件
    write.table(expr_df, paste0("single_cell_", gene, "_", ct, "_expression.txt"), sep = "\t", quote = F, row.names = F)
  }
}

# 10. 保存元数据和 markers -------------------------------------------------
write.table(combined@meta.data, "single_cell_metadata.txt", sep = "\t", quote = F, row.names = T)

# 在 FindAllMarkers 前保存完整对象（该步骤可能较慢）
saveRDS(combined, "single_cell_clustered.rds")
cat("Saved: single_cell_clustered.rds\n")

# 找每个 cluster 的 markers
all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top_markers, "single_cell_marker_genes.txt", sep = "\t", quote = F, row.names = F)

# 保存最终 Seurat 对象（供后续交互式查看或修改）
saveRDS(combined, "single_cell_results.rds")
cat("Saved: single_cell_results.rds\n")

cat("\nDone!\n")
cat("Output files:\n")
cat("- single_cell_UMAP.pdf\n")
cat("- single_cell_dotplot.pdf\n")
cat("- single_cell_violin.pdf\n")
cat("- single_cell_metadata.txt\n")
cat("- single_cell_marker_genes.txt\n")
cat("- single_cell_*_*_expression.txt\n")
