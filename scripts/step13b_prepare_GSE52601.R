# step13b_prepare_GSE52601.R
# 解析 GSE52601 series matrix（Illumina HumanHT-12 V4），提取 postnatal non-failing
# 与 pre-LVAD ischemic cardiomyopathy (ICM) 样本，注释到 Gene Symbol，输出统一格式。

.libPaths(c('Rlib', .libPaths()))

library(data.table)
library(illuminaHumanv4.db)

# ---------------------------- 参数 ----------------------------
series_file <- "GSE52601_series_matrix.txt.gz"
group_file  <- "GSE52601_group.csv"
symbol_file <- "GSE52601_symbol.csv"

# ---------------------------- 读取 series matrix ------------------------
lines <- readLines(gzfile(series_file))

# 提取 sample metadata
get_meta <- function(prefix) {
  idx <- grep(paste0("^!", prefix, "\t"), lines)
  if (length(idx) == 0) return(NULL)
  vals <- strsplit(lines[idx], "\t")[[1]][-1]
  gsub('"', '', vals)
}

geo_ids   <- get_meta("Sample_geo_accession")
samples   <- get_meta("Sample_title")

if (length(geo_ids) != length(samples)) {
  stop("Sample_geo_accession 与 Sample_title 长度不一致")
}

meta <- data.frame(
  geo_accession = geo_ids,
  title  = samples,
  stringsAsFactors = FALSE
)

# 定义分组：只保留 postnatal non-failing 与 postnatal failing ICM（未 VAD 治疗）
meta$group <- NA
meta$group[grep("postnatal-non-failing", meta$title)] <- "Normal"
meta$group[grep("postnatal-failing-ICM-rep", meta$title)] <- "Disease"
# 排除 fetal、DCM、VAD-treated（title 中含 VAD-treatment）
meta <- meta[!is.na(meta$group), ]

cat("GSE52601 可用样本：", nrow(meta), "\n")
print(table(meta$group))

# ---------------------------- 读取表达矩阵 ----------------------------
tbl_start <- grep("!series_matrix_table_begin", lines) + 1
tbl_end   <- grep("!series_matrix_table_end", lines) - 1
tbl_lines <- lines[tbl_start:tbl_end]

tmp <- tempfile()
writeLines(tbl_lines, tmp)
expr <- fread(tmp, check.names = FALSE, data.table = FALSE)
unlink(tmp)

rownames(expr) <- expr[[1]]
expr <- expr[, -1, drop = FALSE]

# 只保留 selected samples
to_keep <- intersect(meta$geo_accession, colnames(expr))
if (length(to_keep) != nrow(meta)) {
  warning("部分样本在表达矩阵中缺失")
  meta <- meta[meta$geo_accession %in% to_keep, ]
}
expr <- expr[, meta$geo_accession, drop = FALSE]

# ---------------------------- 注释到 Symbol ----------------------------
probe_ids <- rownames(expr)
mapped <- select(illuminaHumanv4.db, keys = probe_ids, columns = c("SYMBOL"), keytype = "PROBEID")
# 去重：每个 probe 只保留第一个 Symbol
mapped <- mapped[!duplicated(mapped$PROBEID), ]

expr$PROBEID <- rownames(expr)
expr_annot <- merge(expr, mapped, by = "PROBEID", all.x = TRUE)
expr_annot <- expr_annot[!is.na(expr_annot$SYMBOL) & expr_annot$SYMBOL != "", ]

# 对同一个 Symbol 取中位数
gene_expr <- sapply(split(expr_annot[, meta$geo_accession, drop = FALSE], expr_annot$SYMBOL),
                    function(x) apply(x, 2, median, na.rm = TRUE))
gene_expr <- as.data.frame(t(gene_expr))
gene_expr$symbol <- rownames(gene_expr)
gene_expr <- gene_expr[, c("symbol", meta$geo_accession)]

# ---------------------------- 保存 ----------------------------
write.csv(meta[, c("geo_accession", "group")], group_file, row.names = FALSE)
write.csv(gene_expr, symbol_file, row.names = FALSE)

cat("Saved:", group_file, "\n")
cat("Saved:", symbol_file, "\n")
cat("Genes:", nrow(gene_expr), "\n")
