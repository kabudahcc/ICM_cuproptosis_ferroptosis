# step13c_inspect_GSE141910.R
# 检查 GSE141910（MAGNet RNA-seq）的 etiology 分组，确认是否含有明确 ICM 标签

series_file <- "GSE141910_series_matrix.txt.gz"
lines <- readLines(gzfile(series_file))

get_meta <- function(prefix) {
  idx <- grep(paste0("^!", prefix, "\t"), lines)
  if (length(idx) == 0) return(NULL)
  vals <- strsplit(lines[idx], "\t")[[1]][-1]
  gsub('"', '', vals)
}

geo_ids <- get_meta("Sample_geo_accession")
titles  <- get_meta("Sample_title")

# 提取 characteristics 中 etiology 行
eti_idx <- grep('^!Sample_characteristics_ch1\t"etiology:', lines)
etiology <- gsub('"', '', strsplit(lines[eti_idx], "\t")[[1]][-1])
etiology <- gsub("etiology: ", "", etiology)

meta <- data.frame(
  sample = geo_ids,
  title  = titles,
  etiology = etiology,
  stringsAsFactors = FALSE
)

cat("GSE141910 样本分组统计：\n")
print(table(meta$etiology))

# 是否含 ICM？
has_icm <- any(grepl("ischemic|ICM|IHD|ischemic cardiomyopathy", meta$etiology, ignore.case = TRUE))
cat("含明确 ICM 标签:", has_icm, "\n")

if (!has_icm) {
  cat("结论：GSE141910 仅包含 Non-Failing / PPCM / HCM / DCM，无可用于 ICM 验证的明确分组，跳过作为 ICM 验证集。\n")
}

write.csv(meta, "GSE141910_etiology_check.csv", row.names = FALSE)
cat("Saved: GSE141910_etiology_check.csv\n")
