library(data.table)
files <- c(driver="ferrdb_driver.txt", suppressor="ferrdb_suppressor.txt", marker="ferrdb_marker.txt")
out_list <- list()
for(cat in names(files)){
  df <- fread(files[cat], sep=",", quote='"', header=TRUE)
  # standardize column names
  setnames(df, make.names(names(df)))
  df <- df[Test_organism == "Human" & Confidence == "Validated", .(Symbol, Category)]
  df[, category := cat]
  out_list[[cat]] <- df
}
all <- rbindlist(out_list)
all[, Category := trimws(Category)]
unique_genes <- unique(all[, .(symbol = Symbol, category, confidence="Validated", organism="Human")])
unique_genes <- unique_genes[order(symbol)]
fwrite(unique_genes, "ferroptosis_genes.txt", sep="\t", quote=FALSE)
cat("Ferroptosis validated human genes:", nrow(unique_genes), "\n")
print(head(unique_genes, 20))
