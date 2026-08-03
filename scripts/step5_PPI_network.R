# PPI 网络分析
# 输入：CRG_genes_GSE*.txt / FRG_genes_GSE*.txt, cell_death_intersect_three.txt
# 输出：PPI_edges.txt, PPI_nodes.txt, hub_genes.txt, PPI_network.pdf

rm(list = ls())
library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

# 1. 读取基因 ------------------------------------------------------------
datasets <- c("GSE16499", "GSE5406", "GSE57338")

read_gene_list <- function(file) {
  if (!file.exists(file)) return(character(0))
  g <- read.table(file, header = F, sep = "\t", stringsAsFactors = F)[, 1]
  return(unique(trimws(as.character(g))))
}

crg_union <- unique(unlist(lapply(paste0("CRG_genes_", datasets, ".txt"), read_gene_list)))
frg_union <- unique(unlist(lapply(paste0("FRG_genes_", datasets, ".txt"), read_gene_list)))

# 合并为 PPI 输入基因集
ppi_genes <- unique(c(crg_union, frg_union))
core_genes <- read.table("cell_death_intersect_three.txt", header = F, sep = "\t", stringsAsFactors = F)[, 1]

cat("PPI input genes:", length(ppi_genes), "\n")
cat("Core genes:", length(core_genes), "\n")

# 2. STRINGdb 映射与获取互作 ---------------------------------------------
# 设置 input_directory 为当前目录，缓存 STRING 数据
string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 400,
                          input_directory = getwd())

genes_df <- data.frame(symbol = ppi_genes, stringsAsFactors = F)
mapped <- string_db$map(genes_df, "symbol", removeUnmappedRows = TRUE)
cat("Mapped to STRING:", nrow(mapped), "genes\n")

if (nrow(mapped) < 2) {
  stop("Less than 2 genes mapped to STRING. Cannot build PPI network.")
}

interactions <- string_db$get_interactions(mapped$STRING_id)
cat("Raw interactions (score >= 400):", nrow(interactions), "\n")

# 添加基因 symbol
interactions <- interactions %>%
  left_join(mapped, by = c("from" = "STRING_id")) %>%
  rename(from_symbol = symbol) %>%
  left_join(mapped, by = c("to" = "STRING_id")) %>%
  rename(to_symbol = symbol) %>%
  select(from_symbol, to_symbol, combined_score) %>%
  filter(!is.na(from_symbol), !is.na(to_symbol))

cat("Interactions within input gene set:", nrow(interactions), "\n")

# 3. 构建 igraph 网络 ----------------------------------------------------
if (nrow(interactions) > 0) {
  g <- graph_from_data_frame(interactions, directed = FALSE, vertices = mapped$symbol)
  
  # 计算网络拓扑指标
  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g, normalized = TRUE)
  V(g)$closeness <- closeness(g, normalized = TRUE)
  V(g)$core <- ifelse(V(g)$name %in% core_genes, "Core", "Extended")
  V(g)$category <- ifelse(V(g)$name %in% crg_union,
                          ifelse(V(g)$name %in% frg_union, "Both", "Cuproptosis"),
                          "Ferroptosis")
  
  # 提取节点属性
  node_df <- data.frame(
    symbol = V(g)$name,
    category = V(g)$category,
    core = V(g)$core,
    degree = V(g)$degree,
    betweenness = V(g)$betweenness,
    closeness = V(g)$closeness,
    stringsAsFactors = F
  ) %>% arrange(desc(degree))
  
  write.table(interactions, "PPI_edges.txt", sep = "\t", quote = F, row.names = F)
  write.table(node_df, "PPI_nodes.txt", sep = "\t", quote = F, row.names = F)
  
  # 4. Hub 基因 ------------------------------------------------------------
  top_hub <- head(node_df, 10)
  write.table(top_hub, "hub_genes.txt", sep = "\t", quote = F, row.names = F)
  cat("\nTop 10 hub genes:\n")
  print(top_hub)
  
  # 5. 可视化网络 ----------------------------------------------------------
  # 节点颜色：类别； 节点大小：degree； 边粗细：combined_score
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = combined_score / 200, alpha = combined_score / 1000),
                   color = "grey60", show.legend = FALSE) +
    geom_node_point(aes(color = category, size = degree, shape = core)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 50,
                   family = "sans") +
    scale_color_manual(values = c("Cuproptosis" = "#0073C2",
                                   "Ferroptosis" = "#CD534C",
                                   "Both" = "#7AC36A")) +
    scale_shape_manual(values = c("Core" = 17, "Extended" = 16)) +
    scale_size_continuous(range = c(3, 10)) +
    scale_edge_width(range = c(0.5, 2)) +
    theme_void(base_family = "sans") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, family = "sans")) +
    labs(title = "PPI network of cuproptosis/ferroptosis-related DEGs",
         color = "Category", shape = "Core gene", size = "Degree")
  
  ggsave(p, file = "PPI_network.pdf", width = 10, height = 8)
  cat("\nSaved: PPI_network.pdf\n")
} else {
  cat("No PPI interactions found among input genes.\n")
}

cat("\nDone!\n")
cat("Output files:\n")
cat("- PPI_edges.txt\n")
cat("- PPI_nodes.txt\n")
cat("- hub_genes.txt\n")
cat("- PPI_network.pdf\n")
