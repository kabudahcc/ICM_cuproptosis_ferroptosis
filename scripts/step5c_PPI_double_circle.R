# PPI 双圈图（Double-circle plot）
# 输入：PPI_edges.txt, PPI_nodes.txt
# 输出：PPI_double_circle.pdf

rm(list = ls())
library(ggplot2)
library(dplyr)

# 1. 读取 PPI 结果 --------------------------------------------------------
edges <- read.table("PPI_edges.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nodes <- read.table("PPI_nodes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Nodes:", nrow(nodes), "Edges:", nrow(edges), "\n")

# 2. 确定 Hub 基因（Top 10 by degree）--------------------------------------
top_n_hub <- 10
nodes <- nodes %>% arrange(desc(degree))
hub_genes <- head(nodes$symbol, top_n_hub)
nodes$is_hub <- ifelse(nodes$symbol %in% hub_genes, "Hub", "Other")

# 3. 计算节点在双圈上的坐标 ----------------------------------------------
# 内圈：Hub 基因；外圈：其他基因
# 按 category 分组排序，使同类别基因尽量聚集
nodes <- nodes %>%
  mutate(category = factor(category, levels = c("Cuproptosis", "Ferroptosis", "Both")),
         core = factor(core, levels = c("Core", "Extended"))) %>%
  arrange(category, desc(is_hub), desc(degree))

# 内圈 Hub 角度
hub_df <- nodes %>% filter(is_hub == "Hub") %>%
  mutate(angle = seq(0, 2 * pi * (1 - 1 / n()), length.out = n()),
         radius = 1,
         x = radius * cos(angle),
         y = radius * sin(angle))

# 外圈其他基因角度
other_df <- nodes %>% filter(is_hub == "Other") %>%
  mutate(angle = seq(0, 2 * pi * (1 - 1 / n()), length.out = n()),
         radius = 2,
         x = radius * cos(angle),
         y = radius * sin(angle))

layout_df <- bind_rows(hub_df, other_df) %>%
  select(symbol, x, y, angle, radius, category, core, is_hub, degree)

# 4. 计算边的坐标 --------------------------------------------------------
edge_coords <- edges %>%
  inner_join(layout_df %>% select(symbol, x, y), by = c("from_symbol" = "symbol")) %>%
  rename(x_from = x, y_from = y) %>%
  inner_join(layout_df %>% select(symbol, x, y), by = c("to_symbol" = "symbol")) %>%
  rename(x_to = x, y_to = y)

# 5. 绘制双圈图 ----------------------------------------------------------
p <- ggplot() +
  # 边
  geom_segment(data = edge_coords,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                   alpha = combined_score / 1000, linewidth = combined_score / 400),
               color = "grey70", show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.15, 0.5)) +
  scale_linewidth_continuous(range = c(0.2, 1)) +
  # 节点
  geom_point(data = layout_df,
             aes(x = x, y = y, color = category, shape = core, size = degree),
             stroke = 1) +
  scale_color_manual(values = c("Cuproptosis" = "#0073C2",
                                 "Ferroptosis" = "#CD534C",
                                 "Both" = "#7AC36A")) +
  scale_shape_manual(values = c("Core" = 17, "Extended" = 16)) +
  scale_size_continuous(range = c(3, 9)) +
  # 标签：Hub 基因显示，其他只显示 degree >= 9 的
  geom_text(data = layout_df %>% filter(is_hub == "Hub" | degree >= 9),
            aes(x = x, y = y, label = symbol),
            size = 3, vjust = -1, family = "sans") +
  # 添加内外圈参考圆环
  annotate("path",
           x = cos(seq(0, 2 * pi, length.out = 100)),
           y = sin(seq(0, 2 * pi, length.out = 100)),
           color = "grey80", linetype = "dashed", linewidth = 0.3) +
  annotate("path",
           x = 2 * cos(seq(0, 2 * pi, length.out = 100)),
           y = 2 * sin(seq(0, 2 * pi, length.out = 100)),
           color = "grey80", linetype = "dashed", linewidth = 0.3) +
  coord_fixed() +
  theme_void(base_family = "sans") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, family = "sans")) +
  labs(title = "PPI double-circle network of cuproptosis/ferroptosis-related DEGs",
       color = "Category", shape = "Core gene", size = "Degree")

ggsave(p, file = "PPI_double_circle.pdf", width = 10, height = 8)
cat("Saved: PPI_double_circle.pdf\n")

# 6. 保存双圈图布局 --------------------------------------------------------
write.table(layout_df, "PPI_double_circle_layout.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: PPI_double_circle_layout.txt\n")

cat("\nDone!\n")
