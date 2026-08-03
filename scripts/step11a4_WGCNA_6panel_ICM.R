# WGCNA 6 panel 组合图（用于论文 Figure）
rm(list = ls())

library(ggplot2)
library(cowplot)
library(patchwork)
library(pdftools)
library(png)

panel_files <- c(
  "WGCNA_A_sample_dendrogram_traits.pdf",
  "WGCNA_B_GO_enrichment.pdf",
  "WGCNA_C_module_eigengene_boxplot.pdf",
  "WGCNA_D_module_heatmap.pdf",
  "WGCNA_E_module_trait_correlation.pdf",
  "WGCNA_F_module_network.pdf"
)

tmp_pngs <- sapply(panel_files, function(f) {
  tmp <- tempfile(fileext = ".png")
  bitmap <- pdf_render_page(f, page = 1, dpi = 200)
  writePNG(bitmap, tmp)
  tmp
})

plot_list <- lapply(tmp_pngs, function(png_file) {
  ggdraw() + draw_image(png_file, scale = 0.95)
})

p_all <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("WGCNA_6panel_combined.pdf", p_all, width = 16, height = 20)
ggsave("WGCNA_6panel_combined.png", p_all, width = 16, height = 20, dpi = 300)

cat("Saved: WGCNA_6panel_combined.pdf\n")
cat("Saved: WGCNA_6panel_combined.png\n")
