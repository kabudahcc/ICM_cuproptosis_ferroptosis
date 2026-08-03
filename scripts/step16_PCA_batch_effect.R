# step16: PCA of the three discovery datasets illustrating batch effects (Supplementary Figure S6)
# Panel A: before batch correction, colored by dataset
# Panel B: after ComBat correction, colored by dataset
# Panel C: after ComBat correction, colored by disease group

rm(list = ls())

# ---- paths (workspace = Paper1_铜铁死亡_ICM) ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
proj_dir <- dirname(script_dir)  # Scripts/ -> project root
data_dir <- file.path(proj_dir, "Data")
fig_dir  <- file.path(proj_dir, "Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(sva)
  library(gridExtra)
})

datasets <- c("GSE16499", "GSE5406", "GSE57338")
mats <- list()
groups <- list()

for (g in datasets) {
  e <- new.env()
  load(file.path(data_dir, paste0(g, "_step2_GEO_data.Rdata")), envir = e)
  mats[[g]] <- as.matrix(e$GSE_ann)           # genes x samples
  groups[[g]] <- as.character(e$Group)        # Normal / Disease
  cat(g, dim(mats[[g]]), "\n")
}

# ---- common genes & merge ----
common <- Reduce(intersect, lapply(mats, rownames))
cat("Common genes:", length(common), "\n")
merged <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
batch <- rep(datasets, times = sapply(mats, ncol))
group <- unlist(groups, use.names = FALSE)
stopifnot(ncol(merged) == length(batch))

# ---- ComBat correction (protect disease group) ----
mod <- model.matrix(~ as.factor(group))
merged_cb <- ComBat(dat = merged, batch = batch, mod = mod, par.prior = TRUE)

# ---- PCA helper ----
do_pca <- function(mat) {
  prcomp(t(mat), scale. = TRUE, center = TRUE)
}
pca_df <- function(pca, color_by, label) {
  ve <- round(100 * summary(pca)$importance[2, 1:2], 1)
  data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    color = color_by,
    xlab = paste0("PC1 (", ve[1], "%)"),
    ylab = paste0("PC2 (", ve[2], "%)")
  )
}
plot_pca <- function(df, title, legend_title, pal) {
  ggplot(df, aes(PC1, PC2, color = color)) +
    geom_point(size = 1.8, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5, show.legend = FALSE) +
    scale_color_manual(values = pal) +
    labs(title = title, x = df$xlab[1], y = df$ylab[1], color = legend_title) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 11),
          legend.position = "right")
}

batch_pal <- c("GSE16499" = "#00AFBB", "GSE5406" = "#E7B800", "GSE57338" = "#FC4E07")
group_pal <- c("Normal" = "#00AFBB", "Disease" = "#E7B800")

pA <- plot_pca(pca_df(do_pca(merged), batch),    "Before batch correction", "Dataset", batch_pal)
pB <- plot_pca(pca_df(do_pca(merged_cb), batch), "After ComBat correction", "Dataset", batch_pal)
pC <- plot_pca(pca_df(do_pca(merged_cb), group), "After ComBat correction", "Group",   group_pal)

combined <- arrangeGrob(pA, pB, pC, ncol = 2, layout_matrix = rbind(c(1, 2), c(3, 3)))

ggsave(file.path(fig_dir, "PCA_batch_effect_S6.pdf"), combined, width = 10, height = 8)
ggsave(file.path(fig_dir, "PCA_batch_effect_S6.png"), combined, width = 10, height = 8, dpi = 300)
cat("Saved: Figures/PCA_batch_effect_S6.pdf / .png\n")
