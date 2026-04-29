#!/usr/bin/env Rscript
# =============================================================================
# Domain Heatmap Plotting
# =============================================================================
# Description:
#   Generate publication-ready heatmaps from taxonomy-domain matrices.
#   Row labels include the number of associated sequences.
#
# Input:
#   - *_heatmap_matrix.xlsx
#
# Output:
#   - *_heatmap.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(pheatmap)
  library(grid)
  library(stringi)
})

# =============================================================================
# Input and output directories
# =============================================================================
input_dir  <- "input"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

heatmap_files <- c(
  "host_phylum_heatmap_matrix.xlsx",
  "host_family_top20_heatmap_matrix.xlsx",
  "virus_phylum_heatmap_matrix.xlsx",
  "virus_family_heatmap_matrix.xlsx"
)

# =============================================================================
# Heatmap plotting function
# =============================================================================
plot_heatmap <- function(input_file, output_file, title) {

  data <- read_excel(input_file)
  colnames(data) <- stri_replace_all_fixed(colnames(data), "β", "beta")

  required_cols <- c("taxonomy", "seqs_count")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input file must contain: taxonomy and seqs_count")
  }

  row_labels <- paste0(data$taxonomy, " (", data$seqs_count, ")")

  mat <- data %>%
    select(-taxonomy, -seqs_count) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  rownames(mat) <- row_labels
  mat[is.na(mat)] <- 0

  log_mat <- log10(mat + 1)
  max_value <- max(mat)

  legend_breaks <- if (max_value <= 10) {
    seq(0, max_value, length.out = 5)
  } else if (max_value <= 100) {
    seq(0, max_value, length.out = 6)
  } else {
    unique(c(10^seq(0, floor(log10(max_value))), max_value))
  }

  p <- pheatmap(
    log_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("white", "#0D47A1"))(100),
    main = title,
    fontsize = 10,
    fontsize_row = 6,
    fontsize_col = 7,
    angle_col = 45,
    legend_breaks = log10(legend_breaks + 1),
    legend_labels = round(legend_breaks, 2),
    silent = TRUE
  )

  width  <- 6 + ncol(mat) * 0.25
  height <- 4 + nrow(mat) * 0.20

  pdf(output_file, width = width, height = height)
  grid.draw(p$gtable)
  dev.off()

  message("Saved: ", output_file)
}

# =============================================================================
# Generate heatmaps
# =============================================================================
for (file in heatmap_files) {

  input_file  <- file.path(input_dir, file)
  output_file <- file.path(output_dir, sub("\\.xlsx$", ".pdf", file))
  plot_title  <- gsub("_", " ", sub("_heatmap_matrix\\.xlsx$", "", file))

  plot_heatmap(input_file, output_file, plot_title)
}

message("All heatmaps generated successfully.")
