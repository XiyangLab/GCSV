# ============================================================================
# File: plot_mcp_terl_properties_template.R
# Status: reconstructed template (not the exact original figure script)
# Purpose: plot physicochemical properties of MCP / TerL proteins.
# Inputs:
#   1) calc_protein_properties_6indices.sh output table
#   2) optional metadata table with columns: Sequence_ID, group, protein_type
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_mcp_terl_properties_template.R <properties.tsv> <out_prefix> [metadata.tsv]")
}

prop_file <- args[1]
out_prefix <- args[2]
meta_file <- ifelse(length(args) >= 3, args[3], NA)

x <- read_tsv(prop_file, show_col_types = FALSE)
if (!is.na(meta_file)) {
  meta <- read_tsv(meta_file, show_col_types = FALSE)
  x <- left_join(x, meta, by = "Sequence_ID")
} else {
  x$group <- "All"
}

long <- x %>%
  pivot_longer(cols = c(Hydropathicity, Relative_mutability, Average_flexibility,
                        Refractivity, Polarity, Transmembrane_tendency, Length),
               names_to = "metric", values_to = "value")

p <- ggplot(long, aes(group, value, fill = group)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(paste0(out_prefix, "_violin.pdf"), p, width = 12, height = 8)
