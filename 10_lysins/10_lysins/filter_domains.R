#!/usr/bin/env Rscript
# =============================================================================
# Filter InterProScan Domains (Refined Version)
# =============================================================================
# Description:
#   - Filter lysin-related functional domains from InterProScan TSV
#   - Remove redundant overlapping hits
#   - Keep best hit (lowest E-value) per overlap cluster
#   - Output full filtered set + functional-only subset
#
# Input:
#   args[1] TSV file (InterProScan standard output)
#   args[2] Functional IPR list ((one InterPro(IPR) accession per line)
#   args[3] Output .xlsx file path
#
# Output:
#   - filtered domains (Excel)
#   - functional IPR only (Excel)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(openxlsx)
})

# -------------------------
# Arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript filter_domains_refined.R <input.tsv> <ipr.txt> <output.xlsx>")
}

in_tsv <- args[1]
ipr_file <- args[2]
out_file <- args[3]

# -------------------------
# Load data
# -------------------------
cat("Reading input files...\n")
dat <- read_tsv(in_tsv, col_names = FALSE, show_col_types = FALSE)
ipr_list <- read_lines(ipr_file)

# -------------------------
# Basic validation
# -------------------------
if (ncol(dat) < 12) {
  stop("TSV format error: expected at least 12 columns from InterProScan output")
}

# Standard InterProScan column mapping (position-based)
dat <- dat %>%
  mutate(
    Protein = .data[[1]],
    Start   = as.numeric(.data[[7]]),
    End     = as.numeric(.data[[8]]),
    Evalue  = as.numeric(.data[[9]]),
    IPR     = .data[[12]]
  ) %>%
  filter(!is.na(Evalue))

# -------------------------
# Step 1: functional filtering
# -------------------------
step1 <- dat %>%
  group_by(Protein) %>%
  filter(any(IPR %in% ipr_list & Evalue < 1e-5, na.rm = TRUE)) %>%
  ungroup()

# -------------------------
# Overlap detection
# -------------------------
overlap <- function(s1, e1, s2, e2) {
  !(e1 < s2 | e2 < s1)
}

# -------------------------
# Reduce redundant overlaps
# -------------------------
reduce_group <- function(df) {
  if (nrow(df) <= 1) return(df)
  visited <- rep(FALSE, nrow(df))
  out <- list()
  
  for (i in seq_len(nrow(df))) {
    if (visited[i]) next
    
    cluster <- which(
      overlap(df$Start[i], df$End[i], df$Start, df$End)
    )
    
    visited[cluster] <- TRUE
    sub <- df[cluster, ]
    
    best <- sub %>% slice_min(Evalue, n = 1, with_ties = FALSE)
    out[[length(out) + 1]] <- best
  }
  
  bind_rows(out)
}

# -------------------------
# Apply reduction
# -------------------------
result <- step1 %>%
  group_by(Protein) %>%
  group_modify(~reduce_group(.x)) %>%
  ungroup()

# -------------------------
# Outputs
# -------------------------
write.xlsx(result, out_file)

functional <- result %>%
  filter(IPR %in% ipr_list & Evalue < 1e-5)

func_file <- sub("\\.xlsx$", "_FUNCTIONAL_IPR.xlsx", out_file)
write.xlsx(functional, func_file)

cat("Done:\n")
cat(" - filtered:", out_file, "\n")
cat(" - functional:", func_file, "\n")
