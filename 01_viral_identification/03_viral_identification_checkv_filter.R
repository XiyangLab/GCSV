# =============================================================================
# Viral Identification — Final CheckV-based Quality Filtering
# =============================================================================
# Description: After all predicted viral sequences (from VIBRANT, VirSorter2,
#              and geNomad) are pooled and run through CheckV, this script
#              cross-references the CheckV quality report with the per-tool
#              scores to discard low-confidence sequences.
#
# Inputs:
#   cold_seep.xlsx (sheet 3) : geNomad filtered results  
#   cold_seep.xlsx (sheet 4) : VirSorter2 (score≥0.8)  
#   checkv_cold_seep_all/quality_summary.tsv : CheckV output for all sequences
#
# Output:
#   cold_seep_checkv_discarded_id.txt : sequence IDs to exclude from the
#                                          final viral dataset
# =============================================================================

library(readxl)

INPUT_XLSX <- "cold_seep.xlsx"

# -----------------------------------------------------------------------------
# Step 1: Load tool results and CheckV quality summary
# -----------------------------------------------------------------------------

genomad <- read_xlsx(INPUT_XLSX, sheet = 3)  # geNomad filtered predictions
vs2_0.8 <- read_xlsx(INPUT_XLSX, sheet = 4)  # VirSorter2, max_score >= 0.8

checkv <- read.table(
  "checkv_cold_seep_all/quality_summary.tsv",
  sep    = "\t",
  header = TRUE
)


# -----------------------------------------------------------------------------
# Step 2: Isolate CheckV-ambiguous sequences
# -----------------------------------------------------------------------------
# Flag sequences where CheckV found zero viral genes but at least one host gene

ambiguous_checkv <- subset(checkv, viral_genes == 0 & host_genes != 0)

cat(sprintf("Sequences with 0 viral genes and >=1 host gene: %d\n",
            nrow(ambiguous_checkv)))


# -----------------------------------------------------------------------------
# Step 3: Merge with per-tool scores
# -----------------------------------------------------------------------------

merged_df <- merge(
  ambiguous_checkv, genomad,
  by.x  = names(ambiguous_checkv)[1],
  by.y  = names(genomad)[1],
  all.x = TRUE
)

merged_df <- merge(
  merged_df, vs2_0.8,
  by.x  = names(merged_df)[1],
  by.y  = names(vs2_0.8)[1],
  all.x = TRUE
)

merged_df <- merged_df[, c(1, 3, 6, 7, 8, 10, 20, 22, 28, 31)]


# -----------------------------------------------------------------------------
# Step 4: Define discard criteria
# -----------------------------------------------------------------------------
# A sequence is discarded if it is ambiguous in CheckV AND shows weak evidence in whichever tool(s) detected it.

discarded_1 <- merged_df[
  !is.na(merged_df$max_score)   & merged_df$max_score  < 0.95 &
  !is.na(merged_df$hallmark)    & merged_df$hallmark   <= 2   &
   is.na(merged_df$virus_score) &
   is.na(merged_df$n_hallmarks),
]

discarded_2 <- merged_df[
  !is.na(merged_df$virus_score) & merged_df$virus_score  < 0.95 &
  !is.na(merged_df$n_hallmarks) & merged_df$n_hallmarks  <= 2   &
   is.na(merged_df$max_score)   &
   is.na(merged_df$hallmark),
]

discarded_3 <- merged_df[
  !is.na(merged_df$virus_score) & merged_df$virus_score <= 0.95 &
  !is.na(merged_df$n_hallmarks) & merged_df$n_hallmarks <= 2    &
  !is.na(merged_df$max_score)   & merged_df$max_score   <= 0.95 &
  !is.na(merged_df$hallmark)    & merged_df$hallmark    <= 2,
]

discarded <- rbind(discarded_1, discarded_2, discarded_3)

# -----------------------------------------------------------------------------
# Step 5: Check for duplicate IDs and export
# -----------------------------------------------------------------------------

duplicated_ids <- discarded[duplicated(discarded[, 1]), ]

if (nrow(duplicated_ids) > 0) {
  cat(sprintf("Warning: %d duplicate sequence IDs found in discard list:\n",
              nrow(duplicated_ids)))
  print(duplicated_ids[, 1])
} else {
  cat("No duplicate sequence IDs in discard list.\n")
}

write.table(
  discarded[, 1],
  file      = "cold_seep_checkv_discarded_id.txt",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
