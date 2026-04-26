# =============================================================================
# Viral Identification Results — Merging and Filtering
# =============================================================================
# Description: Merges viral sequence predictions from VIBRANT, VirSorter2,
#              and geNomad, with CheckV quality filtering applied after an
#              initial VIBRANT × VirSorter2 intersection step.
#
# Input Excel file (cold_seep.xlsx) sheet layout:
#   Sheet 1 : VIBRANT results        
#   Sheet 2 : VirSorter2 (score≥0.5)  
#   Sheet 3 : geNomad filtered results 
#   Sheet 4 : VirSorter2 (score≥0.8)  
#   Sheet 5 : CheckV-filtered VIBRANT × VirSorter2 
#
# Output:
#   cold_seep_VIBRANT_vs2_0.5.txt   : intersection of VIBRANT and VirSorter2
#                                        (score≥0.5); used as CheckV input
#   cold_seep_ALL_combined_id.txt   : union of all three tools after
#                                        CheckV filtering
# =============================================================================

library(readxl)

INPUT_XLSX <- "cold_seep.xlsx"

# -----------------------------------------------------------------------------
# Step 1: Intersect VIBRANT and VirSorter2 (score >= 0.5)
# -----------------------------------------------------------------------------
# This conservative intersection is passed to CheckV for quality assessment.

vibrant    <- read_xlsx(INPUT_XLSX, sheet = 1)  # VIBRANT predictions
vs2_0.5    <- read_xlsx(INPUT_XLSX, sheet = 2)  # VirSorter2, max_score >= 0.5

vibrant_vs2_0.5 <- intersect(vibrant$scaffold, vs2_0.5$seqname)
vibrant_vs2_0.5 <- as.data.frame(vibrant_vs2_0.5)
colnames(vibrant_vs2_0.5) <- "seq_id"

write.table(
  vibrant_vs2_0.5,
  file      = "cold_seep_VIBRANT_vs2_0.5.txt",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

# -----------------------------------------------------------------------------
# Step 2: Load CheckV-filtered results
# -----------------------------------------------------------------------------
# After running CheckV on the sequences above, low-quality sequences are
# removed. The filtered set is loaded here for the final merge.
# (CheckV filtering is performed externally; see 01_viral_identification.sh)

vibrant_vs2_0.5_checkv <- read_xlsx(INPUT_XLSX, sheet = 5)  # CheckV-filtered


# -----------------------------------------------------------------------------
# Step 3: Union of all three tools
# -----------------------------------------------------------------------------
# Combine geNomad, VirSorter2 (score>=0.8), and the CheckV-filtered

genomad  <- read_xlsx(INPUT_XLSX, sheet = 3)  # geNomad filtered results
vs2_0.8  <- read_xlsx(INPUT_XLSX, sheet = 4)  # VirSorter2, max_score >= 0.8

all_combined <- union(genomad$seq_name, vs2_0.8$seqname)           
all_combined <- union(all_combined, vibrant_vs2_0.5_checkv$contig_id)  

all_combined <- as.data.frame(all_combined)
colnames(all_combined) <- "seq_id"

write.table(
  all_combined,
  file      = "cold_seep_ALL_combined_id.txt",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
