# =============================================================================
# Provirus Boundary Mapping and Cleanup
# =============================================================================
# Description: Maps vOTUs classified as proviruses back to their original
#              geNomad and VirSorter2 boundary predictions, reconciles
#              boundaries using CheckV and IRanges, and removes rRNA-
#              containing regions from the final provirus sequences.
# =============================================================================

library(tidyverse)
library(readxl)
library(IRanges)
library(purrr)


# =============================================================================
# Part 1: Load vOTU Metadata and Apply Exclusion Filters
# =============================================================================

vOTU_v2 <- read.csv("vOTU_new_info.csv")

# Derive the source contig ID:
#   - For proviruses: strip the trailing numeric suffix (e.g. "_1") from vOTU_id
#   - For non-proviruses: use vOTU_id directly
vOTU_v2$contig_id <- ifelse(
  vOTU_v2$Provirus == "Yes",
  sub("_[0-9]+$", "", vOTU_v2$vOTU_id),
  vOTU_v2$vOTU_id
)

# Remove low-quality vOTUs (plasmid-like, concatemers, or < 1 kb after trimming)
remove_id_1 <- read_xlsx("remove_vOTUs_id.xlsx", sheet = "plasme")
remove_id_2 <- read_xlsx("remove_vOTUs_id.xlsx", sheet = "concatemer")
remove_id_3 <- read_xlsx("remove_vOTUs_id.xlsx", sheet = "<1kb")

remove_id <- data.frame(
  vOTU_id = unique(c(remove_id_1[[1]], remove_id_2[[1]], remove_id_3[[1]]))
)

vOTU_v2 <- vOTU_v2[!(vOTU_v2$vOTU_id %in% remove_id$vOTU_id), ]

# Split into provirus and non-provirus sets
vOTU_v2_nopro <- filter(vOTU_v2, Provirus == "No")
vOTU_v2_pro   <- filter(vOTU_v2, Provirus == "Yes")


# =============================================================================
# Part 2: Match Proviruses to geNomad and VirSorter2 Boundary Predictions
# =============================================================================
# Each vOTU predicted as a provirus is looked up in the original geNomad and
# VirSorter2 results to retrieve the predicted integration boundaries.

vs2_0.5 <- read_xlsx('cold_seep.xlsx', sheet = 2)
genomad <- read_xlsx('cold_seep.xlsx', sheet = 3)

vs2_provirus <- filter(vs2_0.5, TYPE != "full") %>%
  mutate(source = "vs2")
vs2_provirus <- vs2_provirus[,c(1,11)]
genomad_provirus <- filter(genomad, topology == "Provirus") %>%
  mutate(source = "genomad") 
genomad_provirus <- genomad_provirus[,c(1,4,12)]
colnames(genomad_provirus) <- c("contig_id", "coordinates", "source")
colnames(vs2_provirus) <- c("contig_id", "source")

provirus_map_genomad <- merge(vOTU_v2, genomad_provirus, by = "contig_id")
provirus_map_vs2 <- merge(vOTU_v2, vs2_provirus, by = "contig_id")
vs2_provirus_boundary <- read.delim("cold_seep_vs2_provirus_boundary.tsv")
vs2_provirus_boundary <- vs2_provirus_boundary[,c(1,4,5)]
provirus_map_vs2 <- merge(provirus_map_vs2, vs2_provirus_boundary, 
                            by.x = "contig_id", by.y = "seqname", all.x = T)

write.csv(provirus_map_genomad, "provirus_map_genomad.csv", row.names = FALSE)
write.csv(provirus_map_vs2,     "provirus_map_vs2.csv",     row.names = FALSE)


# =============================================================================
# Part 3: Separate Proviruses 
# =============================================================================

vOTU_v2_pro_only_checkv <- vOTU_v2_pro[
  !vOTU_v2_pro$contig_id %in% provirus_map_genomad$contig_id &
  !vOTU_v2_pro$contig_id %in% provirus_map_vs2$contig_id, ]

vOTU_v2_nopro_new <- vOTU_v2_nopro[
  !vOTU_v2_nopro$contig_id %in% provirus_map_genomad$contig_id &
  !vOTU_v2_nopro$contig_id %in% provirus_map_vs2$contig_id, ]

write.csv(vOTU_v2_pro_only_checkv, "vOTU_v2_pro_only_checkv.csv", row.names = FALSE)
write.csv(vOTU_v2_nopro_new,       "vOTU_v2_nopro_new.csv",       row.names = FALSE)

write.table(vOTU_v2_pro_only_checkv[, 1], "vOTU_v2_pro_only_checkv.id",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(vOTU_v2_nopro_new[, 1],       "vOTU_v2_nopro_new.id",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# =============================================================================
# Part 4: Prepare BED File for Initial Boundary Extraction (seqkit)
# =============================================================================
# Combine geNomad and VirSorter2 boundary coordinates, convert to
# 0-based BED format required by seqkit subseq.

provirus_map_genomad <- read.csv("provirus_map_genomad.csv")
provirus_map_vs2     <- read.csv("provirus_map_vs2.csv")

provirus_map_genomad_boundary <- provirus_map_genomad[, c(1, 2, 12, 11)] %>%
  separate(coordinates,
           into   = c("trim_bp_start", "trim_bp_end"),
           sep    = "-",
           remove = TRUE) %>%
  mutate(trim_bp_start = as.numeric(trim_bp_start),
         trim_bp_end   = as.numeric(trim_bp_end))

provirus_map_vs2_boundary <- provirus_map_vs2[, c(1, 2, 11, 12, 13)]

# Combine and convert to seqkit BED format (0-based start)
provirus_boundary_vs2_genomad <- rbind(provirus_map_genomad_boundary,
                                        provirus_map_vs2_boundary) %>%
  mutate(trim_bp_start_seqkit = trim_bp_start - 1,
         length = trim_bp_end - trim_bp_start + 1) %>%
  distinct()

provirus_boundary_vs2_genomad_seqkit <- provirus_boundary_vs2_genomad[, c(1, 6, 5)] %>%
  distinct()

write.table(provirus_boundary_vs2_genomad_seqkit,
            "provirus_boundary_vs2_genomad_seqkit.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.csv(provirus_boundary_vs2_genomad,
          "provirus_boundary_vs2_genomad.csv", row.names = FALSE)


# =============================================================================
# Part 5: Parse CheckV-trimmed Boundaries
# =============================================================================
# After running checkv contamination on the extracted provirus sequences,
# parse the trimmed coordinates from the CheckV output IDs.
# Input Excel sheets:
#   checkV_virus    : sequences CheckV kept as complete viruses (no further trimming)
#   checkV_provirus : sequences CheckV trimmed as proviruses (border_checkv column)

checkv_virus    <- read_xlsx("checkv_provirus_boundary.xlsx", sheet = "checkV_virus")
checkv_provirus <- read_xlsx("checkv_provirus_boundary.xlsx", sheet = "checkV_provirus")

# Parse coordinates embedded in sequence names
checkv_virus_2 <- checkv_virus %>%
  separate(name,
           into   = c("contig", "start_end"),
           sep    = "_(?=[^_]+$)",
           remove = FALSE) %>%
  separate(start_end, into = c("start", "end"), sep = "-", remove = TRUE)

# For CheckV-trimmed proviruses, offset the CheckV boundary relative to the
# extracted sub-sequence coordinates to get absolute positions on the contig
checkv_provirus_2 <- checkv_provirus %>%
  separate(name,
           into   = c("contig", "start_end"),
           sep    = "_(?=[^_]+$)",
           remove = FALSE) %>%
  separate(start_end,    into = c("start", "end"),             sep = "-", remove = TRUE) %>%
  separate(border_checkv, into = c("start_checkv", "end_checkv"), sep = "-", remove = TRUE) %>%
  mutate(start_trim = as.numeric(start) + as.numeric(start_checkv) - 1,
         end_trim   = as.numeric(start) + as.numeric(end_checkv)   - 1)

# Combine into a unified boundary table
checkv_virus_3    <- checkv_virus_2[, 3:5]
checkv_provirus_3 <- checkv_provirus_2[, c(3, 9, 10)]
colnames(checkv_provirus_3) <- c("contig", "start", "end")
checkv_all_border <- rbind(checkv_virus_3, checkv_provirus_3)


# =============================================================================
# Part 6: Merge Overlapping Boundaries with IRanges
# =============================================================================
# When multiple tools predict overlapping regions on the same contig, IRanges
# reduce() merges them into the minimal set of non-overlapping intervals.

result_list <- checkv_all_border %>%
  mutate(start = as.numeric(start),
         end   = as.numeric(end)) %>%
  group_split(contig) %>%
  map(function(df) {
    ir <- IRanges(start = df$start, end = df$end)
    reduced_df <- IRanges::reduce(ir) %>% as.data.frame()
    reduced_df$contig <- unique(df$contig)
    reduced_df$name   <- paste0(reduced_df$contig, "_", seq_len(nrow(reduced_df)))
    return(reduced_df)
  })

result_table <- bind_rows(result_list)

checkv_all_border_f <- result_table[, c(4, 5, 1, 2, 3)]
colnames(checkv_all_border_f) <- c("contig", "virus_name", "start", "end", "length")

# Remove fragments shorter than 1 kb
checkv_all_border_f <- checkv_all_border_f %>%
  filter(length >= 1000) %>%
  mutate(
    name_seqkit   = paste0(contig, "_", start, "-", end, ":."),
    start_seqkit  = start - 1,
    end_seqkit    = end
  )

write.csv(checkv_all_border_f, "checkv_all_border_f.csv", row.names = FALSE)

# Export BED file and name replacement table for seqkit (see shell script Step 5)
checkv_all_border_f_seqkit <- checkv_all_border_f[, c(1, 7, 8)]
write.table(checkv_all_border_f_seqkit,
            "checkv_all_border_f_seqkit.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

seqkit_name_replace <- checkv_all_border_f[, c(6, 2)]
write.table(seqkit_name_replace,
            "provirus_all_boundary_f_name_replace.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# =============================================================================
# Part 7: rRNA Detection and Region Removal
# =============================================================================
# barrnap GFF outputs (bacteria, archaea, eukarya) are merged and filtered.
# Only complete rRNA annotations are kept (partial rRNA = filtering artefact).
# IRanges merges overlapping rRNA hits, then the rRNA-free sub-region of each
# provirus is determined based on rRNA midpoint position.

rrna_checkv_only <- read.delim("rrna/rrna_checkv_only_all.gff",
  header     = FALSE,
  col.names  = c("virus_id", "c2", "c3", "start", "end",
                 "evalue", "strand", "c8", "rRNA"))

rrna_provirus_trimmed <- read.delim("rrna/rrna_provirus_trimmed_boundary_all.gff",
  header     = FALSE,
  col.names  = c("virus_id", "c2", "c3", "start", "end",
                 "evalue", "strand", "c8", "rRNA"))

rrna_all <- rbind(rrna_checkv_only, rrna_provirus_trimmed) %>%
  select(virus_id, start, end, evalue, strand, rRNA) %>%
  arrange(virus_id, start) %>%
  filter(!grepl("partial", rRNA))  # retain only complete rRNA genes

# Merge overlapping rRNA hits per provirus with IRanges
rrna_result_list <- rrna_all %>%
  mutate(start = as.numeric(start),
         end   = as.numeric(end)) %>%
  group_split(virus_id) %>%
  map(function(df) {
    ir <- IRanges(start = df$start, end = df$end)
    reduced_df <- IRanges::reduce(ir) %>% as.data.frame()
    reduced_df$virus_id <- unique(df$virus_id)
    return(reduced_df)
  })

rrna_region_all <- bind_rows(rrna_result_list)[, c(4, 1, 2, 3)] %>%
  group_by(virus_id) %>%
  summarise(
    start        = min(start),
    end          = max(end),
    length_rrna  = end - start + 1,
    .groups      = "drop"
  )

# Attach provirus sequence lengths
checkv_only_length <- read.csv("vOTU_v2_pro_only_checkv.csv") %>%
  mutate(length = Length_kb * 1000) %>%
  select(1, length) %>%
  setNames(c("virus_id", "virus_length"))

provirus_trimmed_length <- read.csv("checkv_all_border_f.csv") %>%
  select(2, 5) %>%
  setNames(c("virus_id", "virus_length"))

provirus_all_length <- rbind(checkv_only_length, provirus_trimmed_length)

# Determine the rRNA-free region to retain:
#   If rRNA is in the first half of the sequence → keep the region after it
#   If rRNA is in the second half              → keep the region before it
rrna_region_all <- merge(rrna_region_all, provirus_all_length,
                          by = "virus_id", all.x = TRUE)

provirus_keep_region <- rrna_region_all %>%
  group_by(virus_id) %>%
  mutate(
    rrna_mid   = (start + end) / 2,

    keep_start = ifelse(rrna_mid < virus_length / 2, end + 1,       1),
    keep_end   = ifelse(rrna_mid < virus_length / 2, virus_length,  start - 1),
    keep_length = keep_end - keep_start + 1,

    seqkit_keep_start = keep_start - 1,
    name_seqkit = paste0(virus_id, "_", keep_start, "-", keep_end, ":.")
  ) %>%
  ungroup() %>%
  filter(keep_length >= 1000)  # discard sub-regions shorter than 1 kb

write.csv(provirus_keep_region,
          "rrna/provirus_all_rrna_trimmed_boundary.csv", row.names = FALSE)

# Export BED and renaming table for seqkit (see shell script Step 6)
provirus_keep_region_bed <- provirus_keep_region[, c(1, 10, 8)]
write.table(provirus_keep_region_bed,
            "rrna/provirus_keep_region.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

provirus_keep_region_seqkit_name <- provirus_keep_region[, c(11, 1)]
write.table(provirus_keep_region_seqkit_name,
            "rrna/provirus_keep_region_seqkit_name.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
