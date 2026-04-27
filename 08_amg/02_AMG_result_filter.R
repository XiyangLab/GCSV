library(tidyverse)
library(readxl)

# VIRBRANT--------
VIBRANT_new_pro <- read.delim("VIBRANT_annotations_all_final.tsv")
VIBRANT_new_pro_AMG_step0 <- VIBRANT_new_pro %>%
  filter(AMG == "AMG")

eggnog_VIBRANT_AMG_step0 <- read_xlsx("VIBRANT_AMGs_step0.emapper.annotations.xlsx")

# ─────────────────────────────────────────────────────────────────
# Preparation: extract position from protein name and build lookup tables
# ─────────────────────────────────────────────────────────────────

get_pos <- function(x) as.integer(str_extract(x, "\\d+$"))

all_proteins <- VIBRANT_new_pro %>%
  mutate(protein_pos = get_pos(protein))

AMG_raw <- VIBRANT_new_pro_AMG_step0 %>%
  mutate(protein_pos = get_pos(protein))

# Full original AMG protein vector (used in Step 1 & Step 3 to identify AMGs)
all_AMG_proteins <- AMG_raw$protein

# Scaffold -> proteins ordered by position (list column)
scf2pros_ordered <- all_proteins %>%
  arrange(scaffold, protein_pos) %>%
  group_by(scaffold) %>%
  summarise(pros_ordered = list(protein), .groups = "drop")


# ══════════════════════════════════════════════════════════════════
# Step 1: Filter tail AMGs
#   Mirrors Python's check_left() / check_right():
#   An AMG is a tail if ALL proteins to its left in scaffold order are AMGs,
#   OR all proteins to its right are AMGs (including the edge case where
#   it sits at position 1 or the last position).
# ══════════════════════════════════════════════════════════════════

check_tail_amg <- function(amg, pros_ordered) {
  n   <- length(pros_ordered)
  idx <- which(pros_ordered == amg)
  
  # Explicit boundary guards to avoid seq(from > to) generating descending sequences
  left_indices  <- if (idx > 1) seq(1, idx - 1) else integer(0)
  right_indices <- if (idx < n) seq(idx + 1, n) else integer(0)
  
  left_all_amg  <- all(pros_ordered[left_indices]  %in% all_AMG_proteins)
  right_all_amg <- all(pros_ordered[right_indices] %in% all_AMG_proteins)
  
  return(left_all_amg | right_all_amg)
}

AMG_step1 <- AMG_raw %>%
  left_join(scf2pros_ordered, by = "scaffold") %>%
  rowwise() %>%
  mutate(is_tail = check_tail_amg(protein, pros_ordered)) %>%
  ungroup() %>%
  filter(!is_tail) %>%
  select(-pros_ordered, -is_tail)          # <-- drop temp columns

cat("Remaining AMGs after Step 1:", nrow(AMG_step1), "\n")


# ══════════════════════════════════════════════════════════════════
# Step 2: Filter AMGs with any v-score (KO or Pfam) >= 1
#   Matches Python: only filter when the value exists AND >= 1
# ══════════════════════════════════════════════════════════════════

AMG_step2 <- AMG_step1 %>%
  filter(
    !( (!is.na(KO.v.score)   & KO.v.score   >= 1) |
         (!is.na(Pfam.v.score) & Pfam.v.score >= 1) )
  )

cat("Remaining AMGs after Step 2:", nrow(AMG_step2), "\n")


# ══════════════════════════════════════════════════════════════════
# Step 3: Filter AMGs whose flanking genes all have KO.v.score < 0.25
#   Mirrors Python exactly:
#   - Each AMG evaluated individually (not as tandem groups)
#   - Up to 2 non-AMG flanking genes collected on each side (flanking_num = 2),
#     skipping over other AMGs using the FULL original AMG list
#   - Only genes WITH a KO.v.score annotation contribute to the low-score count
#   - Filter if: low_score_count == total number of flanking genes
# ══════════════════════════════════════════════════════════════════

# Named vector for fast KO.v.score lookup: protein -> KO.v.score
prot_vscore_lookup <- setNames(all_proteins$KO.v.score, all_proteins$protein)

get_flanking_low_flag <- function(amg, pros_ordered, flanking_num = 2) {
  idx <- which(pros_ordered == amg)
  n   <- length(pros_ordered)
  
  # Collect up to flanking_num non-AMG genes on the left
  left_genes <- c()
  if (idx > 1) {
    for (i in seq(idx - 1, 1, by = -1)) {
      if (!pros_ordered[i] %in% all_AMG_proteins) {
        left_genes <- c(left_genes, pros_ordered[i])
        if (length(left_genes) == flanking_num) break
      }
    }
  }
  
  # Collect up to flanking_num non-AMG genes on the right
  right_genes <- c()
  if (idx < n) {
    for (i in seq(idx + 1, n, by = 1)) {
      if (!pros_ordered[i] %in% all_AMG_proteins) {
        right_genes <- c(right_genes, pros_ordered[i])
        if (length(right_genes) == flanking_num) break
      }
    }
  }
  
  flanking_genes <- c(left_genes, right_genes)
  if (length(flanking_genes) == 0) return(FALSE)  # No flanking genes: do not filter
  
  # Count flanking genes with a non-NA KO.v.score < 0.25
  low_score_count <- 0
  for (gene in flanking_genes) {
    vscore <- prot_vscore_lookup[gene]
    if (!is.na(vscore) && vscore < 0.25) {
      low_score_count <- low_score_count + 1
    }
  }
  
  return(low_score_count == length(flanking_genes))
}

AMG_step3 <- AMG_step2 %>%
  left_join(scf2pros_ordered, by = "scaffold") %>%
  rowwise() %>%
  mutate(all_flanking_low = get_flanking_low_flag(protein, pros_ordered)) %>%
  ungroup() %>%
  filter(!all_flanking_low) %>%
  select(-pros_ordered, -all_flanking_low)

cat("Remaining AMGs after Step 3:", nrow(AMG_step3), "\n")


# ══════════════════════════════════════════════════════════════════
# Step 4: Filter AMGs with COG category "T" or "B"
#   AMGs with no COG assignment (NA) are retained.
#   Note: eggNOG query column is typically named "#query"
# ══════════════════════════════════════════════════════════════════

eggnog_cog <- eggnog_VIBRANT_AMG_step0 %>%
  select(protein = query, COG_category)

AMG_step4 <- AMG_step3 %>%
  left_join(eggnog_cog, by = "protein") %>%
  filter(
    is.na(COG_category) |
      !grepl("[TB]", COG_category)
  )

cat("Remaining AMGs after Step 4:", nrow(AMG_step4), "\n")


# ══════════════════════════════════════════════════════════════════
# Final output
# ══════════════════════════════════════════════════════════════════

VIBRANT_AMG_filtered_final <- AMG_step4

cat("\nFiltering summary:\n",
    "  Input AMGs            :", nrow(VIBRANT_new_pro_AMG_step0), "\n",
    "  After Step 1 (tail)   :", nrow(AMG_step1), "\n",
    "  After Step 2 (v-score):", nrow(AMG_step2), "\n",
    "  After Step 3 (flank)  :", nrow(AMG_step3), "\n",
    "  After Step 4 (COG)    :", nrow(AMG_step4), "\n")

# DRAM-v--------
amg_summary <- read.delim('DRAM-v_amg_summary.tsv')
trna <- read.delim('trnas.tsv')
suspicious_gene <- read.delim('DRAM-v_annotations_all_suspicious-gene.tsv')
ir <- read.table('new_id_ir_no_duplicates.txt')

library(tidyverse)

trna_55 <- trna[trna$Score > 55, ]
#step2
amg_summary$scaffold <- trimws(amg_summary$scaffold)
trna_55$Name <- trimws(trna_55$Name)
amg_summary_trna <- amg_summary %>%
  filter(scaffold %in% trna_55$Name)
amg_summary_trna_ir <- amg_summary_trna %>% filter(scaffold %in% ir$V1)

#step3
amg_summary_trna_ir_filtered <- amg_summary_trna_ir %>%  filter(!scaffold %in% suspicious_gene$scaffold)
amg_summary_trna_ir_filtered <- amg_summary_trna_ir_filtered %>%  filter(!grepl("T", amg_flags))
DRAMv_AMG_filtered_final <- amg_summary_trna_ir_filtered %>%
  filter(gene_id != "")



#merge_VIBRANT_DRAMv-------
VIBRANT_results <- VIBRANT_AMG_filtered_final[,c(2,1,3,5)]
DRAMv_results <- DRAMv_AMG_filtered_final[,c(2,1,3,4)]
colnames(VIBRANT_results) <- colnames(DRAMv_results)
VIBRANT_results_2 <- VIBRANT_results %>%
  mutate(
    scaffold = sub(" \\S+$", "", scaffold),
    gene = sub(" \\S+_", "_", gene)
  )

DRAMv_results_2 <- DRAMv_results %>%
  mutate(
    scaffold = sub("__full-cat_\\d+$", "", scaffold),
    gene = sub("__full-cat_\\d+", "", gene)
  )

both <- inner_join(DRAMv_results_2, VIBRANT_results_2, 
                   by = c("gene", "gene_id", "scaffold")) %>%
  mutate(method = "both") %>%
  select(-gene_description.x)

dram_only <- DRAMv_results_2 %>%
  filter(!gene %in% VIBRANT_results_2$gene) %>%
  mutate(method = "DRAM-v")

vibrant_only <- VIBRANT_results_2 %>%
  filter(!gene %in% DRAMv_results_2$gene) %>%
  mutate(method = "VIBRANT")

colnames(both) <- colnames(dram_only)
merged_results <- bind_rows(both, dram_only, vibrant_only)
conflict <- inner_join(DRAMv_results_2, VIBRANT_results_2,
                       by = "gene", suffix = c("_DRAMv", "_VIBRANT")) %>%
  filter(gene_id_DRAMv != gene_id_VIBRANT)

#viral_hallmark------
vs2 <- read.delim("final-viral-score.tsv")
vs2_2 <- vs2[,c(1,8)] %>%
  group_by(seqname) %>%
  summarise(hallmark_vs2 = sum(hallmark, na.rm = TRUE))

genomad <- read.delim("genomad_final_virus_summary.tsv")
genomad_2 <- genomad[,c(1,9)] %>%
  group_by(seq_name) %>%
  summarise(hallmark_genomad = sum(n_hallmarks, na.rm = TRUE))

AMG_1 <- merged_results
AMG_2 <- merge(AMG_1, vs2_2, by.x = "scaffold", by.y = "seqname", all.x = T)
AMG_3 <- merge(AMG_2, genomad_2, by.x = "scaffold", by.y = "seq_name", all.x = T)
AMG_with_hallmark <- AMG_3 %>%
  filter(
    (hallmark_vs2 != 0 & !is.na(hallmark_vs2)) |
      (hallmark_genomad != 0 & !is.na(hallmark_genomad))
  )

#kegg_pathway------
kegg_pathway <- read.delim("ko_pathway1.tsv")
AMG_pathway <- merge(AMG_with_hallmark[,c(2,3)], kegg_pathway,
                        by.x = "gene_id", by.y = "KO")

level2_counts <- AMG_v2_pathway[,1:6] %>%
  filter(level1_pathway_name == "Metabolism") %>%
  unique() %>%  
  group_by(level2_pathway_id, level2_pathway_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(percentage = count/nrow(AMG_v2)*100)
