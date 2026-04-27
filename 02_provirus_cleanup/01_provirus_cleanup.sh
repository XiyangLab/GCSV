#!/bin/bash
# =============================================================================
# Provirus Boundary Refinement and Cleanup
# =============================================================================
# Description: Refines the genomic boundaries of predicted proviruses by
#              integrating boundary predictions from geNomad and VirSorter2,
#              trimming with CheckV, merging overlapping regions with IRanges
#              (see companion R script), and removing rRNA-containing regions.
#
# Workflow overview:
#   1. Match geNomad and VirSorter2 provirus boundaries in R
#   2. Extract provirus sequences using seqkit + BED file
#   3. Run CheckV contamination mode to trim host flanking regions
#   4. Merge trimmed boundaries using IRanges in R (see IRanges.R)
#   5. Extract final trimmed provirus sequences
#   6. Detect and remove rRNA-containing regions
# =============================================================================


# -----------------------------------------------------------------------------
# Step 2: Extract Initial Provirus Sequences by Boundary
# -----------------------------------------------------------------------------
# BED file prepared in R by matching geNomad (*_summary.tsv) and VirSorter2
# (*-viral-boundary.tsv, columns: trim_bp_start / trim_bp_end) boundary
# predictions. See companion R script: 02a_provirus_boundary_merge.R
#
# Sequences are extracted from each sample's combined viral FASTA and
# appended into a single output file.

BED_INITIAL="provirus_boundary_vs2_genomad_seqkit.bed"
OUTPUT_INITIAL="provirus_boundary_vs2_genomad.fna"

seqkit subseq --bed "${BED_INITIAL}" sample1_combined.fna  > "${OUTPUT_INITIAL}"

# -----------------------------------------------------------------------------
# Step 3: CheckV Contamination Mode
# -----------------------------------------------------------------------------
# Trims remaining host sequence flanking the integrated provirus.
# Output proviruses.fna and viruses.fna contain CheckV-adjusted boundaries.

checkv contamination \
    "${OUTPUT_INITIAL}" \
    checkv_contamination_output \
    -t 4 \
    -d /path/to/checkv-db 

# -----------------------------------------------------------------------------
# Step 4: Merge Boundaries with IRanges (R)
# -----------------------------------------------------------------------------
# IRanges resolves overlapping or conflicting boundary predictions from
# geNomad, VirSorter2, and CheckV into a single non-redundant coordinate
# per provirus.
# Output BED file: checkv_all_border_f_seqkit.bed


# -----------------------------------------------------------------------------
# Step 5: Extract Final Trimmed Provirus Sequences
# -----------------------------------------------------------------------------
# Use the IRanges-reconciled BED file to cut the definitive provirus sequences.

BED_FINAL="checkv_all_border_f_seqkit.bed"
OUTPUT_TRIMMED="provirus_trimmed_boundary.fna"

seqkit subseq --bed "${BED_FINAL}" sample1_combined.fna  > "${OUTPUT_TRIMMED}"

# Rename sequences using a two-column lookup table (old_name <tab> new_name)
cp "${OUTPUT_TRIMMED}" provirus_trimmed_boundary_rename.fna

while read old new; do
    sed -i "s/${old}/${new}/g" provirus_trimmed_boundary_rename.fna
done < provirus_trimmed_boundary_name_replace.txt


# -----------------------------------------------------------------------------
# Step 6: rRNA Detection and Removal
# -----------------------------------------------------------------------------
# Proviruses retaining rRNA genes are likely to be mis-trimmed.
# barrnap is run across all three kingdoms to maximise sensitivity.


mkdir -p rrna

# --- CheckV-only provirus sequences ---
INPUT_CHECKV="vOTU_v2_pro_only_checkv.fna"

barrnap -o rrna/rrna_checkv_bac.fa \
    --kingdom bac --threads 4 --evalue 1e-10 \
    "${INPUT_CHECKV}" > rrna/rrna_checkv_bac.gff

barrnap -o rrna/rrna_checkv_arc.fa \
    --kingdom arc --threads 4 --evalue 1e-10 \
    "${INPUT_CHECKV}" > rrna/rrna_checkv_arc.gff

barrnap -o rrna/rrna_checkv_euk.fa \
    --kingdom euk --threads 4 --evalue 1e-10 \
    "${INPUT_CHECKV}" > rrna/rrna_checkv_euk.gff

# --- Boundary-trimmed provirus sequences ---
INPUT_TRIMMED="provirus_trimmed_boundary_rename.fna"

barrnap -o rrna/rrna_provirus_trimmed_bac.fa \
    --kingdom bac --threads 4 --evalue 1e-10 \
    "${INPUT_TRIMMED}" > rrna/rrna_provirus_trimmed_bac.gff

barrnap -o rrna/rrna_provirus_trimmed_arc.fa \
    --kingdom arc --threads 4 --evalue 1e-10 \
    "${INPUT_TRIMMED}" > rrna/rrna_provirus_trimmed_arc.gff

barrnap -o rrna/rrna_provirus_trimmed_euk.fa \
    --kingdom euk --threads 4 --evalue 1e-10 \
    "${INPUT_TRIMMED}" > rrna/rrna_provirus_trimmed_euk.gff

# Merge all GFF results, remove partial rRNA annotations, and determine
# rRNA-free regions in R 
# Output from R:
#   provirus_keep_region.bed          : sub-regions to retain after rRNA removal
#   provirus_keep_region_seqkit_name.txt : name replacement table (old <tab> new)


# -----------------------------------------------------------------------------
# Step 6 (continued): Extract rRNA-free Provirus Sequences
# -----------------------------------------------------------------------------
# Pool both provirus sets, extract rRNA-trimmed sub-regions, rename them,
# then combine with sequences that had no rRNA (extracted by ID exclusion).

# Pool all proviruses
cat "${INPUT_CHECKV}"   > provirus_all_v2.fna
cat "${INPUT_TRIMMED}" >> provirus_all_v2.fna

# Extract and rename the rRNA-trimmed sub-regions
seqkit subseq --bed provirus_keep_region.bed provirus_all_v2.fna \
    > provirus_rrna_trimmed.fna

dos2unix provirus_keep_region_seqkit_name.txt

cp provirus_rrna_trimmed.fna provirus_rrna_trimmed_rename.fna

while read old new; do
    sed -i "s/${old}/${new}/g" provirus_rrna_trimmed_rename.fna
done < provirus_keep_region_seqkit_name.txt

# Extract sequences that did NOT contain rRNA (complement set)
cut -f 2 provirus_keep_region_seqkit_name.txt > provirus_rrna_trimmed.id

seqkit grep -f provirus_rrna_trimmed.id -v \
    -o provirus_without_rrna.fna \
    provirus_all_v2.fna

# Combine into final provirus dataset
cat provirus_without_rrna.fna        > provirus_all_final.fna
cat provirus_rrna_trimmed_rename.fna >> provirus_all_final.fna
