#!/bin/bash
# =============================================================================
# Viral Protein Functional Annotation Pipeline
# =============================================================================
# Description: Annotates viral proteins in three stages:
#   (1) Protein clustering — dereplicate vOTU proteins into viral Protein
#       Clusters (vPCs) using MMseqs2;
#   (2) Structure-based annotation — search vPCs against structural databases
#       (BFVD, PDB100, AFDB-SwissProt) using Foldseek, and against the PHROG
#       database using Phold;
#   (3) Sequence-based annotation — search vPCs against profile/sequence
#       databases (enVhog, eggNOG) using MMseqs2 and eggNOG-mapper.
#
# Tools required:
#   - Prodigal      : ORF prediction
#   - MMseqs2       : protein clustering and sequence database search
#   - Phold         : structure-based phage protein annotation (PHROG-based)
#   - Foldseek      : protein structure search against BFVD / PDB100 / AFDB-SwissProt
#   - eggNOG-mapper : functional annotation via eggNOG database (diamond mode)
#
# Companion scripts:
#   clean_h.py
#   add_annotation.py
#   enVhog_mmseqs.sh
#   mmseqs_best_hits.py
# =============================================================================


# =============================================================================
# Part 1: Protein Prediction and Clustering — Generate Viral Protein Clusters
# =============================================================================

prodigal -i vOTU_all_final.fna \
-o prodigal_vOTU/coldseeps_vOTU_gene.gff \
-d prodigal_vOTU/coldseeps_vOTU_gene.fna \
-a prodigal_vOTU/coldseeps_vOTU_protein.faa \
-p meta

# Cluster all vOTU proteins at 50% identity / 90% coverage to reduce
# redundancy. The representative sequences (vPCs) are used as input for
# all downstream annotation steps.

mmseqs easy-cluster \
    prodigal_vOTU/coldseeps_vOTU_protein.faa \
    mmseqs2/mmseqs_vOTU_id0.5_c0.9_clusters \
    tmp \
    --min-seq-id 0.5 \
    -c 0.9 \
    --cov-mode 1

# =============================================================================
# Part 2: Structure-based Annotation
# =============================================================================

# -----------------------------------------------------------------------------
# Step 2a: Phold
# -----------------------------------------------------------------------------
# Reference: https://github.com/gbouras13/phold
# Input file should contain fewer than 100k proteins.

# Step 1 — structure prediction
phold proteins-predict \
    -i vPCs.part_001.fasta \
    -o phold_output/phold_proteins_predict_001 \
    -t 16 \
    -d /path/to/phold_database

# Step 2 — compare predicted structures against the database
phold proteins-compare \
    -i vPCs.part_001.fasta \
    --predictions_dir phold_output/phold_proteins_predict_001 \
    -o phold_output/phold_proteins_compare_001 \
    -t 60 \
    -d /path/to/phold_database

# Merge results across all splits
# proteins-predict failures
find phold_proteins_predict_output_part*/ \
    -type f -name "fails.tsv" -exec cat {} + \
    > fails_all.tsv

# Main annotation file (per-CDS predictions)
find phold_proteins_compare_output_part*/ \
    -type f -name "phold_per_cds_predictions.tsv" -exec cat {} + \
    > phold_per_cds_predictions_all.tsv

# Filter: remove unannotated sequences ("No_PHROG") and apply e-value < 1e-10
# Note: results from sub-databases are already integrated into
#       phold_per_cds_predictions; this file is sufficient for downstream use.

awk -F'\t' '$2 != "No_PHROG" && $9 < 1e-10' phold_per_cds_predictions_all.tsv \
    > phold_per_cds_predictions_filtered_1e-10.tsv

# Summarise annotated functional categories
awk -F'\t' '{print $3}' phold_per_cds_predictions_filtered.tsv \
    | sort | uniq -c | sort -nr


# -----------------------------------------------------------------------------
# Step 2b: Foldseek — Search Against Structural Databases
# -----------------------------------------------------------------------------
# Reference: https://github.com/steineggerlab/foldseek
# GPU-accelerated structure search using ProstT5 embeddings.
# Databases used: BFVD, PDB100, AFDB-SwissProt.
# Large input files were split before running.
#
# All searches use:
#   -c 0.5 / --cov-mode 0  : 50% coverage relative to the query

# Prepare padded databases for GPU search (run once per database)
foldseek makepaddedseqdb software_db/bfvd_foldseekdb/bfvd \
    software_db/bfvd_foldseekdb/bfvd_pad
foldseek makepaddedseqdb software_db/pdb100/pdb100 \
    software_db/pdb100/pdb100_pad
foldseek makepaddedseqdb software_db/afdb_swissprot/afdb_swissprot \
    software_db/afdb_swissprot/afdb_swissprot_pad

FOLDSEEK_COLS="query,target,evalue,fident,alnlen,bits,mismatch,gapopen,qstart,qend,tstart,tend"

# Search against BFVD
foldseek easy-search \
    vPCs.part_001.fasta \
    software_db/bfvd_foldseekdb/bfvd_pad \
    results/foldseek_bfvd_001.m8 \
    tmp \
    --gpu 1 \
    --prostt5-model software_db/prostt5-f16.gguf \
    --format-mode 4 \
    -c 0.5 --cov-mode 0 \
    --format-output "${FOLDSEEK_COLS}"

# Search against PDB100
foldseek easy-search \
    vPCs.part_001.fasta \
    software_db/pdb100/pdb100_pad \
    results/foldseek_pdb100_001.m8 \
    tmp \
    --gpu 1 \
    --prostt5-model software_db/prostt5-f16.gguf \
    --format-mode 4 \
    -c 0.5 --cov-mode 0 \
    --format-output "${FOLDSEEK_COLS}"

# Search against AFDB-SwissProt
foldseek easy-search \
    vPCs.part_001.fasta \
    software_db/afdb_swissprot/afdb_swissprot_pad \
    results/foldseek_afdb_swissprot_001.m8 \
    tmp \
    --gpu 1 \
    --prostt5-model software_db/prostt5-f16.gguf \
    --format-mode 4 \
    -c 0.5 --cov-mode 0 \
    --threads 40 \
    --format-output "${FOLDSEEK_COLS}"

# --- Post-processing: extract top hit per query and filter by e-value < 1e-10 ---
# Applied to each database independently (PDB100 shown as example).

# Top hit per query (first occurrence = lowest e-value)
awk -F'\t' '!seen[$1]++' results/foldseek_pdb100_all.m8 \
    > results/tophit_pdb100_all.tsv

head -n 1 results/tophit_pdb100_all.tsv > results/tophit_pdb100_all_1e-10.tsv
awk -F'\t' '$3 < 1e-10' results/tophit_pdb100_all.tsv \
    >> results/tophit_pdb100_all_1e-10.tsv

# --- BFVD: extract UniProt ID from target field ---
# Target names follow the format <UniProtID>_<suffix>;
# extract the UniProt ID as a new column for downstream UniProt ID mapping.
awk 'BEGIN{FS=OFS="\t"}
     NR==1 {print $0, "target_id"}
     NR>1  {split($2,a,"_"); print $0, a[1]}' \
    results/tophit_bfvd_all_1e-10.tsv \
    > results/tophit_bfvd_all_1e-10_new.tsv

awk 'BEGIN{FS="\t"} NR>1{print $NF}' \
    results/tophit_bfvd_all_1e-10_new.tsv \
    | sort -u > add_annotation/bfvd_target_id.txt

# --- PDB100: add annotation from header file ---
cat software_db/pdb100/pdb100_pad_h > pdb_h.tsv
python3 clean_h.py pdb_h.tsv pdb_h_clean.tsv
python3 add_annotation.py \
    add_annotation/pdb_h_clean.tsv \
    results/tophit_pdb100_all_1e-10.tsv \
    results/tophit_pdb100_all_1e-10_annotation.tsv
cut -f1,3,13 results/tophit_pdb100_all_1e-10_annotation.tsv \
    > results/tophit_pdb100_all_1e-10_only_annotation.tsv

# --- AFDB-SwissProt: merge splits, filter, and extract UniProt ID ---
# Target names follow the format AF-<UniProtID>-F1-model_v4.

# Extract UniProt ID as a new column
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1 {print $0, "uniprot_id"; next}
    {split($2,a,"-"); print $0, a[2]}' \
    results/tophit_afdb_swissprot_all_1e-10.tsv \
    > results/tophit_afdb_swissprot_all_1e-10_new.tsv

# Unique UniProt IDs for batch ID mapping on the UniProt website
cut -f2 results/tophit_afdb_swissprot_all_1e-10.tsv \
    > add_annotation/afdb_swissprot_uniprot_id.tsv
awk '!seen[$0]++' add_annotation/afdb_swissprot_uniprot_id.tsv \
    > add_annotation/afdb_swissprot_uniprot_id_nr.tsv
awk -F'-' '{print $2}' add_annotation/afdb_swissprot_uniprot_id_nr.tsv \
    > add_annotation/afdb_swissprot_uniprot_id_clean.txt


# =============================================================================
# Part 3: Sequence-based Annotation
# =============================================================================

# -----------------------------------------------------------------------------
# Step 3a: enVhog — Viral Orthologous Group Annotation
# -----------------------------------------------------------------------------
# See companion script: enVhog_mmseqs.sh
# Searches vPCs against the enVhog profile database using MMseqs2.
# Best hits are filtered by identity (>=30%), e-value (<=1e-5), and
# query coverage (>=80%) using the companion Python script mmseqs_best_hits.py.

bash enVhog_mmseqs.sh

# -----------------------------------------------------------------------------
# Step 3b: eggNOG-mapper — COG / KEGG / GO Annotation
# -----------------------------------------------------------------------------
# Runs diamond-based search against the eggNOG database restricted to
# the Viruses taxonomy scope.

emapper.py \
    -i mmseqs2/mmseqs_vOTU_id0.5_c0.9_clusters_rep_seq.fasta \
    --output vPC_eggnog \
    --data_dir /path/to/eggnog-mapper \
    -m diamond \
    --cpu 35 \
    --tax_scope Viruses
