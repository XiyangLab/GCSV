#!/bin/bash
# =============================================================================
# Viral Sequence Identification Pipeline
# =============================================================================
# Description: Identifies viral sequences from assembled metagenomic contigs
#              using multiple tools (VIBRANT, geNomad, VirSorter2) and
#              validates results with CheckV.
#
# Tools required:
#   - seqkit       : sequence filtering and manipulation
#   - VIBRANT      : virus identification via neural networks
#   - geNomad      : end-to-end virus/plasmid identification
#   - VirSorter2   : virus identification via machine learning
#   - CheckV       : viral sequence quality assessment
#
# =============================================================================


# -----------------------------------------------------------------------------
# Step 1: Sequence Length Filtering
# -----------------------------------------------------------------------------

INPUT_FASTA="contigs.fasta"

# Filter sequences >= 5 kb
seqkit seq -m 5000 -g "${INPUT_FASTA}" \
    > contigs_5000bp.fasta

# Prepare two subsets for geNomad (uses different score thresholds per length):
#   - 5–10 kb sequences
seqkit seq -m 5000 -M 9999 -g contigs_5000bp.fasta \
    > contigs_5000_10000bp.fasta

#   - >= 10 kb sequences
seqkit seq -m 10000 -g contigs_5000bp.fasta \
    > contigs_10000bp.fasta


# -----------------------------------------------------------------------------
# Step 2: VIBRANT
# -----------------------------------------------------------------------------
# Note: Sequences with 'fragment' in their names are predicted proviruses.

python3 VIBRANT_run.py \
    -i contigs_5000bp.fasta \
    -t 30 \
    -l 5000 \
    -folder VIBRANT_output


# -----------------------------------------------------------------------------
# Step 3: geNomad
# -----------------------------------------------------------------------------
# Score thresholds applied after the run (different per length class):
#   >= 10 kb : virus_score >= 0.8 AND (n_hallmarks >= 1 OR virus_enrichment_ratio > 5.0)
#   5–10 kb  : virus_score >= 0.9 AND n_hallmarks >= 1 AND virus_enrichment_ratio > 2.0


# Run geNomad on the >= 10 kb subset
genomad end-to-end --cleanup \
    contigs_10000bp.fasta \
    genomad_output_10000bp \
    genomad_db

# Run geNomad on the 5–10 kb subset
genomad end-to-end --cleanup \
    contigs_5000_10000bp.fasta \
    genomad_output_5000_10000bp \
    genomad_db

# Filter >= 10 kb results
# Columns: $7 = virus_score, $9 = n_hallmarks, $10 = virus_enrichment_ratio
awk -F '\t' '($7 >= 0.8 && ($9 >= 1 || $10 > 5.0))' \
    genomad_output_10000bp/contigs_10000bp_virus_summary.tsv \
    > contigs_10000bp_virus_summary_filtered.tsv

# Filter 5–10 kb results 
awk -F '\t' '($7 >= 0.9 && $9 >= 1 && $10 > 2.0)' \
    genomad_output_5000_10000bp/contigs_5000_10000bp_virus_summary.tsv \
    > contigs_5000_10000bp_virus_summary_filtered.tsv


# -----------------------------------------------------------------------------
# Step 4: VirSorter2
# -----------------------------------------------------------------------------

virsorter run \
    -w virsorter2_output_001 \
    -i contigs_5000bp.fasta \
    --min-length 5000 \
    -j 12 all 


# Filter VirSorter2 results
# Columns: $4 = max_score, $7 = hallmark gene count
#
# High-confidence set (score >= 0.8, >= 1 hallmark gene)
awk -F '\t' '($4 >= 0.8 && $7 >= 1 && $4 != "nan")' \
    7sediments_L-final-viral-score_5000bp.tsv \
    > vs2_0.8_5000bp.tsv

# Relaxed set (score >= 0.5) — used for merging with VIBRANT results
awk -F '\t' '($4 >= 0.5 && $4 != "nan")' \
    7sediments_L-final-viral-score_5000bp.tsv \
    > vs2_0.5_5000bp.tsv

# Notes on VirSorter2 sequence suffixes:
#   ||full      : complete viral contig
#   ||{i}_partial : provirus 
#   ||lt2gene   : short sequences (< 2 genes) with hallmark genes — excluded from downstream analysis


# -----------------------------------------------------------------------------
# Step 5: Merge VIBRANT and VirSorter2 Results
# -----------------------------------------------------------------------------
# Merge is performed in R (see companion R script).
# Key considerations when merging:
#   - VirSorter2 appends suffixes (||full, ||{i}_partial) to sequence IDs
#   - VIBRANT appends 'fragment' to provirus sequence IDs
#   - Harmonize IDs before merging to avoid duplicates

# After R processing, retrieve the merged sequence ID list:
#   Output file: VIBRANT_vs2_0.5.txt

# Extract merged sequences from the original 5 kb-filtered FASTA
seqkit grep \
    -f VIBRANT_vs2_0.5.txt \
    contigs_5000bp.fasta \
    > VIBRANT_vs2_0.5.fna


# -----------------------------------------------------------------------------
# Step 6: CheckV — Quality Assessment
# -----------------------------------------------------------------------------

checkv end_to_end \
    VIBRANT_vs2_0.5.fna \
    checkv_VIBRANT_vs2_0.5 \
    -t 4 \
    -d checkv-db-v1.5 

# Key output files:
#   quality_summary.tsv  : completeness, contamination, quality tier per sequence
#   proviruses.fna       : trimmed provirus sequences
#   viruses.fna          : complete viral sequences
