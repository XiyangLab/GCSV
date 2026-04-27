#!/bin/bash
# =============================================================================
# Viral Abundance and Diversity Estimation Pipeline
# =============================================================================
# Description: Quantifies vOTU abundance and estimates viral diversity using
#              two approaches:
#              (1) MetaPop — read recruitment pipeline that calculates
#                  coverage, microdiversity, and macrodiversity per sample;
#              (2) CoverM — fast RPKM-based abundance estimation with
#                  read identity and alignment filters.
#
# Tools required:
#   - bowtie2    : read alignment and reference index construction
#   - samtools   : BAM file generation and read counting
#   - MetaPop    : population-level diversity and abundance analysis
#   - CoverM     : contig-level RPKM abundance calculation
# =============================================================================


# -----------------------------------------------------------------------------
# Part 1: MetaPop
# -----------------------------------------------------------------------------
# Reference: https://github.com/metaGEM/MetaPop

# --- Step 1a: Build bowtie2 reference index ---
bowtie2-build \
    coldseeps_vOTU.fna \
    index/coldseeps_vOTU \
    --threads 30 

# --- Step 1b: Align clean reads to vOTU reference (per sample) ---
# Repeat for each sample; pipe directly to samtools to produce BAM files.
# Example for a single sample:
bowtie2 \
    -x index/coldseeps_vOTU \
    -1 /path/to/reads/SAMPLE_1.fastq \
    -2 /path/to/reads/SAMPLE_2.fastq \
    -p 30 \
    | samtools view -bS -@ 20 > bam/SAMPLE_vOTU.bam


# --- Step 1c: Generate CT file (mapped read counts per BAM) ---
# The CT file records total read counts for each BAM file, required by MetaPop.
# samtools view -c counts all mapped reads in a BAM file.
# Iterate over all BAM files to produce the CT file:
for bam_file in bam/*.bam; do
    sample_name=$(basename "${bam_file}" _vOTU.bam)
    read_count=$(samtools view -c "${bam_file}")
    echo -e "${sample_name}\t${read_count}"
done > ct_file.tsv

# --- Step 1d: Run MetaPop ---
metapop \
    -i bam \
    -r reference \
    -n ct_file.tsv \
    --threads 30 \
    -o metapop_output \
    --gene coldseeps_vOTU_gene.fna 


# -----------------------------------------------------------------------------
# Part 2: CoverM
# -----------------------------------------------------------------------------
# Reference: https://github.com/wwood/CoverM

# Set a custom TMPDIR to avoid filling the system /tmp partition
export TMPDIR=/path/to/tmp_coverm

mkdir -p vOTU_coverm

# Run CoverM per sample (example for a single-end sample)
coverm contig \
    --single /path/to/reads/SAMPLE.fastq \
    --reference coldseeps_vOTU.fna \
    --trim-min 0.10 \
    --trim-max 0.90 \
    --min-read-percent-identity 0.95 \
    --min-read-aligned-percent  0.75 \
    -m rpkm \
    -o vOTU_coverm/SAMPLE_rpkm.tsv \
    -t 20 

# --- Merge CoverM results into a single matrix ---
# Assumes all per-sample TSV files are in the vOTU_coverm/ directory.
# The first column (contig names) is taken from the first file;
# the second column (RPKM values) is appended from each subsequent file.

cd vOTU_coverm

files=(*.tsv)

# Initialise the merged file with contig IDs from the first sample
cut -f1 "${files[0]}" > ../coverm_coldseeps_vOTU_rpkm_merged.tsv

# Append RPKM columns from all samples
for file in "${files[@]}"; do
    cut -f2 "${file}" > tmp_col.tsv
    paste ../coverm_coldseeps_vOTU_rpkm_merged.tsv tmp_col.tsv > tmp_merged.tsv
    mv tmp_merged.tsv ../coverm_coldseeps_vOTU_rpkm_merged.tsv
done

rm -f tmp_col.tsv

