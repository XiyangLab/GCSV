#!/bin/bash

# ============================================================
# Step 04: Domain annotation and filtering
# Tool: InterProScan (v5.0.39)
# Description:
#   - Annotate lysin candidates with InterProScan
#   - Filter functional domains (E-value ≤ 1e-5)
#   - Perform domain-level refinement in R
#
# Input:
#   - input/lysin_candidates.faa
#
# Output:
#   - interproscan_results.tsv
#   - final_lysin_domains.xlsx
# ============================================================


# ============================================================
# Input preparation:
#   Lysin candidates from Step 02 (DeepMineLys) and Step 03
#   (structure-based identification) were merged and
#   fully deduplicated to generate the final input FASTA file.
# ============================================================


# ========== Configuration ==========
INPUT_FASTA="input/lysin_candidates.faa"
WORKDIR="interproscan_work"
OUT_TSV="${WORKDIR}/interproscan_results.tsv"

THREADS=32

mkdir -p ${WORKDIR}


# ========== Step 1: Run InterProScan ==========
echo "Running InterProScan..."

interproscan.sh \
    -i ${INPUT_FASTA} \
    -f tsv \
    -goterms -iprlookup -pa \
    --cpu ${THREADS} \
    -o ${OUT_TSV}


# ========== Step 2: Domain filtering ==========
echo "Running domain filtering in R..."

Rscript scripts/filter_domains.R \
    ${OUT_TSV} \
    config/functional_ipr.txt \
    ${WORKDIR}/final_lysin_domains.xlsx


echo "Done!"