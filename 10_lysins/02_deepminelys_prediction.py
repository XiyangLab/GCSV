#!/bin/bash

# ============================================================
# Step 02: DeepMineLys prediction
# Description:
#   - Split input FASTA into individual sequences
#   - Run DeepMineLys prediction (via SLURM jobs)
#   - Collect lysin candidates (label = 1 or 2)
#
# Input:
#   - input/candidate_sequences.faa  (from Step 01, filtered ≥20% identity)
#
# Output:
#   - DeepMineLys prediction results (*.csv)
#   - lysin_candidates.faa
# ============================================================cc

# ========== Configuration ==========
INPUT_FASTA="input/candidate_sequences.faa"
WORKDIR="deepminelys_work"
SPLIT_DIR="${WORKDIR}/split_fasta"
LIST_FILE="${WORKDIR}/input_list.txt"

mkdir -p ${WORKDIR} ${SPLIT_DIR}

# ========== Step 1: Split FASTA ==========
echo "Splitting FASTA into individual sequences..."

awk '/^>/{s=sprintf("'"${SPLIT_DIR}"'/%05d.faa", ++d)} {print > s}' ${INPUT_FASTA}

# ========== Step 2: Generate input list ==========
echo "Generating input list..."

ls ${SPLIT_DIR}/*.faa > ${LIST_FILE}

# (Optional: Only take the first N tests)
# head -n 5 ${LIST_FILE} > ${WORKDIR}/test_input.txt

# ========== Step 3: Submit jobs ==========
echo "Submitting DeepMineLys jobs..."

bash scripts/run_deepminelys_jobs.sh ${LIST_FILE}

echo "DeepMineLys jobs submitted."