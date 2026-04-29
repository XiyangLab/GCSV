#!/bin/bash
# =============================================================================
# Sequence similarity search and extraction pipeline (MMseqs2 + seqkit)
# =============================================================================
# Purpose:
#   1. Build query and reference MMseqs2 databases
#   2. Perform sequence similarity search
#   3. Extract query sequences with significant hits
#   4. Remove stop codon "*" from FASTA sequences
#
# Requirements:
#   - MMseqs2
#   - seqkit
# =============================================================================

set -e  # stop if any command fails

# =========================
# Input files (USER INPUT)
# =========================
QUERY_FASTA="query_sequences.fasta"
REFERENCE_CSV="reference_database.csv"

# =========================
# Output directory
# =========================
WORKDIR="mmseqs_work"
mkdir -p ${WORKDIR}

REF_FAA="${WORKDIR}/reference.faa"
QUERY_DB="${WORKDIR}/queryDB"
REF_DB="${WORKDIR}/refDB"
RESULT_DB="${WORKDIR}/resultDB"
RESULT_M8="${WORKDIR}/results.m8"
HIT_IDS="${WORKDIR}/query_hit_ids.txt"
EXTRACT_FASTA="${WORKDIR}/extracted_hits.faa"
CLEAN_FASTA="${WORKDIR}/extracted_hits_clean.faa"

# =========================
# Step 1: Convert reference CSV to FASTA
# =========================
# Assumption:
#   Column 1 = sequence ID
#   Column 2 = protein sequence

awk -F',' 'NR>1 {print ">"$1; print $2}' ${REFERENCE_CSV} > ${REF_FAA}

# =========================
# Step 2: Build MMseqs2 databases
# =========================
mmseqs createdb ${REF_FAA} ${REF_DB}
mmseqs createdb ${QUERY_FASTA} ${QUERY_DB}

# =========================
# Step 3: Sequence similarity search
# =========================
mmseqs search \
    ${QUERY_DB} \
    ${REF_DB} \
    ${RESULT_DB} \
    ${WORKDIR}/tmp \
    --min-seq-id 0.2

# =========================
# Step 4: Convert alignment output
# =========================
mmseqs convertalis \
    ${QUERY_DB} \
    ${REF_DB} \
    ${RESULT_DB} \
    ${RESULT_M8}

# =========================
# Step 5: Extract query IDs with hits
# =========================
awk -F '\t' '{print $1}' ${RESULT_M8} | sort | uniq > ${HIT_IDS}

# =========================
# Step 6: Extract sequences
# =========================
seqkit grep -f ${HIT_IDS} ${QUERY_FASTA} > ${EXTRACT_FASTA}

# =========================
# Step 7: Remove stop codon "*"
# =========================
sed 's/\*//g' ${EXTRACT_FASTA} > ${CLEAN_FASTA}

echo "Pipeline completed successfully!"
echo "Final output: ${CLEAN_FASTA}"