#!/bin/bash

# ============================================================
# Step 03: Structure-based identification of lysins
# Tools: Phold + Foldseek
# Description:
#   - Aggregate Phold annotation results
#   - Perform structure-based similarity search (Foldseek)
#   - Filter high-confidence hits
#   - Extract lysin candidates based on functional keywords
#
# Input:
#   - input/viral_proteins.faa
#
# Output:
#   - structure_lysin_candidates.tsv
# ============================================================



# ========== Configuration ==========
INPUT_FASTA="input/viral_proteins.faa"
WORKDIR="structure_work"
RESULTS_DIR="${WORKDIR}/results"
TMPDIR="${WORKDIR}/tmp"

BFVD_DB="db/bfvd"
PDB_DB="db/pdb100"
AFDB_DB="db/afdb_swissprot"

mkdir -p ${RESULTS_DIR} ${TMPDIR}


# ============================================================
# Step 1: Aggregate Phold results
# ============================================================

echo "Aggregating Phold results..."

find phold_output/ -type f -name "phold_per_cds_predictions.tsv" \
    -exec cat {} + > ${RESULTS_DIR}/phold_all.tsv

# remove "No_PHROG"
awk -F'\t' '$2 != "No_PHROG"' ${RESULTS_DIR}/phold_all.tsv > \
    ${RESULTS_DIR}/phold_filtered.tsv

# remove duplicate headers
awk 'BEGIN {found=0} /bitscore/ && found++ {next} {print}' \
    ${RESULTS_DIR}/phold_filtered.tsv > \
    ${RESULTS_DIR}/phold_clean.tsv

# e-value filtering
awk -F'\t' '$8 < 1e-10' ${RESULTS_DIR}/phold_clean.tsv > \
    ${RESULTS_DIR}/phold_final.tsv


# ============================================================
# Step 2: Foldseek search
# ============================================================

echo "Running Foldseek search..."

foldseek easy-search \
    ${INPUT_FASTA} \
    ${BFVD_DB} \
    ${RESULTS_DIR}/bfvd_results.tsv \
    ${TMPDIR} \
    --format-mode 4 \
    -c 0.5 \
    --cov-mode 0 \
    --format-output "query,target,evalue,fident,alnlen,bits" \
    --threads 16

# extract top hit
awk -F'\t' '!seen[$1]++' ${RESULTS_DIR}/bfvd_results.tsv > \
    ${RESULTS_DIR}/bfvd_tophit.tsv

# e-value filtering
awk -F'\t' '$3 < 1e-10' ${RESULTS_DIR}/bfvd_tophit.tsv > \
    ${RESULTS_DIR}/bfvd_final.tsv


# ============================================================
# Step 3: Extract lysin candidates
# ============================================================

echo "Extracting lysin candidates..."

# from Phold (functional keywords)
awk -F'\t' '
BEGIN{IGNORECASE=1}
$0 ~ /endolysin|amidase|lysozyme|lysM/ {print $1}
' ${RESULTS_DIR}/phold_final.tsv > phold_ids.txt

# from Foldseek (target annotation contains lysin)
awk -F'\t' '
BEGIN{IGNORECASE=1}
$2 ~ /lysin|endolysin/ {print $1}
' ${RESULTS_DIR}/bfvd_final.tsv > bfvd_ids.txt

# merge candidate IDs
cat phold_ids.txt bfvd_ids.txt | sort | uniq > lysin_ids.txt

# extract sequences
seqkit grep -f lysin_ids.txt ${INPUT_FASTA} > \
    ${RESULTS_DIR}/structure_lysin_candidates.faa


echo "Done!"
echo "Results saved to: ${RESULTS_DIR}/structure_lysin_candidates.faa"