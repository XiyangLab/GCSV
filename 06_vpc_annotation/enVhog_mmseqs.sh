#!/bin/bash

# Initialize Conda environment

conda activate enVhog

# Define file paths and database names
INPUT_FILE="./vPCs_new.faa"
QUERY_DB="./query_db/query_db"
PROFILE_DB="enVhog/envhog_mmseqs/envhog_mmseqs_profiles"
RESULTS_DIR="./query_db_envhog_results/query_db_envhog_results"
OUTPUT_FILE="EnVhog_mmseqs.tsv"
TMP_DIR="./tmp"

# Ensure the temporary directory exists
mkdir -p "$TMP_DIR"
mkdir -p ./query_db
mkdir -p ./query_db_envhog_results

# --- Step 1: Create MMseqs2 query database ---
echo "Creating MMseqs query database..."
mmseqs createdb "$INPUT_FILE" "$QUERY_DB"
if [ $? -ne 0 ]; then
  echo "Error: Failed to create query database"
  exit 1
fi

# --- Step 2: Search against enVhog profile database ---
echo "Running MMseqs2 search..."
mmseqs search "$QUERY_DB" "$PROFILE_DB" "$RESULTS_DIR" "$TMP_DIR" \
	--threads 36 \
	--start-sens 2 \
	-s 7 \
	--sens-steps 3 \
	-a
if [ $? -ne 0 ]; then
  echo "Error: MMseqs2 search failed"
  exit 1
fi

# --- Step 3: Convert results to tabular format ---
echo "Generating m8 output..."
mmseqs convertalis "$QUERY_DB" "$PROFILE_DB" "$RESULTS_DIR" resultDB.m8 \
	--format-output "query,target,evalue,pident,fident,qstart,qend,qlen,tstart,tend,tlen,alnlen,qcov,tcov"
if [ $? -ne 0 ]; then
  echo "Error: Failed to create TSV output"
  exit 1
fi

# --- Step 4: Extract best hit per query and merge with enVhog annotations ---
# echo "Cleaning up temporary files..."
# rm -rf "$TMP_DIR"

python mmseqs_best_hits.py \
    -i resultDB.m8 \
    -e enVhog/enVhog_tables/envhog_infos.tsv \
    -o best1.tsv \
    --min_identity 30 \
    --max_evalue 1E-5 \
    --min_qcov 0.8

echo "Process completed successfully!"

