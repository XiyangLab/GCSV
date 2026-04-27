#!/bin/bash
# =============================================================================
# Viral Host Prediction Pipeline
# =============================================================================
# Description: Predicts hosts for vOTUs using two complementary approaches:
#              (1) iPHoP, (2) direct contig-level matching of proviruses to
#              MAGs based on sequence coordinates.
#
# Tools required:
#   - GTDBTk  : phylogenetic classification of MAGs for iPHoP database
#   - iPHoP   : integrated phage host prediction
#   - MArVD2  : machine learning prediction of archaeal virus probability
# =============================================================================


# -----------------------------------------------------------------------------
# Step 1: Build Custom iPHoP Database with Cold Seep MAGs
# -----------------------------------------------------------------------------
# iPHoP's default database can be extended with in-house MAGs to improve
# host prediction for under-represented environments.
# GTDBTk is first used to classify the MAGs, then iPHoP ingests both
# the genomes and their taxonomy to expand the reference database.

# --- GTDBTk: classify MAGs (bacteria and archaea separately) ---
gtdbtk de_novo_wf \
    --genome_dir dereplicated_genomes \
    --bacteria \
    --outgroup_taxon p__Patescibacteria \
    --out_dir MAGs_GTDB-tk_results_bac \
    --cpus 32 \
    --force \
    --extension fa 

gtdbtk de_novo_wf \
    --genome_dir dereplicated_genomes \
    --archaea \
    --outgroup_taxon p__Altiarchaeota \
    --out_dir MAGs_GTDB-tk_results_arc \
    --cpus 32 \
    --force \
    --extension fa 

# --- iPHoP: add MAGs to the reference database ---
# Merge bacterial and archaeal GTDBTk results into a single directory first,
# then run add_to_db.
iphop add_to_db \
    --fna_dir  dereplicated_genomes \
    --gtdb_dir MAGs_GTDB-tk_results_all \
    --out_dir  iphop_custom_db \
    --db_dir   iphop_db \
    -t 20 

# -----------------------------------------------------------------------------
# Step 2: iPHoP Host Prediction
# -----------------------------------------------------------------------------
# The vOTU FASTA is split into subsets to allow parallel execution.
# Results from all subsets are merged in Step 3.

# Split vOTU sequences
seqkit split2 -p 60 coldseeps_vOTU.fna

# Run iPHoP on each subset (repeat for all parts)
for i in $(seq -w 1 60); do
    iphop predict \
        --fa_file  coldseeps_vOTU.part_${i}.fna \
        --db_dir   iphop_custom_db \
        --out_dir  iphop_results/vOTU_${i} \
        -t 30
done


# -----------------------------------------------------------------------------
# Step 3: Merge iPHoP Results
# -----------------------------------------------------------------------------
# Two output files are merged across all subsets:
#   Host_prediction_to_genome_m90.csv  : top host predictions (score >= 90)
#   Detailed_output_by_tool.csv        : per-tool detailed scores
#

RESULT_DIR="iphop_results"

# --- Merge Host_prediction_to_genome_m90.csv ---
OUTPUT_GENOME="Host_prediction_to_genome_m90_all.csv"
files_genome=$(find "${RESULT_DIR}" -type f -name "Host_prediction_to_genome_m90.csv" | sort)

if [ -n "${files_genome}" ]; then
    head -n 1 $(echo "${files_genome}" | head -n 1) > "${OUTPUT_GENOME}"
    for file in ${files_genome}; do
        tail -n +2 "${file}" >> "${OUTPUT_GENOME}"
    done
    echo "Merged: ${OUTPUT_GENOME}"
fi

# --- Merge Detailed_output_by_tool.csv ---
# Note: this file has an 11-line metadata header; data starts at line 12
OUTPUT_DETAIL="Detailed_output_by_tool_all.csv"
files_detail=$(find "${RESULT_DIR}" -type f -name "Detailed_output_by_tool.csv" | sort)

if [ -n "${files_detail}" ]; then
    # Extract column header (line 11) from first file
    sed -n '11p' $(echo "${files_detail}" | head -n 1) > "${OUTPUT_DETAIL}"
    for file in ${files_detail}; do
        tail -n +12 "${file}" >> "${OUTPUT_DETAIL}"
    done
    echo "Merged: ${OUTPUT_DETAIL}"
fi


# -----------------------------------------------------------------------------
# Step 4: Provirus-to-MAG Matching by Contig Coordinates
# -----------------------------------------------------------------------------
# For proviruses, host identity can be inferred directly by finding which MAG
# contains the contig from which the provirus was extracted.
# The matching script compares provirus contig IDs against MAG sequences.

conda activate anvio-8

python match_provirus_to_MAG.py \
    --virus   nr_provirus_all_final.fna \
    --mag_dir dereplicated_genomes \
    --output  nr_provirus_all_final_host_match_results.tsv


# -----------------------------------------------------------------------------
# Step 5: MArVD2 — Archaeal Virus Probability Prediction
# -----------------------------------------------------------------------------
# Reference: https://github.com/plasmidologist/MArVD2
# MArVD2 uses a random forest classifier to estimate the probability that
# each viral sequence is an archaeal virus, based on protein homology to
# pVOGs, NCBI nr, and marine reference databases.

MArVD2.py \
    -i coldseeps_vOTU.fna \
    --load-model marvd2/rf_model.pkl \
    -o MArVD2_output \
    --db-pvog           marvd2/AllvogHMMprofiles.hmm \
    --db-nr             marvd2/nr.faa \
    --marine-jackhmmer-db marvd2/pVOG_prots_ref_marine_pVOG.faa \
    --viral-refseq-txt  marvd2/viruses.txt \
    --pvog-dir          marvd2/ \
    --db-accession2tax  marvd2/prot.accession2taxid.trimmed \
    -c 32 