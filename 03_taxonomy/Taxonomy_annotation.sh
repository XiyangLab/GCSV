#!/bin/bash
# =============================================================================
# Viral Taxonomic Annotation Pipeline
# =============================================================================
# Description: Assigns taxonomy to viral sequences (vOTUs) using three
#              complementary approaches: geNomad, vContact3, and VITAP. 
#              Results are merged in a companion R script.
#
# Tools required:
#   - seqkit    : splitting large FASTA files for parallel processing
#   - geNomad   : taxonomy assignment via viral marker genes
#   - vContact3 : genus-level clustering against reference virus database
#   - VITAP     : ICTV VMR-based taxonomic assignment
# =============================================================================


# -----------------------------------------------------------------------------
# Step 1: geNomad
# -----------------------------------------------------------------------------
# Here it is re-run on sequences that lacked a geNomad annotation in the
# viral identification step (e.g., sequences identified only by VIBRANT
# or VirSorter2).

genomad end-to-end --cleanup \
    coldseeps_viruses_clusters_genomad_NA.fna \
    genomad_taxonomy_output \
    /path/to/genomad_db 


# -----------------------------------------------------------------------------
# Step 2: vContact3
# -----------------------------------------------------------------------------

# Split vOTU sequences into smaller subsets for parallel processing
seqkit split2 -p 5 coldseeps_vOTU.fna

# Run vContact3 on each subset (repeat for all parts)
for i in 001 002 003 004 005; do
    vcontact3 run \
        --db-version 223 \
        --db-path /path/to/vcontact3-db \
        --db-domain "prokaryotes" \
        --nucleotide coldseeps_vOTU.part_${i}.fna \
        --output vcontact3_output_${i}
done


# -----------------------------------------------------------------------------
# Step 3: VITAP
# -----------------------------------------------------------------------------

# (Optional) Update the VITAP database when a new ICTV VMR release is available
VITAP upd \
    --vmr ICTV_VMR-MSL40.csv \
    -o ICTV_VMR-MSL40_reformat.csv \
    -d VMR-MSL40

mkdir -p VITAP_output

# Run VITAP on each subset (repeat for all parts)
for i in 001 002 003 004 005; do
    VITAP assignment \
        -i coldseeps_vOTU.part_${i}.fna \
        -d /path/to/VITAP/DB/DB_VMR-MSL40 \
        -o VITAP_output/coldseeps_vOTU_${i} \
        -p 35
done
