#!/bin/bash
# =============================================================================
# Auxiliary Metabolic Gene (AMG) Annotation Pipeline
# =============================================================================
# Description: Identifies and annotates AMGs in viral sequences using three
#              complementary approaches:
#              (1) DRAM-v — metabolic annotation with downstream filtering
#                  for tRNA, inverted repeats, and suspicious genes;
#              (2) VIBRANT — independent AMG calling with eggNOG validation;
#              (3) Hydrocarbon degradation gene annotation (CANT-HYD + HMDB)
#
#              Final AMG calls from all tools are merged and reconciled in a
#              companion R script.
#
# Tools required:
#   - VirSorter2    : prepares affiliation table required by DRAM-v
#   - DRAM-v        : metabolic annotation and AMG distillation
#   - tRNAscan-SE   : identifies tRNA genes (used to filter false AMGs)
#   - EMBOSS einverted : identifies inverted repeats (used to filter false AMGs)
#   - VIBRANT       : independent AMG annotation
#   - eggNOG-mapper : validates VIBRANT AMG candidates
#   - HMMER hmmsearch : searches CANT-HYD HMM profiles
#   - BLAST blastp  : searches HMDB hydrocarbon degradation database
# =============================================================================


# =============================================================================
# Part 1: DRAM-v
# =============================================================================
# Reference: https://github.com/WrightonLabCSU/DRAM
# DRAM-v requires a VirSorter2-generated affiliation table as input.
# Due to runtime constraints, split the input FASTA into subsets.

# --- Step 1a: Prepare VirSorter2 affiliation table for DRAM-v ---
conda activate vs2

virsorter run \
    --prep-for-dramv \
    -w for_dramv_out \
    -i coldseeps_vOTU.fna \
    -j 30 all

# --- Step 1b: DRAM-v annotation (run per subset) ---
conda activate DRAM

DRAM-v.py annotate \
    -i coldseeps_vOTU.part_001.fna \
    -v for_dramv_out/for-dramv/viral-affi-contigs-for-dramv.tab \
    -o annotation_DRAM-v_001 \
    --threads 40

# --- Step 1c: DRAM-v distillation (per subset) ---
DRAM-v.py distill \
    -i annotation_DRAM-v_001/annotations.tsv \
    -o annotation_DRAM-v_001/distilled

# --- Step 1d: Merge AMG summary tables across all subsets ---
head -n 1 annotation_DRAM-v_001/distilled/amg_summary.tsv > amg_summary_all.tsv
tail -n +2 -q annotation_DRAM-v_*/distilled/amg_summary.tsv >> amg_summary_all.tsv


# =============================================================================
# Part 2: DRAM-v Result Filtering
# =============================================================================
# DRAM-v results require additional filtering following the workflow described
# in Chen et al. (2022) "Prokaryotic-virus-encoded auxiliary metabolic genes
# throughout the global oceans" (DOI: 10.1038/s41467-022-29parallels).

# --- Step 2a: tRNA detection ---
conda activate DRAM

tRNAscan-SE \
    -G \
    -o trna_scaffolds.tsv \
    -m trna_scaffolds.stat \
    scaffolds_all.fna

# --- Step 2b: Inverted repeat detection ---
conda activate emboss

einverted \
    scaffolds_all.fna \
    -gap 12 \
    -threshold 30 \
    -match 3 \
    -mismatch -4 \
    -outfile scaffolds_ir.inv \
    -outseq scaffolds_ir.fasta \
    -maxrepeat 3000

# Extract and clean scaffold IDs from inverted repeat output
grep '>' scaffolds_ir.fasta \
    | sed 's/>//' \
    | sed 's/_[0-9]*_[0-9]*$//' \
    | awk '!seen[$0]++' \
    > new_id_ir_no_duplicates.txt

# --- Step 2c: Filter suspicious genes ---
# Download from:
#   https://bitbucket.org/MAVERICLab/virsorter2-sop/raw/03b8f28bee979e2b7fd99d7375d915c29c938339/resource/suspicious-gene.list
grep -F -f suspicious-gene.list annotations.tsv \
    > annotations_suspicious-gene.tsv


# =============================================================================
# Part 3: VIBRANT AMG Annotation
# =============================================================================
# Reference: https://github.com/AnantharamanLab/VIBRANT
# VIBRANT is re-run on vOTUs to independently identify AMG candidates.
# Filtering follows the three-step approach described in:
#   Aylward et al. (2024) Nature Microbiology (GitHub README):
#   https://github.com/AnantharamanLab/TYMEFLIES_Viral/blob/main/Reconstruct_vMAGs/README.md
# Steps 1–3 of that filtering workflow are implemented in the companion R script.
# Step 4 (eggNOG validation) is run here.

# --- Step 3a: Run VIBRANT on vOTUs ---
VIBRANT_run.py \
    -i coldseeps_vOTU.fna \
    -t 30 \
    -folder VIBRANT_output

# --- Step 3b: Extract AMG protein sequences (after R-based filtering, steps 1–3) ---
awk -F '\t' '{print $1}' VIBRANT_AMG.tsv > VIBRANT_AMG_id.txt

seqkit grep -n \
    -f VIBRANT_AMG_id.txt \
    coldseeps_vOTU.phages_combined.faa \
    -o VIBRANT_AMG.faa

# --- Step 3c: eggNOG-mapper validation (step 4 of filtering) ---
conda activate eggnog-mapper

emapper.py \
    -i VIBRANT_AMG.faa \
    --itype proteins \
    --output VIBRANT_AMG_eggnog \
    --data_dir /path/to/eggnog-mapper \
    -m diamond \
    --cpu 30


# =============================================================================
# Part 4: Hydrocarbon Degradation Gene Annotation
# =============================================================================
# Two databases are searched in parallel; results are merged
# and cross-validated against the DRAM-v filtering criteria.

# --- Step 4a: CANT-HYD HMM search ---
# Reference: https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation
mkdir -p CANT-HYD

hmmsearch \
    --tblout CANT-HYD/hmmsearch_CANT-HYD_001.tblout \
    --cut_nc \
    --cpu 20 \
    -o CANT-HYD/hmmsearch_CANT-HYD_001.out \
    /path/to/CANT-HYD.hmm \
    vPCs.part_001.fasta

# --- Step 4b: HMDB blastp search ---
# Reference: https://github.com/BushmanLab/HMDB

# Build BLAST database (run once)
makeblastdb \
    -in HMDB_v1.1.faa \
    -input_type fasta \
    -dbtype prot \
    -out HMDB_v1.1/HMDB_v1.1

# Search vOTU proteins against HMDB
blastp \
    -db HMDB_v1.1/HMDB_v1.1 \
    -query genes_all.faa \
    -evalue 1e-10 \
    -num_threads 20 \
    -outfmt 6 \
    -out coldseep_vOTU_proteins_vs_HMDB_v1.1.bsp

# Add metadata annotations using the BLAST tabular assistant
perl BLAST_tabular_assistant.pl \
    HMDB_v1.1.faa \
    coldseep_vOTU_proteins_vs_HMDB_v1.1.bsp \
    -meta \
    > coldseep_vOTU_proteins_vs_HMDB_v1.1.txt
