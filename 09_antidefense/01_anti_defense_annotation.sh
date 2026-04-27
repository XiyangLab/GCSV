#!/bin/bash
# =============================================================================
# Anti-Defense Gene Annotation Pipeline
# =============================================================================
# Description: Identifies viral anti-defense genes and host defense systems
#              using three approaches:
#              (1) DefenseFinder — predicts anti-defense and defense systems in 
#                  viral proteins, and defense systems in host MAG proteins;
#              (2) dbAPIS — searches viral proteins against the anti-phage
#                  immune system database using both HMM and diamond blastp.
#
#              Results are merged and filtered in a companion R script.
#
# Tools required:
#   - DefenseFinder  : defense and anti-defense system prediction
#   - Prodigal       : ORF prediction for host genomes
#   - HMMER hmmscan  : profile HMM search against dbAPIS
#   - diamond blastp : sequence similarity search against dbAPIS
# =============================================================================


# =============================================================================
# Part 1: DefenseFinder — Viral Anti-Defense Gene Prediction
# =============================================================================
# Reference: https://github.com/mdmparis/defense-finder

conda activate defensefinder

defense-finder run \
    -a \
    --db-type gembase \
    -o ADF_DF_vOTU \
    -w 35 \
    coldseeps_vOTU_protein.faa


# =============================================================================
# Part 2: DefenseFinder — Host Defense System Prediction
# =============================================================================
# Note: input must be complete host genome sequences.

# --- Step 2a: Collect host genome sequences ---
# Copy dereplicated MAG FASTA files listed in host_id.txt to a local directory.
mkdir -p host_genome

cat host_id.txt | xargs -I {} cp \
    /path/to/dereplicated_genomes/{}.fa \
    host_genome/

# Concatenate all host genomes into a single file for Prodigal
cat host_genome/*.fa > host_genome_all.fna

# --- Step 2b: ORF prediction with Prodigal ---
mkdir -p host_genome_prodigal

prodigal \
    -i host_genome_all.fna \
    -o host_genome_prodigal/host_genome_gene.gff \
    -d host_genome_prodigal/host_genome_gene.fna \
    -a host_genome_prodigal/host_genome_protein.faa \
    -p meta

# --- Step 2c: DefenseFinder on host proteins ---
defense-finder run \
    -a \
    --db-type gembase \
    -o host_defensefinder \
    -w 40 \
    host_genome_prodigal/host_genome_protein.faa


# =============================================================================
# Part 3: dbAPIS — Anti-Phage Immune System Database Search
# =============================================================================
# Reference: https://bcb.unl.edu/dbAPIS/
# Two complementary search methods are used: hmmscan and diamond
# Results from both methods are combined with the provided parse script,
# then further processed in the companion R script.

mkdir -p dbAPIS

# --- Step 3a: Download and press the dbAPIS HMM database ---
wget https://bcb.unl.edu/dbAPIS/downloads/dbAPIS.hmm -P dbAPIS/
hmmpress dbAPIS/dbAPIS.hmm

# --- Step 3b: hmmscan search ---
hmmscan \
    -E 1e-10 \
    --domtblout dbAPIS/hmmscan_dbAPIS.out \
    --noali \
    dbAPIS/dbAPIS.hmm \
    mmseqs2/mmseqs_vOTU_id0.5_c0.9_clusters_rep_seq.fasta

# --- Step 3c: Build diamond database and run blastp search ---
diamond makedb \
    --in anti_defense.pep \
    -d dbAPIS/APIS_db

diamond blastp \
    --db dbAPIS/APIS_db \
    -e 1e-10 \
    --id 30 \
    -q mmseqs2/mmseqs_vOTU_id0.5_c0.9_clusters_rep_seq.fasta \
    -f 6 qseqid sseqid pident length mismatch gapopen \
          qstart qend sstart send evalue bitscore qlen slen \
    -o dbAPIS/diamond_dbAPIS.out \
    --max-target-seqs 10000

# --- Step 3d: Parse and combine HMM and diamond results ---
# The parse script integrates both result files into a unified annotation table.
# Download from: https://bcb.unl.edu/dbAPIS/
bash parse_annotation_result.sh \
    dbAPIS/hmmscan_dbAPIS.out \
    dbAPIS/diamond_dbAPIS.out
