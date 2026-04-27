#!/bin/bash
#SBATCH -J DeepMineLys
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:1
#SBATCH -p normal
#SBATCH -o %j.out
#SBATCH -e %j.err

# ===============================
# Load environment
# ===============================
module purge
module load singularity

INPUT_FILE=$1
BASENAME=$(basename ${INPUT_FILE} .faa)

WORKDIR="deepminelys_work"
OUTDIR="${WORKDIR}/results"
mkdir -p ${OUTDIR}

# ===============================
# Step 1: TAPE embedding
# ===============================
singularity exec <container_path>.sif bash -c "
source /opt/env.sh &&
conda activate DeepMineLys &&
tape-embed unirep ${INPUT_FILE} ${OUTDIR}/${BASENAME}.npz babbler-1900 --tokenizer unirep
"

# ===============================
# Step 2: Read npz
# ===============================
singularity exec <container_path>.sif bash -c "
conda activate DeepMineLys &&
python scripts/read_npz.py ${OUTDIR}/${BASENAME}
"

# ===============================
# Step 3: PHY embedding
# ===============================
singularity exec <container_path>.sif bash -c "
conda activate phy_embedding &&
Rscript scripts/phy_embedding.R ${OUTDIR}/${BASENAME}
"

# ===============================
# Step 4: Combine embeddings
# ===============================
singularity exec <container_path>.sif bash -c "
conda activate DeepMineLys &&
python scripts/combine_embeddings.py ${OUTDIR}/${BASENAME}
"

# ===============================
# Step 5: Prediction
# ===============================
singularity exec <container_path>.sif bash -c "
conda activate DeepMineLys &&
python scripts/predict_model.py \
    --input_path ${OUTDIR}/${BASENAME}_tape_phy.csv \
    --model_path model/ \
    --output_path ${OUTDIR}/${BASENAME}_prediction.csv
"