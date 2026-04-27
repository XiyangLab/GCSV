#!/bin/bash

# Description:
#   Submit DeepMineLys jobs using SLURM

INPUT_LIST=$1

while read file; do
    echo "Submitting job for ${file}"
    sbatch scripts/deepminelys.slurm ${file}
    sleep 2
done < ${INPUT_LIST}