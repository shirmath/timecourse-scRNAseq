#!/bin/bash

#SBATCH --job-name=vi_sim_test
#SBATCH --partition=medium
#SBATCH --array=1
#SBATCH --ntasks=1 
#SBATCH --time=0-12:00
#SBATCH --mem-per-cpu=30000
#SBATCH --output=out/%A_%a.out
#SBATCH --error=error/%A_%a.err

module load R-bundle-CRAN
Rscript cluster_vi_sim_script.R ${SLURM_ARRAY_TASK_ID} ${SLURM_JOB_NAME}
