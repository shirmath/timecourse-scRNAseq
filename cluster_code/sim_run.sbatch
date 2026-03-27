#!/bin/bash

#SBATCH --job-name=mom_sims
#SBATCH --partition=short
#SBATCH --array=1-900
#SBATCH --ntasks=1 
#SBATCH --time=0-12:00
#SBATCH --mem-per-cpu=20000
#SBATCH --output=out/%A_%a.out
#SBATCH --error=error/%A_%a.err

module load R-bundle-CRAN
Rscript cluster_sim_script.R ${SLURM_ARRAY_TASK_ID} ${SLURM_JOB_NAME}
