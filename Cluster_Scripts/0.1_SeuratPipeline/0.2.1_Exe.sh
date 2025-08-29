#!/bin/bash

#SBATCH --job-name=0.2.1            # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=spetrissans@carrerasresearch.org   # Where to send mail
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=12gb                   # Job memory request
#SBATCH --time=00-00:30                # Time limit hrs:min:sec
#SBATCH --output=0.2.1_FilteringQC_%A_%a.log       # Standard output and error log
#SBATCH --array=1-8

# AUTOR: SOFIA PETRISSANS MOLL
# FECHA: 06-06-2025

################################################################################################

module load Anaconda3
source activate R.Seurat

# Llamar al script R con el índice del array
Rscript 0.2_Filtering.R $SLURM_ARRAY_TASK_ID

echo "Pipeline Finished"