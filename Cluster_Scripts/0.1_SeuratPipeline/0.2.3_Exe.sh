#!/bin/bash

#SBATCH --job-name=0.1.1            # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=spetrissans@carrerasresearch.org   # Where to send mail
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=6gb                   # Job memory request
#SBATCH --time=00-00:05                # Time limit hrs:min:sec
#SBATCH --output=0.2.3_%j.log         # Standard output and error log

# AUTOR: SOFIA PETRISSANS MOLL
# FECHA: 06-06-2025

################################################################################################

module load Anaconda3
source activate R.Seurat

Rscript 0.2.2_Summary_Thresholds.R
echo "Pipeline Finished"