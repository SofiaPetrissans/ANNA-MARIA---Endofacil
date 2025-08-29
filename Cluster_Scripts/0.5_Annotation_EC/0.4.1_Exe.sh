#!/bin/bash

#SBATCH --job-name=ZonationScore             # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=spetrissans@carrerasresearch.org   # Where to send mail
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=5gb                   # Job memory request
#SBATCH --time=00:30:00                # Time limit hrs:min:sec
#SBATCH --output=logs/ZonationScore_%j.log         # Standard output and error log

# AUTOR: SOFIA PETRISSANS MOLL
# FECHA: 01-07-2025

################################################################################################

module load Anaconda3
source activate R.Seurat.V2

Rscript 0.4_Zonation_Score.R
echo "Pipeline Finished"