#!/bin/bash

#SBATCH --job-name=Zonation             # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=spetrissans@carrerasresearch.org   # Where to send mail
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=12gb                   # Job memory request
#SBATCH --time=00:30:00                # Time limit hrs:min:sec
#SBATCH --output=logs/Zonation_%j.log         # Standard output and error log

# AUTOR: SOFIA PETRISSANS MOLL
# FECHA: 05-06-2025

################################################################################################

module load Anaconda3
source activate R.Seurat.V2

Rscript 0.2_Zonation_Functions.R
echo "Pipeline Finished"