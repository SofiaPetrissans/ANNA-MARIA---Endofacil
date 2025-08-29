#!/bin/bash
#SBATCH --job-name=Traicing
#SBATCH --output=Traicing_%A_%a.out
#SBATCH --error=Traicing_%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org


module load Anaconda3
source activate R.Seurat

# ---------------------------------------------------------------------

# Ejecuta el script R con los valores definidos
Rscript 0.0_Traicing.R 
#Rscript prueba.R 