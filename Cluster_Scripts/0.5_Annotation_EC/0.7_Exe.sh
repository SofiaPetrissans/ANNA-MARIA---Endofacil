#!/bin/bash
#SBATCH --job-name=DEGs
#SBATCH --array=1
#SBATCH --output=DEGs_%j.log   
#SBATCH --time=00:15:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat

Rscript 0.7_DEGs_EC_phenotype.R 


