#!/bin/bash
#SBATCH --job-name=markers
#SBATCH --array=1
#SBATCH --output=markers_%j.log   
#SBATCH --time=00:15:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat

Rscript 0.9_Annotation_mt15_copy.R 


