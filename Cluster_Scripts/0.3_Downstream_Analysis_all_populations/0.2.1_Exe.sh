#!/bin/bash
#SBATCH --job-name=Senescence
#SBATCH --array=1
#SBATCH --output=Senescence_%j.log   
#SBATCH --time=00:10:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat.V2

Rscript 0.2_Markers_Senescencia.R 


