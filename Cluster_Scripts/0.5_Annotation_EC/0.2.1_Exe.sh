#!/bin/bash
#SBATCH --job-name=Markers
#SBATCH --array=1
#SBATCH --output=logs/Markers_%A_%a.out
#SBATCH --error=logs/Markers_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat

Rscript 0.2_Markers_Endotelio.R
