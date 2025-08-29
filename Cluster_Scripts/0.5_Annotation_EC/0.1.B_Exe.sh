#!/bin/bash
#SBATCH --job-name=Markers_Division
#SBATCH --array=1
#SBATCH --output=Markers_Division_%A_%a.out
#SBATCH --error=Markers_Division_%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat

Rscript 0.1.B_Markers_Division.R 

