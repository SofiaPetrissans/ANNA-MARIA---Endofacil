#!/bin/bash
#SBATCH --job-name=Annotation
#SBATCH --array=1
#SBATCH --output=Annotation_%j.log   
#SBATCH --time=00:10:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat.V2

Rscript 0.0_EvaluateClustering.R 


