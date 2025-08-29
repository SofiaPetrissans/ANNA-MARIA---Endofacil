#!/bin/bash
#SBATCH --job-name=Propeller
#SBATCH --array=1
#SBATCH --output=logs/Propeller_%A_%a.out
#SBATCH --error=logs/Propeller_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat.V2

Rscript 5.0_Propeller.R mt20 Harmony_Log_res.0.5


