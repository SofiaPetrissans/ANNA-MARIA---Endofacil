#!/bin/bash
#SBATCH --job-name=HarmonyPipeline
#SBATCH --output=logs/HarmonyPipeline_%A_%a.out      
#SBATCH --error=logs/HarmonyPipeline_%A_%a.err
#SBATCH --array=1-4        
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# ACTIVAR ENTORNO
module load Anaconda3
source activate R.Seurat

thresholds=(mt25 mt30 mt35 mt40)
umbral=${thresholds[$SLURM_ARRAY_TASK_ID-1]}

Rscript 0.6_Integration.R $umbral


