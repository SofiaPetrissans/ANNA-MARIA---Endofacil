#!/bin/bash
#SBATCH --job-name=SeuratPipeline
#SBATCH --output=logs/SeuratPipeline_%A_%a.out               
#SBATCH --error=logs/SeuratPipeline_%A_%a.err
#SBATCH --array=1-8                                       
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00

# ACTIVAR ENTORNO
module load Anaconda3
source activate R.Seurat

# DEFINIR THRESHOLDS
thresholds=(mt5 mt10 mt15 mt20 mt25 mt30 mt35 mt40)

# SELECCIONAR THRESHOLD SEGÚN JOB ARRAY
umbral=${thresholds[$SLURM_ARRAY_TASK_ID-1]}

# EJECUTAR SCRIPT
Rscript 0.3_Seurat_Pipeline.R $umbral 


