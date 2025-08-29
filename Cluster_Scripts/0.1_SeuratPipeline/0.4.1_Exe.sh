#!/bin/bash

#SBATCH --job-name=DoubletFinder
#SBATCH --output=logs/DoubletFinder_%A_%a.out
#SBATCH --error=logs/DoubletFinder_%A_%a.err
#SBATCH --array=1-8                         # ← Cambia si hay más/menos thresholds
#SBATCH --time=00:01:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org


module load Anaconda3
source activate R.Seurat

# Lista de thresholds
thresholds=("mt5" "mt10" "mt15" "mt20" "mt25" "mt30" "mt35" "mt40")

# Selecciona el umbral correspondiente a esta tarea del array
THRESHOLD=${thresholds[$SLURM_ARRAY_TASK_ID-1]}

# Ruta base de resultados 
RESULTS_BASE="/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/0.1_SeuratPipeline/3.SeuratProcessed"

# Ejecuta el script de DoubletFinder
Rscript 0.4_DoubletFinder.R $THRESHOLD $RESULTS_BASE
