#!/bin/bash
#SBATCH --job-name=Annotate_Int
#SBATCH --array=1-21
#SBATCH --output=logs/Annotate_Int_%A_%a.out
#SBATCH --error=logs/Annotate_Int_%A_%a.err
#SBATCH --time=40:00:00               # Time limit hrs:min:sec
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat

thresholds=("mt10" "mt15" "mt20" "mt25" "mt30" "mt35" "mt40")
resolutions=("0.1" "0.3" "0.5")

# Obtener combinación de threshold y resolución según SLURM_ARRAY_TASK_ID
TH_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) / 3 ))
RES_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) % 3 ))

THRESHOLD=${thresholds[$TH_IDX]}
RESOLUTION=${resolutions[$RES_IDX]}

echo "Procesando THRESHOLD=${THRESHOLD}, RESOLUTION=${RESOLUTION}"

Rscript 0.2.A_Markers_Annotation_Integration.R $THRESHOLD $RESOLUTION
