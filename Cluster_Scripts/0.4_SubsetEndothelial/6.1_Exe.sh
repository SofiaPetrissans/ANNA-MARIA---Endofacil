#!/bin/bash
#SBATCH --job-name=Annotate_EC
#SBATCH --array=1-3
#SBATCH --output=logs/Annotate_EC_%A_%a.out
#SBATCH --error=logs/Annotate_EC_%A_%a.err
#SBATCH --time=00:20:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

module load Anaconda3
source activate R.Seurat

thresholds=("mt15")
resolutions=("0.1" "0.3" "0.5")
types=("Raw" "Integration")

# Cálculo de índices
TH_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) / 6 ))
RES_IDX=$(( ( (SLURM_ARRAY_TASK_ID - 1) / 2 ) % 3 ))
TYPE_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 ))

THRESHOLD=${thresholds[$TH_IDX]}
RESOLUTION=${resolutions[$RES_IDX]}
TYPE=${types[$TYPE_IDX]}

echo "Procesando THRESHOLD=${THRESHOLD}, RESOLUTION=${RESOLUTION}, TYPE=${TYPE}"

Rscript 6.0_Markers_Annotation_EC.R $THRESHOLD $RESOLUTION $TYPE

