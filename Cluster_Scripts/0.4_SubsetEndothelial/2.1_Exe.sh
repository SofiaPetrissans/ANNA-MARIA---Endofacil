#!/bin/bash
#SBATCH --job-name=Annotate
#SBATCH --array=1-6
#SBATCH --output=logs/Annotate_%A_%a.out
#SBATCH --error=logs/Annotate_%A_%a.err
#SBATCH --time=01:00:00
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

Rscript 2.0_Markers_Annotation.R $THRESHOLD $RESOLUTION $TYPE

