#!/bin/bash
#SBATCH --job-name=AUCell_Panglao
#SBATCH --output=logs/AUCell_%A_%a.out
#SBATCH --error=logs/AUCell_%A_%a.err
#SBATCH --array=1-9
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

# --------------------------
# CARGAR ENTORNO
# --------------------------
module load Anaconda3
source activate R.Seurat.V2

# --------------------------
# UMBRALES Y RESOLUCIONES
# --------------------------

thresholds=("mt15" "mt20" "mt25")
resolutions=("0.1" "0.3" "0.5")

# --------------------------
# COMBINACIÓN ARRAY_ID → THRESHOLD + RES
# --------------------------
THRESHOLD_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / 3 ))
RESOLUTION_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) % 3 ))

THRESHOLD=${thresholds[$THRESHOLD_IDX]}
RESOLUTION=${resolutions[$RESOLUTION_IDX]}

echo "Ejecutando AUCell para $THRESHOLD con resolución $RESOLUTION"

# --------------------------
# EJECUTAR SCRIPT R
# --------------------------
Rscript 0.1.C_Annotation_AUCell_Panglao.R $THRESHOLD $RESOLUTION
