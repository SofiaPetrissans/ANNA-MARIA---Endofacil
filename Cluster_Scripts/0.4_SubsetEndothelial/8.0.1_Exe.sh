#!/bin/bash
#SBATCH --job-name=RemoveCluster
#SBATCH --output=logs/RemoveCluster_%A_%a.out
#SBATCH --error=logs/RemoveCluster_%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org


module load Anaconda3
source activate R.Seurat


# ---------------------------------------------------------------------
# defines directamente los valores de los argumentos
#MT_THRESHOLD="mt20"
#RESOLUTION="Harmony_Log_res.0.5"
#CLUSTERS_TO_KEEP="0,1,2,3,4,5,7,8"
#CLUSTERS_TO_KEEP="0,1,2,3,4,5,7,8,9"

MT_THRESHOLD="mt15"
RESOLUTION="Harmony_Log_res.0.3"
CLUSTERS_TO_KEEP="0,1,2,3,4"
# ---------------------------------------------------------------------

# Ejecuta el script R con los valores definidos
Rscript 8.0_RemoveCluster.R $MT_THRESHOLD $RESOLUTION $CLUSTERS_TO_KEEP
