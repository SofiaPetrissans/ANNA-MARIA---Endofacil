#!/bin/bash
#SBATCH --job-name=SubsetEndothelial
#SBATCH --output=logs/SubsetEndothelial_%A_%a.out
#SBATCH --error=logs/SubsetEndothelial_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org


module load Anaconda3
source activate R.Seurat


# ---------------------------------------------------------------------
# defines directamente los valores de los argumentos
#MT_THRESHOLD="mt15"
#RESOLUTION="RNA_snn_res.0.3"
#CLUSTERS_TO_KEEP="1,2,11,12,15,16,17"

MT_THRESHOLD="mt15"
RESOLUTION="Annotation_Layer1"
CLUSTERS_TO_KEEP="Endothelial"

#MT_THRESHOLD="mt20"
#RESOLUTION="RNA_snn_res.0.3"
#CLUSTERS_TO_KEEP="1,2,9,12,13,16,17"

#MT_THRESHOLD="mt25"
#RESOLUTION="RNA_snn_res.0.3"
#CLUSTERS_TO_KEEP="1,2,9,11,12,15,16,18"


# ---------------------------------------------------------------------

# Ejecuta el script R con los valores definidos
Rscript 1.0_Subset_Seurat.R $MT_THRESHOLD $RESOLUTION $CLUSTERS_TO_KEEP
