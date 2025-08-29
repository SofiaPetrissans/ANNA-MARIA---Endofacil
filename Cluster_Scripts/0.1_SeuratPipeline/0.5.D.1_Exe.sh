#!/bin/bash
#SBATCH --job-name=Scrublet
#SBATCH --array=1-8
#SBATCH --output=logs/Scrublet_%A_%a.out
#SBATCH --error=logs/Scrublet_%A_%a.err
#SBATCH --time=00:00:30
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=spetrissans@carrerasresearch.org

# Definir thresholds y rutas
thresholds=("mt5" "mt10" "mt15" "mt20" "mt25" "mt30" "mt35" "mt40")
THRESHOLD=${thresholds[$SLURM_ARRAY_TASK_ID-1]}
RESULTS_BASE="/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.1_SeuratPipeline/4.DoubletFinder"     # MODFICAR !!!
SCRUBLET_BASE="/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA//3.MERGE/0.1_SeuratPipeline/5.Scrublet"        # MODFICAR !!!

# --------------------------
# Paso 1: Exportar matrices desde Seurat
# --------------------------
module load Anaconda3
source activate R.Seurat

echo "Exportando matriz para $THRESHOLD..."
Rscript 0.5.A_Scrublet_ExportMatrix.R $THRESHOLD $RESULTS_BASE

# --------------------------
# Paso 2: Ejecutar Scrublet en Python
# --------------------------
module unload Anaconda3

module load python
source ~/.venvs/scrublet_env/bin/activate

SCRUBLET_FOLDER="${SCRUBLET_BASE}/${THRESHOLD}"
INPUT_FOLDER="${SCRUBLET_FOLDER}/Input"

echo "Ejecutando Scrublet en Python para $THRESHOLD..."

python 0.5.B_Scrublet.py $INPUT_FOLDER

# --------------------------
# Paso 3: Importar resultados en Seurat
# --------------------------
module unload python

module load Anaconda3
source activate R.Seurat

echo "Importando resultados Scrublet a Seurat para $THRESHOLD..."
Rscript 0.5.C_Scrublet_ImportToSeurat.R $THRESHOLD $RESULTS_BASE
