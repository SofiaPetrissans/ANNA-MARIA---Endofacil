################################################################################

# SCRIPT: Scrublet
# AUTHOR: SOFIA PETRISSANS MOLL
# DATE: 10-06-2025

################################################################################


# SUMMARY SCRIPT:
# Qué hace
  # Carga el objeto Seurat procesado
  # Este script exporta por cada Sample:
     # La matriz de conteo .mtx
     # Los genes (genes.tsv)
     # Los barcodes (barcodes.tsv)
     # Un Excel con las dimensiones de PCA estimadas por Sample


# LIBRERÍAS .........................................................................
source('/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(Seurat)
library(Matrix)
library(tidyverse)
library(writexl)


# ...............................................................................
# ARGUMENTOS 
# ...............................................................................

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Uso: Rscript 01_Scrublet_ExportMatrix.R <umbral_mt>")
}

umbral_mt <- args[1]

# ...............................................................................
# PATHS
# ...............................................................................
main_path_objeto <- paths_AnnaMaria$Path_seurat_object
input_path <- file.path(main_path_objeto, "4.DoubletFinder", umbral_mt, "Seurat_DoubletFinder.rds")


# Crear carpeta específica para cada threshold
path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, "5.Scrublet", umbral_mt, "Input")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ...............................................................................
# FUNCIONES AUXILIARES 
# ...............................................................................

# Función para estimar dimensiones de PCA
EstimatePCA.Dims <- function(data){
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  return(min(co1, co2))
}


# ...............................................................................
# CARGAR OBJETO 
# ...............................................................................
seurat_object <- readRDS(input_path)


# ...............................................................................
# EXTRAER LA MATRIZ 
# ...............................................................................

samples <- unique(seurat_object$Sample)
dims <- data.frame(Sample = samples, PCA = NA)

for (s in samples) {
  sub <- subset(seurat_object, subset = Sample == s)
  sub <- NormalizeData(sub) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  dim_val <- EstimatePCA.Dims(sub)
  dims$PCA[dims$Sample == s] <- dim_val

  counts <- GetAssayData(sub, slot = "counts")
  Matrix::writeMM(counts, file = file.path(output_dir, paste0(s, ".mtx")))
  write.table(rownames(counts), file = file.path(output_dir, paste0(s, "_genes.tsv")),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(colnames(counts), file = file.path(output_dir, paste0(s, "_barcodes.tsv")),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}

write_xlsx(dims, file.path(output_dir, "PCA_dimensions.xlsx"))
