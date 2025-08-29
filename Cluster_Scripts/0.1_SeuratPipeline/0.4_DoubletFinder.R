################################################################################

# SCRIPT: Doublet Finder 
# AUTHOR: SOFIA PETRISSANS MOLL
# DATE: 09-06-2025

################################################################################


# SUMMARY SCRIPT:
# Entradas
  # Argumento 1: umbral_mt (por ejemplo mt10, mt20)
  # Argumento 2: res_path, carpeta general que contiene subcarpetas por umbral

# Qué hace
  # Carga el objeto Seurat con PCA y clustering aplicados
  # Divide por Sample y aplica doubletFinder_v3() por muestra
  # Usa una estimación de 7.5% de doublets esperados
  # Añade los resultados al @meta.data principal
  # Guarda el nuevo objeto con dobletes anotados
  # Genera un DimPlot por muestra coloreado por clasificación de DoubletFinder




# LIBRERÍAS .........................................................................
source('/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(RColorBrewer)


# ...............................................................................
# FUNCIONES AUXILIARES 
# ...............................................................................

getPalette <- colorRampPalette(brewer.pal(9, "Set3"))

EstimatePCA.Dims <- function(data) {
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  dim.final <- min(co1, co2)
  return(dim.final)
}

DoubletFinderPipeline <- function(objeto, col_annotation, doublet_rate = 0.076){
  dim.max <- EstimatePCA.Dims(objeto)
  sweep.res <- paramSweep(objeto, PCs = 1:dim.max, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK_val <- bcmvn %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::pull(pK) %>%
    as.numeric()

  annotations <- objeto@meta.data[[col_annotation]]
  homotypic.prop <- modelHomotypic(annotations)
  nExp <- round(doublet_rate * nrow(objeto@meta.data))
  nExp.adj <- round(nExp * (1 - homotypic.prop))

  objeto <- doubletFinder(objeto,
                          PCs = 1:dim.max,
                          pN = 0.25,
                          pK = pK_val,
                          nExp = nExp.adj,
                          reuse.pANN = FALSE,
                          sct = FALSE)

  return(objeto)
}


# ...............................................................................
# ARGUMENTOS 
# ...............................................................................

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Uso: Rscript 02_DoubletFinder.R <umbral_mt>")
}
umbral_mt <- args[1]


# ...............................................................................
# PATHS
# ...............................................................................

main_path_objeto <- paths_AnnaMaria$Path_seurat_object
input_path <- file.path(main_path_objeto, "3.SeuratProcessed", umbral_mt, "Seurat_Processed.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_path <- file.path(path.guardar_original, "4.DoubletFinder", umbral_mt)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)


# ...............................................................................
# CARGAR OBJETO SEURAT
# ...............................................................................
seu.obj <- readRDS(input_path)

# ...............................................................................
# EJECUTAR DOUBLETFINDER 
# ...............................................................................
seu.obj <- DoubletFinderPipeline(seu.obj, col_annotation = "Sample")


# ...............................................................................
# GUARDAR OBJETO 
# ...............................................................................

saveRDS(seu.obj, file = file.path(output_path, "Seurat_DoubletFinder.rds"))


# ...............................................................................
# PLOTS 
# ...............................................................................

# Buscar el nombre de la columna de clasificación de DF
df_col <- colnames(seu.obj@meta.data)[grepl("^DF.classifications", colnames(seu.obj@meta.data))][1]

# UMAP coloreado por clasificación DoubletFinder
if (!is.na(df_col)) {
  p_df <- DimPlot(seu.obj, group.by = df_col, pt.size = 1, raster = FALSE) & NoAxes()
  ggsave(filename = file.path(output_path, "DoubletFinder_Classification.png"), plot = p_df, width = 7, height = 7)
}

# UMAP por Sample
if ("Sample" %in% colnames(seu.obj@meta.data)) {
  p_sample <- DimPlot(seu.obj, group.by = "Sample", pt.size = 1, raster = FALSE) & NoAxes()
  ggsave(filename = file.path(output_path, "Sample.png"), plot = p_sample, width = 7, height = 7)
}

# UMAP por Phenotype 
if ("Phenotype" %in% colnames(seu.obj@meta.data)) {
  p_pheno <- DimPlot(seu.obj, group.by = "Phenotype", pt.size = 1, raster = FALSE) & NoAxes()
  ggsave(filename = file.path(output_path, "Phenotype.png"), plot = p_pheno, width = 7, height = 7)
}
