################################################################################

# SCRIPT: Seurat Pipeline diferentes resoluciones
# AUTHOR: SOFIA PETRISSANS MOLL
# DATE: 07-06-2025

################################################################################

# SUMMARY SCRIPT:

  # Aplica el pipeline completo de Seurat a los objetos ya filtrados.
  # Acepta argumentos desde SLURM (umbral de % mitocondrial y carpeta de salida).
  # Crea carpetas separadas por cada umbral (mt10, mt15, etc.).
  # Aplica múltiples resoluciones (res = 0.1–2), generando los UMAPs correspondientes.
  # Exporta gráficos de calidad (nFeatures, nCounts, %mito) y por Sample.



# LIBRARIES ....................................................................
source('/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(scales)


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

SeuratPipeline <- function(data) {
  data <- NormalizeData(data)
  gene.info <- summary(Matrix::colSums(data@assays$RNA@counts > 0))
  hvg.number <- round(gene.info[4] + 100)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = hvg.number)
  data <- ScaleData(data)
  data <- RunPCA(data)
  dim.final <- EstimatePCA.Dims(data)
  data <- FindNeighbors(data, dims = 1:dim.final)
  data <- RunUMAP(data, dims = 1:dim.final)
  for (res in c(0.1, 0.3, 0.5, 1, 1.5, 2)) {
    data <- FindClusters(data, resolution = res)
  }
  return(data)
}


# ...............................................................................
# ARGUMENTOS DE SLURM 
# ...............................................................................

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Uso: Rscript 01_SeuratPipeline.R <umbral_mt> "
}
umbral_mt <- args[1]

# ...............................................................................
# DIRECTORIOS 
# ...............................................................................
main_path_objeto <- paths_AnnaMaria$Path_seurat_object
data_path <- file.path(main_path_objeto, "2.Filtering", umbral_mt, "QC_Filtered.rds")


# Crear carpeta específica para cada threshold
path.guardar_original <- paths_AnnaMaria$Path_guardar
output_path <- file.path(path.guardar_original, "3.SeuratProcessed", umbral_mt)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)


# ...............................................................................
# CARGAR OBJETO SEURAT 
# ...............................................................................

data <- readRDS(data_path)

# ...............................................................................
# SEURAT PIPELINE
# ...............................................................................

seu.obj <- SeuratPipeline(data)

# ...............................................................................
# GUARDAR OBJETO 
# ...............................................................................

saveRDS(seu.obj, file = file.path(output_path, "Seurat_Processed.rds"))


# ...............................................................................
# GENERAR UMAPS 
# ...............................................................................
resolutions <- c("0.1", "0.3", "0.5", "1", "1.5", "2")
for (res in resolutions) {
  res_id <- paste0("RNA_snn_res.", res)
  if (res_id %in% colnames(seu.obj@meta.data)) {
    n_clusters <- length(unique(seu.obj[[res_id]][, 1]))
    col <- getPalette(n_clusters)
    umap <- DimPlot(seu.obj, reduction = "umap", group.by = res_id,
                    label = TRUE, label.size = 5, cols = col, pt.size = 1, raster = FALSE) &
            NoAxes() & NoLegend()
    ggsave(filename = file.path(output_path, paste0(res_id, "_UMAP.png")), 
           plot = umap, width = 7, height = 7)
  }
}

# ...............................................................................
# OTROS PLOTS QC 
# ...............................................................................
FeaturePlotWrap <- function(feature, name) {
  p <- FeaturePlot(seu.obj, features = feature, pt.size = 1, raster = FALSE) & NoAxes()
  ggsave(filename = file.path(output_path, paste0(name, ".png")), plot = p, width = 7, height = 7)
}

FeaturePlotWrap("percent.mt", "MT")
FeaturePlotWrap("nCount_RNA", "Counts")
FeaturePlotWrap("nFeature_RNA", "Features")


# PLOTS DE METADATOS 

col <- getPalette(length(unique(seu.obj$Sample)))
p <- DimPlot(seu.obj, group.by = "Sample", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(output_path, "Sample.png"), plot = p, width = 10, height = 10)

# PHENOTYPE
col <- getPalette(length(unique(seu.obj$Phenotype)))
p <- DimPlot(seu.obj, group.by = "Phenotype", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(output_path, "Phenotype.png"), plot = p, width = 10, height = 10)
