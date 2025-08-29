# --------------------------------------------------------------------------------
# SCRIPT: Subset Endothelial + Seurat Pipeline + Harmony Integration
# AUTORA: Sofia Petrissans Moll
# FECHA: 20/06/2025
# --------------------------------------------------------------------------------
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.4_SubsetEndothelial/0.0_Paths.R")
library(Seurat)
library(tidyverse)
library(harmony)
library(RColorBrewer)

# ...............................................................................
# ARGUMENTOS
# ...............................................................................
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript 1.0_Subset_Seurat.R <MT_THRESHOLD> <RESOLUTION> <CLUSTERS_TO_KEEP>")
}

MT_THRESHOLD <- args[1]  # Ej: mt5
RESOLUTION <- args[2]    # Ej: RNA_snn_res.0.3
CLUSTERS_TO_KEEP <- strsplit(args[3], ",")[[1]]  # Ej: "2,7,9"

# ...............................................................................
# PATHS
# ...............................................................................
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "3.Annotation_Layer", MT_THRESHOLD, "Seurat_annotated.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, MT_THRESHOLD)
raw_dir <- file.path(output_dir, "1.SeuratProcessed")
harmony_dir <- file.path(output_dir, "2.Integration")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(harmony_dir, recursive = TRUE, showWarnings = FALSE)

# ...............................................................................
# FUNCIONES Y PALETA
# ...............................................................................
getPalette <- colorRampPalette(c(brewer.pal(12, "Set3"), brewer.pal(8, "Pastel1")))

EstimatePCA.Dims <- function(data) {
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
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

getClusterPalette <- function(seu, group_col) {
  n <- length(unique(seu[[group_col]][,1]))
  if (n < 2) return(NULL)
  return(getPalette(n))
}

# ...............................................................................
# CARGA Y SUBSET
# ...............................................................................
seurat_object <- readRDS(input_rds)
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)
subset_object <- subset(seurat_object, idents = CLUSTERS_TO_KEEP)

write.csv(data.frame(CellsBefore = ncol(subset_object)), file = file.path(output_dir, "NumberOfCells_BeforeProcessing.csv"))

# ...............................................................................
# PIPELINE RAW
# ...............................................................................
subset_object <- SeuratPipeline(subset_object)
saveRDS(subset_object, file = file.path(raw_dir, "SubsetEndothelial_Raw.rds"))

for (res in c(0.1, 0.3, 0.5, 1)) {
  group_col <- paste0("RNA_snn_res.", res)
  col <- getClusterPalette(subset_object, group_col)
  
  umap_plot <- DimPlot(subset_object, group.by = group_col, label = TRUE, pt.size = 1, raster = FALSE) & NoAxes()
  if (!is.null(col)) {
    umap_plot <- umap_plot + scale_color_manual(values = col)
  }
  
  ggsave(filename = file.path(raw_dir, paste0("RNA_snn_res.", res, ".png")), plot = umap_plot, width = 7, height = 7)
}

if ("Sample" %in% colnames(subset_object@meta.data)) {
  n_samples <- length(unique(subset_object$Sample))
  umap_sample <- DimPlot(subset_object, group.by = "Sample", pt.size = 1, raster = FALSE) +
    scale_color_manual(values = getPalette(n_samples)) & NoAxes()
  ggsave(filename = file.path(raw_dir, "Sample.png"), plot = umap_sample, width = 7, height = 7)
}

if ("Phenotype" %in% colnames(subset_object@meta.data)) {
  n_pheno <- length(unique(subset_object$Phenotype))
  umap_pheno <- DimPlot(subset_object, group.by = "Phenotype", pt.size = 1, raster = FALSE) +
    scale_color_manual(values = getPalette(n_pheno)) & NoAxes()
  ggsave(filename = file.path(raw_dir, "Phenotype.png"), plot = umap_pheno, width = 7, height = 7)
}

# ...............................................................................
# HARMONY INTEGRATION
# ...............................................................................
dim.final <- EstimatePCA.Dims(subset_object)

ElbowPlot(subset_object, ndims = 30, reduction = "pca")
ggsave(filename = file.path(output_dir, "Harmony_ElbowPlot_Log.png"), width = 5, height = 5)

subset_object <- RunHarmony(subset_object, group.by.vars = "Sample", dims = 1:dim.final, plot_convergence = TRUE, reduction.save = "HarmonyLog", assay = "RNA")
subset_object <- RunUMAP(subset_object, reduction = "HarmonyLog", dims = 1:dim.final)
subset_object <- FindNeighbors(subset_object, reduction = "HarmonyLog", dims = 1:dim.final, graph.name = "Harmony_Log")
subset_object <- FindClusters(subset_object, resolution = c(0.1, 0.3, 0.5, 0.7, 0.9, 1), graph.name = "Harmony_Log")

saveRDS(subset_object, file = file.path(harmony_dir, "SubsetEndothelial_Harmony.rds"))

for (res in c(0.1, 0.3, 0.5, 1)) {
  group_col <- paste0("Harmony_Log_res.", res)
  col <- getClusterPalette(subset_object, group_col)
  
  umap_plot <- DimPlot(subset_object, reduction = "umap", group.by = group_col, label = TRUE, pt.size = 1, raster = FALSE) & NoAxes()
  if (!is.null(col)) {
    umap_plot <- umap_plot + scale_color_manual(values = col)
  }
  
  ggsave(filename = file.path(harmony_dir, paste0("Harmony_Log_res.", res, ".png")), plot = umap_plot, width = 7, height = 7)
}

if ("Sample" %in% colnames(subset_object@meta.data)) {
  n_samples <- length(unique(subset_object$Sample))
  umap_sample <- DimPlot(subset_object, reduction = "umap", group.by = "Sample", pt.size = 1, raster = FALSE) +
    scale_color_manual(values = getPalette(n_samples)) & NoAxes()
  ggsave(filename = file.path(harmony_dir, "Sample.png"), plot = umap_sample, width = 7, height = 7)
}

if ("Phenotype" %in% colnames(subset_object@meta.data)) {
  n_pheno <- length(unique(subset_object$Phenotype))
  umap_pheno <- DimPlot(subset_object, reduction = "umap", group.by = "Phenotype", pt.size = 1, raster = FALSE) +
    scale_color_manual(values = getPalette(n_pheno)) & NoAxes()
  ggsave(filename = file.path(harmony_dir, "Phenotype.png"), plot = umap_pheno, width = 7, height = 7)
}

cat("Análisis completado. Resultados guardados en:\n")
cat("- Raw: ", raw_dir, "\n")
cat("- Harmony: ", harmony_dir, "\n")
