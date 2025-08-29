# --------------------------------------------------------------------------------
# SCRIPT: Subset Endothelial + Remove Clusters
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

MT_THRESHOLD <- args[1]  # Ej: mt15
RESOLUTION <- args[2]    # Ej: RNA_snn_res.0.3
CLUSTERS_TO_KEEP <- as.numeric(strsplit(args[3], ",")[[1]])


# ...............................................................................
# PATHS
# ...............................................................................
base_path <- paths_AnnaMaria$Path_seurat_object
input_rds <- file.path(base_path, MT_THRESHOLD, "2.Integration", "SubsetEndothelial_Harmony.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, MT_THRESHOLD, "4.RemoveCluster")
raw_dir <- file.path(output_dir, "1.SeuratProcessed")
harmony_dir <- file.path(output_dir, "2.Integration")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(harmony_dir, recursive = TRUE, showWarnings = FALSE)

# ...............................................................................
# CARGA Y SUBSET
# ...............................................................................
seurat_object <- readRDS(input_rds)
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)
subset_object <- subset(seurat_object, idents = CLUSTERS_TO_KEEP)

write.csv(data.frame(CellsBefore = ncol(subset_object)), file = file.path(output_dir, "NumberOfCells_BeforeProcessing.csv"))


if (!all(CLUSTERS_TO_KEEP %in% levels(Idents(seurat_object)))) {
  stop("Algunos clusters especificados no existen en el objeto Seurat: ",
       paste(setdiff(CLUSTERS_TO_KEEP, levels(Idents(seurat_object))), collapse = ", "))
}

# ...............................................................................
# FUNCIONES
# ...............................................................................
getPalette <- colorRampPalette(brewer.pal(9, "Set3"))

EstimatePCA.Dims <- function(data) {
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
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


generateClusterColors <- function(object, group_col) {
  cluster_levels <- levels(factor(object@meta.data[[group_col]]))
  n_clusters <- length(cluster_levels)
  colors <- getPalette(n_clusters)
  names(colors) <- cluster_levels
  return(colors)
}
# ...............................................................................
# PIPELINE RAW
# ...............................................................................
subset_object <- SeuratPipeline(subset_object)
saveRDS(subset_object, file = file.path(raw_dir, "SubsetEndothelial_Raw.rds"))

# Plots RAW
for (res in c(0.1, 0.3, 0.5, 1)) {
  group_col <- paste0("RNA_snn_res.", res)
  cluster_colors <- generateClusterColors(subset_object, group_col)
  
  umap_plot <- DimPlot(subset_object, group.by = group_col, label = TRUE,
                       pt.size = 1, raster = FALSE, cols = cluster_colors) & NoAxes()
  
  ggsave(filename = file.path(raw_dir, paste0(group_col, ".png")), plot = umap_plot, width = 7, height = 7)
}


col <- getPalette(length(unique(subset_object$Sample)))
p <- DimPlot(subset_object, group.by = "Sample", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(raw_dir, "Sample.png"), plot = p, width = 10, height = 10)

# PHENOTYPE
col <- getPalette(length(unique(subset_object$Phenotype)))
p <- DimPlot(subset_object, group.by = "Phenotype", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(raw_dir, "Phenotype.png"), plot = p, width = 10, height = 10)

# ...............................................................................
# HARMONY INTEGRATION
# ...............................................................................
integration.var <- "Sample"

dim.final <- EstimatePCA.Dims(subset_object)
print(dim.final)

ElbowPlot(subset_object, ndims = 30, reduction = "pca")
ggsave(filename = file.path(output_dir, "Harmony_ElbowPlot_Log.png"), width = 5, height = 5)

subset_object <- RunHarmony(subset_object, group.by.vars = integration.var, dims = 1:dim.final, plot_convergence = TRUE, reduction.save = "HarmonyLog", assay="RNA")
subset_object <- RunUMAP(subset_object, reduction = "HarmonyLog", dims = 1:dim.final)
subset_object <- FindNeighbors(subset_object, reduction = "HarmonyLog", dims = 1:dim.final, graph.name = "Harmony_Log")
resolutions <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1)
subset_object <- FindClusters(subset_object, resolution = resolutions, graph.name = "Harmony_Log")

saveRDS(subset_object, file = file.path(harmony_dir, "SubsetEndothelial_Harmony.rds"))

# Plots HARMONY
for (res in c(0.1, 0.3, 0.5, 1)) {
  group_col <- paste0("Harmony_Log_res.", res)
  cluster_colors <- generateClusterColors(subset_object, group_col)
  
  umap_plot <- DimPlot(subset_object, reduction = "umap", group.by = group_col, label = TRUE,
                       pt.size = 1, raster = FALSE, cols = cluster_colors) & NoAxes()
  
  ggsave(filename = file.path(harmony_dir, paste0("Harmony_Log_res.", res, ".png")),
         plot = umap_plot, width = 7, height = 7)
}

col <- getPalette(length(unique(subset_object$Sample)))
p <- DimPlot(subset_object, group.by = "Sample", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(harmony_dir, "Sample.png"), plot = p, width = 10, height = 10)

# PHENOTYPE
col <- getPalette(length(unique(subset_object$Phenotype)))
p <- DimPlot(subset_object, group.by = "Phenotype", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(harmony_dir, "Phenotype.png"), plot = p, width = 10, height = 10)


cat("Análisis completado. Resultados guardados en:\n")
cat("- Harmony: ", harmony_dir, "\n")