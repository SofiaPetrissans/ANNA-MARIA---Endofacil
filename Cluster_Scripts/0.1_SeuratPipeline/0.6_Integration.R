# --------------------------------------------------------------------------------
# SCRIPT: Integration with Harmony
# AUTORA: Sofia Petrissans Moll
# FECHA: 16/06/2025
# --------------------------------------------------------------------------------
source('/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(Seurat)
library(harmony)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(scales)

getPalette <- colorRampPalette(brewer.pal(9, "Set3"))

EstimatePCA.Dims <- function(data) {
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  dim.final <- min(co1, co2, na.rm = TRUE)
  return(dim.final)
}

# ARGUMENTOS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Uso: Rscript 0.7_Integration.R <umbral_mt>")
}
umbral_mt <- args[1]


# PATHS
main_path_objeto <- paths_AnnaMaria$Path_seurat_object
input_rds <- file.path(main_path_objeto, "5.Scrublet", umbral_mt, "Seurat_WithScrublet.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original , "6.Integration", umbral_mt)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# CARGAR
seurat_object <- readRDS(input_rds)

# PCA + DIMS
dim.final <- EstimatePCA.Dims(seurat_object)
print(dim.final)

ElbowPlot(seurat_object, ndims = 30, reduction = "pca")
ggsave(filename = file.path(output_dir, "Harmony_ElbowPlot_Log.png"), width = 5, height = 5)

# HARMONY + CLUSTERING
seurat_object.harmony <- RunHarmony(seurat_object, group.by.vars = "Sample", dims = 1:dim.final, plot_convergence = TRUE, reduction.save = "HarmonyLog", assay = "RNA")
seurat_object.harmony <- RunUMAP(seurat_object.harmony, reduction = "HarmonyLog", dims = 1:dim.final)
seurat_object.harmony <- FindNeighbors(seurat_object.harmony, reduction = "HarmonyLog", dims = 1:dim.final, graph.name = "Harmony_Log")
resolutions <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1)
seurat_object.harmony <- FindClusters(seurat_object.harmony, resolution = resolutions, graph.name = "Harmony_Log")

# GUARDAR RDS
saveRDS(seurat_object.harmony, file.path(output_dir, "Harmony.rds"))

# PLOTS
for (res in c("0.1", "0.3", "0.5", "0.7")) {
  col <- getPalette(length(unique(seurat_object.harmony[[paste0("Harmony_Log_res.", res)]])) + 30)
  p <- DimPlot(seurat_object.harmony, group.by = paste0("Harmony_Log_res.", res), pt.size = 1, label = TRUE, cols = col, raster = FALSE) & NoAxes()
  ggsave(filename = file.path(output_dir, paste0("Harmony_Log_res.", res, ".png")), plot = p, width = 10, height = 10)
}

# SAMPLE
col <- getPalette(length(unique(seurat_object.harmony$Sample)))
p <- DimPlot(seurat_object.harmony, group.by = "Sample", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(output_dir, "Sample.png"), plot = p, width = 10, height = 10)

# PHENOTYPE
col <- getPalette(length(unique(seurat_object.harmony$Phenotype)))
p <- DimPlot(seurat_object.harmony, group.by = "Phenotype", pt.size = 1, label = FALSE, cols = col, raster = FALSE) & NoAxes()
ggsave(filename = file.path(output_dir, "Phenotype.png"), plot = p, width = 10, height = 10)



