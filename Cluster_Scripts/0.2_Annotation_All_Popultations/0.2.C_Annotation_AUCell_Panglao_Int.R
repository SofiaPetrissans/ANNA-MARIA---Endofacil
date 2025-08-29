# --------------------------------------------------------------------------------
# SCRIPT: Anotación superficial + estimación de marcadores para múltiples thresholds
# AUTORA: Sofía Petrissans Moll
# FECHA: 11/06/2025
# --------------------------------------------------------------------------------


# ...............................................................................
# LIBRERIAS
# ...............................................................................
source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.2_Annotation_All_Popultations/0.0_Paths.R")
library(Seurat)
library(AUCell)
library(tidyverse)
library(ggplot2)


# ...............................................................................
# ARGUMENTOS
# ...............................................................................

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript 0.8.B_Annotation_AUCell_Panglao.R <THRESHOLD> <RESOLUTION>")
}
THRESHOLD <- args[1]
RESOLUTION <- args[2]


# ...............................................................................
# PATHS
# ...............................................................................

base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "6.Integration", THRESHOLD, "Harmony.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, "2.Annotations_Integrated", THRESHOLD, paste0("res", RESOLUTION))
plot_dir <- file.path(output_dir, "AUCell")
# Crear carpetas
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# ...............................................................................
# Cargar objeto seurat
# ...............................................................................

seurat_object <- readRDS(input_rds)
resolution_col <- paste0("Harmony_Log_res.", RESOLUTION)
seurat_object <- SetIdent(seurat_object, value = resolution_col)


# ...............................................................................
# Cargar markers de PanglaoDB
# ...............................................................................

path_markers <- "/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/markers_panglao"
markers.file <- file.path(path_markers, "PanglaoDB_markers.tsv")
markers_all <- read.csv(markers.file, sep = "\t")


# CONVERTIR -- DE MAYUSCULAS A PRIMERA MAYUSCULA Y RESTO MINUSCULA
mousify <- function(a){
  return(paste0(substr(a,1,1), tolower(substr(a,2,nchar(a)))))
  
}

# FILTRAR POBLACIONES DE INTERES 
list_populations <- c("Cholangiocytes", "Hepatic stellate cells", "Hepatoblasts", "Hepatocytes", "Kupffer cells", "B cells", "T cells", 
  "NK cells", "Neutrophils", "Macrophages", "Monocytes", "Endothelial cells", "Pericytes", "Epithelial cells", "Fibroblasts")

# ...............................................................................
# Conteos y rankings
# ...............................................................................
print(Assays(seurat_object))
counts <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")
cell_rankings <- AUCell_buildRankings(counts, plotStats = FALSE)


# ...............................................................................
# BUCLE POR POBLACIÓN DE INTERÉS 
# ...............................................................................

for (population in list_populations) {
  message("Procesando: ", population)
  
  markers <- markers_all[markers_all$cell.type == population & markers_all$species != "Hs", ]
  genes <- unique(sapply(markers$official.gene.symbol, mousify))
  genes_validos <- genes[genes %in% rownames(seurat_object)]

  if (length(genes_validos) < 5) {
    warning(paste("Pocos genes válidos para", population, "- saltando"))
    next
  }

  # Calcular AUC
  cells_AUC <- AUCell_calcAUC(setNames(list(genes_validos), population), cell_rankings)

  # Guardar AUCs
  write.csv(getAUC(cells_AUC), file = file.path(plot_dir, paste0("AUC_", population, ".csv")))
  saveRDS(cells_AUC, file = file.path(plot_dir, paste0("AUC_", population, ".rds")))

  # Crear y guardar histograma
  hist_file <- file.path(plot_dir, paste0("AUCell_Threshold_", population, ".pdf"))
  pdf(hist_file, width = 8, height = 6)
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign = TRUE)
  dev.off()

  if (is.null(cells_assignment) || length(cells_assignment) == 0) {
    warning(paste("AUCell_exploreThresholds no devolvió umbrales para", population, "- saltando"))
    next
  }

  # Verificar threshold
  gene_set_name <- names(cells_assignment)[1]
  threshold_info <- cells_assignment[[gene_set_name]]$thresholds[[1]]

  if (is.null(threshold_info$aucThr)) {
    warning(paste("No se pudo determinar threshold para", population, "- saltando"))
    next
  }
  threshold_value <- threshold_info$aucThr

  # Extraer células positivas
  new_cells <- names(which(getAUC(cells_AUC)[gene_set_name, ] > threshold_value))

  # Plot
  # Crear columna categórica en meta.data
  meta_colname <- gsub(" ", "_", population)
  meta_column_values <- ifelse(colnames(seurat_object) %in% new_cells, population, paste0("no-", population))
  seurat_object <- AddMetaData(seurat_object, metadata = meta_column_values, col.name = meta_colname)

  # Comprobaciones
  print(table(seurat_object[[population]][,1]))
  print(head(new_cells))
  message("Nº células positivas para ", population, ": ", length(new_cells))
}

# ...............................................................................
# GUARDAR OBJETO FINAL
# ...............................................................................
saveRDS(seurat_object, file = file.path(output_dir, "Seurat_WithAUCell.rds"))






