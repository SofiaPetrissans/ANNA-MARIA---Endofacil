# --------------------------------------------------------------------------------
# SCRIPT: Generar DimPlots con thresholds estadísticos de AUCell
# AUTORA: Sofia Petrissans Moll
# FECHA: 20/06/2025
# --------------------------------------------------------------------------------

source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.2_Annotation_All_Popultations/0.0_Paths.R")
library(Seurat)
library(AUCell)
library(ggplot2)
library(tidyverse)

# ...............................................................................
# ARGUMENTOS
# ...............................................................................
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript 0.6.C_GeneratePlots_AUCell.R <THRESHOLD> <RESOLUTION>")
}
THRESHOLD <- args[1]
RESOLUTION <- args[2]


# ...............................................................................
# PATHS
# ...............................................................................
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "5.Scrublet", THRESHOLD, "Seurat_WithScrublet.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
auc_dir <- file.path(path.guardar_original, "1.Annotations", THRESHOLD, paste0("res", RESOLUTION), "AUCell")
output_dir <- file.path(path.guardar_original, "1.Annotations", THRESHOLD, paste0("res", RESOLUTION), "Plots_AUCell")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ...............................................................................
# Cargar objeto Seurat
# ...............................................................................
seurat_object <- readRDS(input_rds)
seurat_object <- SetIdent(seurat_object, value = paste0("Harmony_Log_res.", RESOLUTION))

# ...............................................................................
# Detectar poblaciones automáticamente
# ...............................................................................
auc_files <- list.files(auc_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(auc_files) == 0) stop("No se encontraron archivos AUC .rds en ", auc_dir)
populations <- gsub("AUC_|\\.rds", "", basename(auc_files))

# ...............................................................................
# Procesar cada población
# ...............................................................................
for (i in seq_along(populations)) {
  population <- populations[i]
  auc_file <- auc_files[i]
  message("Procesando: ", population)

  cells_AUC <- readRDS(auc_file)

  gene_set_name <- rownames(getAUC(cells_AUC))[1]
  message("  Usando gene_set_name: ", gene_set_name)

  auc_values <- getAUC(cells_AUC)[gene_set_name, ]

  thresholds_to_test <- list(
    Percentil_95 = quantile(auc_values, 0.95, na.rm = TRUE),
    Mean = mean(auc_values, na.rm = TRUE),
    Mean_1SD = mean(auc_values, na.rm = TRUE) + sd(auc_values, na.rm = TRUE),
    Mean_2SD = mean(auc_values, na.rm = TRUE) + 2 * sd(auc_values, na.rm = TRUE)
  )

  for (thresh_name in names(thresholds_to_test)) {
    threshold_value <- thresholds_to_test[[thresh_name]]
    message("  Aplicando threshold: ", thresh_name, " (", round(threshold_value, 4), ")")

    new_cells <- names(auc_values)[which(auc_values > threshold_value)]

    col_name <- paste0(gsub(" ", "_", population), "_AUCell_", thresh_name)
    meta_values <- ifelse(colnames(seurat_object) %in% new_cells, "Active", "Inactive")
    seurat_object <- AddMetaData(seurat_object, metadata = meta_values, col.name = col_name)

    p_dim <- DimPlot(seurat_object, group.by = col_name, label = TRUE, cols = c("purple", "grey80")) & NoAxes()
    ggsave(filename = file.path(output_dir, paste0("DimPlot_", gsub(" ", "_", population), "_", thresh_name, ".png")),
           plot = p_dim, width = 8, height = 7)

    pdf(file.path(output_dir, paste0("Histogram_", gsub(" ", "_", population), "_", thresh_name, ".pdf")),
        width = 8, height = 6)
    AUCell_plotHist(cells_AUC[gene_set_name, ], aucThr = threshold_value)
    dev.off()

    message("    Nº células positivas: ", length(new_cells))
  }
}

# ...............................................................................
# Guardar Seurat anotado
# ...............................................................................
saveRDS(seurat_object, file = file.path(output_dir, "Seurat_WithAUCell_Annotations.rds"))
