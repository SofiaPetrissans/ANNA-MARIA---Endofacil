################################################################################

# SCRIPT: Scrublet - Import to seurat
# AUTHOR: SOFIA PETRISSANS MOLL
# DATE: 10-06-2025

################################################################################

# SUMMARY SCRIPT:
# Qué hace
  # Importa los resultados .csv de Scrublet y los añade al objeto Seurat principal:



# LIBRERÍAS .........................................................................
source('/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(Seurat)
library(tidyverse)
library(RColorBrewer)


# ...............................................................................
# ARGUMENTOS 
# ...............................................................................

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Uso: Rscript 03_Scrublet_ImportToSeurat.R <umbral_mt>")
}

umbral_mt <- args[1]


# ...............................................................................
# PATHS
# ...............................................................................
main_path_objeto <- paths_AnnaMaria$Path_seurat_object
rds_path <- file.path(main_path_objeto, "4.DoubletFinder", umbral_mt, "Seurat_DoubletFinder.rds")

scrub_path <- file.path(main_path_objeto, "5.Scrublet", umbral_mt, "Input", "All_Scrublet_Results.csv")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_path <- file.path(path.guardar_original , "5.Scrublet", umbral_mt)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)


# ...............................................................................
# Cargar Objeto Seurat
# ...............................................................................
seurat_object <- readRDS(rds_path)
scrub_data <- read_csv(scrub_path)


# Paleta de colores
getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
colors <- getPalette(length(unique(seurat_object$Sample)))
# ...............................................................................
# MERGE SCRUBLET - SEURAT
# ...............................................................................
# Formatear barcodes si es necesario (verifica con head(scrub_data$Barcode))
rownames(scrub_data) <- scrub_data$Barcode

# Añadir metadatos
seurat_object <- AddMetaData(seurat_object, metadata = scrub_data[colnames(seurat_object), c("Scrublet_Score", "Scrublet_Classification")])
saveRDS(seurat_object, file = file.path(output_path, "Seurat_WithScrublet.rds"))

# UMAP coloreado por clasificación
p <- DimPlot(seurat_object, group.by = "Scrublet_Classification", pt.size = 1, raster = FALSE) & NoAxes()
ggsave(filename = file.path(output_path, "Scrublet_Classification.png"), plot = p, width = 7, height = 7)


# BARPLOT: Proporción de doublets por muestra 

# Asegurarse de tener estas columnas
if (all(c("Scrublet_Classification", "Sample") %in% colnames(seurat_object@meta.data))) {
  
  df_plot <- seurat_object@meta.data %>%
    group_by(Sample, Scrublet_Classification) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(perc = 100 * n / sum(n)) %>%
    filter(Scrublet_Classification == "Doublet")
  write.csv(df_plot, file.path(output_path, "Scrublet_DoubletRate_PerSample.csv"), row.names = FALSE)

  barplot_scrub <- ggplot(df_plot, aes(x = Sample, y = perc, fill = Sample)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    ylab("% Doublets (Scrublet)") +
    xlab("") +
    ggtitle("Proportion of Doublets per sample") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      legend.position = "none"
    )
  ggsave(filename = file.path(output_path, "Scrublet_DoubletRate_Barplot.png"),
         plot = barplot_scrub, width = 8, height = 5)
}


