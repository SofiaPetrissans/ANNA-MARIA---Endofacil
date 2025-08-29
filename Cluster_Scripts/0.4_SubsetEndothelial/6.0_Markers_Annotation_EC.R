# --------------------------------------------------------------------------------
# SCRIPT: Anotación superficial + estimación de marcadores para múltiples thresholds
# AUTORA: Sofía Petrissans Moll
# FECHA: 11/06/2025
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# SCRIPT: Anotación superficial + estimación de marcadores para múltiples thresholds
# AUTORA: Sofia Petrissans Moll
# FECHA: 20/06/2025
# --------------------------------------------------------------------------------

# ...............................................................................
# LIBRERIAS
# ...............................................................................
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.4_SubsetEndothelial/0.0_Paths.R")
library(Seurat)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggplot2)

# ...............................................................................
# ARGUMENTOS
# ...............................................................................
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript 4.0_Markers_Annotation.R <THRESHOLD> <RESOLUTION> <TYPE>")
}

THRESHOLD <- args[1]
RESOLUTION <- args[2]
TYPE <- args[3]  # "Raw" o "Integration"

# ...............................................................................
# PATHS
# ...............................................................................
base_path <- paths_AnnaMaria$Path_seurat_object
input_dir <- if (TYPE == "Raw") "4.RemoveCluster/1.SeuratProcessed" else "4.RemoveCluster/2.Integration"
input_rds <- file.path(base_path, THRESHOLD, input_dir, 
                       if (TYPE == "Raw") "SubsetEndothelial_Raw.rds" else "SubsetEndothelial_Harmony.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, THRESHOLD, "7.MarkersAnnotation_EC", TYPE, paste0("res", RESOLUTION))
plot_dir <- file.path(output_dir, "Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ...............................................................................
# FUNCIONES
# ...............................................................................
SimpleAnnotator <- function(seu, resolution, plot_dir) {
  seu <- SetIdent(seu, value = resolution)
  
  gene_sets <- list(
    Artery = c("Gja5", "Hey1", "Efnb2", "Dll4", "Sox17", "Cxcl12", "Gkn3"),
    Vein = c("Nr2f2", "Nrp2", "Aplnr", "Flt4", "Ephb4", "Cdh5", "Slc38a5"),
    Capillary = c("Car4", "Plvap", "Cd36", "Aqp1", "Gpx3", "Tm4sf1"),
    Capillary_Liver = c("Stab1", "Stab2", "Lyve1", "Clec4g", "Mrc1"),
    Endothelial_Senescent = c("Cdkn1a", "Cdkn2a", "Mmp3", "Mmp12", "Serpine1", "Il6", "Il1b", "Cxcl1", "Cxcl2", "Tnf", "Igfbp3")
  )


  for (celltype in names(gene_sets)) {
    features <- gene_sets[[celltype]]
    features <- features[features %in% rownames(seu)]  # filtra genes existentes
    if (length(features) == 0) next
    p1 <- FeaturePlot(seu, features = features, order = TRUE, raster = FALSE) & NoAxes()
    p2 <- VlnPlot(seu, features = features, pt.size = 0.1, raster = FALSE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    g <- ggarrange(p1, p2, nrow = 1)
    ggsave(filename = file.path(plot_dir, paste0(celltype, ".png")), plot = g, width = 16, height = 7)
  }
}



# ...............................................................................
# EJECUCIÓN
# ...............................................................................
# Cargar el objeto
seurat_object <- readRDS(input_rds)
resolution_col <- if (TYPE == "Raw") {
  paste0("RNA_snn_res.", RESOLUTION)
} else {
  paste0("Harmony_Log_res.", RESOLUTION)
}

if (!resolution_col %in% colnames(seurat_object@meta.data)) {
  stop(paste("Columna de resolución", resolution_col, "no encontrada en el objeto."))
}

SimpleAnnotator(seurat_object, resolution_col, plot_dir)

message("Anotación completada para ", THRESHOLD, " ", TYPE, " resolución ", RESOLUTION)


# DOT PLOTS ADICIONALES (al estilo del script inicial)

  gene_sets <- list(
    Artery = c("Gja5", "Hey1", "Efnb2", "Dll4", "Sox17", "Cxcl12", "Gkn3"),
    Vein = c("Nr2f2", "Nrp2", "Aplnr", "Flt4", "Ephb4", "Cdh5", "Slc38a5"),
    Capillary = c("Car4", "Plvap", "Cd36", "Aqp1", "Gpx3", "Tm4sf1"),
    Capillary_Liver = c("Stab1", "Stab2", "Lyve1", "Clec4g", "Mrc1"),
    Endothelial_Senescent = c("Cdkn1a", "Cdkn2a", "Mmp3", "Mmp12", "Serpine1", "Il6", "Il1b", "Cxcl1", "Cxcl2", "Tnf", "Igfbp3")
  )

features <- unique(unlist(gene_sets))

seurat_object <- SetIdent(seurat_object, value = resolution_col)

# DotPlot conjunto
p_dot <- DotPlot(seurat_object, features = features) &
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(1.1)),
        axis.text.y = element_text(size = rel(1.05)))
ggsave(filename = file.path(plot_dir, "CanonicalMarkers_DotPlot.png"),
       plot = p_dot, width = 18, height = 5)

# DotPlot por categorías
gene_category_df <- bind_rows(lapply(names(gene_sets), function(cat) {
  data.frame(Feature = gene_sets[[cat]], Category = cat)
}))
dp_data <- DotPlot(seurat_object, features = gene_category_df$Feature)$data
dp_data <- dp_data %>%
  left_join(gene_category_df, by = c("features.plot" = "Feature"))
dp <- ggplot(dp_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "grey90", high = "#3C2692") +
  facet_wrap(~Category, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank()) +
  labs(x = "Features", y = "Identity", size = "Percent Expressed", color = "Average Expression")
ggsave(filename = file.path(plot_dir, "CanonicalMarkers_DotPlot_sep.png"),
       plot = dp, width = 18, height = 5)


