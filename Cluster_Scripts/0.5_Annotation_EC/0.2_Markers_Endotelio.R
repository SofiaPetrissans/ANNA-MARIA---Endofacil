# --------------------------------------------------------------------------------
# SCRIPT: Anotación superficial + estimación de marcadores para múltiples thresholds
# AUTORA: Sofia Petrissans Moll
# FECHA: 30/06/2025
# --------------------------------------------------------------------------------

# ...............................................................................
# LIBRERIAS
# ...............................................................................
source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggplot2)

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

# ...............................................................................
# ARGUMENTOS
# ...............................................................................

RESOLUTION <- c("0.1", "0.3", "0.5")
MT_THRESHOLD <- "mt20"
TYPE <- "Integation"

# ...............................................................................
# PATHS
# ...............................................................................

base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path,  MT_THRESHOLD, "8.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")
# Cargar el objeto
seurat_object <- readRDS(input_rds)

path.guardar_original <- paths_AnnaMaria$Path_guardar
path.guardar_original <- fil.path(path.guardar_original, MT_THRESHOLD)

for (res in RESOLUTION) {
  resolution_col <- if (TYPE == "Raw") {
    paste0("RNA_snn_res.", res)
  } else {
    paste0("Harmony_Log_res.", res)
  }

  if (!resolution_col %in% colnames(seurat_object@meta.data)) {
    warning(paste("Columna de resolución", resolution_col, "no encontrada. Saltando..."))
    next
  }

  output_dir <- file.path(path.guardar_original, "1.MarkersAnnotation", paste0("res", res))
  plot_dir <- file.path(output_dir, "Plots")
  marker_dir <- file.path(output_dir, "Markers")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)

  # Anotación y marcadores
  SimpleAnnotator(seurat_object, resolution_col, plot_dir)

  message("Anotación y marcadores completados para resolución ", res)
}


SimpleAnnotator(seurat_object, resolution_col, plot_dir)

message("Anotación y marcadores completados para ", MT_THRESHOLD, " ", TYPE, " resolución ", RESOLUTION)


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
ggsave(filename = file.path(plot_dir, "CanonicalMarkers_DotPlot_Endothelial.png"),
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
ggsave(filename = file.path(plot_dir, "CanonicalMarkers_DotPlot_sep_Endothelial.png"),
       plot = dp, width = 18, height = 5)







seurat_object <- SetIdent(seurat_object,value="Annotation_LSEC_LVEC")
gene_sets <- list(
    LSEC = c("Tiam1", "Gna14", "Stab2", "Fcn2", "Fcn3","Clec4g", "Lyve1", "Ms4a6a", "Cd36"),
    LVEC = c("Pecam1", "Cd34", "Vwf", "Sema3g", "Rgcc", "Cd200", "Esm1")
  )

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
ggsave(filename = file.path(path.guardar, "CanonicalMarkers_DotPlot_sep_LSEC_LVEC.png"),
       plot = dp, width = 18, height = 5)




# Opcional: ordenar las categorías como niveles para el facetado
dp_data$Category <- factor(dp_data$Category, levels = names(gene_sets))

# Crear heatmap
hm <- ggplot(dp_data, aes(x = features.plot, y = id, fill = avg.exp.scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "grey90", mid = "white", high = "#3C2692", midpoint = 0) +
  facet_wrap(~Category, scales = "free_x", nrow = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  labs(x = "Markers", y = "Identity", fill = "Scaled Expression")

# Guardar
ggsave(file.path(path.guardar, "CanonicalMarkers_Heatmap_sep__LSEC_LVE.png"),
       plot = hm, width = 18, height = 5)