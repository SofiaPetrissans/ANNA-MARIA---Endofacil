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
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggplot2)


# ...............................................................................
# ARGUMENTOS
# ...............................................................................
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript 6.1_AnnotateClusters.R <THRESHOLD> <RESOLUTION>")
}

THRESHOLD <- args[1]  # mt5, mt10, etc.
RESOLUTION <- args[2]  # 0.1, 0.3, 0.5


# ...............................................................................
# PATHS
# ...............................................................................
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "5.Scrublet", THRESHOLD, "Seurat_WithScrublet.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, "1.Annotations", THRESHOLD, paste0("res", RESOLUTION))
plot_dir <- file.path(output_dir, "Plots")
marker_dir <- file.path(output_dir, "Markers")

# Crear carpetas
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)

# ...............................................................................
# FUNCIONES
# ...............................................................................
SimpleAnnotator <- function(seu, resolution, plot_dir) {
  seu <- SetIdent(seu, value = resolution)
  
  gene_sets <- list(
    Hepatocytes = c("Alb", "Ttr", "Cyp2e1", "Cyp3a11", "Ass1", "Glul", "Cyp7a1", "Hnf4a"),
    Endothelial = c("Pecam1", "Vwf", "Kdr"),
    LSEC = c("Lyve1", "Stab2", "Clec4g", "Mrc1", "Fcgr2b", "Icam1", "Cd36"),
    LVEC = c("Cdh5", "Cd31",  "Flt1", "Tek", "Nos3"),
    Kupffer = c("Adgre1", "Clec4f", "Cd68", "Marco"),
    Monocytes = c("Cd14", "Ccr2", "Ly6c1"),
    TCells = c("Cd3e", "Cd8a", "Cd4", "Trac"),
    BCells = c("Cd79a", "Ms4a1", "Cd19"),
    NK = c("Nkg7", "Klrb1c", "Gzmb", "Klrf1"),
    Neutrophil = c("S100a8", "S100a9", "Lcn2", "Cxcr2"),
    HSC = c("Rbp1", "Des", "Lrat", "Pdgfrb", "Reln"),
    Cholangiocytes = c("Krt", "Epcam", "Sox9", "Cftr"),
    Fibroblasts = c("Col1a1", "Col1a2", "Dcn", "Pdgfra")
  )

  for (celltype in names(gene_sets)) {
    features <- gene_sets[[celltype]]
    p1 <- FeaturePlot(seu, features = features, order = TRUE, raster = FALSE) & NoAxes()
    p2 <- VlnPlot(seu, features = features, pt.size = 0.1, raster = FALSE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    g <- ggarrange(p1, p2, nrow = 1)
    ggsave(filename = file.path(plot_dir, paste0(celltype, ".png")), plot = g, width = 16, height = 7)
  }
}


ClusterMarkers <- function(seu, resolution, marker_dir, plot_dir) {
  seu <- SetIdent(seu, value = resolution)
  clusters <- levels(seu)

  all_markers <- list()

  for (clus in clusters) {
    subcells <- WhichCells(seu, idents = clus)
    if (length(subcells) < 10) {
      next
    }
    markers <- FindMarkers(seu, ident.1 = clus)
    markers$gene <- rownames(markers)
    markers$cluster <- clus
    all_markers[[paste0("Clus_", clus)]] <- markers
  }

  saveRDS(all_markers, file.path(marker_dir, paste0("Markers_res", RESOLUTION, ".rds")))
  write.xlsx(all_markers, file.path(marker_dir, paste0("Markers_res", RESOLUTION, ".xlsx")))

  # Total markers concatenated
  total <- bind_rows(all_markers)
  write.xlsx(total, file.path(marker_dir, paste0("Markers_res", RESOLUTION, "_Total.xlsx")))

  # NUEVO: DotPlot de top 4 genes por clúster
  top_markers <- total %>%
    group_by(cluster) %>%
    top_n(4, avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

  features_top <- unique(top_markers$gene)

  p_dot <- DotPlot(seu, features = features_top) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(1.1), face = "plain"),
      axis.text.y = element_text(size = rel(1.05), face = "plain")
    )

  ggsave(filename = file.path(plot_dir, paste0("Top4Markers_DotPlot_res", RESOLUTION, ".png")),
         plot = p_dot, width = 16, height = 7)

  message("DotPlot de top 4 genes por cluster guardado en ", plot_dir)
}

# ...............................................................................
# Aplicar las funciones para los diferentes thresholds y resoluciones
# ...............................................................................
seurat_object <- readRDS(input_rds)
resolution_col <- paste0("RNA_snn_res.", RESOLUTION)

if (!resolution_col %in% colnames(seurat_object@meta.data)) {
  stop(paste("Columna de resolución", resolution_col, "no encontrada en el objeto Seurat."))
}


SimpleAnnotator(seurat_object, resolution_col, plot_dir)
ClusterMarkers(seurat_object, resolution_col, marker_dir, plot_dir)

message("Anotación y detección de marcadores completada para ", THRESHOLD, " a resolución ", RESOLUTION)


# ...............................................................................
# DOT PLOT
# ...............................................................................
seurat_object <- SetIdent(seurat_object, value = resolution_col)

gene_sets <- list(
    Hepatocytes = c("Alb", "Ttr", "Cyp2e1", "Cyp3a11", "Ass1", "Glul", "Cyp7a1", "Hnf4a"),
    Endothelial = c("Pecam1", "Vwf", "Kdr"),
    LSEC = c("Lyve1", "Stab2", "Clec4g", "Mrc1", "Fcgr2b", "Icam1", "Cd36"),
    LVEC = c("Cdh5", "Cd31", "Flt1", "Tek", "Nos3"),
    Kupffer = c("Adgre1", "Clec4f", "Cd68", "Marco"),
    Monocytes = c("Cd14", "Ccr2", "Ly6c1"),
    TCells = c("Cd3e", "Cd8a", "Cd4", "Trac"),
    BCells = c("Cd79a", "Ms4a1", "Cd19"),
    NK = c("Nkg7", "Klrb1c", "Gzmb", "Klrf1"),
    Neutrophil = c("S100a8", "S100a9", "Lcn2", "Cxcr2"),
    HSC = c("Rbp1", "Des", "Lrat", "Pdgfrb", "Reln"),
    Cholangiocytes = c("Krt", "Epcam", "Sox9", "Cftr"),
    Fibroblasts = c("Col1a1", "Col1a2", "Dcn", "Pdgfra")
  )

features <- unique(unlist(gene_sets))


# Crea el DotPlot
p_dot <- DotPlot(seurat_object, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))
        
ggsave(filename = file.path(plot_dir, "CanonicalMarkers_DotPlot.png"),
       plot = p_dot, width = 18, height = 5)



# DotPlot separado por categorias ............................
# Crear tabla con las categorías
gene_category_df <- bind_rows(lapply(names(gene_sets), function(cat) {
  data.frame(Feature = gene_sets[[cat]], Category = cat)
}))

# Obtener los datos del DotPlot
dp_data <- DotPlot(seurat_object, features = gene_category_df$Feature)$data

# Añadir la categoría
dp_data <- dp_data %>%
  left_join(gene_category_df, by = c("features.plot" = "Feature"))

# Ahora crear el ggplot directamente con facet
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
