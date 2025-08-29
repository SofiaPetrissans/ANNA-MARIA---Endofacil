# --------------------------------------------------------------------------------
# SCRIPT: Propeller analysis using Seurat clusters
# AUTORA: Sofia Petrissans Moll
# FECHA: 26.06.2025
# --------------------------------------------------------------------------------
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.4_SubsetEndothelial/0.0_Paths.R")
library(Seurat)
library(tidyverse)
library(speckle)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

# --------------------------------------------------------------------------------
# CONFIGURACIÓN
# --------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript Propeller_Clusters.R <MT_THRESHOLD> <RESOLUTION>")
}

MT_THRESHOLD <- args[1]   # Ej: mt15
RESOLUTION <- args[2]     # Ej: Harmony_Log_res.0.3


# Paths
base_path <- paths_AnnaMaria$Path_seurat_object
input_rds <- file.path(base_path, MT_THRESHOLD, "4.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original , MT_THRESHOLD, "6.Propeller", paste0("res", RESOLUTION))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# CARGA Y PREPARACIÓN
# --------------------------------------------------------------------------------
seurat_object <- readRDS(input_rds)
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)

getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
cluster_levels <- levels(seurat_object)
colors_clusters <- getPalette(length(cluster_levels))
names(colors_clusters) <- cluster_levels

# Propeller
cells_info <- data.frame(
  clusters = seurat_object@active.ident,
  samples = seurat_object$Sample,
  group = seurat_object$Phenotype
)
saveRDS(cells_info, file.path(output_dir, "cells_info.rds"))

res_propeller <- propeller(
  clusters = cells_info$clusters,
  sample = cells_info$samples,
  group = cells_info$group
)
saveRDS(res_propeller, file.path(output_dir, "Res_Propeller.rds"))

# --------------------------------------------------------------------------------
# HEATMAP
# --------------------------------------------------------------------------------
#res_selected <- res_propeller %>%
#  select(CellType, starts_with("PropMean.")) %>%
#  column_to_rownames("CellType")
#colnames(res_selected) <- str_replace(colnames(res_selected), "PropMean.", "")


# --------------------------------------------------------------------------------
# HEATMAP CON ANOTACIÓN DE CLUSTERS (color como CHOP/vehicle)
# --------------------------------------------------------------------------------

# Paleta y función de colores por cluster
getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
generateClusterColors <- function(object, group_col) {
  cluster_levels <- levels(factor(object@meta.data[[group_col]]))
  n_clusters <- length(cluster_levels)
  colors <- getPalette(n_clusters)
  names(colors) <- cluster_levels
  return(colors)
}

cluster_colors <- generateClusterColors(seurat_object, RESOLUTION)

# Preparar matriz
res_selected <- res_propeller %>%
  select(starts_with("PropMean."))
rownames(res_selected) <- paste0("Cluster_", rownames(res_propeller))  # añade prefijo opcional
colnames(res_selected) <- str_replace(colnames(res_selected), "PropMean.", "")

# Crear vector de anotación de color por cluster
cluster_ids <- str_replace(rownames(res_selected), "Cluster_", "")  # extraer ID numérico
cluster_colors_named <- setNames(cluster_colors[cluster_ids], rownames(res_selected))

# Anotación de columnas (Phenotype)
col_phenotype <- c("CHOP" = "#99540F", "vehicle" = "#51A3CC")
col_ha <- HeatmapAnnotation(
  Phenotype = colnames(res_selected),
  col = list(Phenotype = col_phenotype)
)

# Anotación de filas (Clusters)
row_ha <- rowAnnotation(
  Cluster = rownames(res_selected),
  col = list(Cluster = cluster_colors_named),
  show_annotation_name = FALSE
)

# Anotación de P-valor
points_ha <- rowAnnotation(
  PVal = anno_points(
    -log10(res_propeller$P.Value),
    pch = 16,
    gp = gpar(col = "#1E88E5"),
    size = unit(2, "mm"),
    ylim = c(0, max(-log10(res_propeller$P.Value), na.rm = TRUE))
  )
)

# Combinar anotaciones de fila
left_ha <- rowAnnotation(
  Cluster = anno_simple(cluster_ids, col = cluster_colors),
  PVal = anno_points(
    -log10(res_propeller$P.Value),
    gp = gpar(col = "#1E88E5"),
    size = unit(2, "mm"),
    ylim = c(0, max(-log10(res_propeller$P.Value), na.rm = TRUE))
  )
)

# Colores para la matriz
col_fun <- colorRamp2(c(0, 0.5), c("white", "#01579B"))

# HEATMAP FINAL
ht <- Heatmap(res_selected,
              name = "Prop.Mean",
              col = col_fun,
              top_annotation = col_ha,
              left_annotation = left_ha,
              show_row_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10))

pdf(file.path(output_dir, "Heatmap_PropMean_ClusterAnnotation.pdf"), width = 5.5, height = 5)
draw(ht, show_annotation_legend = TRUE)
dev.off()


# --------------------------------------------------------------------------------
# BARPLOT
# --------------------------------------------------------------------------------
df_long <- res_selected %>%
  rownames_to_column("Cluster") %>%
  pivot_longer(cols = -Cluster, names_to = "Phenotype", values_to = "Proportion")

ggplot(df_long, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = colors_clusters) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Proportion")

ggsave(file.path(output_dir, "BarPlot_PropMean.png"), width = 7, height = 4)

# --------------------------------------------------------------------------------
# BARPLOT CON COLORES DE CLUSTERS
# --------------------------------------------------------------------------------

# Volver a usar rownames si ya están en formato "Cluster_X"
df_long <- res_selected %>%
  rownames_to_column("Cluster") %>%
  pivot_longer(cols = -Cluster, names_to = "Phenotype", values_to = "Proportion")

# Eliminar "Cluster_" si quieres usar solo los números como clave para los colores
df_long$Cluster_ID <- str_replace(df_long$Cluster, "Cluster_", "")

# Asignar colores usando cluster_colors
cluster_colors_barplot <- cluster_colors[df_long$Cluster_ID]
names(cluster_colors_barplot) <- df_long$Cluster

# Hacer plot con colores personalizados
ggplot(df_long, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = cluster_colors_barplot) +  # aquí se aplica la paleta Set3
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Proportion", x = "")

ggsave(file.path(output_dir, "BarPlot_PropMean_ClusterColors.png"), width = 7, height = 4)


cat("Análisis completado. Resultados en:", output_dir, "\n")