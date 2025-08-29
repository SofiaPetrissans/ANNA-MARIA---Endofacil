# --------------------------------------------------------------------------------
# SCRIPT: IDENTIFICAR CLÚSTERES EN DIVISIÓN CELULAR 
    # Detectar si un clúster representa células en división (o ciclo celular activo) 
    # Células en G1: no están proliferando.
    # Células en S o G2M: están activamente dividiéndose.
    # Si un clúster tiene alta proporción de S/G2M, es probable que sea una población proliferativa


# AUTORA: Sofia Petrissans Moll
# FECHA: 11/07/2025
# --------------------------------------------------------------------------------

# ...............................................................................
# LIBRERIAS
# ...............................................................................
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(stringr)
library(reshape2) 

# =========================
# 1. PATHS
# =========================
MT_THRESHOLD <- "mt15"
RESOLUCION <- "Harmony_Log_res.0.3"

path.guardar_original <- paths_AnnaMaria$Path_guardar
path.guardar <- file.path(path.guardar_original, MT_THRESHOLD, "1.MarkersDivision")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path,  MT_THRESHOLD, "8.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")

# Colores clusteres
getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
generateClusterColors <- function(object, group_col) {
  cluster_levels <- levels(factor(object@meta.data[[group_col]]))
  n_clusters <- length(cluster_levels)
  colors <- getPalette(n_clusters)
  names(colors) <- cluster_levels
  return(colors)
}


# =========================
# 2. CARGAR y ANOTACION
# =========================
seurat_object <- readRDS(input_rds)

# Generar colores para tu agrupación 
cluster_colors <- generateClusterColors(seurat_object, RESOLUCION)

# Genes de ciclo celular
cc.genes <- Seurat::cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convertir a formato de símbolo de ratón (solo primera letra mayúscula)
s.genes.mouse <- str_to_title(tolower(s.genes))
g2m.genes.mouse <- str_to_title(tolower(g2m.genes))

# Asignar scores de ciclo celular
seurat_object <- CellCycleScoring(
  seurat_object,
  s.features = s.genes.mouse,
  g2m.features = g2m.genes.mouse,
  set.ident = TRUE  # Esto pone la fase como identidad activa (opcional)
)

# =========================
# 3. TABLAS DE FASE POR CLÚSTER
# =========================
cluster_vector <- pull(seurat_object@meta.data, RESOLUCION)
phase_vector <- seurat_object$Phase

tab_counts <- table(cluster_vector, phase_vector)
tab_props  <- round(prop.table(tab_counts, margin = 1) * 100, 2)

# Convertir las tablas a data.frame antes de guardar
tab_counts_df <- as.data.frame.matrix(tab_counts)
tab_props_df  <- as.data.frame.matrix(tab_props)

# Guardar
write.csv(tab_counts_df, file = file.path(path.guardar, "Phase_Counts_by_Cluster.csv"))
write.csv(tab_props_df,  file = file.path(path.guardar, "Phase_Percent_by_Cluster.csv"))


# =========================
# 4. VISUALIZACIONES
# =========================
# UMAP coloreado por fase
p1 <- DimPlot(seurat_object, group.by = "Phase") + ggtitle("Ciclo celular") & NoAxes()
ggsave(plot = p1, filename=file.path(path.guardar, "DimPlot_Ciclo_Celular.png"), width=8, height=7)


# Crear por separado los violin plots
vln_s <- VlnPlot(seurat_object, features = "S.Score", group.by = RESOLUCION, pt.size = 0.1) +
  scale_fill_manual(values = cluster_colors) +
  ggtitle("S.Score (fase S)") +
  theme(legend.position = "none")

vln_g2m <- VlnPlot(seurat_object, features = "G2M.Score", group.by = RESOLUCION, pt.size = 0.1) +
  scale_fill_manual(values = cluster_colors) +
  ggtitle("G2M.Score (fase G2/M)") +
  theme(legend.position = "none")

# Unir con patchwork
p2 <- vln_s | vln_g2m

# Guardar
ggsave(plot = p2, filename = file.path(path.guardar, "VlnPlot_Ciclo_Celular.png"),
       width = 12, height = 6)

cat("Número de genes S encontrados:\n")
print(length(intersect(s.genes.mouse, rownames(seurat_object))))

cat("Número de genes G2M encontrados:\n")
print(length(intersect(g2m.genes.mouse, rownames(seurat_object))))

# =========================
# 5. BARPLOT DE PROPORCIÓN DE FASES POR CLÚSTER
# =========================
# Leer la tabla de proporciones
tab_props_df <- read.csv(file.path(path.guardar, "Phase_Percent_by_Cluster.csv"), row.names = 1)

# Agregar columna de clúster
tab_props_df$Cluster <- rownames(tab_props_df)

# Convertir a formato largo
df_long <- melt(tab_props_df, id.vars = "Cluster", variable.name = "Phase", value.name = "Proportion")

# Asegurar orden consistente de fases
df_long$Phase <- factor(df_long$Phase, levels = c("G1", "S", "G2M"))

# Barplot apilado
p_bar <- ggplot(df_long, aes(x = Cluster, y = Proportion, fill = Phase)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Proportion of cells per cell cycle phase and cluster",
       y = "% of cells", x = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Guardar
ggsave(filename = file.path(path.guardar, "Barplot_Ciclo_Celular_por_Cluster.png"),
       plot = p_bar, width = 8, height = 6)


p <- FeaturePlot(seurat_object, features=c("Mki67", "Top2a", "Ccnb1"), order = TRUE, raster = FALSE) & NoAxes()
ggsave(filename = file.path(path.guardar, "FeaturePlot_Prolif.png"), plot = p, width = 8, height = 7)