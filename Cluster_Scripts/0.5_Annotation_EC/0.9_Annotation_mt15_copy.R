# --------------------------------------------------------------------------------
# SCRIPT: Generar anotaciones: nueva capa, plots y propller
# AUTORA: Sofia Petrissans Moll
# FECHA: 27/06/2025
# --------------------------------------------------------------------------------


# =========================
# LIBRERÍAS
# =========================
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(purrr)

RESOLUTION <- "Harmony_Log_res.0.5"


# CARGAR OBJETO
MT_THRESHOLD <- "mt15"
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path,  MT_THRESHOLD, "9.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")
seurat_object <- readRDS(input_rds)


# PATH GUARDAR
path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, MT_THRESHOLD, "2.Annotation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# DOT PLOTS - Harmony_Log_res.0.5 ........................................................................................................
resolution_col <- "Harmony_Log_res.0.5"

gene_sets <- list(
  Mizonal = c("Aass", "Apoc1", "Ccnd1", "Cox7c", "Hint1", "Lyve1", "Ndufb1"),
  Pericentral = c("Alad", "Bmp2", "Clec1b", "Gatm", "Kit", "Lgals1", "Ptgs1", "Rab3b", "Tcim", "Thbd", "Wnt2"),
  Periportal = c("Acer2", "Adgrg6", "Dll4", "Ednrb", "Efnb1", "Efnb2", "Glul", "Jag1", "Lrg1", "Ltbp4"),
  Portal_area = c("Gja5", "Lmo7", "Ntn4"),
  Proliferative = c("Mki67", "Top2a", "Pcna", "Mcm5")
)

features <- unique(unlist(gene_sets))

seurat_object <- SetIdent(seurat_object, value = resolution_col)

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
ggsave(filename = file.path(output_dir, paste0("CanonicalMarkers_DotPlot_", resolution_col, ".png")),
       plot = dp, width = 18, height = 5)






# --------------------------------------------------------------------------------
# ANNOTATION_LAYER_LSEC_LVEC
# --------------------------------------------------------------------------------

RESOLUTION <- "Harmony_Log_res.0.5"
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)

# Annotation_Layer1
cluster_names <- c(
  "0" = "Periportal",
  "1" = "Pericentral",
  "2" = "Periportal",
  "3" = "Midzonal",
  "4" = "Periportal",
  "5" = "Midzonal"
)


Idents(seurat_object) <- RESOLUTION
seurat_object@meta.data$Annotation_Layer2 <- factor(cluster_names[as.character(Idents(seurat_object))])

table(Annotation = seurat_object$Annotation_Layer2, Cluster = Idents(seurat_object))

annotation.colors <- c(
  "Midzonal" = "#D5B3E5", 
  "Pericentral" = "#B9DDF1", 
  "Periportal" = "#FFD5C6",
  "Proliferative" = "#F2D6A0"
)

colors <- annotation.colors
layer_name <- "Annotation_Layer2"

DimPlot(seurat_object, reduction = "umap", group.by = layer_name, 
        label = FALSE, pt.size = 0.3, raster=FALSE) +
  scale_color_manual(values = colors) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right", legend.text = element_text(size = 12)) +  
  guides(colour = guide_legend(override.aes = list(size = 2)))
# Guardar
ggsave(filename = file.path(output_dir, paste0("DimPlot_", layer_name, ".png")), width = 7, height = 6)



# DOT PLOT ANOTACION ........................................................................................................

gene_sets <- list(
  "Mizonal" = c("Hamp", "Apoc1", "Fabp1", "Mt2a", "Lyve1", "Ccnd1"),
  "Pericentral" = c("Rbp1", "Bmp2", "Tcim", "Thbd", "Wnt2"),
  "Periportal" = c("Acer2", "Efnb2", "Meis1", "Zeb2", "Arhgap26", "Mbd5", "Ctnnd1", "Stab2"),
  "Proliferative" = c("Mki67", "Top2a", "Pcna")
)

features <- unique(unlist(gene_sets))
seurat_object <- SetIdent(seurat_object, value = layer_name)

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
ggsave(filename = file.path(output_dir, paste0("CanonicalMarkers_DotPlot_", layer_name, ".png")),
       plot = dp, width = 18, height = 5)




# HEAT MAP ANOTACION ........................................................................................................

# --- 1. Definir los conjuntos de genes  ---
gene_categories <- list(
  "Mizonal" = c("Hamp", "Apoc1", "Fabp1", "Mt2a", "Lyve1", "Ccnd1"),
  "Pericentral" = c("Rbp1", "Bmp2", "Tcim", "Thbd", "Wnt2"),
  "Periportal" = c("Acer2", "Efnb2", "Meis1", "Zeb2", "Arhgap26", "Mbd5", "Ctnnd1", "Stab2"),
  "Proliferative" = c("Mki67", "Top2a", "Pcna")
)

# --- 2. Combinar todos los genes en una lista única ---
marker_genes <- unlist(gene_categories, use.names = FALSE)

# --- 3. Calcular expresión promedio por cluster (Annotation_Layer2) ---
avg.exp <- AverageExpression(seurat_object, group.by = layer_name, return.seurat = FALSE)

# Verificar qué genes están en el objeto
genes_disponibles <- intersect(marker_genes, rownames(avg.exp$RNA))
if (length(genes_disponibles) == 0) {
  stop("Error: Ninguno de los genes de `marker_genes` está presente en `avg.exp$RNA`.")
}

# --- 4. Matriz de expresión ---
expression.matrix <- avg.exp$RNA[genes_disponibles, ]
expression.matrix.scaled <- t(scale(t(expression.matrix)))
expression.matrix.scaled[is.na(expression.matrix.scaled)] <- 0

# --- 5. Crear anotación de filas (genes) ---
gene.groups <- data.frame(
  Tipos_Celulares = rep(names(gene_categories), times = map_int(gene_categories, length))
)
rownames(gene.groups) <- make.unique(marker_genes)
gene.groups <- gene.groups[rownames(expression.matrix.scaled), , drop = FALSE]

# --- 6. Anotación de filas (colores) ---
left_annotation <- rowAnnotation(
  Tipos_Celulares = gene.groups$Tipos_Celulares,
  col = list(Tipos_Celulares = annotation.colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = "Categorias",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm")
  )
)

# --- 7. Colores para el heatmap ---
col_fun <- colorRamp2(
  c(min(expression.matrix.scaled), 0, max(expression.matrix.scaled)),
  c("#1874CD", "#EEEEE0", "#CD2626")
)

# Crear y guardar el heatmap
heatmap <- Heatmap(
  expression.matrix.scaled,
  name = "Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "Scaled Expression"),
  border = TRUE,
  column_title = paste("Canonical Markers by", layer_name),
  left_annotation = left_annotation,
  split = gene.groups$Zonation,
  gap = unit(2, "mm")
)

# Guardar
png(file.path(output_dir, paste0("Markers_Heatmap_", layer_name,".png")), width = 1600, height = 2000, res = 300)
draw(heatmap, heatmap_legend_side = "right")
dev.off()
