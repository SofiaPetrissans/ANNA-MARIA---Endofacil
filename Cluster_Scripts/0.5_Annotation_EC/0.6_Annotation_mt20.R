# --------------------------------------------------------------------------------
# SCRIPT: Generar anotaciones: nueva capa, plots y propller
# AUTORA: Sofia Petrissans Moll
# FECHA: 27/06/2025
# --------------------------------------------------------------------------------


# =========================
# LIBRERÍAS
# =========================
source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(speckle)
library(RColorBrewer)


RESOLUTION <- "Harmony_Log_res.0.5"


# CARGAR OBJETO
MT_THRESHOLD <- "mt20"
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path,  MT_THRESHOLD, "8.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")
seurat_object <- readRDS(input_rds)

# Reordenar niveles de Phenotype
seurat_object$Phenotype <- factor(seurat_object$Phenotype, levels = c("vehicle", "CHOP"))

# PATH GUARDAR
path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, MT_THRESHOLD, "5.Annotation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# CARGA Y PREPARACIÓN
# --------------------------------------------------------------------------------
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)


# --------------------------------------------------------------------------------
# ANNOTATION_LAYER_LSEC_LVEC
# --------------------------------------------------------------------------------
# Annotation_Layer1
cluster_names <- c(
  "0" = "LSEC",
  "1" = "LSEC",
  "2" = "LSEC",
  "3" = "LSEC",
  "4" = "LSEC",
  "5" = "LSEC",
  "6" = "LVEC"
)

Idents(seurat_object) <- RESOLUTION
seurat_object@meta.data$Annotation_LSEC_LVEC <- factor(cluster_names[as.character(Idents(seurat_object))])

table(Annotation = seurat_object$Annotation_LSEC_LVEC, Cluster = Idents(seurat_object))


Annotation_Colors_LSEC_LVEC <- c(
  "LSEC" = "#B2DE81",
  "LVEC" = "#D8B2C6"
)

colors <- Annotation_Colors_LSEC_LVEC
layer_name <- "Annotation_LSEC_LVEC"

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



DimPlot(seurat_object, reduction = "umap", group.by = layer_name, 
        split.by = "Phenotype",
        label = FALSE, pt.size = 0.3, raster=FALSE) +
  scale_color_manual(values = colors) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right", legend.text = element_text(size = 12)) +  
  guides(colour = guide_legend(override.aes = list(size = 2)))
# Guardar
ggsave(filename = file.path(output_dir, paste0("DimPlot_", layer_name, "_SplitBy_Pheno.png")), width = 11, height = 6)


# BARRAS APILADAS
# Crear tabla resumen: número de células por condición y tipo celular
df_barplot <- seurat_object@meta.data %>%
  count(Phenotype, Annotation_LSEC_LVEC) %>%
  group_by(Phenotype) %>%
  mutate(Proportion = n / sum(n) * 100)

# Barplot apilado
ggplot(df_barplot, aes(x = Phenotype, y = Proportion, fill = Annotation_LSEC_LVEC)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Annotation_Colors_LSEC_LVEC) +
  labs(y = "Percentage of Cells", x = "Condition", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# Guardar
ggsave(file.path(output_dir, "Barplot_Stacked_Percentages_By_Phenotype.png"), width = 6, height = 8)


# --------------------------------------------------------------------------------
# PROPELLER ANALYSIS
# --------------------------------------------------------------------------------

# Crear tabla de metadatos
cells_info <- data.frame(
  clusters = seurat_object$Annotation_LSEC_LVEC,
  samples = seurat_object$Sample,
  group = seurat_object$Phenotype
)
saveRDS(cells_info, file.path(output_dir, "cells_info.rds"))

# Ejecutar propeller
res_propeller <- propeller(
  clusters = cells_info$clusters,
  sample = cells_info$samples,
  group = cells_info$group
)

# --------------------------------------------------------------------------------
# HEATMAP CON ANOTACIÓN DE CELULAS
# --------------------------------------------------------------------------------

# Preparar matriz
res_selected <- res_propeller %>%
  select(starts_with("PropMean."))
rownames(res_selected) <- rownames(res_propeller)  # usa los nombres de las poblaciones
colnames(res_selected) <- str_replace(colnames(res_selected), "PropMean.", "")

# Anotaciones
col_phenotype <- c("CHOP" = "#99540F", "vehicle" = "#51A3CC")

# Columnas
col_ha <- HeatmapAnnotation(
  Phenotype = colnames(res_selected),
  col = list(Phenotype = col_phenotype)
)

# Filas
left_ha <- rowAnnotation(
  CellType = anno_simple(rownames(res_selected), col = Annotation_Colors_LSEC_LVEC),
  PVal = anno_points(
    -log10(res_propeller$P.Value),
    gp = gpar(col = "#1E88E5"),
    size = unit(2, "mm"),
    ylim = c(0, max(-log10(res_propeller$P.Value), na.rm = TRUE))
  )
)

# Colores para la matriz
col_fun <- colorRamp2(c(0, 0.5), c("white", "#01579B"))

# Heatmap
ht <- Heatmap(res_selected,
              name = "Prop.Mean",
              col = col_fun,
              top_annotation = col_ha,
              left_annotation = left_ha,
              show_row_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10))

pdf(file.path(output_dir, "Heatmap_PropMean_CellType_LSEC_LVEC.pdf"), width = 5.5, height = 5)
draw(ht, show_annotation_legend = TRUE)
dev.off()

# --------------------------------------------------------------------------------
# BARPLOT CON COLORES DE CELULAS
# --------------------------------------------------------------------------------

df_long <- res_selected %>%
  rownames_to_column("CellType") %>%
  pivot_longer(cols = -CellType, names_to = "Phenotype", values_to = "Proportion")

ggplot(df_long, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = Annotation_Colors_LSEC_LVEC) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Proportion", x = "")

ggsave(file.path(output_dir, "BarPlot_PropMean_CellTypes_LVEC_LSEC.png"), width = 7, height = 4)





# --------------------------------------------------------------------------------
# ANNOTATION_LAYER2
# --------------------------------------------------------------------------------
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)

# Annotation_Layer2
cluster_names <- c(
  "0" = "Midzonal",
  "1" = "Pericentral",
  "2" = "Pericentral",
  "3" = "Periportal",
  "4" = "Midzonal",
  "5" = "Periportal",
  "6" = "Portal area"
)

Idents(seurat_object) <- RESOLUTION
seurat_object@meta.data$Annotation_Layer2 <- factor(cluster_names[as.character(Idents(seurat_object))])

table(Annotation = seurat_object$Annotation_Layer2, Cluster = Idents(seurat_object))


Annotation_Colors_Layer2 <- c(
  "Pericentral" = "#EDF8BC",
  "Midzonal" = "#B6E4B3",
  "Periportal" = "#63C3BF",
  "Portal area" = "#3A95B1"
)

colors <- Annotation_Colors_Layer2
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



DimPlot(seurat_object, reduction = "umap", group.by = layer_name, 
        split.by = "Phenotype",
        label = FALSE, pt.size = 0.3, raster=FALSE) +
  scale_color_manual(values = colors) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right", legend.text = element_text(size = 12)) +  
  guides(colour = guide_legend(override.aes = list(size = 2)))
# Guardar
ggsave(filename = file.path(output_dir, paste0("DimPlot_", layer_name, "_SplitBy_Pheno.png")), width = 11, height = 6)


# BARRAS APILADAS
# Crear tabla resumen: número de células por condición y tipo celular
df_barplot <- seurat_object@meta.data %>%
  count(Phenotype, Annotation_Layer2) %>%
  group_by(Phenotype) %>%
  mutate(Proportion = n / sum(n) * 100)

# Barplot apilado
ggplot(df_barplot, aes(x = Phenotype, y = Proportion, fill = Annotation_Layer2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Annotation_Colors_Layer2) +
  labs(y = "Percentage of Cells", x = "Condition", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# Guardar
ggsave(file.path(output_dir, "Barplot_Stacked_Percentages_By_Phenotype_Endo.png"), width = 6, height = 8)

# --------------------------------------------------------------------------------
# PROPELLER ANALYSIS
# --------------------------------------------------------------------------------

# Crear tabla de metadatos
cells_info <- data.frame(
  clusters = seurat_object$Annotation_Layer2,
  samples = seurat_object$Sample,
  group = seurat_object$Phenotype
)
saveRDS(cells_info, file.path(output_dir, "cells_info.rds"))

# Ejecutar propeller
res_propeller <- propeller(
  clusters = cells_info$clusters,
  sample = cells_info$samples,
  group = cells_info$group
)

# --------------------------------------------------------------------------------
# HEATMAP CON ANOTACIÓN DE CELULAS
# --------------------------------------------------------------------------------

# Preparar matriz
res_selected <- res_propeller %>%
  select(starts_with("PropMean."))
rownames(res_selected) <- rownames(res_propeller)  # usa los nombres de las poblaciones
colnames(res_selected) <- str_replace(colnames(res_selected), "PropMean.", "")

# Anotaciones
col_phenotype <- c("CHOP" = "#99540F", "vehicle" = "#51A3CC")

# Columnas
col_ha <- HeatmapAnnotation(
  Phenotype = colnames(res_selected),
  col = list(Phenotype = col_phenotype)
)

# Filas
left_ha <- rowAnnotation(
  CellType = anno_simple(rownames(res_selected), col = Annotation_Colors_Layer2),
  PVal = anno_points(
    -log10(res_propeller$P.Value),
    gp = gpar(col = "#1E88E5"),
    size = unit(2, "mm"),
    ylim = c(0, max(-log10(res_propeller$P.Value), na.rm = TRUE))
  )
)

# Colores para la matriz
col_fun <- colorRamp2(c(0, 0.5), c("white", "#01579B"))

# Heatmap
ht <- Heatmap(res_selected,
              name = "Prop.Mean",
              col = col_fun,
              top_annotation = col_ha,
              left_annotation = left_ha,
              show_row_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10))

pdf(file.path(output_dir, "Heatmap_PropMean_CellType_Endo.pdf"), width = 5.5, height = 5)
draw(ht, show_annotation_legend = TRUE)
dev.off()

# --------------------------------------------------------------------------------
# BARPLOT CON COLORES DE CELULAS
# --------------------------------------------------------------------------------

df_long <- res_selected %>%
  rownames_to_column("CellType") %>%
  pivot_longer(cols = -CellType, names_to = "Phenotype", values_to = "Proportion")

ggplot(df_long, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = Annotation_Colors_Layer2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Proportion", x = "")

ggsave(file.path(output_dir, "BarPlot_PropMean_CellTypes_Endo.png"), width = 7, height = 4)






# --------------------------------------------------------------------------------
# ANNOTATION_LAYER_LSEC_LEVELS 
# --------------------------------------------------------------------------------
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)

# Annotation_Layer1
cluster_names <- c(
  "0" = "LSEC_2",
  "1" = "LSEC_3",
  "2" = "LSEC_3",
  "3" = "LSEC_1",
  "4" = "LSEC_2",
  "5" = "LSEC_1",
  "6" = "LVEC"
)


Idents(seurat_object) <- RESOLUTION
seurat_object@meta.data$Annotation_Layer3 <- factor(cluster_names[as.character(Idents(seurat_object))])

table(Annotation = seurat_object$Annotation_Layer3, Cluster = Idents(seurat_object))


Annotation_Colors_Layer3 <- c(
  "LSEC_1" = "#66A61E",
  "LSEC_2" = "#B2DE81",
  "LSEC_3" = "#1B9E77" ,
  "LVEC" = "#D8B2C6"
)

colors <- Annotation_Colors_Layer3
layer_name <- "Annotation_Layer3"

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



DimPlot(seurat_object, reduction = "umap", group.by = layer_name, 
        split.by = "Phenotype",
        label = FALSE, pt.size = 0.3, raster=FALSE) +
  scale_color_manual(values = colors) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right", legend.text = element_text(size = 12)) +  
  guides(colour = guide_legend(override.aes = list(size = 2)))
# Guardar
ggsave(filename = file.path(output_dir, paste0("DimPlot_", layer_name, "_SplitBy_Pheno.png")), width = 11, height = 6)


# BARRAS APILADAS
# Crear tabla resumen: número de células por condición y tipo celular
df_barplot <- seurat_object@meta.data %>%
  count(Phenotype, Annotation_Layer3) %>%
  group_by(Phenotype) %>%
  mutate(Proportion = n / sum(n) * 100)

# Barplot apilado
ggplot(df_barplot, aes(x = Phenotype, y = Proportion, fill = Annotation_Layer3)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Annotation_Colors_Layer3) +
  labs(y = "Percentage of Cells", x = "Condition", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# Guardar
ggsave(file.path(output_dir, "Barplot_Stacked_Percentages_By_Phenotype_L3.png"), width = 6, height = 4.5)



# --------------------------------------------------------------------------------
# PROPELLER ANALYSIS
# --------------------------------------------------------------------------------

# Crear tabla de metadatos
cells_info <- data.frame(
  clusters = seurat_object$Annotation_Layer3,
  samples = seurat_object$Sample,
  group = seurat_object$Phenotype
)
saveRDS(cells_info, file.path(output_dir, "cells_info.rds"))

# Ejecutar propeller
res_propeller <- propeller(
  clusters = cells_info$clusters,
  sample = cells_info$samples,
  group = cells_info$group
)

# --------------------------------------------------------------------------------
# HEATMAP CON ANOTACIÓN DE CELULAS
# --------------------------------------------------------------------------------

# Preparar matriz
res_selected <- res_propeller %>%
  select(starts_with("PropMean."))
rownames(res_selected) <- rownames(res_propeller)  # usa los nombres de las poblaciones
colnames(res_selected) <- str_replace(colnames(res_selected), "PropMean.", "")

# Anotaciones
col_phenotype <- c("CHOP" = "#99540F", "vehicle" = "#51A3CC")

# Columnas
col_ha <- HeatmapAnnotation(
  Phenotype = colnames(res_selected),
  col = list(Phenotype = col_phenotype)
)

# Filas
left_ha <- rowAnnotation(
  CellType = anno_simple(rownames(res_selected), col = Annotation_Colors_Layer3),
  PVal = anno_points(
    -log10(res_propeller$P.Value),
    gp = gpar(col = "#1E88E5"),
    size = unit(2, "mm"),
    ylim = c(0, max(-log10(res_propeller$P.Value), na.rm = TRUE))
  )
)

# Colores para la matriz
col_fun <- colorRamp2(c(0, 0.5), c("white", "#01579B"))

# Heatmap
ht <- Heatmap(res_selected,
              name = "Prop.Mean",
              col = col_fun,
              top_annotation = col_ha,
              left_annotation = left_ha,
              show_row_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10))

pdf(file.path(output_dir, "Heatmap_PropMean_CellType_L3.pdf"), width = 5.5, height = 5)
draw(ht, show_annotation_legend = TRUE)
dev.off()

# --------------------------------------------------------------------------------
# BARPLOT CON COLORES DE CELULAS
# --------------------------------------------------------------------------------

df_long <- res_selected %>%
  rownames_to_column("CellType") %>%
  pivot_longer(cols = -CellType, names_to = "Phenotype", values_to = "Proportion")

ggplot(df_long, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = Annotation_Colors_Layer3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Proportion", x = "")

ggsave(file.path(output_dir, "BarPlot_PropMean_CellTypes_L3.png"), width = 7, height = 4)




# --------------------------------------------------------------------------------
# ANNOTATION_LAYER_LSEC_LEVELS 
# --------------------------------------------------------------------------------
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)

# Annotation_Layer1
cluster_names <- c(
  "0" = "LSEC_4",
  "1" = "LSEC_5",
  "2" = "LSEC_6",
  "3" = "LSEC_1",
  "4" = "LSEC_3",
  "5" = "LSEC_2",
  "6" = "LVEC"
)


Idents(seurat_object) <- RESOLUTION
seurat_object@meta.data$Annotation_Layer4 <- factor(cluster_names[as.character(Idents(seurat_object))])

table(Annotation = seurat_object$Annotation_Layer4, Cluster = Idents(seurat_object))


Annotation_Colors_Layer4 <- c(
  "LSEC_1" = "#A6CEE3",   
  "LSEC_2" = "#B2DF8A",   
  "LSEC_3" = "#FDBF6F",   
  "LSEC_4" = "#CAB2D6",   
  "LSEC_5" = "#FFCCBC",   
  "LSEC_6" = "#41B7C4",   
  "LVEC"   = "#D8B2C6"
)

colors <- Annotation_Colors_Layer4
layer_name <- "Annotation_Layer4"

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



DimPlot(seurat_object, reduction = "umap", group.by = layer_name, 
        split.by = "Phenotype",
        label = FALSE, pt.size = 0.3, raster=FALSE) +
  scale_color_manual(values = colors) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right", legend.text = element_text(size = 12)) +  
  guides(colour = guide_legend(override.aes = list(size = 2)))
# Guardar
ggsave(filename = file.path(output_dir, paste0("DimPlot_", layer_name, "_SplitBy_Pheno.png")), width = 11, height = 6)


# BARRAS APILADAS
# Crear tabla resumen: número de células por condición y tipo celular
df_barplot <- seurat_object@meta.data %>%
  count(Phenotype, Annotation_Layer4) %>%
  group_by(Phenotype) %>%
  mutate(Proportion = n / sum(n) * 100)

# Barplot apilado
ggplot(df_barplot, aes(x = Phenotype, y = Proportion, fill = Annotation_Layer4)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Annotation_Colors_Layer4) +
  labs(y = "Percentage of Cells", x = "Condition", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# Guardar
ggsave(file.path(output_dir, "Barplot_Stacked_Percentages_By_Phenotype_L4.png"), width = 6, height = 8)


# --------------------------------------------------------------------------------
# PROPELLER ANALYSIS
# --------------------------------------------------------------------------------

# Crear tabla de metadatos
cells_info <- data.frame(
  clusters = seurat_object$Annotation_Layer4,
  samples = seurat_object$Sample,
  group = seurat_object$Phenotype
)
saveRDS(cells_info, file.path(output_dir, "cells_info.rds"))

# Ejecutar propeller
res_propeller <- propeller(
  clusters = cells_info$clusters,
  sample = cells_info$samples,
  group = cells_info$group
)

# --------------------------------------------------------------------------------
# HEATMAP CON ANOTACIÓN DE CELULAS
# --------------------------------------------------------------------------------

# Preparar matriz
res_selected <- res_propeller %>%
  select(starts_with("PropMean."))
rownames(res_selected) <- rownames(res_propeller)  # usa los nombres de las poblaciones
colnames(res_selected) <- str_replace(colnames(res_selected), "PropMean.", "")

# Anotaciones
col_phenotype <- c("CHOP" = "#99540F", "vehicle" = "#51A3CC")

# Columnas
col_ha <- HeatmapAnnotation(
  Phenotype = colnames(res_selected),
  col = list(Phenotype = col_phenotype)
)

# Filas
left_ha <- rowAnnotation(
  CellType = anno_simple(rownames(res_selected), col = Annotation_Colors_Layer4),
  PVal = anno_points(
    -log10(res_propeller$P.Value),
    gp = gpar(col = "#1E88E5"),
    size = unit(2, "mm"),
    ylim = c(0, max(-log10(res_propeller$P.Value), na.rm = TRUE))
  )
)

# Colores para la matriz
col_fun <- colorRamp2(c(0, 0.5), c("white", "#01579B"))

# Heatmap
ht <- Heatmap(res_selected,
              name = "Prop.Mean",
              col = col_fun,
              top_annotation = col_ha,
              left_annotation = left_ha,
              show_row_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10))

pdf(file.path(output_dir, "Heatmap_PropMean_CellType_L4.pdf"), width = 5.5, height = 5)
draw(ht, show_annotation_legend = TRUE)
dev.off()

# --------------------------------------------------------------------------------
# BARPLOT CON COLORES DE CELULAS
# --------------------------------------------------------------------------------

df_long <- res_selected %>%
  rownames_to_column("CellType") %>%
  pivot_longer(cols = -CellType, names_to = "Phenotype", values_to = "Proportion")

ggplot(df_long, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = Annotation_Colors_Layer4) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Proportion", x = "")

ggsave(file.path(output_dir, "BarPlot_PropMean_CellTypes_L4.png"), width = 7, height = 4)



# --------------------------------------------------------------------------------
# SAVE SEURAT_OBJECT
# --------------------------------------------------------------------------------
saveRDS(seurat_object, file.path(output_dir, "Seurat_annotated.rds"))

