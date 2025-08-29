# --------------------------------------------------------------------------------
# SCRIPT: Generar anotaciones: nueva capa, plots y propller
# AUTORA: Sofia Petrissans Moll
# FECHA: 27/06/2025
# --------------------------------------------------------------------------------

# LIBRARY
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.2_Annotation_All_Popultations/0.0_Paths.R")
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



MT_THRESHOLD <- "mt20"
RESOLUTION <- "Harmony_Log_res.0.3"


# Paths
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "6.Integration", MT_THRESHOLD, "Harmony.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, "3.Annotation_Layer", MT_THRESHOLD)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# CARGA Y PREPARACIÓN
# --------------------------------------------------------------------------------
seurat_object <- readRDS(input_rds)
seurat_object <- SetIdent(seurat_object, value = RESOLUTION)


# Annotation_Layer1
cluster_names <- c(
  "0" = "NK",
  "1" = "Endothelial",
  "2" = "Neutrophils",
  "3" = "Kupffer",
  "4" = "B Cells",
  "5" = "Endothelial",
  "6" = "Neutrophils",
  "7" = "T Cells",
  "8" = "HSC",
  "9" = "NK",
  "10" = "Endothelial",
  "11" = "Hepatocytes",
  "12" =  "Kupffer",
  "13" = "Monocytes",
  "14" = "Neutrophils",
  "15" = "B Cells",
  "16" = "Monocytes"
)

Idents(seurat_object) <- RESOLUTION
seurat_object@meta.data$Annotation_Layer1 <- factor(cluster_names[as.character(Idents(seurat_object))])

table(Annotation = seurat_object$Annotation_Layer1, Cluster = Idents(seurat_object))


Annotation_Colors_L1 <- c(
  "NK" = "#B29CA6",
  "Endothelial" = "#B5DDC3",
  "Neutrophils" = "#CEDD9A", 
  "Kupffer" = "#ff9f9b",
  "B Cells" = "#BDBAD7",
  "T Cells" = "#debb9b",
  "HSC" = "#fab0e4",
  "Hepatocytes" = "#BADCDE",
  "Monocytes" = "#E0A0A7"
)

colors <- Annotation_Colors_L1
layer_name <- "Annotation_Layer1"

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
  count(Phenotype, Annotation_Layer1) %>%
  group_by(Phenotype) %>%
  mutate(Proportion = n / sum(n) * 100)

# Barplot apilado
ggplot(df_barplot, aes(x = Phenotype, y = Proportion, fill = Annotation_Layer1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Annotation_Colors_L1) +
  labs(y = "Percentage of Cells", x = "Condition", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# Guardar
ggsave(file.path(output_dir, "Barplot_Stacked_Percentages_By_Phenotype.png"), width = 6, height = 4.5)



# --------------------------------------------------------------------------------
# PROPELLER ANALYSIS
# --------------------------------------------------------------------------------

# Crear tabla de metadatos
cells_info <- data.frame(
  clusters = seurat_object$Annotation_Layer1,
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
  CellType = anno_simple(rownames(res_selected), col = Annotation_Colors_L1),
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

pdf(file.path(output_dir, "Heatmap_PropMean_CellType.pdf"), width = 5.5, height = 5)
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
  scale_fill_manual(values = Annotation_Colors_L1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Proportion", x = "")

ggsave(file.path(output_dir, "BarPlot_PropMean_CellTypes.png"), width = 7, height = 4)



# ...............................................................................
# DOT PLOT
# ...............................................................................
seurat_object <- SetIdent(seurat_object, value = layer_name)

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
ggsave(filename = file.path(output_dir, "CanonicalMarkers_DotPlot_sep.png"),
       plot = dp, width = 18, height = 5)



# --------------------------------------------------------------------------------
# SAVE SEURAT_OBJECT
# --------------------------------------------------------------------------------
saveRDS(seurat_object, file.path(output_dir, "Seurat_annotated.rds"))



