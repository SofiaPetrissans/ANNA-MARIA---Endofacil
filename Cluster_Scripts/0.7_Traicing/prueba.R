# --------------------------------------------------------------------------------
# SCRIPT: Traicing:
#           - Ane tenia anotadas las subpoblaciones de EC en: (/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/1.Unenriched_samples/0.0_Previously_processed_by_ane/mt15/Data_Annotated_mt15_EC.rds)
#           - Sin embargo, hemos estado haciendo las anotaciones el objeto filtrado con un threshold de 20 ARN mitocondrial
#           - Por eso, estamos haciendo un treacing de estas células anotadas en el nuevo subset de EC del objeto filtrado con mt20
# AUTORA: Sofia Petrissans Moll
# FECHA: 08/07/2025
# --------------------------------------------------------------------------------
source("/ijc/USERS/spetrissans/ANNA_MARIA/1.Unenriched_samples/0.1_Scripts/0.0_Paths.R")
library(Seurat)
library(tidyverse)
library(harmony)
library(RColorBrewer)

# ...............................................................................
# PATHS
# ...............................................................................
path_mt15_EC <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.4_SubsetEndothelial/mt15/9.RemoveCluster/2.Integration/SubsetEndothelial_Harmony.rds"
path_mt20_EC <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.5_Annotation_EC/mt20/5.Annotation/Seurat_annotated.rds"

output_dir <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.7_Traicing/mt15/prueba"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ...............................................................................
# CARGA Y SUBSET
# ...............................................................................
seurat_mt15_EC <- readRDS(path_mt15_EC)
seurat_mt20_EC <- readRDS(path_mt20_EC)

#Obtener los niveles (ordenados como quieras)
layer_levels <- levels(seurat_mt20_EC$Annotation_Layer2)

#Crear una paleta de colores fija
n <- length(layer_levels)
palette_fixed <- setNames(brewer.pal(n, "Set2")[1:n], layer_levels)

# guardar el dimplot con las anotaciones del objeto anterior que queremos hacer el traicing
p1 <- DimPlot(seurat_mt20_EC, group.by = "Annotation_Layer2", label = FALSE, pt.size = 1.5, cols = palette_fixed) & NoAxes()
ggsave(plot = p1, filename = file.path(output_dir, "mt15_Layers_EC.png"), width = 10, height =8)


# 1. Extraer anotación original
tracing_annotation <- seurat_mt20_EC@meta.data[, "Annotation_Layer2", drop = FALSE]

# Crear versiones limpias de los nombres (sin sufijos)
colnames_mt20_clean <- sub("_\\d+$", "", colnames(seurat_mt15_EC))
rownames_mt15_clean <- sub("_\\d+$", "", rownames(tracing_annotation))

# Crear tabla de mapping para intersección
names(colnames_mt20_clean) <- colnames(seurat_mt15_EC)
names(rownames_mt15_clean) <- rownames(tracing_annotation)

# Intersección de nombres "limpios"
common_clean <- intersect(rownames_mt15_clean, colnames_mt20_clean)

# ¿Qué células originales corresponden?
cells_mt20_matched <- names(colnames_mt20_clean)[colnames_mt20_clean %in% common_clean]
cells_mt15_matched <- names(rownames_mt15_clean)[rownames_mt15_clean %in% common_clean]

# Subset y asignación segura
matched_annotation <- tracing_annotation[cells_mt15_matched, , drop = FALSE]
seurat_mt15_EC$Layer_EC_1_traced <- NA
seurat_mt15_EC$Layer_EC_1_traced[cells_mt20_matched] <- as.character(matched_annotation$Annotation_Layer2)

# Convertir a factor
seurat_mt15_EC$Layer_EC_1_traced <- factor(
  seurat_mt15_EC$Layer_EC_1_traced,
  levels = levels(tracing_annotation$Annotation_Layer2)
)


# 3. Visualizar
p2 <- DimPlot(seurat_mt15_EC, group.by = "Layer_EC_1_traced", label = FALSE, pt.size = 1.5, cols = palette_fixed) & NoAxes()
ggsave(plot = p2, filename = file.path(output_dir, "Tracing_Layer1.png"), width = 10, height =8)



# 4. Generar DimPlot individual para cada subpoblación trazada

# Crear un vector con las categorías
categorias <- levels(seurat_mt15_EC$Layer_EC_1_traced)

# Reemplazo para color gris (fondo)
color_gris <- "#d9d9d9"  # gris claro agradable

for (cat in categorias) {
  # Crear un vector de colores: todo gris excepto la categoría destacada
  cols_destacadas <- setNames(rep(color_gris, length(categorias)), categorias)
  cols_destacadas[cat] <- palette_fixed[cat]

  # Crear el plot
  p <- DimPlot(seurat_mt15_EC, group.by = "Layer_EC_1_traced", cols = cols_destacadas, pt.size = 1.5) & NoAxes() & NoLegend()

  # Guardar
  filename_plot <- paste0("Tracing_", gsub(" ", "_", cat), ".png")
  ggsave(plot = p, filename = file.path(output_dir, filename_plot), width = 9, height = 8)
}
