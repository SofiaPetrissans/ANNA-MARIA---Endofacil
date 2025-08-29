# -----------------------------------------------------
# SCRIPT: DEGs entre CHOP y vehicle en todo el endotelio refinado
  # Volver a hacer el análisis de DEGs entre CHOP vs vehicle, pero esta vez sobre todas las células del subset refinado de endotelio 
  # (sin dividir por subpoblaciones), para confirmar que las diferencias 
  # observadas originalmente no estaban sesgadas por las células ruidosas que has eliminado.
# -----------------------------------------------------

# LIBRERÍAS
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(dplyr)
library(tibble)
library(openxlsx)
library(ggplot2)


MT_THRESHOLD <- "mt15"
RESOLUCION <- "Harmony_Log_res.0.3"


# Paths
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, MT_THRESHOLD, "8.RemoveCluster", "2.Integration", "SubsetEndothelial_Harmony.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, MT_THRESHOLD, "3.DEGs_phenotype")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# CARGA Y PREPARACIÓN
# --------------------------------------------------------------------------------
seurat_object <- readRDS(input_rds)


# SET IDENT Y VERIFICACIÓN
table(seurat_object$Phenotype)  # Asegúrate que existan ambas condiciones
Idents(seurat_object) <- "Phenotype"


# HACER EL ANÁLISIS
deg_results <- FindMarkers(seurat_object, ident.1 = "CHOP", ident.2 = "vehicle",
                           logfc.threshold = 0.25, min.pct = 0.1)

deg_results <- deg_results %>% rownames_to_column("gene")

# GUARDAR RESULTADOS
write.xlsx(deg_results, file = file.path(output_dir, "DEGs_CHOP_vs_vehicle.xlsx"))

# BARPLOT: Número de genes up/down
deg_sig <- deg_results %>% filter(p_val_adj < 0.05)
n_up <- sum(deg_sig$avg_log2FC > 0)
n_down <- sum(deg_sig$avg_log2FC < 0)

deg_bar <- tibble(Direction = c("Up", "Down"), DEGs = c(n_up, n_down))

p <- ggplot(deg_bar, aes(x = Direction, y = DEGs, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Up" = "#B3E0A6", "Down" = "#F78D80")) +
  labs(title = "DEGs in Endothelial cells (CHOP vs vehicle)",
       y = "Number of significant DEGs (adj.p < 0.05)", x = "") +
  theme_classic(base_size = 13)

ggsave(filename = file.path(output_dir, "Barplot_DEGs_Endothelial_Total.png"),
       plot = p, width = 5, height = 5)
