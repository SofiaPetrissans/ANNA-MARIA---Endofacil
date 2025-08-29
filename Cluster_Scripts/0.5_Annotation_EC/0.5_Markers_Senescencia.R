# LIBRARY
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(AUCell)
library(openxlsx)
library(patchwork)

MT_THRESHOLD <- "mt15"
RESOLUCION <- "Harmony_Log_res.0.3"


# Paths
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, MT_THRESHOLD, "8.RemoveCluster", "2.Integration", "SubsetEndothelial_Harmony.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
path.guardar <- file.path(path.guardar_original, MT_THRESHOLD, "2.Senescence")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# CARGA Y PREPARACIÓN
# --------------------------------------------------------------------------------
seurat_object <- readRDS(input_rds)


senescence_core <- c("Cdkn1a","Ctnnb1","Cxcl16", "Hmgb1",
                     "Igfbp4", "Igfbp7", "Il1a", "Mmp14", "Kitl",  "Tnfrsf1a",  
                     "Ppp1ca", "Smurf2")
fp <- FeaturePlot(seurat_object, features = senescence_core) & NoAxes()
ggsave(filename = file.path(path.guardar, "Senescence_mix.png"), plot = fp, width = 14, height = 12)

seurat_object <- SetIdent(seurat_object, value = RESOLUCION)
vlp <- VlnPlot(seurat_object, features = senescence_core, split.by = "Phenotype") + theme(legend.position = "right")
ggsave(filename = file.path(path.guardar, "VlnPlot_Senescence_mix.png"), plot = vlp, width = 18, height = 16)



# Asegurar que los clusters estén como identidades
seurat_object <- SetIdent(seurat_object, value = RESOLUCION)

# Genes de interés
senescence_core <- c("Cdkn2a", "Bmi1", "Trp53", "Hmga1", "Chek1", "Chek2", "Prodh", "Tnfrsf10b", "Cdkn1a", "Dao")

# Extraer expresión + metadatos
df <- FetchData(seurat_object, vars = c(senescence_core, RESOLUCION, "Phenotype")) %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = all_of(senescence_core), names_to = "Gene", values_to = "Expression") %>%
  mutate(Cluster = factor(RESOLUCION))

# Crear el plot
vln_plot <- ggplot(df, aes(x = Cluster, y = Expression, fill = Phenotype)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = TRUE) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(x = "Cluster", y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Guardar
ggsave(filename = file.path(path.guardar, "VlnPlot_Senescence_mix2.png"),
       plot = vln_plot, width = 18, height = 16)

gene_sets <- list(
  senescence_core <- c("Cdkn2a", "Bmi1", "Trp53", "Hmga1", "Chek1", "Chek2", "Prodh", "Tnfrsf10b", "Cdkn1a", "Dao")
)



# SIGNATURAS ......................................................................
# 1. Lista de firmas con nombres únicos
gene_sets <- list(
  senescence_markers = c("Cdkn1a","Ctnnb1","Cxcl16", "Hmgb1",
                     "Igfbp4", "Igfbp7", "Il1a", "Mmp14", "Kitl",  "Tnfrsf1a",  
                     "Ppp1ca", "Smurf2")
)

# 2. Añadir los module scores
# ¡Aquí está el truco!: usa `name = "GO_Score"`
seurat_object <- AddModuleScore(seurat_object, features = gene_sets, name = "GO_Score")

# Esto genera columnas en seurat_object@meta.data como:
# GO_Score1, GO_Score2, GO_Score3, ..., en el mismo orden que la lista

# 3. Mapear nombres manualmente
score_names <- paste0("GO_Score", seq_along(gene_sets))
names(score_names) <- names(gene_sets)  # así puedes acceder por nombre original

for (sig in names(score_names)) {
  score_col <- score_names[[sig]]
  
  p <- FeaturePlot(seurat_object, features = score_col) &
       NoAxes() &
       ggtitle(sig)  # <--- aquí se asigna el nombre original como título
  
  ggsave(filename = file.path(path.guardar, paste0("Score_", sig, ".png")), 
         plot = p, width = 8, height = 8)
}




# AUCell - SENSCENCE MARKERS

# ...............................................................................
# FILTRAR GENES PRESENTES
# ...............................................................................
counts <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")
valid_genes <- gene_sets$senescence_markers[gene_sets$senescence_markers %in% rownames(counts)]
if (length(valid_genes) < 5) stop("Muy pocos genes válidos encontrados.")

# ...............................................................................
# AUCell
# ...............................................................................
rankings <- AUCell_buildRankings(counts, plotStats = FALSE)
auc_result <- AUCell_calcAUC(list(senescence_markers = valid_genes), rankings)


# Calcular valores de AUC
auc_values <- getAUC(auc_result)[1, ]

# Definir thresholds manuales
thresholds <- list(
  Percentil_95 = quantile(auc_values, 0.95, na.rm = TRUE),
  Mean = mean(auc_values, na.rm = TRUE),
  Mean_1SD = mean(auc_values, na.rm = TRUE) + sd(auc_values, na.rm = TRUE),
  Mean_2SD = mean(auc_values, na.rm = TRUE) + 2 * sd(auc_values, na.rm = TRUE)
)

# Aplicar cada uno
for (thresh_name in names(thresholds)) {
  threshold_value <- thresholds[[thresh_name]]

  active_cells <- names(which(auc_values > threshold_value))

  col_name <- paste0("AUCell_senescence_", thresh_name)
  seurat_object[[col_name]] <- ifelse(colnames(seurat_object) %in% active_cells, "Active", "Inactive")

  # DimPlot
  p <- DimPlot(seurat_object, group.by = col_name, cols = c("firebrick", "grey80")) & NoAxes()
  ggsave(
    filename = file.path(path.guardar, paste0("DimPlot_", col_name, ".png")),
    plot = p, width = 8, height = 7
  )

  # Histograma con línea del threshold
  pdf(file.path(path.guardar, paste0("Histogram_", col_name, ".pdf")), width = 8, height = 6)
  hist(auc_values, breaks = 100, col = "lightgrey", main = paste0("AUC - ", col_name),
       xlab = "AUC score", xlim = c(min(auc_values), max(auc_values)))
  abline(v = threshold_value, col = "red", lwd = 2)
  dev.off()
}

# Aplicar cada uno split.by phenotype
for (thresh_name in names(thresholds)) {
  threshold_value <- thresholds[[thresh_name]]

  active_cells <- names(which(auc_values > threshold_value))

  col_name <- paste0("AUCell_senescence_", thresh_name)
  seurat_object[[col_name]] <- ifelse(colnames(seurat_object) %in% active_cells, "Active", "Inactive")

  # DimPlot
  p <- DimPlot(seurat_object, group.by = col_name, cols = c("firebrick", "grey80"), split.by="Phenotype") & NoAxes()
  ggsave(
    filename = file.path(path.guardar, paste0("DimPlot_", col_name, "Split.By.Phenotype.png")),
    plot = p, width = 14, height = 7
  )
}



# SUMMARY NÚMEROS ----------------------------------------------
# Inicializar workbook
wb <- createWorkbook()

# Tabla resumen total de activas por threshold
summary_table <- data.frame(
  Threshold = names(thresholds),
  Threshold_value = sapply(thresholds, round, 4),
  N_active_cells = sapply(thresholds, function(thr) {
    sum(auc_values > thr)
  })
)
addWorksheet(wb, sheetName = "Summary")
writeData(wb, "Summary", summary_table)

# Una hoja por cada threshold
for (th in names(thresholds)) {
  thr_val <- thresholds[[th]]
  active_cells <- names(which(auc_values > thr_val))
  status <- ifelse(colnames(seurat_object) %in% active_cells, "Active", "Inactive")

  seurat_object[[paste0("AUCell_", th)]] <- status

  tab_counts <- table(Status = status, Phenotype = seurat_object$Phenotype)
  tab_prop <- round(prop.table(tab_counts, margin = 2) * 100, 2)

  # Preparar hoja
  sheet <- paste0("Threshold_", th)
  addWorksheet(wb, sheet)

  # Escribir texto explicativo
  writeData(wb, sheet, paste0("Threshold used: ", round(thr_val, 4)), startRow = 1)
  writeData(wb, sheet, "Number of cells (count):", startRow = 3)
  writeData(wb, sheet, as.data.frame.matrix(tab_counts), startRow = 4)

  writeData(wb, sheet, "Proportion of Active/Inactive (%):", startRow = 8)
  writeData(wb, sheet, as.data.frame.matrix(tab_prop), startRow = 9)
}

#Guardar archivo
saveWorkbook(wb, file = file.path(path.guardar, "AUCell_Proportions_by_Threshold.xlsx"), overwrite = TRUE)


# Guardar resultados
saveRDS(auc_result, file = file.path(path.guardar, "AUC_senescence.rds"))
write.csv(getAUC(auc_result), file = file.path(path.guardar, "AUC_senescence.csv"))

# Histograma + umbral automático
pdf(file.path(path.guardar, "AUCell_Threshold_senescence.pdf"), width = 8, height = 6)
cells_assignment <- AUCell_exploreThresholds(auc_result, plotHist = TRUE, assign = TRUE)
dev.off()

# Extraer umbral
gene_set_name <- names(cells_assignment)[1]
aucThr <- cells_assignment[[gene_set_name]]$thresholds[[1]]$aucThr
if (is.null(aucThr)) stop("No se pudo calcular un threshold AUC automáticamente.")

# Células positivas
active_cells <- names(which(getAUC(auc_result)[1, ] > aucThr))

# Anotar en Seurat
col_name <- "AUCell_senescence"
seurat_object[[col_name]] <- ifelse(colnames(seurat_object) %in% active_cells, "Active", "Inactive")

# Plot
p <- DimPlot(seurat_object, group.by = col_name, cols = c("firebrick", "grey80")) & NoAxes()
ggsave(filename = file.path(path.guardar, "DimPlot_AUCell_senescence.png"), plot = p, width = 8, height = 7)

# Guardar objeto anotado
saveRDS(seurat_object, file = file.path(path.guardar, "Seurat_With_AUCell_senescence.rds"))





# AJUSTAR POR SI NO APARECE UN THRESHOLD AUTOMÁTICO
threshold_value <- 0.5  # o thresholds$Percentil_95

active_cells <- names(which(auc_values > threshold_value))
seurat_object$AUCell_senescence_manual <- ifelse(colnames(seurat_object) %in% active_cells, "Active", "Inactive")


P <- DimPlot(seurat_object, group.by = "AUCell_senescence_manual", split.by = "Phenotype", cols = c("firebrick", "grey80")) & NoAxes()  & 
     ggtitle("AUCell_senescence_automatic")
ggsave(file.path(path.guardar, "DimPlot_threshold_automatico_Phenotype.png"),
       plot = P, width = 14, height = 7)

P <- DimPlot(seurat_object, group.by = "AUCell_senescence_manual", cols = c("firebrick", "grey80")) & NoAxes()  & 
     ggtitle("AUCell_senescence_automatic")
ggsave(file.path(path.guardar, "DimPlot_threshold_automatico.png"),
       plot = P, width = 8, height = 7)