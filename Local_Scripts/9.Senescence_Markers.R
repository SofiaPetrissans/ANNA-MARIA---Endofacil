# ================================================================
# 15.AUCell_Senescencia_SASP_P53_AL5.R
# ================================================================
# OVERVIEW
# - Evalúa actividad de senescencia (panel manual), SASP, SENESCENCE (Reactome) y P53
#   en SubsetEndothelial_Harmony_Annotated.rds usando AUCell.
# - Anota AUC continuo + versión binaria (umbral "Mean" por defecto).
# - Plots: DimPlot binario, FeaturePlot continuo, violines, barras de proporciones.
#
# INPUTS
# - SubsetEndothelial_Harmony_Annotated.rds
#   meta.data debe contener: Phenotype, Sample, Annotation_Layer3, Annotation_Layer5
#
# OUTPUTS (se crean en /9.Senescence_AUCell)
# - AUC_por_celda.csv, thresholds.csv, resumenes_por_cluster.csv
# - PNG/PDF de visualizaciones (UMAP, violines, proporciones)
# ================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(viridis)
  library(msigdbr)
  library(openxlsx)
})

set.seed(1234)

# -----------------------------
# 0) PARÁMETROS / RUTAS
# -----------------------------
directory <- "/Users/graupera_lab/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/ANNA_MARIA"
setwd(directory)

input_rds <- file.path("SubsetEndothelial_Harmony_Annotated.rds")
output_dir <- file.path("9.Senescence_Markers")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(input_rds))
#seurat_object <- readRDS(input_rds)

# Chequeo de metadatos clave
need_cols <- c("Phenotype","Sample","Annotation_Layer3","Annotation_Layer5")
miss <- setdiff(need_cols, colnames(seurat_object@meta.data))
if (length(miss)) stop("Faltan columnas en meta.data: ", paste(miss, collapse = ", "))

# Orden Phenotype
seurat_object$Phenotype <- factor(seurat_object$Phenotype, levels = c("vehicle","CHOP"))

# Paleta AL5 (ajústala si quieres)
pal_AL5 <- c(
  "LVEC"   = "#86D0B9",
  "LSEC_1" = "#F4CDA5",
  "LSEC_2" = "#EBAFC3",
  "LSEC_3" = "#A8C8E5",
  "LSEC_4" = "#dce6b6",
  "LSEC_5" = "#C5A7E1"
)

# -----------------------------
# 1) UTILIDADES
# -----------------------------
get_counts_matrix <- function(obj, assay = "RNA") {
  # Intenta slot="counts" (Seurat clásico) y si falla prueba layer="counts" (Seurat v5 layered)
  m <- tryCatch(GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  if (is.null(m) || length(m) == 0) {
    m <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
  }
  if (is.null(m) || nrow(m) == 0) {
    # Último recurso: datos por defecto del assay
    m <- GetAssayData(obj, assay = assay)
  }
  if (nrow(m) == 0) stop("No se pudo obtener una matriz de conteos válida.")
  m
}

compute_thresholds <- function(v) {
  list(
    Percentile95 = as.numeric(quantile(v, 0.95, na.rm = TRUE)),
    Mean         = mean(v, na.rm = TRUE),
    Mean_1SD     = mean(v, na.rm = TRUE) + sd(v, na.rm = TRUE),
    Mean_2SD     = mean(v, na.rm = TRUE) + 2 * sd(v, na.rm = TRUE)
  )
}

add_auc_to_meta <- function(obj, auc_mat, set_name, thr_name = "Mean") {
  v <- as.numeric(auc_mat[set_name, colnames(obj)])
  obj[[paste0("AUC_", set_name)]] <- v
  thrs <- compute_thresholds(v)
  thr  <- thrs[[thr_name]]
  obj[[paste0("AUC_", set_name, "_bin")]] <- ifelse(v > thr, "Active", "Inactive")
  list(object = obj, thresholds = data.frame(Set = set_name, Threshold = thr_name, Value = thr))
}

# Plots helpers
save_dimplot_binary <- function(obj, col_name, file, active_col = "firebrick", inactive_col = "grey80") {
  p <- DimPlot(obj, group.by = col_name, cols = c(Active = active_col, Inactive = inactive_col),
               pt.size = 0.4) + NoAxes() + ggtitle(col_name)
  ggsave(file, p, width = 7, height = 6, dpi = 300)
}

save_featureplot <- function(obj, feat, file) {
  p <- FeaturePlot(obj, features = feat, reduction = "umap", order = TRUE,
                   cols = c("grey90","navy"), min.cutoff = "q5", max.cutoff = "q95",
                   pt.size = 0.3)
  ggsave(file, p, width = 7.5, height = 6, dpi = 300)
}

# -----------------------------
# 2) GENE SETS
# -----------------------------
# Panel manual (tu lista)
gs_manual <- list(
  SENESC_PANEL = c("Cdkn1a","Ctnnb1","Cxcl16","Hmgb1","Igfbp4","Igfbp7",
                   "Il1a","Mmp14","Kitl","Tnfrsf1a","Ppp1ca","Smurf2")
)

# MSigDB (Mus musculus)
mm <- msigdbr(species = "Mus musculus")

pick_one <- function(mm, name_vec) {
  nm <- name_vec[name_vec %in% unique(mm$gs_name)][1]
  if (is.na(nm)) return(list(name = NA_character_, genes = character()))
  genes <- mm %>% filter(gs_name == nm) %>% pull(gene_symbol) %>% unique()
  list(name = nm, genes = genes)
}

SASP <- pick_one(mm, c("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                       "WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                       "FRIDMAN_SENESCENCE_UP"))
SENESC <- pick_one(mm, "REACTOME_CELLULAR_SENESCENCE")
P53   <- pick_one(mm, c("HALLMARK_P53_PATHWAY","TANG_SENESCENCE_TP53_TARGETS_UP"))

# -----------------------------
# 3) PREPARAR MATRIZ & FILTRAR GENES PRESENTES
# -----------------------------
counts <- get_counts_matrix(seurat_object, assay = "RNA")
present <- rownames(counts)

# Intersección con matriz
gs_list <- list(
  SENESC_PANEL = intersect(gs_manual$SENESC_PANEL, present),
  SASP         = intersect(SASP$genes, present),
  SENESCENCE   = intersect(SENESC$genes, present),
  P53          = intersect(P53$genes, present)
)

# Fallbacks suaves para SASP si quedó muy corta
if (length(gs_list$SASP) < 8 && SASP$name != "WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP") {
  SASP <- pick_one(mm, "WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP")
  gs_list$SASP <- intersect(SASP$genes, present)
}
if (length(gs_list$SASP) < 8 && SASP$name != "FRIDMAN_SENESCENCE_UP") {
  SASP <- pick_one(mm, "FRIDMAN_SENESCENCE_UP")
  gs_list$SASP <- intersect(SASP$genes, present)
}

sizes <- vapply(gs_list, length, integer(1))
message("Tamaños gene sets tras intersección -> ",
        paste(names(sizes), sizes, sep = ":", collapse = " | "))
if (any(sizes < 5)) warning("Algún gene set quedó muy corto (<5). Revisa el log y/o paneles.")

# -----------------------------
# 4) AUCell
# -----------------------------
rankings <- AUCell_buildRankings(counts, plotStats = FALSE)  # nCores por defecto
auc_res  <- AUCell_calcAUC(gs_list, rankings)
auc_mat  <- getAUC(auc_res)  # filas = gene sets, cols = células

# Guardar matriz AUC
write.csv(as.data.frame(t(auc_mat)),
          file.path(output_dir, "AUC_por_celda.csv"),
          row.names = TRUE)

# Añadir AUCs + binarios (umbral = "Mean") al objeto
thr_table <- list()
for (set_nm in rownames(auc_mat)) {
  res <- add_auc_to_meta(seurat_object, auc_mat, set_nm, thr_name = "Mean")
  seurat_object <- res$object
  thr_table[[set_nm]] <- res$thresholds
}
thr_table <- bind_rows(thr_table)
write.csv(thr_table, file.path(output_dir, "thresholds.csv"), row.names = FALSE)

# -----------------------------
# 5) PLOTS BÁSICOS (UMAP)
# -----------------------------
# Continuos
for (set_nm in rownames(auc_mat)) {
  save_featureplot(seurat_object, paste0("AUC_", set_nm),
                   file.path(output_dir, paste0("UMAP_AUC_", set_nm, ".png")))
}

# Binarios (Active/Inactive)
for (set_nm in rownames(auc_mat)) {
  save_dimplot_binary(seurat_object,
                      paste0("AUC_", set_nm, "_bin"),
                      file.path(output_dir, paste0("UMAP_BIN_", set_nm, ".png")))
}

# -----------------------------
# 6) VIOLINES por AL3 (split por condición)
# -----------------------------
vln_feats <- paste0("AUC_", rownames(auc_mat))
p_vln <- VlnPlot(seurat_object, features = vln_feats,
                 group.by = "Annotation_Layer3", split.by = "Phenotype", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "VlnPlot_AUC_by_AL3_splitPhenotype.png"),
       p_vln, width = 12, height = 6, dpi = 300)

# -----------------------------
# 7) BARRAS DE PROPORCIONES (AL5)
# -----------------------------
for (set_nm in rownames(auc_mat)) {
  bin_col <- paste0("AUC_", set_nm, "_bin")
  # asegurar niveles
  seurat_object@meta.data[[bin_col]] <- factor(
    seurat_object@meta.data[[bin_col]],
    levels = c("Inactive","Active")
  )
  
  df <- seurat_object@meta.data %>%
    filter(!is.na(.data[[bin_col]])) %>%
    count(Annotation_Layer5, .data[[bin_col]], name = "cells") %>%
    group_by(Annotation_Layer5) %>%
    mutate(prop = cells / sum(cells)) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = Annotation_Layer5, y = prop, fill = .data[[bin_col]])) +
    geom_col(width = 0.8) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(Inactive = "grey80", Active = "firebrick"),
                      name = paste0(set_nm, " (bin)")) +
    labs(x = NULL, y = "Proporción de células",
         title = paste0("Actividad ", set_nm, " por Annotation_Layer5")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, paste0("Barplot_", set_nm, "_by_AL5.png")),
         p, width = 7, height = 4, dpi = 300)
}

# -----------------------------
# 8) BARRAS POR AL3 × Phenotype (facet)
# -----------------------------
for (set_nm in rownames(auc_mat)) {
  bin_col <- paste0("AUC_", set_nm, "_bin")
  df2 <- seurat_object@meta.data %>%
    filter(!is.na(.data[[bin_col]])) %>%
    count(Phenotype, Annotation_Layer3, .data[[bin_col]], name = "cells") %>%
    group_by(Phenotype, Annotation_Layer3) %>%
    mutate(prop = cells / sum(cells)) %>%
    ungroup()
  
  p2 <- ggplot(df2, aes(x = Phenotype, y = prop, fill = .data[[bin_col]])) +
    geom_col(width = 0.8) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(Inactive = "grey80", Active = "firebrick"),
                      name = paste0(set_nm, " (bin)")) +
    labs(x = NULL, y = "Proporción de células",
         title = paste0("Actividad ", set_nm, " por AL3 y condición")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    facet_wrap(~ Annotation_Layer3, nrow = 1)
  
  ggsave(file.path(output_dir, paste0("Barplot_", set_nm, "_by_AL3_byPhenotype.png")),
         p2, width = 14, height = 4, dpi = 300)
}

# -----------------------------
# 9) RESÚMENES Y GUARDADO
# -----------------------------
# Resumen por AL5 de proporciones Active/Inactive para cada set
resumen_list <- lapply(rownames(auc_mat), function(set_nm){
  bin_col <- paste0("AUC_", set_nm, "_bin")
  seurat_object@meta.data %>%
    filter(!is.na(.data[[bin_col]])) %>%
    count(Annotation_Layer5, .data[[bin_col]], name = "cells") %>%
    group_by(Annotation_Layer5) %>%
    mutate(prop = cells / sum(cells)) %>%
    mutate(GeneSet = set_nm) %>%
    ungroup()
})
resumen_df <- bind_rows(resumen_list)
write.csv(resumen_df, file.path(output_dir, "resumenes_por_cluster.csv"), row.names = FALSE)

# Guardar objeto con anotaciones AUC
saveRDS(seurat_object, file.path(output_dir, "SubsetEndothelial_Harmony_Annotated_AUCell.rds"))

message("Listo. Resultados en: ", output_dir)
