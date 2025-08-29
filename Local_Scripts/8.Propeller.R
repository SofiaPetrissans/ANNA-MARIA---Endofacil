# ================================================================
# 14.Propeller_ComparacionProporciones_AL1-AL5.R
# ================================================================
# OVERVIEW
# Compara las proporciones celulares entre CHOP vs vehicle usando speckle::propeller,
# para varias capas de anotación (AL1–AL5). Genera:
#  - Heatmap de proporciones medias por fenotipo con anotación de -log10(p)
#  - Barras de PropMean por fenotipo
#  - Barras de log2(CHOP/vehicle) con estrellas de significación
#  - Curva puntos+línea de conteos por fenotipo
#
# INPUTS
#  - Seurat RDS anotado local: SubsetEndothelial_Harmony_Annotated.rds
#  - Metadatos requeridos en seurat_object@meta.data: Sample, Phenotype
#
# OUTPUTS
#  - 8.Propeller/ALx/ ... PNG/PDF + RDS + tablas (xlsx/csv)
#
# ================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
})

# ---------- Comprobar speckle ----------
if (!requireNamespace("speckle", quietly = TRUE)) {
  stop("El paquete 'speckle' no está instalado. Instálalo con: remotes::install_github('Oshlack/speckle')")
}
library(speckle)

# -----------------------------
# 0) PARÁMETROS / RUTAS
# -----------------------------
# Directorio
directory <- "/Users/graupera_lab/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/ANNA_MARIA"
setwd(directory)

input_rds  <- file.path("SubsetEndothelial_Harmony_Annotated.rds")
out_root   <- file.path("0.Propeller")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Capas a evaluar (incluye AL4 y AL5 como querías)
layers <- c("Annotation_Layer1", "Annotation_Layer2", "Annotation_Layer3",
            "Annotation_Layer4", "Annotation_Layer5")

# Paletas por capa (ajusta si lo deseas)
pal_AL1 <- c("Central Vein"="#86D0B9", "LSEC"="#CCEDB1")
pal_AL2 <- c("Central Vein"="#86D0B9", "Periportal"="#5A8BB7", "Midzonal"="#99C5E3")
pal_AL3 <- c("Central Vein"="#86D0B9","LSEC_1"="#F4CDA5","LSEC_2"="#EBAFC3",
             "LSEC_3"="#A8C8E5","LSEC_4"="#dce6b6","LSEC_5"="#C5A7E1")
pal_AL4 <- c("LVEC"="#86D0B9","LSEC"="#CCEDB1")
pal_AL5 <- c("LVEC"="#86D0B9","LSEC_1"="#F4CDA5","LSEC_2"="#EBAFC3",
             "LSEC_3"="#A8C8E5","LSEC_4"="#dce6b6","LSEC_5"="#C5A7E1")

paletas <- list(
  Annotation_Layer1 = pal_AL1,
  Annotation_Layer2 = pal_AL2,
  Annotation_Layer3 = pal_AL3,
  Annotation_Layer4 = pal_AL4,
  Annotation_Layer5 = pal_AL5
)

# Colores condición
pal_cond <- c(vehicle = "#E4F6DD", CHOP = "#FFD7D9")

# -----------------------------
# 1) UTILIDADES
# -----------------------------

# Asegura metadatos requeridos
check_metadata <- function(obj) {
  req <- c("Sample","Phenotype")
  miss <- setdiff(req, colnames(obj@meta.data))
  if (length(miss)) stop("Faltan columnas en meta.data: ", paste(miss, collapse=", "))
  # factorizamos (evita problemas)
  obj$Sample    <- factor(obj$Sample)
  obj$Phenotype <- factor(obj$Phenotype, levels = c("vehicle","CHOP"))
  obj
}

# Detecta columna de clúster en salida de propeller (robusto)
detect_cluster_col <- function(df) {
  cands <- c("BaselineProp.clusters", "clusters", "Level", "Cluster")
  hit <- intersect(cands, names(df))
  if (length(hit)) return(hit[1])
  # fallback: primera columna
  names(df)[1]
}

# Ejecuta propeller con seurat ident actual
run_propeller <- function(obj, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cells_info <- data.frame(
    clusters = obj@active.ident,
    samples  = obj$Sample,
    group    = obj$Phenotype
  )
  saveRDS(cells_info, file.path(out_dir, "cells_info.rds"))
  
  res <- speckle::propeller(
    clusters = cells_info$clusters,
    sample   = cells_info$samples,
    group    = cells_info$group
  )
  saveRDS(res, file.path(out_dir, "Res_Propeller.rds"))
  openxlsx::write.xlsx(res, file.path(out_dir, "Res_Propeller.xlsx"), overwrite = TRUE)
  res
}

# Plot puntos+línea de conteos por fenotipo
plot_counts_curve <- function(obj, layer, pal_ct, out_dir) {
  df <- as.data.frame(table(obj$Phenotype, obj@active.ident))
  colnames(df) <- c("Phenotype","Cluster","Freq")
  p <- ggplot(df, aes(x = Phenotype, y = Freq, col = Cluster, group = Cluster)) +
    geom_point(size = 4) + geom_line(linewidth = 0.6) +
    scale_color_manual(values = pal_ct) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Conteos por condición -", layer), x = NULL, y = "Nº células")
  ggsave(file.path(out_dir, "Counts_DotLine_ByCondition.png"), p, width = 7, height = 5, dpi = 300)
}

# Heatmap de PropMean por fenotipo con -log10(p) a la derecha
plot_propmean_heatmap <- function(res, pal_ct, out_dir) {
  ct_col <- detect_cluster_col(res)
  
  mat <- res |>
    dplyr::select(all_of(ct_col), dplyr::starts_with("PropMean.")) |>
    as.data.frame()
  rownames(mat) <- mat[[ct_col]]
  mat[[ct_col]] <- NULL
  mat <- as.matrix(mat)
  
  # columnas limpias: "PropMean.vehicle" -> "vehicle"
  colnames(mat) <- stringr::str_replace(colnames(mat), "^PropMean\\.", "")
  
  # ordenar filas con paleta
  ord <- intersect(names(pal_ct), rownames(mat))
  mat <- mat[ord, , drop = FALSE]
  
  # -log10(p)
  pvals <- setNames(res$P.Value, res[[ct_col]])[rownames(mat)]
  
  col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("white", "#5A8BB7", "#2A5783"))
  
  ha_top <- HeatmapAnnotation(
    Phenotype = factor(colnames(mat), levels = c("vehicle","CHOP")),
    col = list(Phenotype = pal_cond),
    show_annotation_name = FALSE
  )
  ha_right_ct <- rowAnnotation(
    CellType = factor(rownames(mat), levels = names(pal_ct)),
    col = list(CellType = pal_ct),
    show_annotation_name = TRUE
  )
  ha_right_p <- rowAnnotation(
    `-log10(p)` = anno_points(-log10(pvals), pch = 16, gp = gpar(col = "#1E88E5"),
                              size = unit(2, "mm"),
                              ylim = c(0, max(-log10(pvals), na.rm = TRUE)),
                              axis_param = list(labels_rot = 90))
  )
  
  ht <- Heatmap(mat, name = "Prop.Mean", col = col_fun,
                top_annotation = ha_top,
                right_annotation = c(ha_right_ct, ha_right_p),
                show_row_names = TRUE, show_column_names = TRUE,
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_gp = gpar(fontsize = 9))
  
  pdf(file.path(out_dir, "Heatmap_PropMean_ByPhenotype.pdf"), width = 4.2, height = 5)
  draw(ht, show_annotation_legend = FALSE)
  dev.off()
}

# Barras PropMean por fenotipo (lado a lado)
plot_propmean_bars <- function(res, pal_ct, out_dir) {
  ct_col <- detect_cluster_col(res)
  df <- res %>%
    dplyr::select(all_of(ct_col), dplyr::starts_with("PropMean.")) %>%
    dplyr::rename_with(~ sub("^PropMean\\.", "", .x), dplyr::starts_with("PropMean.")) %>%
    mutate(across(-all_of(ct_col), as.numeric)) %>%
    pivot_longer(cols = c("vehicle","CHOP"),
                 names_to = "Phenotype", values_to = "Proportion") %>%
    dplyr::rename(CellType = !!ct_col)
  
  df$CellType  <- factor(df$CellType, levels = names(pal_ct))
  df$Phenotype <- factor(df$Phenotype, levels = c("vehicle","CHOP"))
  
  p <- ggplot(df, aes(CellType, Proportion, fill = Phenotype)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = pal_cond) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Proportion")
  ggsave(file.path(out_dir, "BarPlot_PropMean_byPhenotype.png"), p, width = 7, height = 4, dpi = 300)
}

# Barras log2(CHOP/vehicle) + estrellas de significación
plot_log2fc_bars <- function(res, pal_ct, out_dir) {
  ct_col <- detect_cluster_col(res)
  eps <- 1e-6
  
  # columnas de PropMean por regex (robusto a nombres)
  veh_col  <- grep("^PropMean\\.veh", names(res), value = TRUE)
  chop_col <- grep("^PropMean\\.CHOP", names(res), value = TRUE)
  if (!length(veh_col) || !length(chop_col)) {
    # fallback nombres exactos
    veh_col  <- "PropMean.vehicle"
    chop_col <- "PropMean.CHOP"
  }
  
  df <- res %>%
    transmute(
      CellType = .data[[ct_col]],
      vehicle  = .data[[veh_col]],
      CHOP     = .data[[chop_col]],
      pval     = P.Value,
      FDR      = FDR
    ) %>%
    mutate(log2FC = log2((CHOP + eps) / (vehicle + eps)))
  
  df$CellType <- factor(df$CellType, levels = names(pal_ct))
  df$stars <- cut(df$FDR,
                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                  labels = c("***","**","*",""))
  
  p <- ggplot(df, aes(x = CellType, y = log2FC, fill = CellType)) +
    geom_col(width = 0.75) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_text(aes(label = stars,
                  y = ifelse(log2FC >= 0, log2FC + 0.05, log2FC - 0.05)),
              vjust = ifelse(df$log2FC >= 0, 0, 1), size = 3) +
    scale_fill_manual(values = pal_ct) +
    labs(y = "log2(CHOP / vehicle)", x = NULL,
         title = "Cambios de proporción por subpoblación") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  ggsave(file.path(out_dir, "BarLog2FC_propeller.png"), p, width = 4.5, height = 3.2, dpi = 300)
}

# -----------------------------
# 2) MAIN
# -----------------------------
stopifnot(file.exists(input_rds))
seurat_object <- readRDS(input_rds)
seurat_object <- check_metadata(seurat_object)

for (layer in layers) {
  message("\n==> Procesando ", layer)
  if (!layer %in% names(paletas)) {
    warning("No hay paleta definida para ", layer, ". Se omite.")
    next
  }
  pal_ct <- paletas[[layer]]
  
  # Subcarpeta de salida por capa
  out_dir <- file.path(out_root, layer)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Fijar ident a la capa y plots básicos de conteos
  if (!layer %in% colnames(seurat_object@meta.data)) {
    warning("La columna ", layer, " no existe en meta.data. Se omite.")
    next
  }
  Idents(seurat_object) <- seurat_object[[layer, drop = TRUE]]
  
  # Curva de conteos por condición
  plot_counts_curve(seurat_object, layer, pal_ct, out_dir)
  
  # Ejecutar propeller
  res <- run_propeller(seurat_object, out_dir)
  
  # Guardar también CSV rápido
  write.csv(res, file.path(out_dir, "Res_Propeller.csv"), row.names = FALSE)
  
  # Heatmap de PropMean con -log10(p)
  plot_propmean_heatmap(res, pal_ct, out_dir)
  
  # Barras PropMean por fenotipo
  plot_propmean_bars(res, pal_ct, out_dir)
  
  # Barras log2(CHOP/vehicle) + estrellas
  plot_log2fc_bars(res, pal_ct, out_dir)
}

message("\nListo. Resultados en: ", out_root)

