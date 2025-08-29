# --------------------------------------------------------------------------------
# SCRIPT: Generar DEGs entre tratamientos - CHOP y vehicle 
# AUTORA: Sofia Petrissans Moll
# FECHA: 30/06/2025
# --------------------------------------------------------------------------------

# LIBRARY
source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.6_Downstream_Analysis_EC/0.0_Paths.R")
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(openxlsx)

MT_THRESHOLD <- "mt20"
# Paths
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "Seurat_annotated.rds")


path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, MT_THRESHOLD, "1.DEGs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# CARGA Y PREPARACIÓN
# --------------------------------------------------------------------------------
seurat_object <- readRDS(input_rds)

# --------------------------------------------------------------------------------
# SCRIPT: Generar DEGs entre tratamientos - CHOP y vehicle para múltiples capas
# AUTORA: Sofia Petrissans Moll
# FECHA: 01/07/2025
# --------------------------------------------------------------------------------

# LIBRARIES
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(openxlsx)

# ------------------------------------------------------------------------------
# 1. PATHS Y OBJETO SEURAT
# ------------------------------------------------------------------------------
input_rds <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/0.3_Annotation_EC_mt20/5.Annotation/Seurat_annotated.rds"
base_output_dir <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/0.3_Annotation_EC_mt20/6.DEGs"
seurat_object <- readRDS(input_rds)

# ------------------------------------------------------------------------------
# 2. DEFINIR CAPAS Y COLORES
# ------------------------------------------------------------------------------

annotation_settings <- list(

  Annotation_LSEC_LVEC = c(
    "LSEC" = "#B2DE81",
    "LVEC" = "#D8B2C6"
  ),

  Annotation_Layer2 = c(
    "Pericentral" = "#EDF8BC",
    "Midzonal" = "#B6E4B3",
    "Periportal" = "#63C3BF",
    "Portal area" = "#3A95B1"
  ),

  Annotation_Layer3 = c(
    "LSEC_1" = "#66A61E",
    "LSEC_2" = "#B2DE81",
    "LSEC_3" = "#1B9E77",
    "LVEC"    = "#D8B2C6"
  ),

  Annotation_Layer4 = c(
    "LSEC_1" = "#A6CEE3",
    "LSEC_2" = "#B3DE69",
    "LSEC_3" = "#FDBF6F",
    "LSEC_4" = "#CAB2D6",
    "LSEC_5" = "#FFCCBC",
    "LSEC_6" = "#B2DF8A",
    "LVEC"   = "#D8B2C6"
  )
)

# ------------------------------------------------------------------------------
# 3. BUCLE PRINCIPAL PARA CADA CAPA
# ------------------------------------------------------------------------------

for (annotation_col in names(annotation_settings)) {
  message("Procesando capa: ", annotation_col)

  output_dir <- file.path(base_output_dir, annotation_col)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  Annotation_Colors <- annotation_settings[[annotation_col]]
  seurat_object <- SetIdent(seurat_object, value = annotation_col)

  cell_types <- unique(seurat_object[[annotation_col]][, 1])
  deg_summary <- tibble()
  deg_list <- list()

  # Loop por tipo celular
  for (cell_type in cell_types) {
    message("  Procesando tipo celular: ", cell_type)

    cells_to_keep <- colnames(seurat_object)[seurat_object[[annotation_col]][,1] == cell_type]
    subset_obj <- subset(seurat_object, cells = cells_to_keep)

    cond_counts <- table(subset_obj$Phenotype)

    if (all(c("CHOP", "vehicle") %in% names(cond_counts)) && all(cond_counts >= 10)) {
      Idents(subset_obj) <- "Phenotype"
      markers <- FindMarkers(subset_obj, ident.1 = "CHOP", ident.2 = "vehicle",
                             logfc.threshold = 0.25, min.pct = 0.1)
      markers <- markers %>% rownames_to_column("gene")
      deg_list[[cell_type]] <- markers
      n_degs <- sum(markers$p_val_adj < 0.05)
    } else {
      deg_list[[cell_type]] <- data.frame(gene = NA, warning = "Insuficientes células para comparar")
      n_degs <- NA
    }

    deg_summary <- bind_rows(deg_summary, tibble(CellType = cell_type, DEGs = n_degs))
  }

  # Añadir colores al resumen
  deg_summary$Color <- Annotation_Colors[as.character(deg_summary$CellType)]

  # Guardar resumen
  write.xlsx(deg_summary, file = file.path(output_dir, "DEGs_per_CellType_summary.xlsx"), rowNames = FALSE)

  # Guardar DEGs completos en hojas separadas
  wb <- createWorkbook()
  for (cell_type in names(deg_list)) {
    sheet_name <- gsub(" ", "_", cell_type)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = deg_list[[cell_type]])
  }
  saveWorkbook(wb, file = file.path(output_dir, "DEGs_per_CellType_all_genes.xlsx"), overwrite = TRUE)

  # Barplot total
  p <- ggplot(deg_summary, aes(x = reorder(CellType, DEGs), y = DEGs, fill = CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Annotation_Colors) +
    coord_flip() +
    labs(x = "Cell Type", y = "Number of DEGs (CHOP vs vehicle)",
         title = paste("DEGs per Cell Type -", annotation_col)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")

  ggsave(filename = file.path(output_dir, "Barplot_DEGs_per_CellType.png"), plot = p, width = 6.5, height = 5)

  # Up/Downregulated DEGs
  deg_updown <- tibble()
  for (cell_type in names(deg_list)) {
    df <- deg_list[[cell_type]]
    if (!"gene" %in% colnames(df)) next

    df_sig <- df %>% filter(p_val_adj < 0.05)
    n_up <- sum(df_sig$avg_log2FC > 0)
    n_down <- sum(df_sig$avg_log2FC < 0)

    deg_updown <- bind_rows(deg_updown, tibble(
      CellType = cell_type,
      Direction = c("Up", "Down"),
      DEGs = c(n_up, n_down)
    ))
  }

  deg_totals <- deg_updown %>%
    group_by(CellType) %>%
    summarise(Total = sum(DEGs)) %>%
    arrange(desc(Total))

  deg_updown$CellType <- factor(deg_updown$CellType, levels = rev(deg_totals$CellType))

  # Stacked barplot
  p2 <- ggplot(deg_updown, aes(x = CellType, y = DEGs, fill = Direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Up" = "#B3E0A6", "Down" = "#F78D80")) +
    coord_flip() +
    labs(x = "Cell Type", y = "Number of DEGs", title = "Up/Downregulated DEGs per Cell Type") +
    theme_classic(base_size = 13)

  ggsave(filename = file.path(output_dir, "Barplot_DEGs_UpDown_Stacked.png"), plot = p2, width = 7, height = 5)

  # Split facets por dirección
  p3 <- ggplot(deg_updown, aes(x = CellType, y = DEGs, fill = CellType)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Direction, scales = "free_y") +
    scale_fill_manual(values = Annotation_Colors) +
    coord_flip() +
    labs(x = "Cell Type", y = "Number of DEGs", title = "Up vs Downregulated DEGs") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")

  ggsave(filename = file.path(output_dir, "Barplot_DEGs_UpDown_SplitFacets.png"), plot = p3, width = 8, height = 5)
}
