# ============================================================
# 12.DEGs_CHOP_vs_vehicle.R
# ------------------------------------------------------------
# OVERVIEW
# - Calcula DEGs CHOP vs vehicle por subpoblación para múltiples
#   capas de anotación (Layer1..Layer5).
# - Resume nº de DEGs up/down y los visualiza (conteo).
# - Enriquecimiento GO:BP en genes up-regulados (por subpoblación).
# - UpSet plots de intersección de up-DEGs entre subpoblaciones.
# - Correlación de log2FC subpoblacionales frente a DEGs globales.
# - Descubre GO:BP "exclusivos" por subpoblación y los grafica.
#
# INPUTS (local):
# - Seurat RDS anotado a varios layers:
#     SubsetEndothelial_Harmony_Annotated.rds
# -  Excel con DEGs endoteliales: `DEGs_CHOP_vs_vehicle.xlsx` (una hoja; debe incluir columnas gene, avg_log2FC, p_val_adj).
#     ("/ijc/LABS/GRAUPERA/RAW/.../ANNA_MARIA/3.MERGE/0.6_Downstream_Analysis_EC/mt15/1.DEGs_phenotype/DEGs_CHOP_vs_vehicle.xlsx")
# 
#
# OUTPUTS (en ./12.DEGs_CHOP_vs_vehicle):
# - DEGs_<Layer>_CHOP_vs_vehicle.xlsx (una hoja por subpoblación)
# - Barplot_DEGs_<Layer>_*.png (conteo)
# - GO_Upregulated/<Layer>/*.xlsx + barplots
# - UpSet_DEGs_<Layer>_CHOP_vs_vehicle.png
# - Correlation_Global_vs_Subpop/Correlation_<Layer>.csv + heatmap
# - Exclusive_GO_*.xlsx + barplots exclusivos por subpoblación
# ============================================================

# ================================
# 1) LIBRERÍAS
# ================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(tidyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ComplexHeatmap)
  library(circlize)          # <- colorRamp2 está aquí
  library(ComplexUpset)
  library(stringr)
  library(purrr)
})

# ================================
# 2) PARÁMETROS / RUTAS
# ================================
# Directorio base (ajusta si lo necesitas)
directory <- "/Users/graupera_lab/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/ANNA_MARIA"
setwd(directory)

# Rutas de entrada
input_rds <- "SubsetEndothelial_Harmony_Annotated.rds"
global_deg_file <- file.path("DEGs_CHOP_vs_vehicle.xlsx")  # una hoja con gene / avg_log2FC / p_val_adj

# Carpeta de salida
output_dir <- "6.subpop_DEGs_CHOPvsVehicle_layers1-5_GO_UpSet_corr"
dir.create(output_dir, showWarnings = FALSE)

# Parámetros FindMarkers
LOGFC_THRESHOLD <- 0.25
MIN_PCT         <- 0.10
TEST_USE        <- "wilcox"
MIN_CELLS_PER_GROUP <- 20  # evita comparaciones con muy pocas células por condición

# Capas de anotación a procesar (incluye 4 y 5)
annotations <- c("Annotation_Layer1", "Annotation_Layer2", "Annotation_Layer3",
                 "Annotation_Layer4", "Annotation_Layer5")

# Paletas por Layer (coherentes con tus scripts previos)
annotation.colors_1 <- c("Central Vein"="#86D0B9","LSEC"="#CCEDB1")
annotation.colors_2 <- c("Central Vein"="#86D0B9","Periportal"="#5A8BB7","Midzonal"="#99C5E3")
annotation.colors_3 <- c("Central Vein"="#86D0B9","LSEC_1"="#F4CDA5","LSEC_2"="#EBAFC3",
                         "LSEC_3"="#A8C8E5","LSEC_4"="#dce6b6","LSEC_5"="#C5A7E1")
annotation.colors_4 <- c("LVEC"="#86D0B9","LSEC"="#CCEDB1")
annotation.colors_5 <- c("LVEC"="#86D0B9","LSEC_1"="#F4CDA5","LSEC_2"="#EBAFC3","LSEC_3"="#A8C8E5",
                         "LSEC_4"="#dce6b6","LSEC_5"="#C5A7E1")

annotation_colors_list <- list(
  Annotation_Layer1 = annotation.colors_1,
  Annotation_Layer2 = annotation.colors_2,
  Annotation_Layer3 = annotation.colors_3,
  Annotation_Layer4 = annotation.colors_4,
  Annotation_Layer5 = annotation.colors_5
)

# ================================
# 3) LECTURA DEL SEURAT OBJECT 
# ================================
seurat_object <- readRDS(input_rds)


# ================================
# 4) FUNCIONES AUXILIARES
# ================================

standardize_deg_cols <- function(df) {
  # fuerza columnas gene, avg_log2FC, p_val_adj
  if (!("gene" %in% colnames(df))) {
    if (!is.null(rownames(df)) && all(rownames(df) != "")) {
      df$gene <- rownames(df)
    } else {
      colnames(df)[1] <- "gene"
    }
  }
  if (!("avg_log2FC" %in% colnames(df))) {
    cand <- grep("log2", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(cand) > 0) colnames(df)[match(cand[1], colnames(df))] <- "avg_log2FC"
  }
  if (!("p_val_adj" %in% colnames(df))) {
    cand <- grep("adj", colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(cand) > 0) colnames(df)[match(cand[1], colnames(df))] <- "p_val_adj"
  }
  return(df)
}

plot_go_terms_inside_bars <- function(go_xlsx, output_dir, annotation_colors) {
  dir.create(output_dir, showWarnings = FALSE)
  sheets <- getSheetNames(go_xlsx)
  for (sheet in sheets) {
    df <- read.xlsx(go_xlsx, sheet = sheet)
    if (nrow(df) == 0) next
    top10 <- df %>%
      arrange(p.adjust) %>%
      slice_head(n = 10) %>%
      mutate(
        Term = Description,
        NegativeLog10P = -log10(p.adjust),
        ID = factor(Term, levels = rev(Term)),
        Group = sheet
      )
    p <- ggplot(top10, aes(x = ID, y = NegativeLog10P, fill = Group)) +
      geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
      geom_text(aes(label = Term), y = 0, hjust = 0, color = "black", size = 4, fontface = "bold") +
      scale_fill_manual(values = annotation_colors, guide = "none") +
      coord_flip() +
      labs(title = paste("Top 10 GO:BP -", gsub("__", ": ", sheet)),
           x = NULL, y = expression(-log[10]~"(adj. p-value)")) +
      theme_classic(base_size = 13) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
    ggsave(file.path(output_dir, paste0("GOBarplot_", sheet, "_terms_in_bars.png")), p, width = 8, height = 6)
  }
}

safe_FindMarkers <- function(obj, ident.1, ident.2, logfc, minpct, test) {
  # evita errores si hay muy pocas células en algún grupo
  n1 <- sum(Idents(obj) == ident.1)
  n2 <- sum(Idents(obj) == ident.2)
  if (n1 < MIN_CELLS_PER_GROUP || n2 < MIN_CELLS_PER_GROUP) {
    message(sprintf("  - SKIP: '%s' (%d) vs '%s' (%d) < %d cells",
                    ident.1, n1, ident.2, n2, MIN_CELLS_PER_GROUP))
    return(NULL)
  }
  FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2,
              logfc.threshold = logfc, min.pct = minpct, test.use = test)
}

# ================================
# 5) DEGs CHOP vs vehicle POR LAYER
# ================================

for (layer in annotations) {
  message("Procesando: ", layer)
  Idents(seurat_object) <- seurat_object[[layer, drop = TRUE]]
  subpopulations <- levels(seurat_object)
  wb <- createWorkbook()
  
  for (pop in subpopulations) {
    message("  Subpoblación: ", pop)
    sub_obj <- subset(seurat_object, idents = pop)
    Idents(sub_obj) <- sub_obj$Phenotype
    
    degs <- safe_FindMarkers(sub_obj, "CHOP", "vehicle",
                             LOGFC_THRESHOLD, MIN_PCT, TEST_USE)
    if (is.null(degs)) {
      addWorksheet(wb, pop)
      writeData(wb, pop, data.frame())  # hoja vacía -> indica skip
      next
    }
    degs$gene <- rownames(degs)
    addWorksheet(wb, pop)
    writeData(wb, pop, degs)
  }
  
  saveWorkbook(wb, file.path(output_dir, paste0("DEGs_", layer, "_CHOP_vs_vehicle.xlsx")), overwrite = TRUE)
}

# ================================
# 6) BARPLOTS: Nº DEGs (UP/DOWN) y PROPORCIONES
# ================================

deg_summary <- list()
for (layer in annotations) {
  file <- file.path(output_dir, paste0("DEGs_", layer, "_CHOP_vs_vehicle.xlsx"))
  sheets <- getSheetNames(file)
  for (sheet in sheets) {
    df <- read.xlsx(file, sheet = sheet)
    if (nrow(df) == 0) {
      up_count <- 0; down_count <- 0
    } else {
      df <- standardize_deg_cols(df)
      up_count   <- sum(df$avg_log2FC > 0 & df$p_val_adj < 0.05, na.rm = TRUE)
      down_count <- sum(df$avg_log2FC < 0 & df$p_val_adj < 0.05, na.rm = TRUE)
    }
    deg_summary[[length(deg_summary) + 1]] <- data.frame(
      Annotation_Layer = layer, Subpopulation = sheet,
      Up_DEGs = up_count, Down_DEGs = down_count
    )
  }
}
deg_summary <- bind_rows(deg_summary)
write.xlsx(deg_summary, file.path(output_dir, "DEG_counts_summary.xlsx"), overwrite = TRUE)

deg_summary_long <- deg_summary %>%
  pivot_longer(c(Up_DEGs, Down_DEGs), names_to = "Regulation", values_to = "Count") %>%
  mutate(Regulation = factor(Regulation, levels = c("Up_DEGs","Down_DEGs")))

for (layer in annotations) {
  df_plot <- deg_summary_long %>% filter(Annotation_Layer == layer)
  p <- ggplot(df_plot, aes(x = Subpopulation, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Up_DEGs" = "#B3E0A6", "Down_DEGs" = "#F78D80")) +
    theme_classic(base_size = 12) +
    labs(title = paste("Número de DEGs CHOP vs vehicle -", layer),
         x = "Subpoblación", y = "DEGs significativos (padj < 0.05)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
  ggsave(file.path(output_dir, paste0("Barplot_DEGs_", layer, "_CHOP_vs_vehicle.png")), p, width = 6, height = 5)
}


# ================================
# 7) GO:BP EN UP-REGULADOS POR SUBPOBLACIÓN
# ================================
output_dir_go <- file.path(output_dir, "GO_Upregulated")
dir.create(output_dir_go, showWarnings = FALSE)

for (layer in annotations) {
  message("GO para: ", layer)
  file <- file.path(output_dir, paste0("DEGs_", layer, "_CHOP_vs_vehicle.xlsx"))
  sheets <- getSheetNames(file)
  wb <- createWorkbook()
  
  for (sheet in sheets) {
    df <- read.xlsx(file, sheet = sheet)
    if (nrow(df) == 0) next
    df <- standardize_deg_cols(df)
    up_genes <- df %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>% pull(gene)
    if (length(up_genes) == 0) next
    
    entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    if (is.null(entrez) || nrow(entrez) == 0) next
    
    ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05, readable = TRUE)
    addWorksheet(wb, sheet)
    writeData(wb, sheet, as.data.frame(ego))
  }
  
  go_xlsx <- file.path(output_dir_go, paste0("GO_Upregulated_", layer, ".xlsx"))
  saveWorkbook(wb, go_xlsx, overwrite = TRUE)
  
  plot_go_terms_inside_bars(
    go_xlsx = go_xlsx,
    output_dir = file.path(output_dir_go, paste0("Plots_", layer)),
    annotation_colors = annotation_colors_list[[layer]]
  )
}

# ================================
# 8) UPSET PLOTS: intersección de up-DEGs
# ================================
for (layer in annotations) {
  message("UpSet: ", layer)
  file <- file.path(output_dir, paste0("DEGs_", layer, "_CHOP_vs_vehicle.xlsx"))
  sheets <- getSheetNames(file)
  
  deg_up_list <- lapply(sheets, function(sh) {
    df <- read.xlsx(file, sheet = sh)
    if (nrow(df) == 0) return(character(0))
    df <- standardize_deg_cols(df)
    df %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>% pull(gene) %>% unique()
  })
  names(deg_up_list) <- sheets
  
  all_genes <- unique(unlist(deg_up_list))
  if (length(all_genes) == 0) next
  
  gene_df <- as.data.frame(sapply(deg_up_list, function(g) all_genes %in% g))
  gene_df$gene <- all_genes
  
  p <- ComplexUpset::upset(
    gene_df,
    intersect = names(deg_up_list),
    name = "Subpopulations",
    min_size = 5,
    width_ratio = 0.2
  ) + ggtitle(paste("UpSet plot of Upregulated DEGs -", layer)) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  ggsave(file.path(output_dir, paste0("UpSet_DEGs_", layer, "_CHOP_vs_vehicle.png")),
         p, width = 8, height = 5)
}

# ================================
# 9) CORRELACIÓN: global vs subpoblaciones
# ================================
output_dir_corr <- file.path(output_dir, "Correlation_Global_vs_Subpop")
dir.create(output_dir_corr, showWarnings = FALSE)

# Leer DEGs globales
global_deg <- read.xlsx(global_deg_file, sheet = 1)
global_deg <- standardize_deg_cols(global_deg)
global_deg_up <- global_deg %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  dplyr::select(gene, avg_log2FC)
colnames(global_deg_up)[2] <- "Global"

for (layer in annotations) {
  message("Correlación: ", layer)
  file <- file.path(output_dir, paste0("DEGs_", layer, "_CHOP_vs_vehicle.xlsx"))
  sheets <- getSheetNames(file)
  
  deg_sub_list <- list(Global = global_deg_up)
  for (sheet in sheets) {
    df <- read.xlsx(file, sheet = sheet)
    if (nrow(df) == 0) next
    df <- standardize_deg_cols(df)
    df_up <- df %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>% dplyr::select(gene, avg_log2FC)
    colnames(df_up)[2] <- sheet
    deg_sub_list[[sheet]] <- df_up
  }
  
  if (length(deg_sub_list) < 2) next
  all_degs <- Reduce(function(x, y) full_join(x, y, by = "gene"), deg_sub_list)
  all_degs[is.na(all_degs)] <- 0
  rownames(all_degs) <- all_degs$gene
  mat <- as.matrix(all_degs[,-1, drop = FALSE])
  if (ncol(mat) < 2) next
  
  cor_matrix <- cor(mat, method = "pearson")
  write.csv(cor_matrix, file.path(output_dir_corr, paste0("Correlation_", layer, ".csv")))
  
  # Triángulo inferior+diagonal
  cor_plot <- cor_matrix
  cor_plot[upper.tri(cor_plot)] <- NA
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("#1874CD", "#EEEEE0", "#CD2626"))
  hp <- Heatmap(
    cor_plot, name = "Correlation", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(col = "white"),
    row_names_side = "left",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(cor_plot[i, j])) {
        grid.text(sprintf("%.2f", cor_plot[i, j]), x, y, gp = gpar(fontsize = 10))
      }
    },
    na_col = "white"
  )
  png(file.path(output_dir_corr, paste0("Correlation_Heatmap_", layer,".png")), width = 2200, height = 2000, res = 300)
  draw(hp)
  dev.off()
}

# ================================
# 10) GO EXCLUSIVOS por subpoblación (ejemplo generalizado)
# ================================
# Toma el GO_Upregulated de cada layer y busca términos exclusivos en cada subpoblación;
# luego guarda un Excel y un barplot Top10 por subpoblación.

exclusive_dir <- file.path(output_dir, "Exclusive_GO")
dir.create(exclusive_dir, showWarnings = FALSE)

for (layer in annotations) {
  go_file <- file.path(output_dir_go, paste0("GO_Upregulated_", layer, ".xlsx"))
  if (!file.exists(go_file)) next
  sheets <- getSheetNames(go_file)
  if (length(sheets) == 0) next
  
  go_list <- lapply(sheets, function(sh) {
    df <- read.xlsx(go_file, sheet = sh)
    df$Subpopulation <- sh
    df
  })
  names(go_list) <- sheets
  
  # Términos exclusivos por subpoblación (Description único)
  exclusives <- lapply(names(go_list), function(sp) {
    this_df <- go_list[[sp]]
    others  <- unique(unlist(lapply(go_list[names(go_list) != sp], function(x) x$Description)))
    this_df %>% filter(!(Description %in% others))
  })
  names(exclusives) <- names(go_list)
  
  # Guarda Excel
  wb <- createWorkbook()
  for (sp in names(exclusives)) {
    addWorksheet(wb, paste0(sp, "_Exclusive"))
    writeData(wb, paste0(sp, "_Exclusive"), exclusives[[sp]])
  }
  excl_xlsx <- file.path(exclusive_dir, paste0("Exclusive_GO_", layer, ".xlsx"))
  saveWorkbook(wb, excl_xlsx, overwrite = TRUE)
  
  # Barplots Top10 exclusivos por subpoblación
  pal <- annotation_colors_list[[layer]]
  for (sp in names(exclusives)) {
    df <- exclusives[[sp]]
    if (is.null(df) || nrow(df) == 0) next
    top10 <- df %>%
      arrange(p.adjust) %>%
      slice_head(n = 10) %>%
      mutate(
        Term = Description,
        NegativeLog10P = -log10(p.adjust),
        ID = factor(Term, levels = rev(Term))
      )
    fill_col <- pal[ifelse(sp %in% names(pal), sp, 1)]
    p <- ggplot(top10, aes(x = ID, y = NegativeLog10P)) +
      geom_bar(stat = "identity", width = 0.6, alpha = 0.85, fill = fill_col) +
      geom_text(aes(label = Term), y = 0, hjust = 0, color = "black", size = 3.5, fontface = "bold") +
      coord_flip() +
      labs(title = paste("Top 10 Exclusive GO:BP -", sp, "(", layer, ")"),
           x = NULL, y = expression(-log[10]~"(adj. p-value)")) +
      theme_classic(base_size = 13) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
    ggsave(file.path(exclusive_dir, paste0("Exclusive_GO_", layer, "_", sp, ".png")), p, width = 8, height = 6)
  }
}

message("DONE ✅")
