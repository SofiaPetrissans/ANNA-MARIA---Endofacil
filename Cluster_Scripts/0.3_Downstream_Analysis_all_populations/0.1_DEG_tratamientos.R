# --------------------------------------------------------------------------------
# SCRIPT: Generar DEGs entre tratamientos - CHOP y vehicle 
# AUTORA: Sofia Petrissans Moll
# FECHA: 30/06/2025
# --------------------------------------------------------------------------------

# LIBRARY
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.3_Downstream_Analysis_all_populations/0.0_Paths.R")
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(openxlsx)

MT_THRESHOLD <- "mt15"


# Paths
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "3.Annotation_Layer", MT_THRESHOLD, "Seurat_annotated.rds")

path.guardar_original <- paths_AnnaMaria$Path_guardar
output_dir <- file.path(path.guardar_original, MT_THRESHOLD, "1.DEGs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------------
# CARGA Y PREPARACIÓN
# --------------------------------------------------------------------------------
seurat_object <- readRDS(input_rds)
seurat_object <- SetIdent(seurat_object, value = "Annotation_Layer1")


# Define colores según el MT_THRESHOLD
if (MT_THRESHOLD == "mt20") {
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
} else if (MT_THRESHOLD == "mt15") {
  Annotation_Colors_L1 <- c(
    "NK" = "#B29CA6",
    "Endothelial" = "#B5DDC3",
    "Neutrophils" = "#CEDD9A", 
    "Kupffer" = "#ff9f9b",
    "B Cells" = "#BDBAD7",
    "T Cells" = "#debb9b",
    "Fibroblast" = "#fab0e4",
    "Cholangiocytes" = "#BADCDE"
  )
} else {
  stop("No se han definido colores para este valor de MT_THRESHOLD.")
}


# Crear resumen de número de DEGs
deg_summary <- tibble()
deg_list <- list()  # aquí guardaremos los DEGs completos

cell_types <- unique(seurat_object$Annotation_Layer1)

for (cell_type in cell_types) {
  message("Procesando: ", cell_type)
  
  subset_obj <- subset(seurat_object, subset = Annotation_Layer1 == cell_type)
  cond_counts <- table(subset_obj$Phenotype)
  
  if (all(c("CHOP", "vehicle") %in% names(cond_counts)) && all(cond_counts >= 10)) {
    Idents(subset_obj) <- "Phenotype"
    markers <- FindMarkers(subset_obj, ident.1 = "CHOP", ident.2 = "vehicle", 
                           logfc.threshold = 0.25, min.pct = 0.1)
    
    # Añadir columna con gene name
    markers <- markers %>% rownames_to_column("gene")
    
    # Guardar tabla en lista
    deg_list[[cell_type]] <- markers
    
    # Contar DEGs significativos
    n_degs <- sum(markers$p_val_adj < 0.05)
    
  } else {
    n_degs <- NA
    deg_list[[cell_type]] <- data.frame(gene = NA, warning = "Insuficientes células para comparar")
  }
  
  deg_summary <- bind_rows(deg_summary, tibble(CellType = cell_type, DEGs = n_degs))
}

# Añadir colores
deg_summary$Color <- Annotation_Colors_L1[as.character(deg_summary$CellType)]

# Guardar resumen .xlsx
write.xlsx(deg_summary, file = file.path(output_dir, "DEGs_per_CellType_summary.xlsx"), rowNames = FALSE)

# Guardar DEGs por hoja en un solo .xlsx
wb <- createWorkbook()
for (cell_type in names(deg_list)) {
  sheet_name <- gsub(" ", "_", cell_type)  # <-- CONSISTENTE
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = deg_list[[cell_type]])  # <-- MISMO NOMBRE
}
saveWorkbook(wb, file = file.path(output_dir, "DEGs_per_CellType_all_genes.xlsx"), overwrite = TRUE)

# Barplot coloreado
p <- ggplot(deg_summary, aes(x = reorder(CellType, DEGs), y = DEGs, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Annotation_Colors_L1) +
  coord_flip() +
  labs(x = "Cell Type", y = "Number of DEGs (CHOP vs vehicle)", 
       title = "DEGs per Cell Type") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave(filename = file.path(output_dir, "Barplot_DEGs_per_CellType.png"),
       plot = p, width = 6.5, height = 5)



# Crear tabla para up/down DEGs por tipo celular
deg_updown <- tibble()

for (cell_type in names(deg_list)) {
  df <- deg_list[[cell_type]]
  
  if (!"gene" %in% colnames(df)) next  # si no hay DEGs, continuar
  
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

# p2: Stacked barplot
p2 <- ggplot(deg_updown, aes(x = CellType, y = DEGs, fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Up" = "#B3E0A6", "Down" = "#F78D80")) +
  coord_flip() +
  labs(x = "Cell Type", y = "Number of DEGs", title = "Up/Downregulated DEGs per Cell Type") +
  theme_classic(base_size = 13)

ggsave(filename = file.path(output_dir, "Barplot_DEGs_UpDown_Stacked.png"),
       plot = p2, width = 7, height = 5)

# p3: Split facets, con colores por población
p3 <- ggplot(deg_updown, aes(x = CellType, y = DEGs, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Direction, scales = "free_y") +
  scale_fill_manual(values = Annotation_Colors_L1) +
  coord_flip() +
  labs(x = "Cell Type", y = "Number of DEGs", title = "Up vs Downregulated DEGs") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave(filename = file.path(output_dir, "Barplot_DEGs_UpDown_SplitFacets.png"),
       plot = p3, width = 8, height = 5)




# Asegúrate de estar trabajando con la matriz de counts
DefaultAssay(seurat_object) <- "RNA"

gene_count_summary <- tibble()

for (cell_type in cell_types) {
  for (cond in c("CHOP", "vehicle")) {
    subset_cells <- subset(seurat_object, subset = Annotation_Layer1 == cell_type & Phenotype == cond)
    
    if (ncol(subset_cells) >= 1) {
      counts_mat <- GetAssayData(subset_cells, slot = "counts")
      # Número de genes con al menos 10 counts sumados en todas las células
      genes_detected <- Matrix::rowSums(counts_mat) >= 10
      n_genes <- sum(genes_detected)
    } else {
      n_genes <- NA
    }

    gene_count_summary <- bind_rows(gene_count_summary, tibble(
      CellType = cell_type,
      Condition = cond,
      GenesDetected = n_genes
    ))
  }
}

# Guardar tabla
write.xlsx(gene_count_summary, file = file.path(output_dir, "GeneCounts_Above10_PerCondition.xlsx"))

# Asegurar el orden de los niveles en la variable Condition
gene_count_summary$Condition <- factor(gene_count_summary$Condition, levels = c("vehicle", "CHOP"))

# Gráfico con los nuevos colores y orden
p4 <- ggplot(gene_count_summary, aes(x = CellType, y = GenesDetected, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("vehicle" = "#D9D9D9", "CHOP" = "#9BD0C7")) +
  coord_flip() +
  labs(x = "Cell Type", y = "Genes with more than 10 counts", title = "Detected genes per cell type and condition") +
  theme_classic(base_size = 13)

ggsave(filename = file.path(output_dir, "Barplot_GeneCounts_Above10.png"),
       plot = p4, width = 8, height = 5)



# Calcular proporciones
deg_props <- deg_updown %>%
  group_by(CellType) %>%
  mutate(Proportion = DEGs / sum(DEGs))

# Plot de proporciones
p5 <- ggplot(deg_props, aes(x = CellType, y = Proportion, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#B3E0A6", "Down" = "#F78D80")) +
  labs(x = "Cell Type", y = "Proportion of DEGs", title = "Proportion of Up vs Down DEGs per Cell Type") +
  theme_classic(base_size = 13)

ggsave(filename = file.path(output_dir, "Barplot_DEGs_Proportions.png"),
       plot = p5, width = 7, height = 5)



library(scales)  # Asegúrate de tener esta librería

# Calcular proporciones
deg_props <- deg_updown %>%
  group_by(CellType) %>%
  mutate(Proportion = DEGs / sum(DEGs))

# Plot de proporciones con eje en porcentaje
p5 <- ggplot(deg_props, aes(x = CellType, y = Proportion, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # mostrar como 0–100%
  scale_fill_manual(values = c("Up" = "#B3E0A6", "Down" = "#F78D80")) +
  labs(x = "Cell Type", y = "Proportion of DEGs (%)", title = "Proportion of Up vs Down DEGs per Cell Type") +
  theme_classic(base_size = 13)

ggsave(filename = file.path(output_dir, "Barplot_DEGs_Proportions_Percent.png"),
       plot = p5, width = 7, height = 5)

