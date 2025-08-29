
############################################
# Propósito
# Tomar los DEGs por tipo celular (obtenidos en el clúster) y realizar enriquecimiento GO:BP con clusterProfiler:
#   -  Para todos los DEGs significativos (adj. p < 0.05),
#   -  Para UP (avg_log2FC > 0),
#   -  Para DOWN (avg_log2FC < 0).
# 
# Entrada:
#   Excel con una hoja por población:
#   "/ijc/LABS/GRAUPERA/RAW/.../ANNA_MARIA/3.MERGE/0.3_Downstream_Analysis_all_populations/mt15/1.DEGs/DEGs_per_CellType_all_genes.xlsx"
# (copia local en: .../ANNA_MARIA/DEGs_per_CellType_all_genes.xlsx)
# 
# Salidas:
#   -  3 ficheros Excel (una hoja por población):
#           -  GO_BP_All_CellTypes.xlsx
#           -  GO_BP_Upregulated_CellTypes.xlsx
#           -  GO_BP_Downregulated_CellTypes.xlsx
#  -  Barplots Top-10 GO:BP por población:
#           -  Plots_GO_BP_All/ · Plots_GO_BP_Up/ · Plots_GO_BP_Down/
#           -  Versión “terms-in-bar” con el texto del término dentro de la barra:Plots_GO_BP_All_InBar/, Plots_GO_BP_Up_InBar/, Plots_GO_BP_Down_InBar/
#  -  (Opcional) Barplot para GO terms seleccionados por endotelio --> **selected_go_terms_endothelial.xlsx**.
############################################


# ================================
# 1) LIBRERÍAS
# ================================
library(openxlsx)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)


# ================================
# 2) PARÁMETROS / RUTAS
# ================================
# Directorio
directory <- "/Users/graupera_lab/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/ANNA_MARIA"
setwd(directory)

# Excel con una hoja por población (descargado del clúster)
xlsx_file <- "DEGs_per_CellType_all_genes.xlsx"

# Colores por población (ajustar si cambian las hojas)
Annotation_Colors_L1 <- c(
  "NK"          = "#B29CA6",
  "Endothelial" = "#B5DDC3",
  "Neutrophils" = "#CEDD9A",
  "Kupffer"     = "#FF9F9B",
  "B_Cells"     = "#BDBAD7",
  "T_Cells"     = "#DEBB9B",
  "HSC"         = "#FAB0E4",
  "Hepatocytes" = "#BADCDE",
  "Monocytes"   = "#E0A0A7"
)


# ================================
# 3) LECTURA DE DEGs (todas las hojas)
# ================================
sheet_names <- getSheetNames(xlsx_file)
deg_list <- lapply(sheet_names, function(sheet) read.xlsx(xlsx_file, sheet = sheet))
names(deg_list) <- sheet_names


# ================================
# 4) FUNCIONES
# ================================
write_go_list <- function(go_list, file_path) {
  wb <- createWorkbook()
  for (ct in names(go_list)) {
    addWorksheet(wb, ct)
    writeData(wb, sheet = ct, go_list[[ct]])
  }
  saveWorkbook(wb, file_path, overwrite = TRUE)
}



plot_go_terms_inside_bars <- function(go_list, name_prefix, output_dir, annotation_colors) {
  dir.create(output_dir, showWarnings = FALSE)
  
  for (ct in names(go_list)) {
    df <- go_list[[ct]]
    if (nrow(df) == 0) next
    
    top10 <- df %>%
      arrange(p.adjust) %>%
      slice_head(n = 10) %>%
      mutate(
        Term = Description,
        NegativeLog10P = -log10(p.adjust),
        ID = factor(Term, levels = rev(Term))  # <- aquí el cambio
      )
    
    
    p <- ggplot(top10, aes(x = ID, y = NegativeLog10P, fill = ct)) +
      geom_bar(stat = "identity", width = 0.6, alpha = 0.7) +
      geom_text(aes(label = Term), 
                y = 0,
                hjust = 0,
                color = "black",
                size = 4,
                fontface = "bold") +
      scale_fill_manual(values = annotation_colors, guide = "none") +
      coord_flip() +
      labs(
        title = paste("Top 10 GO:BP -", ct),
        x = NULL,
        y = expression(-log[10]~"(adj. p-value)")
      ) +
      theme_classic(base_size = 13) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
      )
    
    ggsave(filename = file.path(output_dir, paste0(name_prefix, "_", ct, "_terms_in_bars.png")),
           plot = p, width = 7, height = 6)
  }
}


# ================================
# 5) PIPELINE: construir listas ALL / UP / DOWN y correr GO
# ================================
# Crear lista donde guardar resultados
go_results_all <- list()
go_results_up <- list()
go_results_down <- list()

for (cell_type in names(deg_list)) {
  df <- deg_list[[cell_type]]
  
  if (!"avg_log2FC" %in% colnames(df)) next
  
  df_sig <- df %>% filter(p_val_adj < 0.05)
  genes_all <- df_sig$gene
  genes_up <- df_sig %>% filter(avg_log2FC > 0) %>% pull(gene)
  genes_down <- df_sig %>% filter(avg_log2FC < 0) %>% pull(gene)
  
  # Convertir a ENTREZ IDs
  convert_genes <- function(gene_symbols) {
    bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
  }
  
  genes_all_entrez <- convert_genes(genes_all)
  genes_up_entrez <- convert_genes(genes_up)
  genes_down_entrez <- convert_genes(genes_down)
  
  # GO:BP con clusterProfiler
  go_all <- enrichGO(gene = genes_all_entrez, OrgDb = org.Mm.eg.db,
                     ont = "BP", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  go_up <- enrichGO(gene = genes_up_entrez, OrgDb = org.Mm.eg.db,
                    ont = "BP", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  go_down <- enrichGO(gene = genes_down_entrez, OrgDb = org.Mm.eg.db,
                      ont = "BP", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  go_results_all[[cell_type]] <- as.data.frame(go_all)
  go_results_up[[cell_type]] <- as.data.frame(go_up)
  go_results_down[[cell_type]] <- as.data.frame(go_down)
}

# ================================
# 6) GUARDAR RESULTADOS (Excel)
# ================================
write_go_list(go_results_all, "GO_BP_All_CellTypes.xlsx")
write_go_list(go_results_up, "GO_BP_Upregulated_CellTypes.xlsx")
write_go_list(go_results_down, "GO_BP_Downregulated_CellTypes.xlsx")


# ================================
# 7) PLOTS: Top-10 por población
# ================================

plot_go_terms_inside_bars(go_results_all, "GO_BP_All", "Plots_GO_BP_All_InBar", Annotation_Colors_L1)
plot_go_terms_inside_bars(go_results_up,  "GO_BP_Up",  "Plots_GO_BP_Up_InBar",  Annotation_Colors_L1)
plot_go_terms_inside_bars(go_results_down,"GO_BP_Down","Plots_GO_BP_Down_InBar",Annotation_Colors_L1)


# ================================
# 8) (OPCIONAL) GO TERMS SELECCIONADOS PARA ENDOTELIO
# ================================
# Si preparaste un Excel con los términos elegidos para endotelio:
#   columnas típicas: Description, p.adjust (u otras)
selected_terms_path <- "selected_go_terms_endothelial.xlsx"
selected_goTerms_df <- read_xlsx(selected_terms_path)

# Crear la lista con un solo elemento llamado "Endothelial"
selected_goTerms_list <- list("Endothelial" = selected_goTerms_df)

Annotation_Colors_EC <- c("Endothelial" = "#B5DDC3")

plot_go_terms_inside_bars(
  go_list = selected_goTerms_list,
  name_prefix = "GO_BP_All",
  output_dir = "Plot_GO_BP_Endothelial",
  annotation_colors = Annotation_Colors_EC
)




