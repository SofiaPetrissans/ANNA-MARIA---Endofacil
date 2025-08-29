
############################################
# Propósito
# A partir de una tabla de **DEGs (CHOP vs vehicle)** para el subset de células endoteliales, realizar enriquecimiento GO:BP (clusterProfiler) para:
#   -  todos los genes significativos (adj. p < 0.05),
#   -  solo UP (avg_log2FC > 0),
#   -  solo DOWN (avg_log2FC < 0),
# 
# y generar barplots Top-10 (estilo “terms-in-bar”).
# 
# Entrada:
#   -  Excel con DEGs endoteliales: `DEGs_CHOP_vs_vehicle.xlsx` (una hoja; debe incluir columnas gene, avg_log2FC, p_val_adj).
#     ("/ijc/LABS/GRAUPERA/RAW/.../ANNA_MARIA/3.MERGE/0.6_Downstream_Analysis_EC/mt15/1.DEGs_phenotype/DEGs_CHOP_vs_vehicle.xlsx")
# Salidas
#   -  Excel:
#           -  GO_BP_All_Subset.xlsx
#           -  GO_BP_Up_Subset.xlsx
#           -  GO_BP_Down_Subset.xlsx
# 
# Barplots (Top 10 términos):
#   -  Plots_Subset_GO_BP_All_InBar/GO_BP_All_terms_in_bars.png
#   -  Plots_Subset_GO_BP_Up_InBar/GO_BP_Up_terms_in_bars.png
#   -  Plots_Subset_GO_BP_Down_InBar/GO_BP_Down_terms_in_bars.png
#   -  (Opcional) Plot de términos seleccionados desde un Excel externo --> **selected_go_terms_endothelial_subset.xlsx**
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
library(scales)

# ================================
# 2) PARÁMETROS / RUTAS
# ================================
# Directorio
directory <- "/Users/graupera_lab/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/ANNA_MARIA"
setwd(directory)

# Excel con una hoja por población (descargado del clúster)
xlsx_file <- "DEGs_CHOP_vs_vehicle.xlsx"

outdir     <- "2.GOTerms_mt15_Subset_EC_merge"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Color para plots (endotelio)
color_plots<- "#B5DDC3"


# ================================
# 3) LECTURA DE DEGs 
# ================================
# Leer directamente el Excel (una sola hoja)
deg_df <- read.xlsx(xlsx_file)


# ================================
# 4) FUNCIONES
# ================================
# Convertir a ENTREZ IDs
convert_genes <- function(gene_symbols) {
  bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
}

# BarPlots esteticos go terms 
plot_go_terms_inside_bars <- function(df, name_prefix, output_dir, color = color_plots) {
  dir.create(output_dir, showWarnings = FALSE)
  
  if (nrow(df) == 0) return(NULL)
  
  top10 <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(
      Term = Description,
      NegativeLog10P = -log10(p.adjust),
      ID = factor(Term, levels = rev(Term))
    )
  
  p <- ggplot(top10, aes(x = ID, y = NegativeLog10P)) +
    geom_bar(stat = "identity", fill = color, width = 0.6, alpha = 0.8) +
    geom_text(aes(label = Term), y = 0, hjust = 0, color = "black",
              size = 4, fontface = "bold") +
    coord_flip() +
    labs(
      title = paste("Top 10 GO:BP -", name_prefix),
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
  
  ggsave(filename = file.path(output_dir, paste0(name_prefix, "_terms_in_bars.png")),
         plot = p, width = 7, height = 6)
}



# ================================
# 5) PIPELINE
# ================================

# Filtrar genes significativos
deg_df_sig <- deg_df %>% filter(p_val_adj < 0.05)
genes_all  <- deg_df_sig$gene
genes_up   <- deg_df_sig %>% filter(avg_log2FC > 0) %>% pull(gene)
genes_down <- deg_df_sig %>% filter(avg_log2FC < 0) %>% pull(gene)

# Convertir a ENTREZ
genes_all_entrez  <- convert_genes(genes_all)
genes_up_entrez   <- convert_genes(genes_up)
genes_down_entrez <- convert_genes(genes_down)


# Enrichment
go_all <- enrichGO(gene = genes_all_entrez, OrgDb = org.Mm.eg.db,
                   ont = "BP", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)

go_up <- enrichGO(gene = genes_up_entrez, OrgDb = org.Mm.eg.db,
                  ont = "BP", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)

go_down <- enrichGO(gene = genes_down_entrez, OrgDb = org.Mm.eg.db,
                    ont = "BP", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)


# ================================
# 6) GUARDAR TABLAS GO
# ================================
# Guardar resultados
write.xlsx(as.data.frame(go_all),  "GO_BP_All_Subset.xlsx",  rowNames = FALSE)
write.xlsx(as.data.frame(go_up),   "GO_BP_Up_Subset.xlsx",   rowNames = FALSE)
write.xlsx(as.data.frame(go_down), "GO_BP_Down_Subset.xlsx", rowNames = FALSE)


# ================================
# 7) PLOTS Top-10 (terms-in-bar)
# ================================
plot_go_terms_inside_bars(as.data.frame(go_all),  "GO_BP_All",  "Plots_Subset_GO_BP_All_InBar")
plot_go_terms_inside_bars(as.data.frame(go_up),   "GO_BP_Up",   "Plots_Subset_GO_BP_Up_InBar")
plot_go_terms_inside_bars(as.data.frame(go_down), "GO_BP_Down", "Plots_Subset_GO_BP_Down_InBar")



# ================================
# 8) (OPCIONAL) TÉRMINOS SELECCIONADOS
# ================================
# GO TERMS SELECCIONADOS .......................................................
# Leer la hoja (por defecto si es la única, no necesitas especificarla)
selected_terms_path <- "selected_go_terms_endothelial_subset.xlsx"
selected_goTerms_df <- read.xlsx(selected_terms_path, sheet = 1)

# Crear la lista con un solo elemento llamado "Endothelial"
selected_goTerms_list <- list("Endothelial" = selected_goTerms_df) 

plot_go_terms_inside_bars(
  df = selected_goTerms_list$Endothelial,
  name_prefix = "GO_BP",
  output_dir = "Plot_GO_BP_Endothelial_Refinado"
)
