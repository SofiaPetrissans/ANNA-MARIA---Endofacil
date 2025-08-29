############################################
# Objetivo. Identificar funciones biológicas (GO:BP) que definen cada subpoblación endotelial en Annotation_Layer5 (LVEC y LSEC_1…LSEC_5), partiendo de DEGs 1-vs-rest (solo up-regulated), y resumirlas mediante:
#   
#   -  Barplots de términos top por subpoblación.
#   -  Especificidad (términos compartidos vs específicos).
#
# Entrada.
#   -  `SubsetEndothelial_Harmony_Annotated.rds` (con Annotation_Layer5, Phenotype, etc.)
# 
# Salidas (se crean si no existen):
#   -  10.Results_GO_AnnotationLayer/GO_Results_AnnotationLayer5.xlsx (una hoja por subpoblación).
#   -  10.Results_GO_AnnotationLayer/GOBarplot_* (barplots por subpoblación).
############################################


# ================================
# 1) LIBRERÍAS
# ================================
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)  # Suponiendo que es mouse
library(openxlsx)
library(dplyr)
library(ggplot2)
library(rrvgo)
library(GOSemSim)
library(ComplexHeatmap)



# ================================
# 2) PARÁMETROS / RUTAS
# ================================
# Directorio
directory <- "/Users/graupera_lab/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/ANNA_MARIA"
setwd(directory)

input_rds <- file.path("SubsetEndothelial_Harmony_Annotated.rds")

out_go   <- file.path("4.GO_Terms_AnnotationLayer5")
out_sem  <- file.path("5.RRVGO_AnnotationLayer")
dir.create(out_go,  showWarnings = FALSE)
dir.create(out_sem, showWarnings = FALSE)



# ================================
# 3) LECTURA DEL SEURAT OBJECT 
# ================================
seurat_object <- readRDS(input_rds)

annotation.colors <- c(
  "LVEC" = "#86D0B9", 
  "LSEC_1" = "#F4CDA5",
  "LSEC_2" = "#EBAFC3",
  "LSEC_3" = "#A8C8E5",
  "LSEC_4" = "#dce6b6",
  "LSEC_5" = "#C5A7E1"
)



# ================================
# 4) FUNCIONES
# ================================
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
    
    # Debug opcional
    print(paste("Usando grupo:", unique(top10$Group)))
    
    p <- ggplot(top10, aes(x = ID, y = NegativeLog10P, fill = Group)) +
      geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
      geom_text(aes(label = Term), y = 0, hjust = 0, color = "black", size = 4, fontface = "bold") +
      scale_fill_manual(values = annotation_colors, guide = "none") +
      coord_flip() +
      labs(
        title = paste("Top 10 GO:BP -", gsub("__", ": ", sheet)),
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
    
    filename <- paste0("GOBarplot_", sheet, "_terms_in_bars.png")
    ggsave(filename = file.path(output_dir, filename), plot = p, width = 8, height = 6)
  }
}

