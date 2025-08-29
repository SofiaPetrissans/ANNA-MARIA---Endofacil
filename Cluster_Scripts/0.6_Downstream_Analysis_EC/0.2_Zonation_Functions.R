# --------------------------------------------------------------------------------
# SCRIPT: Zonación hepática con GSVA y UCell (Mus musculus)
# OBJETO: seurat_object
# --------------------------------------------------------------------------------

# =========================
# 1. LIBRERÍAS
# =========================
source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.6_Downstream_Analysis_EC/0.0_Paths.R")
library(Seurat)
library(GSVA)
library(msigdbr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(UCell)
library(openxlsx)
library(ggplot2)


MT_THRESHOLD <- "mt20"

path.guardar_original <- paths_AnnaMaria$Path_guardar
path.guardar <- file.path(path.guardar_original, MT_THRESHOLD, "2.Zonation_functions")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path, "Seurat_annotated.rds")
seurat_object <- readRDS(input_rds)
seurat_object <- SetIdent(seurat_object, value = "Annotation_Layer2")

# =========================
# 2. HALLMARK GENE SETS (ratón)
# =========================
hallmark_mouse <- msigdbr(species = "Mus musculus", category = "H")

relevant_sets <- c(
  "HALLMARK_GLYCOLYSIS", "HALLMARK_LIPID_METABOLISM", "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_P53_PATHWAY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_COMPLEMENT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

hallmark_subset <- hallmark_mouse %>% 
  filter(gs_name %in% relevant_sets) %>% 
  distinct(gs_name, gene_symbol)

gene_sets <- split(hallmark_subset$gene_symbol, hallmark_subset$gs_name)

# =========================
# 3. GSVA por grupo de anotación
# =========================
avg_expr <- AverageExpression(seurat_object, assays = "RNA", slot = "data", group.by = "Annotation_Layer2")$RNA

gsva_mat <- gsva(as.matrix(avg_expr), gene_sets, method = "gsva", kcdf = "Gaussian")

# Diagnóstico (opcional)
missing <- setdiff(relevant_sets, rownames(gsva_mat))
if (length(missing) > 0) {
  warning("Los siguientes pathways no están en el GSVA matrix (posiblemente por falta de genes expresados):\n",
          paste(missing, collapse = ", "))
}

# Filtrar y reordenar según el orden original
relevant_sets_present <- relevant_sets[relevant_sets %in% rownames(gsva_mat)]
gsva_mat <- gsva_mat[relevant_sets_present, ]

# Anotación para destacar vías de interés
highlight_terms <- relevant_sets_present  # o solo un subconjunto si prefieres
row_anno <- rowAnnotation(
  Highlight = anno_simple(
    rownames(gsva_mat) %in% highlight_terms,
    col = c("TRUE" = "black", "FALSE" = NA)
  )
)

# Guardar resultado
write.xlsx(as.data.frame(gsva_mat), file = file.path(path.guardar, "GSVA_Zonation_Scores_Layer2.xlsx"), rowNames = TRUE)

# Heatmap
png(file.path(path.guardar, "GSVA_Zonation_Heatmap_Layer2.png"), width = 2000, height = 1500, res = 300)
Heatmap(gsva_mat,
        name = "GSVA score",
        cluster_columns = FALSE,  # Desactiva si quieres mantener el orden de Annotation_Layer2
        cluster_rows = FALSE,
        right_annotation = row_anno,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_title = "Annotation Layer 2",
        row_title = "Zonation pathways")
dev.off()

# =========================
# 4. UCell por célula (opcional)
# =========================
seurat_object <- AddModuleScore_UCell(seurat_object, features = gene_sets)

fp <- FeaturePlot(seurat_object, features = paste0(names(gene_sets), "_UCell"), ncol = 3) +
  ggtitle("Functional zonation UCell scores") & NoAxes()
ggsave(filename = file.path(path.guardar, "Zonation_UCell.png"), plot = fp, width = 12, height = 10)
