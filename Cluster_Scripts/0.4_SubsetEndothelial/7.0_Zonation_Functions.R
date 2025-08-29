# --------------------------------------------------------------------------------
# SCRIPT: Zonación hepática con GSVA y UCell (Mus musculus)
# OBJETO: seurat_object
# --------------------------------------------------------------------------------

# =========================
# 1. LIBRERÍAS
# =========================
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.4_SubsetEndothelial/0.0_Paths.R")
library(Seurat)
library(GSVA)
library(msigdbr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(UCell)
library(openxlsx)
library(ggplot2)

path.guardar <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/0.2_SubsetEndothelial/mt20/0.Pruebas"
seurat_object <- readRDS("/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/0.2_SubsetEndothelial/mt20/4.RemoveCluster/2.Integration/SubsetEndothelial_Harmony.rds")
seurat_object <- SetIdent(seurat_object, value="Harmony_Log_res.0.5")
# =========================
# 2. HALLMARK GENE SETS (ratón)
# =========================
hallmark_mouse <- msigdbr(species = "Mus musculus", category = "H")

# Firmas funcionales periportal/pericentral
relevant_sets <- c(
  "HALLMARK_GLYCOLYSIS",                    # pericentral
  "HALLMARK_LIPID_METABOLISM",              # pericentral
  "HALLMARK_XENOBIOTIC_METABOLISM",         # pericentral
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",    # pericentral
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",     # periportal
  "HALLMARK_FATTY_ACID_METABOLISM",         # periportal
  "HALLMARK_BILE_ACID_METABOLISM"           # periportal
)

hallmark_subset <- hallmark_mouse %>% 
  filter(gs_name %in% relevant_sets) %>% 
  distinct(gs_name, gene_symbol)

# Lista usable
gene_sets <- split(hallmark_subset$gene_symbol, hallmark_subset$gs_name)

# =========================
# 3. GSVA por CLUSTER
# =========================
avg_expr <- AverageExpression(seurat_object, assays = "RNA", slot = "data")$RNA
gsva_mat <- gsva(as.matrix(avg_expr), gene_sets, method = "gsva", kcdf = "Gaussian")

# HEATMAP estilo CAF
highlight_terms <- c("HALLMARK_GLYCOLYSIS", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", 
                     "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM")

row_anno <- rowAnnotation(
  Highlight = anno_simple(
    rownames(gsva_mat) %in% highlight_terms,
    col = c("TRUE" = "black", "FALSE" = NA)
  )
)

write.xlsx(as.data.frame(gsva_mat), file = file.path(path.guardar, "GSVA_Zonation_Scores.xlsx"), rowNames = TRUE)

png(file.path(path.guardar, "GSVA_Zonation_Heatmap.png"), width = 2000, height = 1500, res = 300)
Heatmap(gsva_mat,
        name = "GSVA score",
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        right_annotation = row_anno,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_title = "Seurat clusters",
        row_title = "Zonation pathways")
dev.off()
# =========================
# 4. UCELL por CÉLULA
# =========================
seurat_object <- AddModuleScore_UCell(seurat_object, features = gene_sets)

# Plot UMAP para cada pathway
fp <- FeaturePlot(seurat_object, features = paste0(names(gene_sets), "_UCell"), ncol = 3) +
  ggtitle("Functional zonation UCell scores") & NoAxes()
ggsave(filename = file.path(path.guardar, "Zonation_UCell.png"), plot = fp, width = 24, height = 10)