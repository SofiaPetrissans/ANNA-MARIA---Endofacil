
# =========================
# LIBRERÍAS
# =========================
source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.5_Annotation_EC_mt20/0.0_Paths.R")
library(Seurat)
library(ggplot2)
library(tidyverse)

# =========================
# 1. PATHS
# =========================
MT_THRESHOLD <- "mt15"

path.guardar_original <- paths_AnnaMaria$Path_guardar
path.guardar <- file.path(path.guardar_original, "4.Senescence_Markers")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path,  MT_THRESHOLD, "8.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")
seurat_object <- readRDS(input_rds)


# Genes centrales de la senescencia
senescence_core <- c("Cdkn2a", "Bmi1", "Trp53", "Hmga1", "Chek1", "Chek2", 
                     "Prodh", "Tnfrsf10b", "Cdkn1a", "Dao")

# Genes efectores de la senescencia
senescence_effectors <- c("Ppp1ca", "Ahcy", "Brf1", "Mapa2k3", "Mapa2k6", 
                          "Smurf2", "Tgfb1i1", "Srsf1", "Angptl2")

# Genes SASP
sasp_genes <- c("Ccl2", "Ccl24", "Ccl3", "Ccl5", "Ctnnb1", "Cxcl1", "Cxcl10", 
                "Cxcl12", "Cxcl2", "Cxcl16", "Hgf", "Hmgb1", "Icam1", "Igfbp2", 
                "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6", "Igfbp7", "Il15", 
                "Il18", "Il1a", "Il1b", "Il2", "Il6", "Mif", "Mmp12", "Mmp13", 
                "Mmp14", "Pgf", "Plat", "Timp2", "Serpine1", "Ccl4", "Ang", 
                "Csf2", "Kitl", "Serpine2", "Tnfrsf1a", "Nrg1", "Ereg", "Areg")

fp <- FeaturePlot(seurat_object, features = senescence_core) & NoAxes()
ggsave(filename = file.path(path.guardar, "Senescence_core.png"), plot = fp, width = 18, height = 12)

fp <- FeaturePlot(seurat_object, features = senescence_effectors) & NoAxes()
ggsave(filename = file.path(path.guardar, "Senescence_effectors.png"), plot = fp, width = 18, height = 12)

fp <- FeaturePlot(seurat_object, features = sasp_genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "Sasp_genes.png"), plot = fp, width = 20, height = 30)


senescence_core <- c("Cdkn1a","Ctnnb1","Cxcl16", "Hmgb1",
                     "Igfbp4", "Igfbp7", "Il1a", "Mmp14", "Kitl",  "Tnfrsf1a",  
                     "Ppp1ca", "Smurf2")
fp <- FeaturePlot(seurat_object, features = senescence_core) & NoAxes()
ggsave(filename = file.path(path.guardar, "Senescence_mix.png"), plot = fp, width = 14, height = 12)

seurat_object <- SetIdent(seurat_object, value = "Annotation_Layer4")
vlp <- VlnPlot(seurat_object, features = senescence_core, split.by = "Phenotype") + theme(legend.position = "right")
ggsave(filename = file.path(path.guardar, "VlnPlot_Senescence_mix_L4.png"), plot = vlp, width = 18, height = 16)









# Asegurar que los clusters estén como identidades
seurat_object <- SetIdent(seurat_object, value = "Annotation_Layer1")

# Genes de interés
senescence_core <- c("Cdkn2a", "Bmi1", "Trp53", "Hmga1", "Chek1", "Chek2", "Prodh", "Tnfrsf10b", "Cdkn1a", "Dao")

# Extraer expresión + metadatos
df <- FetchData(seurat_object, vars = c(senescence_core, "Annotation_Layer1", "Phenotype")) %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = all_of(senescence_core), names_to = "Gene", values_to = "Expression") %>%
  mutate(Cluster = factor(Annotation_Layer1))

# Crear el plot
vln_plot <- ggplot(df, aes(x = Cluster, y = Expression, fill = Phenotype)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = TRUE) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(x = "Cluster", y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Guardar
ggsave(filename = file.path(path.guardar, "VlnPlot_Senescence_mix2.png"),
       plot = vln_plot, width = 18, height = 16)







gene_sets <- list(
  p53_signal_transduction = c("Cdkn1a", "Phlda3", "Eda2r", "Bax", "Rps27l"),
  oxidative_stress = c("Ndufb4", "Epas1", "Psap", "Hgf", "Apoe"),
  leukocyte_adhesion = c("Itpkb", "Il6st", "Chst2", "Cd59a", "Slc39a8", "Clec4g", "Hes1", "Cd81"),
  apoptotic = c("Bax", "Gata4"),
  leukocyte_migration = c("Chst2", "Itga9", "Itga1", "Calr", "Cd81"),
  leukocyte_differentiation = c("Itpkb", "Bax", "Sh3rf1", "Clec4g", "Rora", "Cyp26b1"),
  Tcells_activation = c("Itpkb", "Il6st", "Cd59a", "Sh3rf1", "Clec4g", "Cd81", "Cyp26b1"),
  leukocyte_differentiation_2 = c("Itpkb", "Sh3rf1", "Clec4g", "Cyp26b1"),
  chemotaxis = c("Itga9", "Hgf", "Itga1", "Calr")
)

# 2. Añadir los module scores
seurat_object <- AddModuleScore(seurat_object, features = gene_sets, name = "GO_Score")

# Esto genera columnas en seurat_object@meta.data como:
# GO_Score1, GO_Score2, GO_Score3, ..., en el mismo orden que la lista

# 3. Mapear nombres manualmente
score_names <- paste0("GO_Score", seq_along(gene_sets))
names(score_names) <- names(gene_sets)  # así puedes acceder por nombre original

for (sig in names(score_names)) {
  score_col <- score_names[[sig]]
  
  p <- FeaturePlot(seurat_object, features = score_col) &
       NoAxes() &
       ggtitle(sig)  # <--- aquí se asigna el nombre original como título
  
  ggsave(filename = file.path(path.guardar, paste0("Score_", sig, ".png")), 
         plot = p, width = 8, height = 8)
}







# Establecer identidad por capa
seurat_object <- SetIdent(seurat_object, value = "Annotation_Layer2")

# ViolinPlot con orden correcto y colores
library(patchwork)

vlp_list <- VlnPlot(
  seurat_object,
  features = senescence_core,
  split.by = "Phenotype",
  combine = FALSE
)

# Aplicar colores y ajustes a cada panel
vlp_list <- lapply(vlp_list, function(p) {
  p +
    scale_fill_manual(values = c("vehicle" = "#D9D9D9", "CHOP" = "#9BD0C7")) +
    theme(legend.position = "right")
})

# Combinar y guardar
vlp <- wrap_plots(vlp_list, ncol = 4)
ggsave(filename = file.path(path.guardar, "VlnPlot_Senescence_mix_L2.png"),
       plot = vlp, width = 20, height = 16)