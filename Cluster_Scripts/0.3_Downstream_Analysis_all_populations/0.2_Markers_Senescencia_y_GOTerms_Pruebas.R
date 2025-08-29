# LIBRERIAS
library(Seurat)
library(ggplot2)
library(tidyverse)

# CARGAR OBJETO
seurat_object <- readRDS("/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/0.1_SeuratPipeline/9.Annotation_Layer/mt20/Seurat_annotated.rds")
# PATH GUARDAR
path.guardar <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/0.1_SeuratPipeline/10.Pruebas_Markers"


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

seurat_object <- SetIdent(seurat_object, value = "Annotation_Layer1")
vlp <- VlnPlot(seurat_object, features = senescence_core, split.by = "Phenotype") + theme(legend.position = "right")
ggsave(filename = file.path(path.guardar, "VlnPlot_Senescence_mix.png"), plot = vlp, width = 18, height = 16)



# Subset a endotelio ..........................................
endothelial_obj <- subset(seurat_object, subset = Annotation_Layer1 == "Endothelial")

# Reordenar niveles de Phenotype
endothelial_obj$Phenotype <- factor(endothelial_obj$Phenotype, levels = c("vehicle", "CHOP"))
Idents(endothelial_obj) <- "Phenotype"
library(patchwork)

plots <- VlnPlot(endothelial_obj, features = senescence_core, combine = FALSE)
plots <- lapply(plots, function(p) {
  p + scale_fill_manual(values = c("vehicle" = "#D9D9D9", "CHOP" = "#9BD0C7")) +
    theme(legend.position = "right")
})

vlp_endothelial <- wrap_plots(plots, ncol = 4)
ggsave(filename = file.path(path.guardar, "VlnPlot_Senescence_Endothelial_colored.png"),
       plot = vlp_endothelial, width = 20, height = 16)
# .............................................................................


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






# GENES LISTAS TOP 10 GO TERMS ......................................................................
genes <- c("Cdkn1a", "Phlda3", "Eda2r", "Bax", "Rps27l", "Rpl23", "Atrx", "Rpl11", "Ubb", "Mdm2", "Bbc3", "Aen", "Rps20", "Sesn2", "Rps7", "Rpl5")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_p53.png"), plot = fp, width = 18, height = 20)

genes <- c("Ndufb4", "Epas1", "Psap", "Hspa8", "Hgf", "Apoe", "Mdm2", "Trp53inp1", "Ppia", "Sesn2", "Anxa1", "Txn1", "Rps3", "Rack1", "Rbx1")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_oxidativeStress.png"), plot = fp, width = 18, height = 20)

genes <- c("S100a9", "Itpkb", "Il6st", "S100a8", "Chst2", "Cd59a", "Slc39a8", "Hsp90aa1", "Lgals1", "Clec4g", "Hes1", "Cd81", "Anxa1", "Rps3", "Lrg1")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_adhesion.png"), plot = fp, width = 18, height = 20)

genes <- c("S100a9", "Bax", "S100a8", "Tpt1", "Serinc3", "Gata4", "Ubb", "Mdm2", "Bbc3", "Ppia", "Ackr3", "Rps7", "Rps3", "Rack1")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_apoptotic.png"), plot = fp, width = 18, height = 20)

genes <- c("S100a9", "Camk1d", "S100a8", "Chst2", "Rps19", "Itga9", "Itga1", "Retnlg", "Calr", "Cd81", "Ppia", "Anxa1", "Cd9", "Cd34")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_migration.png"), plot = fp, width = 18, height = 20)

genes <- c("Igkc", "Itpkb", "Bax", "Psap", "Hsp90aa1", "Lgals1", "Id2", "Sh3rf1", "Clec4g", "Rora", "Ctla2a", "Anxa1", "Cyp26b1", "Rps6")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_diff.png"), plot = fp, width = 18, height = 20)

genes <- c("Itpkb", "Il6st", "Cd59a", "Hsp90aa1", "Lgals1", "Sh3rf1", "Clec4g", "Hes1", "Cd81", "Ctla2a", "Anxa1", "Cyp26b1", "Rps3")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_Tcells_activation.png"), plot = fp, width = 18, height = 20)

genes <- c("Itpkb", "Hsp90aa1", "Lgals1", "Id2", "Sh3rf1", "Clec4g", "Ctla2a", "Anxa1", "Tmem176a", "Tmem176b", "Cyp26b1", "Nme2")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_diff2.png"), plot = fp, width = 18, height = 20)

genes <- c("S100a9", "Camk1d", "S100a8", "Rps19", "Itga9", "Hgf", "Itga1", "Retnlg", "Calr", "Ppia", "Ackr3", "Anxa1")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_chemo.png"), plot = fp, width = 18, height = 20)




# SELECCION DE CADA UNA DE LA LISTA DE LOS TOP 10 GO TERMS ......................................................................
genes <- c("Cdkn1a", "Phlda3", "Eda2r", "Bax", "Rps27l", "Rpl23", "Atrx", "Rpl11", "Ubb", "Rps20", "Rps7", "Rpl5")
genes <- c("Cdkn1a", "Phlda3", "Eda2r", "Bax", "Rps27l")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_p53.png"), plot = fp, width = 12, height = 16)

genes <- c("Ndufb4", "Epas1", "Psap", "Hgf", "Apoe", "Ppia", "Rps3", "Rack1", "Rbx1")
genes <- c("Ndufb4", "Epas1", "Psap", "Hgf", "Apoe")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_oxidativeStress.png"), plot = fp, width = 12, height = 16)

genes <- c("Itpkb", "Il6st", "Chst2", "Cd59a", "Slc39a8", "Lgals1", "Clec4g", "Hes1", "Cd81", "Rps3")
genes <- c("Itpkb", "Il6st", "Chst2", "Cd59a", "Slc39a8", "Clec4g", "Hes1", "Cd81")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_adhesion.png"), plot = fp, width = 18, height = 20)

genes <- c("Bax", "Tpt1", "Serinc3", "Gata4", "Ubb", "Ppia", "Rps7", "Rps3", "Rack1")
genes <- c("Bax", "Gata4")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_apoptotic.png"), plot = fp, width = 10, height = 6)

genes <- c("Chst2", "Rps19", "Itga9", "Itga1", "Calr", "Cd81", "Ppia")
genes <- c("Chst2", "Itga9", "Itga1", "Calr", "Cd81")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_migration.png"), plot = fp, width = 12, height = 16)

genes <- c("Itpkb", "Bax", "Psap", "Lgals1", "Sh3rf1", "Clec4g", "Rora", "Cyp26b1", "Rps6")
genes <- c("Itpkb", "Bax", "Sh3rf1", "Clec4g", "Rora", "Cyp26b1")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_diff.png"), plot = fp, width = 12, height = 16)

genes <- c("Itpkb", "Il6st", "Cd59a", "Lgals1", "Sh3rf1", "Clec4g", "Cd81", "Cyp26b1", "Rps3")
genes <- c("Itpkb", "Il6st", "Cd59a", "Sh3rf1", "Clec4g", "Cd81", "Cyp26b1")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_Tcells_activation.png"), plot = fp, width = 18, height = 20)

genes <- c("Itpkb", "Lgals1", "Id2", "Sh3rf1", "Clec4g", "Cyp26b1", "Nme2")
genes <- c("Itpkb", "Sh3rf1", "Clec4g", "Cyp26b1")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_leuco_diff2.png"), plot = fp, width = 10, height = 12)

genes <- c("Rps19", "Itga9", "Hgf", "Itga1", "Calr", "Ppia")
genes <- c("Itga9", "Hgf", "Itga1", "Calr")
fp <- FeaturePlot(seurat_object, features = genes) & NoAxes()
ggsave(filename = file.path(path.guardar, "GoTerms_chemo.png"), plot = fp, width = 10, height = 12)



# SIGNATURAS ......................................................................
# 1. Lista de firmas con nombres únicos
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
# ¡Aquí está el truco!: usa `name = "GO_Score"`
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


