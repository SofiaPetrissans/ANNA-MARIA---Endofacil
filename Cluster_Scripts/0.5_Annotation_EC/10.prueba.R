
seurat_object <- readRDS("/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.4_SubsetEndothelial/mt15/8.RemoveCluster/2.Integration/SubsetEndothelial_Harmony.rds")
plot_dir <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.5_Annotation_EC/mt15/prueba"
output_dir <- plot_dir


Proliferative <-  c("Rhoc", "Tmsb10", "Dll4", "Sox4", "Col4a1", "Col4a2", "Kdr", "Flt1", "Acvrl1", "Itgb1")
fp <- FeaturePlot(seurat_object, features = Proliferative)
ggsave(filename = file.path(plot_dir, "Proliferative.png"),
       plot = fp, width = 18, height = 10)


library(openxlsx)
# INTENTAR DEFINIR UNA SIGNATURE DE VEHICLE A CHOP....................................................................................................
path <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.6_Downstream_Analysis_EC/mt15/subset1/1.DEGs_phenotype/DEGs_CHOP_vs_vehicle.xlsx"

DEGs_CHOP_vs_vehicle <- read.xlsx(path)

# Filtrar genes upregulados significativos
DEGs_up <- DEGs_CHOP_vs_vehicle %>%
  dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))


# Seleccionar los 30 genes top
signature_CHOP <- head(DEGs_up$gene, 30)

signature_CHOP

signature_CHOP <- c("Cdkn1a",    "Pvt1"  ,    "Bax"   ,    "Rps27l"   , "Exoc4"    , "Ltc4s"   ,  "S100a9" ,   "Hgf" ,      "Cyp26b1"  ,
"Eda2r" ,    "Phlda3"    ,"Ces2e"   ,  "Xist"    ,  "Rps19"  ,   "S100a8"   , "Lgals1"   , "Ly6e"    ,  "Gm42047"  ,
"Trp53inp1" ,"Mdm2"   ,   "Ephx1"   , "Ifitm3"  ,  "Map3k20" ,  "Rnf169"  ,  "Dglucy"  ,  "Gngt2"   ,  "Rpl15",   "Ccng1"     ,"Rpl12"    , "Id1")

seurat_object <- AddModuleScore(seurat_object,
                                features = list(signature_CHOP),
                                name = "CHOP_signature")
p <- VlnPlot(seurat_object, features = "CHOP_signature1", split.by = "Phenotype")
ggsave(filename = file.path(plot_dir, "Signature1_VlnPlot.png"),
       plot = p, width = 10, height = 10)

p1 <- FeaturePlot(seurat_object, features = "CHOP_signature1", split.by="Phenotype") &NoAxes()
ggsave(filename = file.path(plot_dir, "Signature1_FeaturePlot.png"),
       plot = p1, width = 16, height = 10)



# Genes extraidos de los go terms seleccionados
genes_vector <- c(
  "Cdkn1a", "Phlda3", "Eda2r", "Bax", "Rps27l", "Rpl23", "Ubb", "Rpl11", 
  "Mdm2", "Bbc3", "Aen", "Sesn2", "Rps20", "Rps7", "Tpt1", "Ackr3", "S100a9", 
  "S100a8", "Ephx1", "Trp53inp1", "Txn1", "Inhbb", "Rps19", "Cd81", "Ppia", 
  "Cd9", "Anxa1", "Rpl13a", "Kit", "Hgf", "Rack1", "Lgals1", "Ctla2a"
)
seurat_object <- AddModuleScore(seurat_object,
                                features = list(genes_vector),
                                name = "CHOP_signature2")
p <- VlnPlot(seurat_object, features = "CHOP_signature21", split.by = "Phenotype")
ggsave(filename = file.path(plot_dir, "Signature2_VlnPlot.png"),
       plot = p, width = 10, height = 10)

p1 <- FeaturePlot(seurat_object, features = "CHOP_signature21", split.by="Phenotype") &NoAxes()
ggsave(filename = file.path(plot_dir, "Signature2_FeaturePlot.png"),
       plot = p1, width = 16, height = 10)



# Genes comunes entre las dos firmas propuestas
common_genes <- intersect(signature_CHOP, genes_vector) # Hay 13 genes que se comparten con los top 30 DEGs
seurat_object <- AddModuleScore(seurat_object,
                                features = list(common_genes),
                                name = "CHOP_signature3")
p <- VlnPlot(seurat_object, features = "CHOP_signature31", split.by = "Phenotype")
ggsave(filename = file.path(plot_dir, "Signature3_VlnPlot.png"),
       plot = p, width = 10, height = 10)

p1 <- FeaturePlot(seurat_object, features = "CHOP_signature31", split.by="Phenotype") &NoAxes()
ggsave(filename = file.path(plot_dir, "Signature3_FeaturePlot.png"),
       plot = p1, width = 16, height = 10)



library(RColorBrewer)
getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
generateClusterColors <- function(object, group_col) {
  cluster_levels <- levels(factor(object@meta.data[[group_col]]))
  n_clusters <- length(cluster_levels)
  colors <- getPalette(n_clusters)
  names(colors) <- cluster_levels
  return(colors)
}

col <- getPalette(length(unique(seurat_object$Harmony_Log_res.0.7)))

umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "Harmony_Log_res.0.7", label = TRUE,
                       pt.size = 1, raster = FALSE, cols = col) & NoAxes()
  
ggsave(filename = file.path(output_dir, "Harmony_Log_res.0.7.png"), plot = umap_plot, width = 7, height = 7)

p <- DimPlot(seurat_object, group.by = "Harmony_Log_res.0.9", split.by="Phenotype", label = FALSE, pt.size = 1, raster = FALSE) +
  scale_color_manual(values = col) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right", legend.text = element_text(size = 12)) +  
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave(filename = file.path(output_dir, "DimPlot_Harmony_Log_res.0.9_Pheno.png"), plot = p, width = 12, height = 6)


col <- getPalette(length(unique(seurat_object$Harmony_Log_res.0.9)))
# BARRAS APILADAS
# Crear tabla resumen: número de células por condición y tipo celular
df_barplot <- seurat_object@meta.data %>%
  count(Phenotype, Harmony_Log_res.0.9) %>%
  group_by(Phenotype) %>%
  mutate(Proportion = n / sum(n) * 100)

# Barplot apilado
ggplot(df_barplot, aes(x = Phenotype, y = Proportion, fill = Harmony_Log_res.0.9)) +
  geom_bar(stat = "identity") +
  labs(y = "Percentage of Cells", x = "Condition", fill = "Cell Type") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# Guardar
ggsave(file.path(output_dir, "Barplot_Stacked_Percentages_By_Phenotype_res.0.9.png"), width = 4, height = 6)


