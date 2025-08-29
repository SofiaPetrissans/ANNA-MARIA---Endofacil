# --------------------------------------------------------------------------------
# SCRIPT: Anotación superficial + estimación de marcadores para múltiples thresholds
# AUTORA: Sofia Petrissans Moll
# FECHA: 30/06/2025
# --------------------------------------------------------------------------------

# ...............................................................................
# LIBRERIAS
# ...............................................................................
source("/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggplot2)



# ...............................................................................
# FUNCIONES
# ...............................................................................
SimpleAnnotator <- function(seu, resolution, plot_dir) {
  seu <- SetIdent(seu, value = resolution)
  
  gene_sets <- list(
    Hepatocytes = c("Alb", "Ttr", "Cyp2e1", "Cyp3a11", "Ass1", "Glul", "Cyp7a1", "Hnf4a"),
    Endothelial = c("Pecam1", "Vwf", "Kdr"),
    LSEC = c("Lyve1", "Stab2", "Clec4g", "Mrc1", "Fcgr2b", "Icam1", "Cd36"),
    LVEC = c("Cdh5", "Cd31",  "Flt1", "Tek", "Nos3"),
    Kupffer = c("Adgre1", "Clec4f", "Cd68", "Marco"),
    Monocytes = c("Cd14", "Ccr2", "Ly6c1"),
    TCells = c("Cd3e", "Cd8a", "Cd4", "Trac"),
    BCells = c("Cd79a", "Ms4a1", "Cd19"),
    NK = c("Nkg7", "Klrb1c", "Gzmb", "Klrf1"),
    Neutrophil = c("S100a8", "S100a9", "Lcn2", "Cxcr2"),
    HSC = c("Rbp1", "Des", "Lrat", "Pdgfrb", "Reln"),
    Cholangiocytes = c("Krt", "Epcam", "Sox9", "Cftr"),
    Fibroblasts = c("Col1a1", "Col1a2", "Dcn", "Pdgfra"),
    Mizonal = c("Igfbp2", "Hamp", "Hamp2", "Cyp8b1", "Mt2a", "Mt1g", "Ndufb1", "Hint1", "Cox7c", "Apoc1", "Fabp1", "Aldob", "Adh1b"),
    Pericentral = c("Igfbp1", "Nt5e", "Cyp3a4", "Adh4", "Clul", "Bche", "Rspo3", "Wnt9b", "Wnt2", "Thdb", "Cdh13", 
    "Clec1b", "Clec4m", "Foxn2", "Egr", "Nfatc1", "Maf", "Nr2f1", "Gata4"),
    Periportal = c("Dll4", "Efnb2", "Ltpp4", "Ednrb", "Jag1", "Lrg1", "Efnb1", "Ltbp4", "Adgrg6"),
    Periportal_2 = c("Cyp2f2", "Hal", "Hsd17b13", "Sds", "Ctsc", "A1bg", "Aldh1b1", "Pck1"),
    Pericentral_2 = c("Cyp4a14", "Cyp2d9", "Gstm3", "Cyp4a10", "Mup17", "Slc1a2", "Slc22a1", "Cyp1a2", "Aldh1a1", "Cyp2a5", "Gulo", "Cyp2c37", "Lect2", "Oat")
  )


  for (celltype in names(gene_sets)) {
    features <- gene_sets[[celltype]]
    features <- features[features %in% rownames(seu)]  # filtra genes existentes
    if (length(features) == 0) next
    p1 <- FeaturePlot(seu, features = features, order = TRUE, raster = FALSE) & NoAxes()
    p2 <- VlnPlot(seu, features = features, pt.size = 0.1, raster = FALSE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    g <- ggarrange(p1, p2, nrow = 1)
    ggsave(filename = file.path(plot_dir, paste0(celltype, ".png")), plot = g, width = 16, height = 7)
  }
}

ClusterMarkers <- function(seu, resolution, marker_dir, plot_dir) {
  seu <- SetIdent(seu, value = resolution)
  clusters <- levels(seu)
  all_markers <- list()

  for (clus in clusters) {
    subcells <- WhichCells(seu, idents = clus)
    if (length(subcells) < 10) next
    markers <- FindMarkers(seu, ident.1 = clus)
    markers$gene <- rownames(markers)
    markers$cluster <- clus
    all_markers[[paste0("Clus_", clus)]] <- markers
  }

  saveRDS(all_markers, file.path(marker_dir, paste0("Markers_res", resolution, ".rds")))
  write.xlsx(all_markers, file.path(marker_dir, paste0("Markers_res", resolution, ".xlsx")))

  total <- bind_rows(all_markers)
  write.xlsx(total, file.path(marker_dir, paste0("Markers_res", resolution, "_Total.xlsx")))

  top_markers <- total %>%
    group_by(cluster) %>%
    top_n(4, avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

  features_top <- unique(top_markers$gene)

  p_dot <- DotPlot(seu, features = features_top) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = rel(1.05)))
  ggsave(filename = file.path(plot_dir, paste0("Top4Markers_DotPlot_res", resolution, ".png")), plot = p_dot, width = 16, height = 7)
}

# ...............................................................................
# EJECUCIÓN
# ...............................................................................

# ...............................................................................
# ARGUMENTOS
# ...............................................................................

RESOLUTION <- c("0.1", "0.3", "0.5", "0.7", "0.9", "1")
MT_THRESHOLD <- "mt15"
TYPE <- "Integation"

# ...............................................................................
# PATHS
# ...............................................................................
base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path,  MT_THRESHOLD, "9.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")
# Cargar el objeto
seurat_object <- readRDS(input_rds)

path.guardar_original <- paths_AnnaMaria$Path_guardar
path.guardar_original <- file.path(path.guardar_original, MT_THRESHOLD)

for (res in RESOLUTION) {
  resolution_col <- if (TYPE == "Raw") {
    paste0("RNA_snn_res.", res)
  } else {
    paste0("Harmony_Log_res.", res)
  }

  if (!resolution_col %in% colnames(seurat_object@meta.data)) {
    warning(paste("Columna de resolución", resolution_col, "no encontrada. Saltando..."))
    next
  }

  output_dir <- file.path(path.guardar_original, "1.B.MarkersAnnotation", paste0("res", res))
  plot_dir <- file.path(output_dir, "Plots")
  marker_dir <- file.path(output_dir, "Markers")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)

  # Anotación y marcadores
  SimpleAnnotator(seurat_object, resolution_col, plot_dir)
  ClusterMarkers(seurat_object, resolution_col, marker_dir, plot_dir)

  message("Anotación y marcadores completados para resolución ", res)
}


SimpleAnnotator(seurat_object, resolution_col, plot_dir)
ClusterMarkers(seurat_object, resolution_col, marker_dir, plot_dir)

message("Anotación y marcadores completados para ", MT_THRESHOLD, " ", TYPE, " resolución ", RESOLUTION)


# DOT PLOTS ADICIONALES (al estilo del script inicial)
gene_sets <- list(
  Hepatocytes = c("Alb", "Ttr", "Cyp2e1", "Cyp3a11", "Ass1", "Glul", "Cyp7a1", "Hnf4a"),
  Endothelial = c("Pecam1", "Vwf", "Kdr"),
  LSEC = c("Lyve1", "Stab2", "Clec4g", "Mrc1", "Fcgr2b", "Icam1", "Cd36"),
  LVEC = c("Cdh5", "Cd31",  "Flt1", "Tek", "Nos3"),
  Kupffer = c("Adgre1", "Clec4f", "Cd68", "Marco"),
  Monocytes = c("Cd14", "Ccr2", "Ly6c1"),
  TCells = c("Cd3e", "Cd8a", "Cd4", "Trac"),
  BCells = c("Cd79a", "Ms4a1", "Cd19"),
  NK = c("Nkg7", "Klrb1c", "Gzmb", "Klrf1"),
  Neutrophil = c("S100a8", "S100a9", "Lcn2", "Cxcr2"),
  HSC = c("Rbp1", "Des", "Lrat", "Pdgfrb", "Reln"),
  Cholangiocytes = c("Krt", "Epcam", "Sox9", "Cftr"),
  Fibroblasts = c("Col1a1", "Col1a2", "Dcn", "Pdgfra"),
  Mizonal = c("Igfbp2", "Hamp", "Hamp2", "Cyp8b1", "Mt2a", "Mt1g", "Ndufb1", "Hint1", "Cox7c", "Apoc1", "Fabp1", "Aldob", "Adh1b"),
  Pericentral = c("Igfbp1", "Nt5e", "Cyp3a4", "Adh4", "Clul", "Bche", "Rspo3", "Wnt9b", "Wnt2", "Thdb", "Cdh13", 
                  "Clec1b", "Clec4m", "Foxn2", "Egr", "Nfatc1", "Maf", "Nr2f1", "Gata4"),
  Periportal = c("Dll4", "Efnb2", "Ltpp4", "Ednrb", "Jag1", "Lrg1", "Efnb1", "Ltbp4", "Adgrg6"),
  Periportal_2 = c("Cyp2f2", "Hal", "Hsd17b13", "Sds", "Ctsc", "A1bg", "Aldh1b1", "Pck1"),
  Pericentral_2 = c("Cyp4a14", "Cyp2d9", "Gstm3", "Cyp4a10", "Mup17", "Slc1a2", "Slc22a1", "Cyp1a2", "Aldh1a1", 
                    "Cyp2a5", "Gulo", "Cyp2c37", "Lect2", "Oat")
)

features <- unique(unlist(gene_sets))

seurat_object <- SetIdent(seurat_object, value = resolution_col)

# DotPlot por categorías
gene_category_df <- bind_rows(lapply(names(gene_sets), function(cat) {
  data.frame(Feature = gene_sets[[cat]], Category = cat)
}))
dp_data <- DotPlot(seurat_object, features = gene_category_df$Feature)$data
dp_data <- dp_data %>%
  left_join(gene_category_df, by = c("features.plot" = "Feature"))
dp <- ggplot(dp_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "grey90", high = "#3C2692") +
  facet_wrap(~Category, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank()) +
  labs(x = "Features", y = "Identity", size = "Percent Expressed", color = "Average Expression")
ggsave(filename = file.path(plot_dir, "CanonicalMarkers_DotPlot_sep.png"),
       plot = dp, width = 18, height = 5)
