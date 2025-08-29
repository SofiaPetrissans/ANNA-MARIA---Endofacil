# ================================================================
# 13_GO_and_marker_signatures_AL5_visualization.R
# ================================================================
# OVERVIEW
# - Toma los resultados de GO:BP por subpoblación endotelial (Annotation_Layer5)
# - Selecciona y resume GO terms (únicos / compartidos, top N por clúster)
# - Genera barplots de GO seleccionados
# - Define firmas de genes por subpoblación (LVEC, LSEC_1–5), calcula module scores
#   y visualiza (UMAP/violin)
# - Construye DotPlots y Heatmaps (ComplexHeatmap) con “canonical markers”
#
# INPUTS (prioridad):
#   1) Objeto Seurat anotado localmente:
#       <path_base>/SubsetEndothelial_Harmony_Annotated.rds
#   2) Lista en memoria `go_results_list` (data.frames por subpoblación con columnas
#      al estilo clusterProfiler: ID, Description, p.adjust, Count, geneID, ...)
#      —SI NO EXISTE— intenta cargar desde un Excel:
#       <path_base>/10.Results_GO_AnnotationLayer/GO_Results_AnnotationLayer5.xlsx
#
# OUTPUTS:
#   - Carpeta <path_base>/13.Markers_GO_AL5 con:
#     * Excel con Top35 GO únicos por subpoblación (+summary)
#     * Barplots de GO seleccionados por subpoblación
#     * UMAP/Violin de module scores por firmas AL5
#     * DotPlot de marcadores “canónicos” por AL5
#     * Heatmaps (vertical/horizontal) de marcadores canónicos por AL5
#
# REQUISITOS:
#   R >= 4.2, Seurat >= 4, clusterProfiler, ComplexHeatmap, etc.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(openxlsx)
  library(readxl)
  library(janitor)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
})

# -----------------------------
# 0) PARÁMETROS / RUTAS
# -----------------------------
path_base <- "/Users/graupera_lab/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/ANNA_MARIA"
input_rds <- file.path(path_base, "SubsetEndothelial_Harmony_Annotated.rds")

# Fallback para GO si no está `go_results_list` en memoria:
go_xlsx_fallback <- "4.GO_Terms_AnnotationLayer5/GO_Results_AnnotationLayer5.xlsx"

# Salidas
output_dir <- file.path(path_base, "7.Markers_GO_HeatMaps_AL5")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Semilla reproducible para sampling/ordenes aleatorios (si hubiera)
set.seed(42)

# -----------------------------
# 1) UTILIDADES
# -----------------------------

# Asegura nombres de columnas esperadas (p.adjust, Count, Description, geneID)
normalize_go_cols <- function(df) {
  df <- janitor::clean_names(df)
  # Posibles alias
  desc_alias <- c("description", "term", "go_term", "go_bp", "go_description")
  padj_alias <- c("p_adjust", "p_adjusted", "p_adjustment", "p_adjusted_value", 
                  "p_adjust_bh", "p_adjust_fdr", "p_adjusted_fdr", "p_adjustment_fdr",
                  "p_adjusted_p_value", "qvalue", "p_adjusted_qvalue", "fdr", "p_adjust_bh",
                  "p_adjusted_bh", "p_adjust_fdr")
  padj_alias <- c(padj_alias, "p_adjusted", "p_adjust", "padj", "p_adjusted_p", "p_adjusted_value", "p_value_adjusted", "p_adjusted_pvalue","p_adjusted_p_value", "p_adjustment")
  padj_alias <- unique(c(padj_alias, "p_adjust", "p_adjusted", "p_adjusted_value", "p_adjusted_pvalue", "p_adjusted_p_value", "p_adjust_fdr", "p_adjust_bh", "qvalue", "fdr", "p_adjusted_fdr"))
  # compact
  desc_col <- intersect(desc_alias, names(df))[1]
  padj_col <- intersect(padj_alias, names(df))[1]
  
  # Count / geneID
  count_col <- if ("count" %in% names(df)) "count" else NA
  geneid_col <- if ("geneid" %in% names(df)) "geneid" else NA
  
  # Si falta alguna clave mínima, devuelve tal cual
  if (is.na(desc_col) || is.na(padj_col)) return(df)
  
  out <- df %>%
    dplyr::rename(
      Description = !!desc_col,
      p.adjust    = !!padj_col
    )
  
  if (!is.na(count_col))  out <- dplyr::rename(out, Count = !!count_col)
  if (!is.na(geneid_col)) out <- dplyr::rename(out, geneID = !!geneid_col)
  
  # Convertir p.adjust (si es char con coma)
  if (is.character(out$p.adjust)) {
    out <- out %>% mutate(p.adjust = as.numeric(gsub(",", ".", p.adjust)))
  }
  out
}

# Cargar Excel de GO -> lista por hoja (subpoblación)
load_go_xlsx_as_list <- function(xlsx_path) {
  if (!file.exists(xlsx_path)) stop("No se encontró el Excel de GO: ", xlsx_path)
  sheets <- readxl::excel_sheets(xlsx_path)
  go_list <- setNames(vector("list", length(sheets)), sheets)
  for (ct in sheets) {
    df <- readxl::read_xlsx(xlsx_path, sheet = ct)
    df <- normalize_go_cols(df)
    go_list[[ct]] <- df
  }
  go_list
}

# Barplots de GO (Top10 por p.adjust) con el nombre de la subpoblación en el fill
plot_go_terms_inside_bars <- function(go_list, name_prefix, output_dir, annotation_colors) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (ct in names(go_list)) {
    df <- go_list[[ct]]
    if (is.null(df) || nrow(df) == 0 || !all(c("Description","p.adjust") %in% colnames(df))) next
    
    top10 <- df %>% arrange(p.adjust) %>% slice_head(n = 10) %>%
      mutate(
        Term = Description,
        NegativeLog10P = -log10(p.adjust),
        ID = factor(Term, levels = rev(Term))
      )
    
    p <- ggplot(top10, aes(x = ID, y = NegativeLog10P, fill = ct)) +
      geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
      geom_text(aes(label = Term), y = 0, hjust = 0, color = "black",
                size = 3.8, fontface = "bold") +
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
    
    ggsave(file.path(output_dir, paste0(name_prefix, "_", ct, "_terms_in_bars.png")),
           plot = p, width = 8, height = 6, dpi = 300)
  }
}

# Helper heatmap de expresión media (genes x grupos)
scale_rows <- function(x) {
  m <- t(scale(t(x)))
  m[is.na(m)] <- 0
  m
}

# -----------------------------
# 2) COLORES Y CAPA
# -----------------------------
layer_name <- "Annotation_Layer5"
annotation.colors_5 <- c(
  "LVEC"   = "#86D0B9",
  "LSEC_1" = "#F4CDA5",
  "LSEC_2" = "#EBAFC3",
  "LSEC_3" = "#A8C8E5",
  "LSEC_4" = "#dce6b6",
  "LSEC_5" = "#C5A7E1"
)

# -----------------------------
# 3) CARGAS: Seurat y GO
# -----------------------------
stopifnot(file.exists(input_rds))
seurat_object <- readRDS(input_rds)

# Usa `go_results_list` si ya existe; si no, intenta desde Excel
if (!exists("go_results_list")) {
  message("`go_results_list` no existe en memoria. Intentando cargar desde Excel fallback...")
  go_results_list <- load_go_xlsx_as_list(go_xlsx_fallback)
}

# Normalizar columnas en todos los data.frames de la lista
go_results_list <- lapply(go_results_list, normalize_go_cols)

# -----------------------------
# 4) TOP GO TERMS (ÚNICOS / FRECUENCIA)
# -----------------------------
# a) Unir y anotar subpoblación (únicos por ID y Subpopulation)
go_all <- imap_dfr(go_results_list, ~ mutate(.x, Subpopulation = .y)) %>%
  distinct(Subpopulation, id, .keep_all = TRUE)

# b) Filtro base: significativos + tamaño mínimo
go_filt <- go_all %>%
  filter(!is.na(p.adjust), p.adjust < 0.05, !is.na(Count), Count >= 10)

# c) Contar en cuántos clústeres aparece cada GO tras el filtro
ncl <- go_filt %>%
  group_by(id) %>%
  summarise(n_clusters = n_distinct(Subpopulation), .groups = "drop")

# d) Términos específicos (n_clusters == 1)
go_unique <- go_filt %>% left_join(ncl, by = "id") %>% filter(n_clusters == 1)

# e) Top-35 por subpoblación (prioriza p.adjust, desempata por Count)
top35_unique <- go_unique %>%
  group_by(Subpopulation) %>%
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  slice_head(n = 35) %>%
  ungroup() %>%
  dplyr::select(Subpopulation, id, Description, Count, gene_id, p.adjust)

# f) Resumen “cuántos específicos disponibles” antes de cortar a 35
summary_unique <- go_unique %>% count(Subpopulation, name = "n_unique_available")

# g) Exportar
openxlsx::write.xlsx(
  list(
    "Top35_Unique_perCluster" = top35_unique,
    "Summary_unique_avail"    = summary_unique
  ),
  file = file.path(output_dir, "Top35_Unique_perCluster_sig_min10.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# 5) BARPLOTS de GO seleccionados (si tienes un Excel curado)
#    (opcional) Cambia `selected_go_xlsx` si ya hiciste una curación manual
# -----------------------------
selected_go_xlsx <- file.path(path_base, "11.Analysis_GO_AnnotationLayer", "Top10_go_terms_selected.xlsx")
if (file.exists(selected_go_xlsx)) {
  go_selected_list <- load_go_xlsx_as_list(selected_go_xlsx)
  plot_go_terms_inside_bars(
    go_list = go_selected_list,
    name_prefix = "Top10_Selected",
    output_dir = file.path(output_dir, "GO_Selected_Barplots"),
    annotation_colors = annotation.colors_5
  )
} else {
  message("No se encontró Excel de GO seleccionados: ", selected_go_xlsx)
}

# -----------------------------
# 6) FIRMAS DE GENES (Module Scores) PARA AL5
# -----------------------------
Idents(seurat_object) <- seurat_object[[layer_name, drop = TRUE]]

# Firmas (curadas a partir del análisis previo)
sig_LVEC <- c("Rspo3","Wnt2","Wnt9b","Tek","Flt1","Ephb4","Cdh5","Ramp2","Rras","Prkd1",
              "Selp","Jam2","Plvap","Cd47","Cd200","Gas6","Il1r1","Lbp",
              "Thbd","Entpd1","F2r","Tfpi2","Anxa2",
              "Col4a1","Col5a2","Tgfbi","Serpinh1","Flrt2","Cdh13","Plcb1","Gja4","Calm1")

sig_LSEC1 <- c("Cdkn1a","Bax","Mdm2","Phlda3","Trp53inp1","Ei24","Aen","Rps27l","Egr1",
               "Gpx1","Gpx4","Prdx1","Sod1","Txn1","Selenow","Sesn2","Cd36","Apoe",
               "Klf2","Klf4","Cldn5","Ecscr","Sparc","Adamts1","Hoxa5","Plk2",
               "Zfp36","Socs3","Ctla2a","Ada","Ifitm2","Ifitm3","Bst2","Irf1","Cxcl9","Cxcl10")

sig_LSEC2 <- c("Rpn1","Rpn2","Dad1","Ost4","Ostc","Magt1","Krtcap2","Dpm3","Dpagt1","Slc39a8","Tmem59",
               "Calr","Pdia3","Ctsl","Lgmn","B2m","H2-K1","H2-D1","H2-T23","Fcgr2b","Fcer1g","Clec4g",
               "Ap2m1","Calm1","Rac1","Arpc3","Tmsb4x","Cfl1","Lamp1","Map1lc3b","Gabarap",
               "Kdr","Flt4","Sparc","Ifitm2")

sig_LSEC3 <- c("Prickle1","Lrp6","Apc","Tcf7l1","Rora","Pard3b","Phldb2","Shroom2",
               "Nrp1","Nrp2","Robo1","Plxna2","Plxna4","Sema6a",
               "Dock1","Dock4","Dock7","Rapgef2","Vav3","Ralgapa1","Ralgapa2","Rap1gds1","Arfgef1",
               "Rabgap1","Asap1","Cdc42bpa","Ptk2","Macf1","Myo6","Exoc4","Exoc6b","Vti1a","Stxbp5","Eps15","Picalm","Vps13b")

sig_LSEC4 <- c("Btg2","Patl2","Tob1","Zfp36","Zfp36l1","Zc3h12a","Arid5a","Ythdf3","Ythdc1","Qki","Smg6","Srsf7","Nup98","Ddx19b",
               "Ern1","Hsp90aa1","Hspa1a","Hspa1b","Hsph1","Bag3","Ppp1r15a","Dnajb9","Hspb1","Atf3",
               "Map3k3","Gadd45b","Gadd45g","Ripk1","Lck","Camk2d","Ulk2","Clk1","Clk4","Fgfr2",
               "Nfkbiz","Stat2","Irf1")

sig_LSEC5 <- c("Lars2","Mrps21","Mrps14","Mrps33","Mrpl42","Mrpl18","Mrpl30","Mrpl51","Uqcc2","Coa3",
               "Vdac1","Vdac2","Slc25a3","Hspd1","Grpel1","Bnip3l","Higd1a","Egln1","Fam162a","Chchd2",
               "Eif4ebp1","Aqp1","Stub1","Prelid1","Tspo","Gpihbp1","Apoe","Lrp10","Acsl5","Glud1","Ech1","Park7","Vcp","Ifitm2","Ifitm3")

signatures <- list(
  LVEC   = intersect(sig_LVEC,   rownames(seurat_object)),
  LSEC_1 = intersect(sig_LSEC1,  rownames(seurat_object)),
  LSEC_2 = intersect(sig_LSEC2,  rownames(seurat_object)),
  LSEC_3 = intersect(sig_LSEC3,  rownames(seurat_object)),
  LSEC_4 = intersect(sig_LSEC4,  rownames(seurat_object)),
  LSEC_5 = intersect(sig_LSEC5,  rownames(seurat_object))
)

# Calcula ModuleScores (un AddModuleScore por firma)
for (nm in names(signatures)) {
  if (length(signatures[[nm]]) < 3) next
  seurat_object <- AddModuleScore(seurat_object, features = list(signatures[[nm]]), name = paste0(nm, "_Signature"))
}

# UMAP y Violin por firma
score_cols <- grep("_Signature1$", colnames(seurat_object@meta.data), value = TRUE)
if (length(score_cols)) {
  # UMAP split por Phenotype (si existe)
  if ("Phenotype" %in% colnames(seurat_object@meta.data)) {
    up <- lapply(score_cols, function(sc) {
      FeaturePlot(seurat_object, features = sc, split.by = "Phenotype", raster = FALSE) +
        ggtitle(sc) + theme(plot.title = element_text(hjust = 0.5)) & NoAxes()
    })
    gg <- patchwork::wrap_plots(up, ncol = 1)
    ggsave(file.path(output_dir, "UMAP_ModuleScores_AL5_splitPhenotype.png"), gg, width = 6, height = 16, dpi = 300)
  }
  # Violin por AL5
  vp <- lapply(score_cols, function(sc) {
    VlnPlot(seurat_object, features = sc, group.by = layer_name, pt.size = 0) +
      scale_fill_manual(values = annotation.colors_5) + ggtitle(sc)
  })
  gg2 <- patchwork::wrap_plots(vp, ncol = 1)
  ggsave(file.path(output_dir, "Violin_ModuleScores_AL5.png"), gg2, width = 8, height = 17, dpi = 300)
}

# -----------------------------
# 7) DOTPLOT Y HEATMAP DE “CANONICAL MARKERS” (AL5)
# -----------------------------
gene_sets <- list(
  "LVEC"   = c("Rspo3","Wnt2","Wnt9b","Tek","Flt1","Ephb4","Cdh5","Ramp2","Prkd1",
               "Selp","Jam2","Plvap","Cd47","Cd200","Gas6","Il1r1","Lbp","Thbd","Entpd1","F2r","Tfpi2",
               "Col4a1","Col5a2","Tgfbi","Serpinh1","Flrt2","Cdh13","Plcb1","Gja4"),
  "LSEC_1" = c("Cdkn1a","Bax","Mdm2","Phlda3","Trp53inp1","Ei24","Aen","Rps27l",
               "Cldn5","Ecscr","Ada","Bst2","Cxcl9","Cxcl10","Sesn2"),
  "LSEC_2" = c("Rpn1","Rpn2","Ostc","Magt1","Dpm3","Dpagt1","Slc39a8","Calr","Ctsl","Lgmn","B2m","Fcgr2b","Clec4g",
               "Ap2m1","Rac1","Lamp1","Kdr","Flt4"),
  "LSEC_3" = c("Prickle1","Lrp6","Apc","Tcf7l1","Rora","Phldb2","Shroom2",
               "Nrp1","Nrp2","Robo1","Plxna2","Plxna4","Sema6a",
               "Dock1","Dock4","Dock7","Rapgef2","Vav3","Ralgapa1","Rap1gds1","Arfgef1","Rabgap1","Asap1","Cdc42bpa","Ptk2","Macf1",
               "Exoc4","Exoc6b","Vti1a","Stxbp5","Eps15"),
  "LSEC_4" = c("Btg2","Tob1","Zfp36l1","Zc3h12a","Arid5a","Ythdc1","Srsf7","Nup98","Ddx19b",
               "Ern1","Hsp90aa1","Hspa1a","Hspa1b","Hsph1","Bag3","Ppp1r15a","Dnajb9","Atf3",
               "Map3k3","Gadd45b","Gadd45g","Ripk1","Lck","Ulk2","Clk1","Clk4","Fgfr2","Nfkbiz","Stat2"),
  "LSEC_5" = c("Mrps21","Mrps14","Mrps33","Mrpl42","Mrpl18","Mrpl30","Mrpl51",
               "Uqcc2","Coa3","Vdac1","Vdac2","Slc25a3","Grpel1","Bnip3l","Higd1a","Egln1","Aqp1",
               "Gpihbp1","Acsl5","Glud1","Ech1","Park7","Vcp")
)

# Asegurar que todas las features están en el objeto
gene_sets <- lapply(gene_sets, function(v) intersect(v, rownames(seurat_object)))
features <- unique(unlist(gene_sets))
seurat_object <- SetIdent(seurat_object, value = layer_name)

# DotPlot por categorías
gene_category_df <- bind_rows(lapply(names(gene_sets), function(cat) {
  data.frame(Feature = gene_sets[[cat]], Category = cat)
}))

dp_data <- DotPlot(seurat_object, features = gene_category_df$Feature)$data %>%
  left_join(gene_category_df, by = c("features.plot" = "Feature"))

dp <- ggplot(dp_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "grey90", high = "#3C2692") +
  facet_wrap(~Category, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10, face = "bold")) +
  labs(x = "Markers", y = "AL5", size = "% Expr.", color = "Avg. Expr (scaled)")
ggsave(file.path(output_dir, paste0("DotPlot_CanonicalMarkers_", layer_name, ".png")),
       dp, width = 18, height = 6, dpi = 300)

# HEATMAP (genes x subpoblación AL5)
avg.exp <- AverageExpression(seurat_object, group.by = layer_name, return.seurat = FALSE)$RNA
marker_genes <- unlist(gene_sets, use.names = FALSE)
genes_disponibles <- intersect(marker_genes, rownames(avg.exp))
stopifnot(length(genes_disponibles) > 0)

mat <- avg.exp[genes_disponibles, , drop = FALSE]
mat.z <- scale_rows(mat)

gene.groups <- data.frame(
  Grupo = rep(names(gene_sets), times = purrr::map_int(gene_sets, length))
)
rownames(gene.groups) <- make.unique(marker_genes)
gene.groups <- gene.groups[rownames(mat.z), , drop = FALSE]

left_annotation <- rowAnnotation(
  Grupo = gene.groups$Grupo,
  col = list(Grupo = annotation.colors_5),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = "Subpoblación",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7)
  )
)

col_fun <- colorRamp2(c(min(mat.z), 0, max(mat.z)), c("#1874CD", "#EEEEE0", "#CD2626"))

hp <- Heatmap(
  mat.z,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  left_annotation = left_annotation,
  border = TRUE,
  column_title = paste("Canonical markers by", layer_name),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_split = gene.groups$Grupo,
  row_title_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title = "Z-score", title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 6))
)
png(file.path(output_dir, paste0("Markers_Heatmap_", layer_name, ".png")),
    width = 1100, height = 2800, res = 300)
draw(hp, heatmap_legend_side = "right", padding = unit(c(10, 10, 2, 2), "mm"))
dev.off()

# HEATMAP horizontal (clusters x genes)
hp_h <- Heatmap(
  t(mat.z),
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 5),3
  top_annotation = HeatmapAnnotation(
    Grupo = gene.groups$Grupo,
    col = list(Grupo = annotation.colors_5),
    show_annotation_name = FALSE
  ),
  column_split = gene.groups$Grupo,
  border = TRUE,
  column_title = paste("Canonical markers by", layer_name),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(title = "Z-score", title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 6))
)
png(file.path(output_dir, paste0("Markers_Heatmap_Horizontal_", layer_name, ".png")),
    width = 3000, height = 650, res = 300)
draw(hp_h, heatmap_legend_side = "right", padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# -----------------------------
# FIN DEL SCRIPT
# -----------------------------
message("Listo. Resultados guardados en: ", output_dir)
