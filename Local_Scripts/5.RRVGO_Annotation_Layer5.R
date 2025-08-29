############################################
# Objetivo. Identificar funciones biológicas (GO:BP) que definen cada subpoblación endotelial en Annotation_Layer5 (LVEC y LSEC_1…LSEC_5), partiendo de DEGs 1-vs-rest (solo up-regulated), y resumirlas mediante:
#
#   -  Reducción semántica (rrvgo) y macro-categorías (agrupación léxica manual).
#   -  Heatmap presencia/ausencia de términos específicos por clúster.
# 
# Entrada.
#   -  `SubsetEndothelial_Harmony_Annotated.rds` (con Annotation_Layer5, Phenotype, etc.)
# 
# Salidas (se crean si no existen):
#   -  11.Analysis_GO_AnnotationLayer/ (figuras y tablas de especificidad, reducción semántica y macro-categorías).
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

out_sem  <- file.path("5.RRVGO_AnnotationLayer")
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
# 4) PIPELINE
# ================================
# 0. Configuración inicial
# ------------------------------------------------
Idents(seurat_object) <- seurat_object$Annotation_Layer5
populations <- levels(seurat_object$Annotation_Layer5)

# 1. Identificación de DEGs para cada población
# ------------------------------------------------
deg_list <- list()

for (pop in populations) {
  # Comparar población vs resto
  deg <- FindMarkers(
    seurat_object, 
    ident.1 = pop, 
    only.pos = TRUE,      # Solo genes upregulados
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  deg$gene <- rownames(deg)
  deg_list[[pop]] <- deg
}

# 2. GO:BP para cada población
# -------------------------------------------------
go_results_list <- list()

for (pop in names(deg_list)) {
  # Filtrar genes significativos
  genes_sig <- deg_list[[pop]] %>%
    filter(p_val_adj < 0.05) %>%
    pull(gene)
  
  if (length(genes_sig) == 0) next
  
  # Mapear a ENTREZ IDs
  entrez <- bitr(genes_sig, fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # Enriquecimiento GO:BP
  ego <- enrichGO(
    gene          = entrez$ENTREZID,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # Guardar resultado
  go_results_list[[pop]] <- as.data.frame(ego)
}


# ================================
# 5) Reducción semántica (rrvgo) y conteos por categoría (parentTerm)
# ================================

# Crear tabla con subpoblación, min y max
count_summary <- lapply(go_results_list, function(df) {
  if (nrow(df) == 0) return(NULL)
  data.frame(
    MinCount = min(df$Count, na.rm = TRUE),
    MeanCount  = mean(df$Count, na.rm = TRUE),
    MaxCount = max(df$Count, na.rm = TRUE)
  )
}) %>%
  bind_rows(.id = "Cluster")

# Mostrar
count_summary


# Identificar GO Terms específicos de cada subpoblación
# --------------------------------
# 0) Unir y dejar únicos por ID+Cluster

go_all <- bind_rows(lapply(names(go_results_list), function(pop) {
  df <- go_results_list[[pop]]
  df$Cluster <- pop
  df
})) %>%
  distinct(Cluster, ID, .keep_all = TRUE)

# -----------------------------
# 1) Contar en cuántos clusters aparece cada GO (solo si es significativo y Count>=10)

go_counts <- go_all %>%
  filter(p.adjust < 0.05, Count >= 10) %>%
  group_by(ID) %>%
  summarise(n_clusters = n_distinct(Cluster), .groups = "drop")

# -----------------------------
# 2) Anotar estado: Count<10, FALSE (compartido), TRUE (específico)

go_annot <- go_all %>%
  left_join(go_counts, by = "ID") %>%
  mutate(
    Status = dplyr::case_when(
      Count < 10                        ~ "Count<10",
      !is.na(n_clusters) & n_clusters==1 ~ "TRUE",   # específicos
      !is.na(n_clusters) & n_clusters>=2 ~ "FALSE",  # no específicos
      TRUE                               ~ "Count<10"  # p.adj≥0.05 o sin n_clusters -> fuera del conteo principal
    ),
    Status = factor(Status, levels = c("Count<10","FALSE","TRUE"))
  )

# -----------------------------
# 3) Barplot apilado con el orden: Count<10 (abajo) -> FALSE -> TRUE

p <- go_annot %>%
  group_by(Cluster, Status) %>%
  summarise(n_terms = n(), .groups = "drop") %>%
  ggplot(aes(x = Cluster, y = n_terms, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c("Count<10" = "#EEEEEE",  # gris claro
               "FALSE"    = "#BDBDBD",  # gris medio
               "TRUE"     = "#8AC8CC")  # tu color de específicos
  ) +
  labs(title = "Número de GO terms por subpoblación",
       y = "Nº de GO terms", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
ggsave(filename = file.path(out_sem, "BarPlot_num_GO_Count10_FALSE_TRUE.png"),
       plot = p, width = 9, height = 7, dpi = 300)



# Agrupar automáticamente GO terms en categorías funcionales mediante similitud semántica,
# Obtener términos GO significativos con p.adjust < 0.05 y Count >= 10
go_filtered <- bind_rows(go_results_list) %>%
  filter(p.adjust < 0.05, Count >= 10)

# Extraer IDs y scores
go_terms <- unique(go_filtered$ID)
scores <- setNames(-log10(go_filtered$p.adjust), go_filtered$ID)


# Preparar ontología y base de datos para Mus musculus
simMatrix <- calculateSimMatrix(go_terms,
                                orgdb = "org.Mm.eg.db",
                                ont = "BP",
                                method = "Rel")  # Puedes usar también "Resnik" o "Lin"


# Reducir y agrupar
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores = scores,
                                threshold = 0.7,  # ajusta para mayor o menor agrupación
                                orgdb = "org.Mm.eg.db")

treemapPlot(reducedTerms)
scatterPlot(simMatrix, reducedTerms)

# Primero, combinamos los resultados con la información de clúster
go_all <- bind_rows(lapply(names(go_results_list), function(pop) {
  df <- go_results_list[[pop]]
  df$Cluster <- pop
  df
}))

# Filtrar por p.adjust y Count
go_filtered <- go_all %>%
  filter(p.adjust < 0.05, Count >= 10)

# Marcar términos específicos (solo aparecen en 1 clúster)
go_counts <- go_filtered %>%
  group_by(ID) %>%
  summarise(n_clusters = n_distinct(Cluster), .groups = "drop")

go_filtered <- go_filtered %>%
  left_join(go_counts, by = "ID") %>%
  mutate(Specific = ifelse(n_clusters == 1, TRUE, FALSE))


go_category_map <- reducedTerms %>%
  transmute(ID = go, FunctionalCategory = parentTerm) %>%
  distinct()

# 3) Ahora sí, unir solo los términos específicos y añadir su categoría
top_go_specific <- go_filtered %>%
  filter(Specific) %>%                       # antes pedías: ID %in% go_category_map$ID & Specific==TRUE
  left_join(go_category_map, by = "ID")

# Unirlo con los datos de GO terms específicos
top_go_specific <- go_filtered %>%
  filter(ID %in% go_category_map$ID, Specific == TRUE) %>%
  left_join(go_category_map, by = "ID")

# Crear matriz binaria actualizada
binary_mat <- top_go_specific %>%
  dplyr::select(Cluster, Description) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Cluster, values_from = value, values_fill = 0) %>%
  column_to_rownames("Description")

# Crear anotación de fila por categoría funcional
row_anno <- rowAnnotation(
  Category = top_go_specific$FunctionalCategory[match(rownames(binary_mat), top_go_specific$Description)],
  col = list(Category = setNames(rainbow(length(unique(top_go_specific$FunctionalCategory))), 
                                 unique(top_go_specific$FunctionalCategory)))
)

# Dibujar heatmap con anotación funcional
hp <- Heatmap(
  as.matrix(binary_mat),
  name = "Específico",
  col = c("0" = "white", "1" = "#D55E00"),
  right_annotation = row_anno,
  show_column_names = TRUE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 10)
)

# Guardar
png(file.path(out_sem, paste0("Hp.png")), width = 2600, height = 3400, res = 300)
draw(hp, heatmap_legend_side = "right")
dev.off()


# Unir cluster funcional (de rrvgo) con cluster celular (de go_filtered)
reduced_with_clusters <- reducedTerms %>%
  left_join(go_filtered[, c("ID", "Cluster")], by = c("go" = "ID"))

# Gurdar la tabla entera
write.xlsx(reduced_with_clusters, file.path(out_sem, "GOTerms_agrupados_funciones.xlsx"), rowNames = TRUE)



# Contar cuántos GO terms hay por categoría funcional y por subpoblación
summary_table <- reduced_with_clusters %>%
  group_by(parentTerm, Cluster) %>%
  summarise(n_terms = n(), .groups = "drop")

# Barplot apilado
ggplot(summary_table, aes(x = Cluster, y = n_terms, fill = parentTerm)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "GO terms agrupados por subpoblación y categoría funcional",
       x = "Subpoblación endotelial",
       y = "Número de GO terms",
       fill = "Categoría funcional (rrvgo)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
ggsave(filename = file.path(out_sem, "0.BarPlot.png"), width =18, height = 7)


macro_category_map <- c(
  # Metabolism
  "amide metabolic process" = "Metabolism",
  "macromolecule glycosylation" = "Metabolism",
  "oxidative phosphorylation" = "Metabolism",
  "purine ribonucleotide biosynthetic process" = "Metabolism",
  "purine ribonucleotide metabolic process" = "Metabolism",
  "proton motive force-driven mitochondrial ATP synthesis" = "Metabolism",
  "negative regulation of phosphorus metabolic process" = "Metabolism",
  
  # RNA Processing / Translation
  "cytoplasmic translation" = "RNA processing / Translation",
  "rRNA processing" = "RNA processing / Translation",
  "ribosomal small subunit biogenesis" = "RNA processing / Translation",
  "protein-RNA complex assembly" = "RNA processing / Translation",
  
  # Cell migration / Motility
  "chemotaxis" = "Cell migration / Motility",
  "leukocyte migration" = "Cell migration / Motility",
  "locomotory behavior" = "Cell migration / Motility",
  
  # Development / Morphogenesis
  "developmental cell growth" = "Development / Morphogenesis",
  "developmental maturation" = "Development / Morphogenesis",
  "dendrite morphogenesis" = "Development / Morphogenesis",
  "morphogenesis of a branching epithelium" = "Development / Morphogenesis",
  "establishment or maintenance of cell polarity" = "Development / Morphogenesis",
  "muscle cell differentiation" = "Development / Morphogenesis",
  "male gonad development" = "Development / Morphogenesis",
  "heart morphogenesis" = "Development / Morphogenesis",
  "regulation of vasculature development" = "Development / Morphogenesis",
  "sprouting angiogenesis" = "Development / Morphogenesis",
  
  # Hemostasis / Coagulation
  "blood coagulation" = "Hemostasis / Coagulation",
  "coagulation" = "Hemostasis / Coagulation",
  
  # Cell Cycle / Proliferation / Growth
  "cell cycle G1/S phase transition" = "Cell Cycle / Proliferation / Growth",
  "homeostasis of number of cells" = "Cell Cycle / Proliferation / Growth",
  "negative regulation of epithelial cell proliferation" = "Cell Cycle / Proliferation / Growth",
  
  # Apoptosis / Cell Death
  "positive regulation of intrinsic apoptotic signaling pathway" = "Apoptosis / Cell Death",
  "regulation of cell killing" = "Apoptosis / Cell Death",
  "signal transduction by p53 class mediator" = "Apoptosis / Cell Death",
  
  # Cell Adhesion / Junctions
  "cell-substrate adhesion" = "Cell adhesion / Junctions",
  
  # Signaling / MAPK / TGF-beta
  "canonical Wnt signaling pathway" = "Signaling / MAPK / TGF-beta",
  "signal transduction by p53 class mediator" = "Signaling / MAPK / TGF-beta",
  
  # Stress / Hypoxia / ROS
  "cellular response to reactive oxygen species" = "Stress / Hypoxia / ROS",
  "cellular oxidant detoxification" = "Stress / Hypoxia / ROS",
  "response to radiation" = "Stress / Hypoxia / ROS",
  "cellular response to radiation" = "Stress / Hypoxia / ROS",
  "regulation of reactive oxygen species metabolic process" = "Stress / Hypoxia / ROS",
  "response to steroid hormone" = "Stress / Hypoxia / ROS",
  "response to starvation" = "Stress / Hypoxia / ROS",
  
  # Immune Response / Antigen Presentation
  "antigen processing and presentation of peptide antigen" = "Immune response / Antigen presentation",
  "biological process involved in interaction with host" = "Immune response / Antigen presentation",
  "leukocyte migration" = "Immune response / Antigen presentation",
  "negative regulation of inflammatory response" = "Immune response / Antigen presentation",
  "negative regulation of cytokine production" = "Immune response / Antigen presentation",
  "response to molecule of bacterial origin" = "Immune response / Antigen presentation",
  "regulation of T cell activation" = "Immune response / Antigen presentation",
  "viral life cycle" = "Immune response / Antigen presentation",
  
  # Homeostasis / Transport
  "calcium ion homeostasis" = "Homeostasis / Transport",
  "mitochondrial transport" = "Homeostasis / Transport",
  "protein localization to plasma membrane" = "Homeostasis / Transport",
  
  # Differentiation / Specialization
  "central nervous system neuron differentiation" = "Differentiation / Specialization",
  "regulation of muscle system process" = "Differentiation / Specialization",
  "regulation of osteoblast differentiation" = "Differentiation / Specialization",
  "regulation of epithelial cell differentiation" = "Differentiation / Specialization",
  
  # Protein / DNA Organization
  "regulation of protein catabolic process" = "Protein / DNA Organization",
  
  # Physiological Rhythms
  "rhythmic process" = "Physiological rhythms",
  
  # Mitochondrial Function
  "mitochondrial translation" = "Mitochondrial Function",
  "mitochondrial respirasome assembly" = "Mitochondrial Function",
  
  # Actin / Cytoskeleton
  "regulation of actin filament-based process" = "Actin / Cytoskeleton",
  
  # GTPase signaling
  "regulation of GTPase activity" = "GTPase signaling",
  
  # Multi-organism process
  "multi-multicellular organism process" = "Multi-organism process",
  
  # Cellular maintenance
  "cellular component maintenance" = "Cellular maintenance"
)

# Asignar la macro categoría
summary_table$MacroCategory <- recode(summary_table$parentTerm, !!!macro_category_map)
summary_table %>% filter(is.na(MacroCategory))  # para revisar los no asignados

ggplot(summary_table, aes(x = Cluster, y = n_terms, fill = MacroCategory)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "GO terms por macro-categoría funcional",
       x = "Subpoblación",
       y = "Número de términos GO",
       fill = "Macro-categoría funcional") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


macro_colors <- c(
  # Colores ya definidos
  "Metabolism"                             = "#B2DF8A",  # verde claro
  "RNA processing / Translation"           = "#984EA3",  # púrpura fuerte
  "Cell migration / Motility"              = "#1F78B4",  # azul intermedio
  "Development / Morphogenesis"            = "#A6611A",  # marrón
  "Hemostasis / Coagulation"               = "#FFD92F",  # amarillo brillante
  "Cell Cycle / Proliferation / Growth"    = "#41B7C4",  # azul cielo
  "Apoptosis / Cell Death"                 = "#F781BF",  # rosa
  "Cell adhesion / Junctions"              = "#A6CEE3",  # azul muy claro
  "Signaling / MAPK / TGF-beta"            = "#FF7F00",  # naranja
  "Stress / Hypoxia / ROS"                 = "#542788",  # violeta intenso
  "Immune response / Antigen presentation" = "#FB8072",  # rojo
  "Homeostasis / Transport"                = "#CCEBC5",  # verde menta
  "Differentiation / Specialization"       = "#BEBADA",  # gris medio
  "Protein / DNA Organization"             = "#6A3D9A",  # morado
  "Physiological rhythms"                  = "#F4DBB5",  # lavanda
  
  # Nuevas categorías añadidas
  "Mitochondrial Function"                 = "#BE9F77",  # verde azulado
  "Actin / Cytoskeleton"                   = "#D66982",  # mostaza
  "GTPase signaling"                       = "#8DD3C7",  # verde turquesa claro
  "Multi-organism process"                = "#FDB462",  # salmón claro
  "Cellular maintenance"                   = "#BC80BD"   # lila
)



ggplot(summary_table, aes(x = Cluster, y = n_terms, fill = MacroCategory)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = macro_colors) +
  labs(
    title = "GO terms por macro-categoría funcional",
    x = "Subpoblación",
    y = "Número de términos GO",
    fill = "Macro-categoría funcional"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = file.path(out_sem, "1.BarPlot_Categ.png"), width =12, height = 7)



# Calcular PROPORCIONES por clúster
summary_table_prop <- summary_table %>%
  group_by(Cluster) %>%
  mutate(Proportion = n_terms / sum(n_terms)) %>%
  ungroup()

# Plot de proporciones
ggplot(summary_table_prop, aes(x = Cluster, y = Proportion, fill = MacroCategory)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = macro_colors) +
  labs(
    title = "Proporción de GO terms por macro-categoría funcional",
    x = "Subpoblación",
    y = "Proporción",
    fill = "Macro-categoría funcional"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Guardar
ggsave(filename = file.path(out_sem, "1.BarPlot_Categ_Proporcion.png"), width = 12, height = 7)




# -------------------------------------------------------
# 1. Combinar resultados GO de todas las subpoblaciones
# -------------------------------------------------------
go_all <- bind_rows(lapply(names(go_results_list), function(pop) {
  df <- go_results_list[[pop]]
  df$Cluster <- pop
  df
}))

# -------------------------------------------------------
# 2. Filtrar por significancia y número de genes
# -------------------------------------------------------
go_filtered <- go_all %>%
  filter(p.adjust < 0.05, Count >= 10)

# -------------------------------------------------------
# 3. Marcar términos específicos (solo en una subpoblación)
# -------------------------------------------------------
go_counts <- go_filtered %>%
  group_by(ID) %>%
  summarise(n_clusters = n_distinct(Cluster), .groups = "drop")

go_filtered <- go_filtered %>%
  left_join(go_counts, by = "ID") %>%
  mutate(Specific = ifelse(n_clusters == 1, TRUE, FALSE))

# -------------------------------------------------------
# 4. Filtrar solo términos específicos
# -------------------------------------------------------
go_specific <- go_filtered %>% filter(Specific == TRUE)

# -------------------------------------------------------
# 5. Agrupar semánticamente con rrvgo
# -------------------------------------------------------

go_terms <- unique(go_specific$ID)
scores <- setNames(-log10(go_specific$p.adjust), go_specific$ID)

simMatrix <- calculateSimMatrix(go_terms,
                                orgdb = "org.Mm.eg.db",
                                ont = "BP",
                                method = "Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores = scores,
                                threshold = 0.7,
                                orgdb = "org.Mm.eg.db")

# -------------------------------------------------------
# 6. Añadir clúster (subpoblación) a cada término reducido
# -------------------------------------------------------
reduced_with_clusters <- reducedTerms %>%
  left_join(go_specific[, c("ID", "Cluster")], by = c("go" = "ID"))

# -------------------------------------------------------
# 7. Crear tabla resumen y asignar macro-categorías
# -------------------------------------------------------
summary_table <- reduced_with_clusters %>%
  group_by(parentTerm, Cluster) %>%
  summarise(n_terms = n(), .groups = "drop")

# -------------------------------------------------------
# 8. Barplot absoluto: términos específicos por categoría
# -------------------------------------------------------
library(Polychrome)
library(ggplot2)

# Niveles fijos para que el mapping sea estable
levels_parent <- sort(unique(summary_table$parentTerm))
summary_table$parentTerm <- factor(summary_table$parentTerm, levels = levels_parent)

set.seed(42)
pal <- Polychrome::createPalette(
  length(levels_parent),
  seedcolors = c("#FFFFFF", "#000000"),  # ayuda a maximizar contraste
  M = 50000                               # muestreo alto → mejor separación
)
pal <- setNames(pal, levels_parent)

ggplot(summary_table, aes(x = Cluster, y = n_terms, fill = parentTerm)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal, drop = FALSE) +
  theme_minimal() +
  labs(
    title = "GO terms específicos por macro-categoría funcional",
    x = "Subpoblación",
    y = "Número de términos GO específicos",
    fill = "Macro-categoría funcional"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(output_dir, "2.BarPlot_SpecificOnly.png"), width = 18, height = 7)

# -------------------------------------------------------
# 9. Barplot proporcional (porcentaje dentro del clúster)
# -------------------------------------------------------
summary_table_prop <- summary_table %>%
  group_by(Cluster) %>%
  mutate(Proportion = n_terms / sum(n_terms)) %>%
  ungroup()

ggplot(summary_table_prop, aes(x = Cluster, y = Proportion, fill = parentTerm)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  #scale_fill_manual(values = macro_colors) +
  labs(
    title = "Proporción de GO terms específicos por macro-categoría funcional",
    x = "Subpoblación",
    y = "Proporción",
    fill = "Macro-categoría funcional"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(output_dir, "2.BarPlot_SpecificOnly_Proporcion.png"), width = 18, height = 7)



macro_category_map <- c(
  # Metabolism
  "amide metabolic process" = "Metabolism",
  "macromolecule glycosylation" = "Metabolism",
  "small molecule catabolic process" = "Metabolism",
  "mitochondrial translation" = "Metabolism",
  
  # RNA Processing / Translation
  "translational initiation" = "RNA processing / Translation",
  "mRNA destabilization" = "RNA processing / Translation",
  "positive regulation of mRNA metabolic process" = "RNA processing / Translation",
  "ribosomal large subunit biogenesis" = "RNA processing / Translation",
  
  # Cell migration / Motility
  "locomotory behavior" =  "Cell migration / Motility",
  "positive regulation of endothelial cell migration" = "Cell migration / Motility",
  
  # Development / Morphogenesis
  "developmental cell growth" = "Development / Morphogenesis",
  "developmental maturation" = "Development / Morphogenesis",
  "morphogenesis of a branching epithelium" = "Development / Morphogenesis",
  "establishment or maintenance of cell polarity" = "Development / Morphogenesis",
  "mesenchymal cell differentiation" = "Development / Morphogenesis",
  "bone development" = "Development / Morphogenesis",
  "male gonad development" = "Development / Morphogenesis",
  "dendrite morphogenesis" = "Development / Morphogenesis",
  "neural tube development" = "Development / Morphogenesis",
  
  # Hemostasis / Coagulation
  "blood coagulation" = "Hemostasis / Coagulation",
  "coagulation" = "Hemostasis / Coagulation",
  
  # Cell Cycle / Proliferation / Growth
  "positive regulation of mitotic cell cycle" = "Cell Cycle / Proliferation / Growth",
  "negative regulation of epithelial cell proliferation" = "Cell Cycle / Proliferation / Growth",
  "cellular component maintenance" = "Cell Cycle / Proliferation / Growth",
  
  # Apoptosis / Cell Death
  "intrinsic apoptotic signaling pathway by p53 class mediator" = "Apoptosis / Cell Death",
  "regulation of cell killing" = "Apoptosis / Cell Death",
  
  # Cell Adhesion / Junctions
  "homotypic cell-cell adhesion" = "Cell adhesion / Junctions",
  "positive regulation of supramolecular fiber organization" = "Cell adhesion / Junctions",
  "synapse assembly" = "Cell adhesion / Junctions",
  
  # Signaling / MAPK / TGF-beta
  "p38MAPK cascade" = "Signaling / MAPK / TGF-beta",
  "cellular response to transforming growth factor beta stimulus" = "Signaling / MAPK / TGF-beta",
  "canonical Wnt signaling pathway" = "Signaling / MAPK / TGF-beta",
  
  # Oxidative Stress / ROS / Radiation
  "cellular response to reactive oxygen species" = "Stress / Hypoxia / ROS",
  "regulation of reactive oxygen species metabolic process" = "Stress / Hypoxia / ROS",
  "response to radiation" = "Stress / Hypoxia / ROS",
  
  # Immune Response / Antigen Presentation
  "antigen processing and presentation of peptide antigen" = "Immune response / Antigen presentation",
  "regulation of leukocyte migration" = "Immune response / Antigen presentation",
  "negative regulation of inflammatory response" = "Immune response / Antigen presentation",
  "regulation of interleukin-6 production" = "Immune response / Antigen presentation",
  "response to virus" = "Immune response / Antigen presentation",
  "regulation of viral life cycle" = "Immune response / Antigen presentation",
  "positive regulation of peptidase activity" = "Immune response / Antigen presentation",
  
  # Homeostasis / Transport
  "calcium ion homeostasis" = "Homeostasis / Transport",
  "mitochondrial transport" = "Homeostasis / Transport",
  "protein localization to endoplasmic reticulum" = "Homeostasis / Transport",
  
  # Differentiation / Specialization
  "central nervous system neuron differentiation" = "Differentiation / Specialization",
  "regulation of muscle system process" = "Differentiation / Specialization",
  "regulation of osteoblast differentiation" = "Differentiation / Specialization",
  
  # Protein / DNA Organization
  "protein-DNA complex organization" = "Protein / DNA Organization",
  
  # Physiological Rhythms
  "rhythmic process" = "Physiological rhythms"
)

# Asignar la macro categoría
summary_table$MacroCategory <- recode(summary_table$parentTerm, !!!macro_category_map)
summary_table %>% filter(is.na(MacroCategory))  # para revisar los no asignados

macro_colors <- c(
  "Metabolism"                             = "#B2DF8A",  # verde claro
  "RNA processing / Translation"           = "#984EA3",  # púrpura fuerte
  "Cell migration / Motility"              = "#1F78B4",  # azul intermedio
  "Development / Morphogenesis"            = "#A6611A",  # marrón
  "Hemostasis / Coagulation"               = "#FFD92F",  # amarillo brillante
  "Cell Cycle / Proliferation / Growth"    = "#41B7C4",  # azul cielo
  "Apoptosis / Cell Death"                 = "#F781BF",  # gris oscuro
  "Cell adhesion / Junctions"              = "#A6CEE3",  # azul muy claro
  "Signaling / MAPK / TGF-beta"            = "#FF7F00",  # naranja
  "Stress / Hypoxia / ROS"                 = "#542788",  # violeta intenso
  "Immune response / Antigen presentation" = "#FB8072",  # rojo
  "Homeostasis / Transport"                = "#CCEBC5",  # verde menta
  "Differentiation / Specialization"       = "#BEBADA",  # gris medio
  "Protein / DNA Organization"             = "#6A3D9A",  # morado
  "Physiological rhythms"                  = "#F4DBB5"   # lavanda
)


ggplot(summary_table, aes(x = Cluster, y = n_terms, fill = MacroCategory)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_manual(values = macro_colors) +
  labs(
    title = "GO terms agrupados en macro-categorías funcionales",
    x = "Subpoblación",
    y = "Número de términos GO",
    fill = "Macro-categoría"
  )

ggsave(filename = file.path(output_dir, "3.BarPlot_Categ_SpecificOnly.png"), width = 12, height = 7)



summary_table_prop <- summary_table %>%
  group_by(Cluster) %>%
  mutate(Proportion = n_terms / sum(n_terms)) %>%
  ungroup()

ggplot(summary_table_prop, aes(x = Cluster, y = Proportion, fill = MacroCategory)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = macro_colors) +
  labs(
    title = "Proporción de GO terms específicos por macro-categoría funcional",
    x = "Subpoblación",
    y = "Proporción",
    fill = "Macro-categoría funcional"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(output_dir, "3.BarPlot_Categ_SpecificOnly_Proporcion.png"), width = 12, height = 7)

