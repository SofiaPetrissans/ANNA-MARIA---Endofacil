
# =========================
# LIBRERÍAS
# =========================
source("/ijc/USERS/spetrissans/ANNA_MARIA/0.1_Scripts/0.5_Annotation_EC/0.0_Paths.R")
library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)
library(paletteer)


# =========================
# 1. PATHS
# =========================
MT_THRESHOLD <- "mt20"

path.guardar_original <- paths_AnnaMaria$Path_guardar
path.guardar <- file.path(path.guardar_original, MT_THRESHOLD, "3.Zonation_Score")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

base_path <- paths_AnnaMaria$Path_seurat_object_first
input_rds <- file.path(base_path,  MT_THRESHOLD, "8.RemoveCluster/2.Integration", "SubsetEndothelial_Harmony.rds")
seurat_object <- readRDS(input_rds)



# Vectores de genes para cada zona
periportal <- c("Adam23", "Msr1", "Dll4", "Efnb2", "Glul")
pericentral <- c("Rspo3", "Thbd", "Wnt9b", "Fabp4", "Alad", "Lgals1")
midzonal <- c("Lyve1", "Ccnd1", "Bok")

# Agrupar en lista con nombres
zonation_genes <- list(
  Periportal = periportal,
  Midzonal = midzonal,
  Pericentral = pericentral
)

# Suponiendo que tu objeto Seurat se llama `seurat_object`
seurat_object <- AddModuleScore_UCell(seurat_object, features = zonation_genes)

fp <- FeaturePlot(seurat_object, features = c("Periportal_UCell", "Midzonal_UCell", "Pericentral_UCell"), ncol = 3)
ggsave(filename = file.path(path.guardar, "Zonation_Score.png"), plot = fp, width = 12, height = 6)



# ..................................................................................................................

periportal <- c("Acer2", "Adam23", "Alb", "Anpep", "Apof", "Apom", "Asgr2", "Ass", "Atp5a1", "Btnl9", "C1s", "C8b", "Ca36", "Cpt2", "Cyp2f2", "Dak", 
                "Dll4", "Eef1b2", "Efnb2", "Elovl2", "Fads1", "Fbp1", "Gc", "Glul", "Gnmt", "Hsd17b13", "Ifgfals", "Ifitm3", "Igf1", "Insr", 
                "Itga9", "Lest", "Ltbp4", "Ltpb4", "Msr1", "Ndufb10", "Nena", "Pck1", "Pigr", "S100a1", "Serpina1c", "Serpina1e", "Serpind1", 
                "Serpinf1", "Trf", "Uqcrh", "Vtn")
midzonal<- c("Aass", "Ac126280.1", "Apoe", "Bok", "Ccnd1", "Clectb", "Cts1", "Cyp1a2", "Cyp2e2", "Cyp4bt", "Eng", "Fcn3", "Fita", "Icam1", "Lyve1")
pericentral <- c("Alad", "Aldh1a1", "Bmp2", "C6", "Cdh13", "Cml2", "Cpox", "Csad", "Cyb5", "Cyp1a2", "Cyp2c37", "Cyp2c50", "Cyp2e1", "Cyp3a11", 
                 "Fabp4", "Gas6", "Gatm", "Gstm1", "Hpd", "Kit", "Lect2", "Lgals1", "Mgst1", "Oat", "Pippt", "Pixnct", "Pon1", "Prodh", "Ptgs1", 
                 "Rab3b", "Rgn", "Rspo3", "Slc16a10", "Stab1", "Tcim", "Thbd", "Wnt2", "Wnt9b")
portal_area <- c("Adgrg6", "Atp13a3", "Cavin3", "Cmkirt", "Gja5", "Lmo7", "Nrgt", "Ntn4", "Plac8", "Sde1", "Vw")
central_vein <- c("Akap13", "Entpdt", "Jptt", "Lhx6", "Ndufa8", "Ramp3", "Rnft65", "Thbd", "Wnt9b")


# Agrupar en lista con nombres
zonation_genes <- list(
  Periportal = periportal,
  Midzonal = midzonal,
  Pericentral = pericentral,
  PortalArea = portal_area,
  CentralVein = central_vein
)

# Suponiendo que tu objeto Seurat se llama `seurat_object`
seurat_object <- AddModuleScore_UCell(seurat_object, features = zonation_genes)

fp <- FeaturePlot(seurat_object, features = c("Periportal_UCell", "Midzonal_UCell", "Pericentral_UCell", "PortalArea_UCell", "CentralVein_UCell"), ncol = 3) & NoAxes()
ggsave(filename = file.path(path.guardar, "UCell_2.png"), plot = fp, width = 18, height = 12)


# ..........................................
# Después de ejecutar AddModuleScore_UCell, puedes usar apply() para extraer la zona con mayor score para cada célula 
# y guardarlo como una nueva columna en @meta.data:
# Vectores de genes para cada zona
periportal <- c("Adam23", "Msr1", "Dll4", "Efnb2", "Glul")
pericentral <- c("Rspo3", "Thbd", "Wnt9b", "Fabp4", "Alad", "Lgals1")
midzonal <- c("Lyve1", "Ccnd1", "Bok")
portal_area <- c("Adgrg6", "Atp13a3", "Cavin3", "Cmkirt", "Gja5", "Lmo7", "Nrgt", "Ntn4", "Plac8", "Sde1", "Vw")

# Agrupar en lista con nombres
zonation_genes <- list(
  Periportal = periportal,
  Midzonal = midzonal,
  Pericentral = pericentral,
  PortalArea = portal_area
)

# Suponiendo que tu objeto Seurat se llama `seurat_object`
seurat_object <- AddModuleScore_UCell(seurat_object, features = zonation_genes)

# Obtener solo las columnas de scores
score_matrix <- seurat_object@meta.data[, c("Periportal_UCell", "Midzonal_UCell", "Pericentral_UCell", "PortalArea_UCell")]
# Calcular cuál zona tiene el mayor score
seurat_object$Dominant_Zonation <- colnames(score_matrix)[apply(score_matrix, 1, which.max)]
dp <- DimPlot(seurat_object, group.by = "Dominant_Zonation", label = TRUE) +
  ggtitle("Dominant Zonation per Cell") & NoAxes()
ggsave(filename = file.path(path.guardar, "Which_Max.png"), plot = dp, width = 7, height = 6)


# Crear un Zonation Score continuo con esa lógica. La idea sería calcular la diferencia entre los scores periportales y pericentrales, de forma que:
  # Valores positivos → perfil más periportal
  # Valores negativos → perfil más pericentral
  # Cercano a cero → expresión intermedia o "midzonal"
library(dplyr)
library(paletteer)
# Crear el Zonation Score como diferencia
seurat_object$Zonation_Score <- seurat_object$Periportal_UCell - seurat_object$Pericentral_UCell

seurat_object$Zonation_Category <- case_when(
  seurat_object$Zonation_Score > 0.3 ~ "Periportal",
  seurat_object$Zonation_Score < -0.3 ~ "Pericentral",
  TRUE ~ "Midzonal"
)

fp <- FeaturePlot(seurat_object, features = "Zonation_Score", slot = "data") +
  scale_colour_gradient2(
    low = "#4575b4",      # azul pericentral
    mid = "#f7f7f7",      # gris claro visible sobre blanco
    high = "#d73027",     # rojo periportal
    midpoint = 0,
    name = "Zonation"
  ) +
  ggtitle("Zonation Score (Periportal - Pericentral)") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave(filename = file.path(path.guardar, "Zonation_Score.png"), plot = fp, width = 6, height = 6)


dp <- DimPlot(seurat_object, group.by = "Zonation_Category", label = TRUE) & NoAxes()
ggsave(filename = file.path(path.guardar, "Zonation_Category.png"), plot = dp, width = 6, height = 6)




