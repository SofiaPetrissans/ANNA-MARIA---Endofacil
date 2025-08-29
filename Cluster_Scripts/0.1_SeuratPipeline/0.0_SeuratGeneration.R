################################################################################

# SCRIPT: Generate Seurat Object 
# AUTOR: SOFIA PETRISSANS MOLL
# FECHA: 05.06.2025

########################################################################################################################

# LIBRERIAS 
source('/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(Seurat)
library(tidyverse)
library(readxl)
library(writexl)




# Crear la carpeta de resultados
path.guardar_original<-paths_AnnaMaria$Path_guardar
path.guardar <- file.path(path.guardar_original, "0.Seurat_Object")
dir.create(path.guardar, showWarnings = FALSE, recursive = TRUE)


########################################################################################################################
# Ruta base de los datos
base_path <- paths_AnnaMaria$Path_datos

# Leer Excel con info de las muestras
samples_info <- read_excel(paths_AnnaMaria$Path_sample_info)
samples_name <- unique(samples_info$Sample)

# QC definition for future QC control

pattern.mito<-"^Mt-"
pattern.ribo<-"^Rps"

# Inicializar
Lista_Seurat_Object <- vector("list", length(samples_name))
names(Lista_Seurat_Object) <- samples_name
Samples_Info <- data.frame(Sample=samples_name, NumberOfGenes=0, NumberOfCells=0, MeanOfGenesPerSample=0)

# Crear objetos Seurat
for (i in seq_along(samples_name)) {
  sn <- samples_name[i]
  cat("Procesando:", sn, "\n")
  
  path_data <- file.path(base_path, paste0(sn, "_output"))
  print(paste("Intentando leer:", path_data))
  print(list.files(path_data))

  counts <- Read10X(data.dir = path_data)
  
  seu <- CreateSeuratObject(counts = counts, project = sn)
  info <- samples_info[samples_info$Sample == sn, ]
  
  # Añadir metadatos
  seu$Sample <- sn
  seu$Phenotype <- info$Phenotype
  seu$Tissue <- info$Tissue
  seu$percent.mt <- PercentageFeatureSet(seu, pattern = "(?i)^mt-")
  seu$percent.rb <- PercentageFeatureSet(seu, pattern = "(?i)^rps")

  # Guardar resumen
  Samples_Info$NumberOfGenes[i] <- nrow(seu@assays$RNA@counts)
  Samples_Info$NumberOfCells[i] <- ncol(seu@assays$RNA@counts)
  Samples_Info$MeanOfGenesPerSample[i] <- mean(seu$nFeature_RNA)

  # Guardar objeto individual
  saveRDS(seu, file = file.path(path.guardar, paste0(sn, ".rds")))
  Lista_Seurat_Object[[i]] <- seu
}

# Guardar resumen en Excel
write_xlsx(Samples_Info, file.path(path.guardar, "Samples_Info.xlsx"))

# Combinar todos
Seurat_Object_Total <- merge(x = Lista_Seurat_Object[[1]], y = Lista_Seurat_Object[-1])
saveRDS(Seurat_Object_Total, file.path(path.guardar, "RawObject_Total.rds"))






