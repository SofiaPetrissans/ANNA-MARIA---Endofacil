################################################################################

#SCRIPT: Funciones para la pipeline de Seurat
#AUTHOR: Sofía Petrissans Moll
#DATE: 20.feb.2025

# INDEX:
  # 1) Seurat Pipeline
  # 2) EstimatePCA
  # 3) DoubletFinder 
  # 4) Harmony_Integration

################################################################################


# 1) Seurat Pipeline  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# Argumentos de entrada: 
  # obj: un seurat object individual 
  # resolution: un vector con resolucion/es
# Objetivo: Procesamiento de los datos: 
  # Nomalización de los datos
  # FindVariablesFeatures
  # Escalado
  # RunPCA
  # FindNeighbors
  # FindCluster
  # RunUmap
# La salida es: el objeto Seurat procesado

SeuratPipeline <- function(objeto, resolution) {
  # Normalizar los datos
  objeto <- NormalizeData(object = objeto)
  # Encontrar características variables iniciales (HVG)
  objeto <- FindVariableFeatures(object = objeto, 
                                 selection.method = "vst", 
                                 nfeatures = 5000)  # Se inicia con un número alto
  variable_genes <- objeto[["RNA"]]@meta.features[["vst.variance.standardized"]]
  # Determinar el punto de codo usando el porcentaje acumulado de varianza
  cumulative_var <- cumsum(sort(variable_genes, decreasing = TRUE))
  optimal_hvg <- which(cumulative_var / max(cumulative_var) > 0.90)[1]  # Ajuste al 90% de la varianza
  # Volver a calcular HVG con el número óptimo encontrado
  objeto <- FindVariableFeatures(object = objeto, 
                                 selection.method = "vst", 
                                 nfeatures = optimal_hvg)
  # Escalar los datos
  objeto <- ScaleData(object = objeto)
  # Ejecutar PCA utilizando las características variables
  objeto <- RunPCA(object = objeto, features = VariableFeatures(object = objeto))
  ## Elección de PCA correcto:
  # Calcular el porcentaje de varianza explicado por cada PC
  pct <- objeto[["pca"]]@stdev / sum(objeto[["pca"]]@stdev) * 100
  # Calcular el porcentaje acumulado
  cumu <- cumsum(pct)
  # Encontrar el primer PC en el que el acumulado supera el 90% y el porcentaje de varianza es menor a 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Encontrar la posición donde la diferencia entre PCs consecutivas es mayor a 0.1
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), 
                decreasing = TRUE)[1] + 1
  # Seleccionar el mínimo de ambos valores
  dim.max <- min(co1, co2)
  # Encontrar vecinos basándose en las dimensiones calculadas
  objeto <- FindNeighbors(object = objeto, dims = 1:dim.max)
  # Realizar el clustering con la resolución especificada
  objeto <- FindClusters(object = objeto, resolution = resolution)
  # Ejecutar UMAP para la reducción de dimensionalidad
  objeto <- RunUMAP(object = objeto, dims = 1:dim.max)
  
  # Retornar el objeto Seurat procesado
  return(objeto)
}

# 2) Estiamte PCA dimensions  ||||||||||||||||||||||||||||||||||||||||||||||||||

# Argumento de entrada: Objeto Seurat
# Objetivo: estimar el número óptimo de dimensiones del PCA
# Salida: número 

EstimatePCA.Dims<-function(objeto){
  # Estimate the PCA Dimmensions.
  pct <- objeto[["pca"]]@stdev / sum(objeto[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  dim.max <- min(co1, co2)
  return(dim.max)
}

# 3) Pipeline de DoubletFinder  ||||||||||||||||||||||||||||||||||||||||||||||||

# Nombre de la función: DoubletFinderPipeline
# Argumentos
  # objeto: Debe ser un objeto Seurat ya procesado (SeuratPipeline)
  # col_annotation: nombre de la columna donde están anotados los clusters
  # doublet_rate: 0.076 por defecto
# Objetivo: determinar los dobletes
# Salida: objeto seurat 

DoubletFinderPipeline <- function(objeto, col_annotation, doublet_rate = 0.076){
  # Utilizar la función para estimar el número óptimo de dimensiones del PCA
  dim.max <- EstimatePCA.Dims(objeto)
  
  ## IDENTIFICACIÓN DE pK CON paramSweep_v3 ------------------------------------
  sweep.res <- paramSweep(objeto, PCs = 1:dim.max, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Seleccionar el pK que maximiza la métrica BCmetric
  pK_val <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>% 
    dplyr::select(pK)
  pK_val <- as.numeric(as.character(pK_val[[1]]))
  
  ## ESTIMACIÓN DE DOUBLETS HOMOTÍPICOS ----------------------------------------
  # Extraer la columna de anotación especificada
  annotations <- objeto@meta.data[[col_annotation]]
  # Calcular la proporción de doublets homotípicos
  homotypic.prop <- modelHomotypic(annotations)
  
  # Calcular el número esperado de doublets según la tasa proporcionada
  nExp <- round(doublet_rate * nrow(objeto@meta.data))
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  ## EJECUCIÓN DE DoubletFinder_v3 ##
  objeto <- doubletFinder(objeto, 
                          PCs = 1:dim.max, 
                          pN = 0.25, 
                          pK = pK_val, 
                          nExp = nExp.adj,
                          reuse.pANN = FALSE, 
                          sct = FALSE)
  
  return(objeto)
}

# 4) Integración con Harmony  ||||||||||||||||||||||||||||||||||||||||||||||||||
# Argumentos de entrada:
  # objeto: Un objeto Seurat ya preprocesado
  # resoluciones: Un vector con las resoluciones a aplicar en el clustering.
  # integration.var: La variable de `meta.data` que se usará para corregir batch effects con Harmony.
# Objetivo:
  # Aplicar integración con Harmony para corregir batch effects.
  # Generar una reducción UMAP basada en la representación corregida de Harmony.
  # Construir la matriz de vecinos en el espacio de Harmony.
  # Aplicar clustering en la representación integrada con diferentes resoluciones.
# Salida:
  # Devuelve el objeto Seurat procesado con la integración de Harmony, UMAP y clustering aplicados.

Harmony_Integration <- function(objeto, resoluciones, integration.var){
  # Utilizar la función para estimar el número óptimo de dimensiones del PCA
  dim.max <- EstimatePCA.Dims(objeto)
  objeto <- RunHarmony(objeto,
                       group.by.vars = integration.var, 
                       dims = 1:dim.max,
                       plot_convergence = TRUE,
                       reduction.save = "HarmonyLog",
                       assay = "RNA")
  if (!integration.var %in% colnames(objeto@meta.data)) {
    stop("⚠️ Error: La variable de integración '", integration.var, "' no existe en `meta.data`. Verifica el nombre.")
  }
  objeto <- RunUMAP(objeto,
                    reduction = "HarmonyLog",
                    dims = 1:dim.max)
  objeto <- FindNeighbors(objeto,
                          reduction = "HarmonyLog", 
                          dims = 1:dim.max,
                          graph.name ="Harmony_Log")
  objeto <- FindClusters(objeto,
                         resolution = resoluciones,
                         graph.name = "Harmony_Log")
  return(objeto)
}

