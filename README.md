# ANNA-MARIA - Endofacil

## 💡 Introducción

Este repositorio documenta el análisis de un conjunto de scRNA-seq de hígado de ratón (proyecto Anna Maria) con dos condiciones (vehicle y CHOP). Las muestras son:

| Condición |	Enriquecimiento |	ID de muestra |
|---|---|---|
| CHOP | EC-enriched | L1019 | 
| CHOP	| no-enriched	| L569 |
| Vehicle |	EC-enriched |	L1071 |
| Vehicle	| no-enriched |	L566 |

El hígado presenta alto contenido mitocondrial basal, lo que convierte el **filtrado por %MT en un punto crítico**: si el umbral es muy estricto, se pierden ECs válidas; si es laxo, se conserva ruido. Para abordar esto:

Ejecutamos **toda la pipeline de Seurat y la anotación inmune/endotelial/global** con una rejilla de umbrales de MT: **mt5, mt10, mt15, …, mt40**.

A partir de los resultados de calidad, estructura y anotación, **focalizamos el análisis en mt15 y mt20**.

Para los análisis endoteliales finales **seleccionamos mt15**, realizamos **integración con Harmony**, **subseteamos ECs**, anotamos subpoblaciones endoteliales y ejecutamos el downstream analysis.

## 🗂️ Dónde se ejecutó cada parte

### En el clúster
Descarga/preparación de datos, QC, SeuratPipeline (normalización → HVGs → escala → PCA → UMAP → clustering), detección de dobletes (Scrublet/DoubletFinder), anotación de poblaciones globales (EC, NK, B cells, etc.), pruebas multiumbral %MT, subset endotelial, integración (Harmony) y generación del objeto Seurat de ECs listo para descargar.

### En local (RStudio)
Carga del subset endotelial desde el clúster, anotación fina de subpoblaciones EC y análisis downstream específicos (p. ej., DEGs, programas, enriquecimientos, visualizaciones).

## ☁️ Clúster
### 🗄️ Estructura de carpetas 

Ruta base de scripts:

```
/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts
```

-  `0.1_SeuratPipeline/`  →  Generación del objeto, QC y filtrado por rejilla %MT; pipeline de Seurat; doublets (Scrublet export/import + DoubletFinder); integración (Harmony).
        -  `0.0_SeuratGeneration.R` (crea objetos)
        -  `0.1_Estimating_QC.R` (métricas QC, decisión de %MT)
        -  `0.2_Filtering.R` (aplica umbrales %MT)
        -  `0.4_DoubletFinder.R` + `0.5.A/B/C/D` (Scrublet & DF)
        -  `0.6_Integration.R` (Harmony)

-  `0.2_Annotation_All_Popultations/`  →   Anotación global por marcadores y apoyo con AUCell (Panglao), tanto antes como después de integración; paneles de porcentaje por tipo celular. Ejemplos:
        -  `0.1.A_Markers_Annotation.R`, `0.1.C_Annotation_AUCell_Panglao.R`, `0.4_Plots_Porcetages.R`.

-  `0.3_Downstream_Analysis_all_populations/`  →  Análisis global (no restringido a EC):
        -  DEGs por tratamiento, firmas de senescencia y GO de prueba.

-  `0.4_SubsetEndothelial/`  →  Subsetting de ECs, limpieza fina, anotación por marcadores, tests de propeller (composición), funciones de zonación.
  
-  `0.5_Annotation_EC/`  →  Evaluación de granularidad de clustering, anotación endotelial detallada, firmas de zonación/senescencia, versiones mt15 y mt20 para comparación.

-  `0.6_Downstream_Analysis_EC/`  →  Downstream específico de EC: DEGs entre tratamientos, utilidades de zonación para EC.
  
-  `0.7_Traicing/`   →  Trazabilidad y union con análisis de las muestras sin enriquecer

Cada subcarpeta incluye `*_Exe.sh` para facilitar la ejecución por lotes en el clúster y `0.0_Paths.R` para centralizar rutas.


## 💻 Local (RStudio)

### ⬇️ Qué debes descargar del cluster antes de empezar

-  **Objeto Seurat** (subset EC, mt15 + Harmony)
    -  Fichero: `SubsetEndothelial_Harmony.rds`  →  Generado en cluster tras 0.4_SubsetEndothelial y la integración con Harmony.
    -  Ruta orientativa: `/ijc/LABS/GRAUPERA/RAW/.../ANNA_MARIA/3.MERGE/0.4_SubsetEndothelial/mt15/9.RemoveCluster/2.Integration/SubsetEndothelial_Harmony.rds`

-  **DEGs por tipo celular** (mt15, “all populations”)
    -  Excel: `DEGs_per_CellType_all_genes.xlsx`
    -  Producido por: `/ijc/.../ANNA_MARIA/3.MERGE/0.1_Scripts/0.3_Downstream_Analysis_all_populations/0.1_DEG_tratamientos.R`
    -  Ruta orientativa: `/ijc/LABS/GRAUPERA/RAW/.../ANNA_MARIA/3.MERGE/0.3_Downstream_Analysis_all_populations/mt15/1.DEGs/DEGs_per_CellType_all_genes.xlsx`

-  **DEGs CHOP vs vehicle** (global o subset EC, una sola hoja)
    -  Excel: `DEGs_CHOP_vs_vehicle.xlsx`
    -  Ruta orientativa: `/ijc/LABS/GRAUPERA/RAW/.../ANNA_MARIA/3.MERGE/0.6_Downstream_Analysis_EC/mt15/1.DEGs_phenotype/DEGs_CHOP_vs_vehicle.xlsx`

### 🗄️ Estructura de carpetas 

-  `1.GOTerms_mt15_merge.R`
  
        -  **Qué hace**: enriquece GO:BP con clusterProfiler usando `DEGs_per_CellType_all_genes.xlsx` (todas las poblaciones, mt15). Exporta tres Excels (All/Up/Down) y barplots por población.
        -  **Input**: `DEGs_per_CellType_all_genes.xlsx`
        -  **Output**: `GO_BP_All_CellTypes.xlsx`, `GO_BP_Upregulated_CellTypes.xlsx`, `GO_BP_Downregulated_CellTypes.xlsx` + /Plots_…

-  `2.GOTerms_mt15_Subset_EC_merge.R`
  
        -  **Qué hace**: idem pero usando la lista de DEGs CHOP vs vehicle (una sola hoja) enfocada a EC. Exporta GO (All/Up/Down) y barplots.
        -  **Input**: `DEGs_CHOP_vs_vehicle.xlsx`
        -  **Output**: `GO_BP_All_Subset.xlsx`, `GO_BP_Up_Subset.xlsx`, `GO_BP_Down_Subset.xlsx` + /Plots_…

-  `3.Annotations_Layer.Rmd`
  
        -  **Qué hace**: carga `SubsetEndothelial_Harmony.rds`, define capas de anotación (Layer1–Layer5), genera DimPlots, DotPlots y Heatmaps de canonical markers; explora subclustering de un cluster concreto y compara proporciones CHOP vs vehicle por Layer.
        -  **Input**: `SubsetEndothelial_Harmony.rds`
        -  **Output**: `SubsetEndothelial_Harmony_Annotated.rds` + PNGs en `9.Annotations_Layer/`

-  `4.GO_Terms_AnnotationLayer5.R`
  
        -  **Qué hace**: para cada subpoblación de Layer5: DEGs “vs resto” y GO:BP (genes up, padj<0.05). Exporta Excel por subpoblación y barplots.
        -  **Input**: `SubsetEndothelial_Harmony_Annotated.rds`
        -  **Output**: `10.Results_GO_AnnotationLayer/GO_Results_AnnotationLayer5.xlsx` + figuras

-  `5.RRVGO_AnnotationLayer.R`
        -  **Qué hace**: agrupa GO por similitud semántica (rrvgo/GOSemSim), crea categorías funcionales y macro-categorías (mapeo manual), y resume por subpoblación (barras absolutas y proporciones).
        -  **Input**: resultados del paso 4
        -  **Output**: `11.Analysis_GO_AnnotationLayer/` (treemaps/scatter, tablas, barplots)

-  `6.subpop_DEGs_CHOPvsVehicle_layers1-5_GO_UpSet_corr.R`
        -  **Qué hace**: para cada Layer (1–5) y sus subpoblaciones, calcula DEGs CHOP vs vehicle; resume #Up/#Down, proporciones, GO Up por subpoblación, UpSet de DEGs Up compartidos y correlaciones entre DEGs globales y DEGs de subpoblación (identifica qué subpoblación “conduce” el efecto).
        -  **Input**: `SubsetEndothelial_Harmony_Annotated.rds` + `DEGs_CHOP_vs_vehicle.xlsx`
        -  **Output**: `12.DEGs_CHOP_vs_vehicle/` (Excels de DEGs por Layer, barplots, GO Up, UpSet, matrices de correlación)

-  `7.Markers_GO_HeatMaps_AL5.R`
        -  **Qué hace**: a partir de los GO más informativos (del paso 4–5), define firmas génicas por subpoblación de Layer5 (LVEC, LSEC_1…LSEC_5) y genera ModuleScores, DotPlots y Heatmaps (vertical y horizontal) con AverageExpression.
        -  **Input**: `SubsetEndothelial_Harmony_Annotated.rds`
        -  **Output**: PNGs de DotPlots y Heatmaps por AL5

-  `8.Propeller.R`
        -  **Qué hace**: con `speckle::propeller` estima cambios de composición CHOP vs vehicle por Layer5; guarda heatmap de PropMean (con p-values/FDR) y barplots (medias y log2FC).
        -  **Input**: `SubsetEndothelial_Harmony_Annotated.rds`
        -  **Output**: `Res_Propeller.rds`, `Heatmap_PropMean_ByPheno_ByCt.pdf`, barras PNG

-  `9.Senescence_Markers.R`
        -  **Qué hace**: puntúa senescencia/SASP/p53 con AUCell (y firmas de msigdbr), genera FeaturePlots, DimPlots binarios, violines y barras de proporciones por Layer/condición.
        -  **Input**: `SubsetEndothelial_Harmony_Annotated.rds`
        -  **Output**: PNGs (UMAP, violines, barras), columnas AUC_* añadidas al objeto



