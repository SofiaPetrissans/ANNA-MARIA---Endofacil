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

## 🗄️ Estructura de carpetas - clúster

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
