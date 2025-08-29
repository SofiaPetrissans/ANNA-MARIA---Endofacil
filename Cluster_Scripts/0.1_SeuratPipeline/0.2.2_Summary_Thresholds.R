################################################################################

# SCRIPT: Resumen Filtering Seurat Object threshold ARN mitocondrial
# AUTHOR: SOFIA PETRISSANS MOLL
# DATE: 06-06-2025

################################################################################


# LIBRARIES --------------------------------------------------------------------
source('/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(readxl)
library(writexl)
library(dplyr)

# THRESHOLDS USADOS -----------------------------------------------------------
thresholds <- c(5, 10, 15, 20, 25, 30, 35, 40)

# Ruta base donde están las carpetas mt5, mt10, etc.
base_path <- paths_AnnaMaria$Path_seurat_object
data_path <- file.path(base_path,"2.Filtering")  

# 1) RESUMEN FILTRADO ---------------------------------------------------------

summary_list <- lapply(thresholds, function(th) {
  path <- file.path(data_path, paste0("mt", th), "Results_Filtered.xlsx")
  df <- read_xlsx(path)
  df$Mito_Threshold <- th
  return(df)
})

resumen_filtrado <- bind_rows(summary_list)
write_xlsx(resumen_filtrado, file.path(data_path, "Resumen_Filtering_Thresholds.xlsx"))
