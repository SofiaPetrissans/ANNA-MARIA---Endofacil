# SCRIPT: Remove the clusters that make noise
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 16.07.2024

##############################################################################################################################

source("/ijc/USERS/amartinezl/EndoTension/Paths.R")
directory<-setwd(Lista_Paths_Main$MainPath_Scripts)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(SoupX)
library(harmony)
library(Rcpp)
library(optparse)
library(AUCell)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
col <-  getPalette(30)
# Path guardar 

# ........................................................................................................
# Selection of Arguments 
# ........................................................................................................

option_list <- list(make_option(c('-p', '--phenotype'), type="character"),make_option(c('-f', '--foldername'), type="character")) 
opt <-parse_args(OptionParser(option_list = option_list))

if (is.null(opt$p)) opt$p <- "PSVD"
if (is.null(opt$f)) opt$f <- "Remove_ALB_RemoveRes05_Clus9"

# ........................................................................................................
# Stablized the path depending on the pheno 
# ........................................................................................................

# Obtener path de anotaciones endoteliales
path.guardar_original <- paste(Lista_Paths_Main$MainPath_Res, "0.2_EndothelialComparment",sep="/")
path.guardar<-paste(path.guardar_original,opt$p,opt$f,"Integration/AUCell",sep="/")
dir.create(path.guardar)

source(paste(Lista_Paths_Main$Path_Utils,"Util_Annotations.R",sep="/"))
##############################################################################################################################

path.obj<-paste(Lista_Paths_Main$MainPath_Res,"0.2_EndothelialComparment",opt$p,opt$f,"Integration",sep="/")
data<-readRDS(paste(path.obj,"Harmony.rds",sep="/"))

signature<-c("ALB","TF","TTR","HNF4A","CYP2A6")

# ..........................................................................................................................................................
# AUCell Analysis
# ..........................................................................................................................................................

# FUNCIONES
# AUCELL Estimate Signature .......................................................................................................................................

Funcion_AUCell_Signature <- function(data, geneSets, path.guardar, assay = "RNA") {
  require(AUCell)
  # Obtener matriz de expresión del assay indicado
  expMatrix <- GetAssayData(data, assay = assay, slot = "counts")

  # Calcular scores AUC
  cells_AUC <- AUCell_run(expMatrix, geneSets)
  saveRDS(cells_AUC, file = file.path(path.guardar, "cells_AUC.rds"))

  # Determinar thresholds automáticamente
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign = TRUE) 
  saveRDS(cells_assignment, file = file.path(path.guardar, "cells_assignment.rds"))

  # Guardar histogramas claramente
  pdf(file.path(path.guardar, "Histogram_GeneSetActivation_CriteriaSelection.pdf"),
      width = 10, height = 10)
  AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign = TRUE) 
  dev.off()

  return(list(cells_AUC = cells_AUC, cells_assignment = cells_assignment))
}

Select_ActiveCells_AUCell <- function(data, cells_AUC, geneSet_name, 
                                      threshold_selected = NULL, 
                                      threshold_name=NULL,
                                      path.guardar, assay = "RNA") {

  # Si no se da threshold, parar claramente la ejecución
  if (is.null(threshold_selected)) {
    stop("Por favor, especifica un threshold seleccionado manualmente o utiliza el automático de AUCell_exploreThresholds.")
  }

  # Validar existencia del geneSet
  if (!(geneSet_name %in% rownames(getAUC(cells_AUC)))) {
    stop("El nombre del geneSet no existe en cells_AUC.")
  }

  # Crear nombre de archivo para guardar plot con threshold claramente indicado
  file_name <- paste0("ThresholdSelected_", geneSet_name, "_", threshold_name, ".pdf")
  
  # Guardar gráfico claramente indicando el threshold usado
  pdf(file.path(path.guardar, file_name), width = 10, height = 6)
  AUCell_plotHist(cells_AUC[geneSet_name, ], aucThr = threshold_selected)
  abline(v = threshold_selected, col = "red", lwd = 2)
  dev.off()

  # Clasificación binaria células activas/inactivas usando threshold seleccionado
  active_cells <- colnames(cells_AUC)[
    getAUC(cells_AUC)[geneSet_name, ] >= threshold_selected
  ]

  # Añadir resultado a objeto Seurat incluyendo claramente el threshold usado
  new_col_name <- paste0(geneSet_name, "_AUCell_Active_", threshold_name)
  
  data[[new_col_name]] <- ifelse(
    colnames(data) %in% active_cells, "Active", "Inactive"
  )
  return(data)
}



# Estimate the Values for AUCell 
res_aucell_list<-Funcion_AUCell_Signature(data, signature, path.guardar, assay = "RNA")
cells_AUC<-res_aucell_list$cells_AUC
cells_assignment <- res_aucell_list$cells_assignment

# The reference threshold: 
# Extraer threshold recomendado para "PericyteSignature"
threshold_selected <- cells_assignment$geneSet$aucThr$selected
threshold_selected<-round(threshold_selected,2)

# ..........................................................................................................................................................
# Plot different thresholds: 
# ..........................................................................................................................................................

# Obtener valores AUC para tu firma
auc_values <- getAUC(cells_AUC)["geneSet", ]

# Base Thresholds
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_selected,threshold_name="Default",path.guardar, assay = "RNA")

# Percentil 95%
# Calcular el percentil 95
threshold_95 <- quantile(auc_values, probs = 0.95, na.rm = TRUE)
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_95,threshold_name="Percentil_95",path.guardar, assay = "RNA")

# Media 
threshold_mean<-mean(auc_values)
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_mean,threshold_name="Mean",path.guardar, assay = "RNA")

# SD 
threshold_sd<-sd(auc_values)
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_sd,threshold_name="SD",path.guardar, assay = "RNA")

# mean + 2sd 
threshold_2sd_mean<-mean(auc_values) + 2 * sd(auc_values)
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_2sd_mean,threshold_name="Mean_2SD",path.guardar, assay = "RNA")

# mean + 2sd 
threshold_1sd_mean<-mean(auc_values) + 1 * sd(auc_values)
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_1sd_mean,threshold_name="Mean_1SD",path.guardar, assay = "RNA")

# manual
threshold_manual<-0.15
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_manual,threshold_name="Manual",path.guardar, assay = "RNA")

# manual
threshold_manual<-0.2
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_manual,threshold_name="Manual_0.2",path.guardar, assay = "RNA")

# manual
threshold_manual<-0.6
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_manual,threshold_name="Manual_0.6",path.guardar, assay = "RNA")

# manual
threshold_manual<-0.4
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_manual,threshold_name="Manual_0.4",path.guardar, assay = "RNA")

# manual
threshold_manual<-0.5
data<-Select_ActiveCells_AUCell(data, cells_AUC, "geneSet", threshold_selected = threshold_manual,threshold_name="Manual_0.5",path.guardar, assay = "RNA")

# ..........................................................................................................................................................
# Feature Plot with the different option
# Select Thresholds
# ..........................................................................................................................................................

metadata<-data@meta.data
umaps_dim<-data@reductions$umap@cell.embeddings
metadata<-cbind(metadata,umaps_dim)

names_threshold<-c("Default","Percentil_95","Mean","SD","Mean_2SD","Mean_1SD","Manual","Manual_0.2")
names_threshold<-paste("geneSet_AUCell_Active_",names_threshold,sep="")

Lista_Threshold<-list()

for(i in seq_along(names_threshold)){
    thres<-names_threshold[i]

    Fp<-ggplot(metadata, aes(x = umap_1, y = umap_2, col = .data[[thres]])) +
    geom_jitter(size = 1.5) +
    theme_bw() + 
    scale_color_manual(values = c("#7B1FA2","grey")) +  # Ensure 'col_final_perissubtypes' is defined
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = "none"
    ) +
    xlab("UMAP_1")+
    ylab("UMAP_2")+
    labs(title=thres)
    file_name<-paste("FP",thres,".png",sep="_")
    ggsave(plot=Fp,filename=paste(path.guardar,file_name,sep="/"),width=7,height=7)
}

saveRDS(data,paste(path.guardar,"Data_AUCell.rds",sep="/"))