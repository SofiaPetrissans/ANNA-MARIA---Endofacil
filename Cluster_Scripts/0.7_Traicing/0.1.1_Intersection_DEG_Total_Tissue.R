# SCRIPT: Check from the tissues DEG how many of the genes intersect with the TOP genes from the total 
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 16.12.2024

########################################################################################################################

directory<-setwd("/mnt/beegfs/amartinezl/AnalisisByTissue/")

# Setting working parameters

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(ComplexHeatmap)
library(VennDiagram)
library(ggVennDiagram)
library(UpSetR)
library(dplyr)
library(reshape2)
library(gridExtra)

library(gprofiler2)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


library(org.Hs.eg.db)
library(clusterProfiler)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
col <-  getPalette(40)

# Load the script with the functions 
source("/mnt/beegfs/amartinezl/AnalisisByTissue/0.1_PERICYTES/Utils/0.3_Util_SeuratPipeline.R")
source("/mnt/beegfs/amartinezl/AnalisisByTissue/0.1_PERICYTES/Utils/0.7_Util_Integration.R")
source("/mnt/beegfs/amartinezl/AnalisisByTissue/0.1_PERICYTES/Utils/0.5_Util_Annotations.R")

########################################################################################################################

path.guardar<-"/mnt/beegfs/amartinezl/AnalisisByTissue/AnalysisSetUp/DEG_Analysis/Peris_Tissue/Intersection_Markers"
dir.create(path.guardar)

markers_total<-readRDS("/mnt/beegfs/amartinezl/AnalisisByTissue/0.1_PERICYTES/0.2_Analysis/0.2_Pericytes_Characterization_Phenotype_Peris_VS_Rest_Clusters/Markers/DEG_Peris_Total/Markers_AnnotLayer1_TumorPeris.rds")
markers_total<-markers_total[order(markers_total$avg_log2FC,decreasing=T),]
markers_total_top<-markers_total[1:250,]

#Read the markers for each of the tissues ........................................................................................
path_files<- "0.1_PERICYTES/0.2_Analysis/0.2_Pericytes_Characterization_Phenotype_Peris_VS_Rest_Clusters/Markers/DEG_Peris_Tissue"
list_files_markers<-list.files(path_files)
file_names<-"Markers_TumorPeris.rds"

Lista_Markers_Peris<-vector(mode="list",length=length(list_files_markers))
names(Lista_Markers_Peris)<-list_files_markers

Lista_Markers_Peris_Filtered<-vector(mode="list",length=length(list_files_markers))
names(Lista_Markers_Peris_Filtered)<-list_files_markers

for(i in seq_along(list_files_markers)){
    folder<-list_files_markers[i]
    print(folder)

    file<-readRDS(paste(path_files,folder,file_names,sep="/"))
    file<-file[order(file$avg_log2FC,decreasing=T),]
    file$Tissue<-folder
    file$Gene<-rownames(file)
    Lista_Markers_Peris[[i]]<-file

    # Filter the genes 
    file_filtered<-file[1:250,]
    Lista_Markers_Peris_Filtered[[i]]<-file_filtered

}

# UpSetPlot 

markers_total_vector<-list(TotalMarkers=markers_total_top$Gene)
gene_list_tissue <- lapply(Lista_Markers_Peris_Filtered, function(x) x$Gene)
Lista_total<-c(markers_total_vector,gene_list_tissue)

lt <- Lista_total
m <- make_comb_mat(lt)
m <- m[comb_size(m) >= 3]
cs <- comb_size(m)
  
# Generar el gráfico UpSet y guardarlo en un PDF
file_name<-paste("UpSetPlot_Top_250_UpRegulated_Markers.pdf",sep="_")
pdf(paste(path.guardar, file_name, sep="/"), width=12, height=5)
UpSet(m,
    comb_order = order(comb_size(m), decreasing = TRUE),
    comb_col = c("#B0DFDD","#8ECFA3","#F9BF97","#F3E79A","#EB9DA5","#E0BFF0","#FBCCCC","#B9DDF1","#D3CFE9"),
    top_annotation = upset_top_annotation(m, add_numbers = TRUE),
    bg_col = c("white", "white"), 
    bg_pt_col = "grey",
    pt_size = unit(4, "mm"),
    lwd = unit(2, "mm"),
    right_annotation = NULL,
    row_names_gp = gpar(fontsize = 20))
dev.off()