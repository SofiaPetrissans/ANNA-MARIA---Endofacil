# SCRIPT: Evaluation Clustering
# AUTOR: ANE AMRTINEZ LARRINAGA
# FEHCA: 30.04.2024

##############################################################################################################################


library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cluster)
library(mclust)
library(writexl)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
cols <-  getPalette(12)

path.guardar <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.5_Annotation_EC_mt20/1.EvaluateClustering"
dir.create(path.guardar)

set.seed(12345)

Calculo.Dimensiones.PCA <- function(data) {
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  dim.final <- min(co1, co2)
  return(dim.final)
}

##############################################################################################################################

data<-readRDS("/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.4_SubsetEndothelial/mt15/8.B.RemoveCluster/2.Integration/SubsetEndothelial_Harmony.rds")
data<-SetIdent(data,value="Harmony_Log_res.0.5")

# Define Paramaters: 
resolutions_params<- seq(0.01, 2.00, by = 0.01)


# Step1: Cluster the cells in the different resolutions 
Resultados_Res<-data.frame(resolution=resolutions_params,mean=0,sd=0,ari=0)

for(i in seq_along(resolutions_params)){
    # Step 1: Cluster the cell in a specific clustering algorithm
    res<-resolutions_params[i]
    res_name_complete<- paste("Harmony_Log_res.",res,sep="")
    dim.final<-Calculo.Dimensiones.PCA(data)
    data <- FindNeighbors(data,reduction = "HarmonyLog",dims = 1:dim.final,graph.name = "Harmony_Log")
    data <- FindClusters(data, resolution = res,graph.name = "Harmony_Log")

    data<-SetIdent(data,value=res_name_complete)
    L<-length(unique(Idents(data)))

    if(L==1){
        print(res)
        print("At least two cluster should be found in the resolution used")
        Resultados_Res$mean[i]<-0
        Resultados_Res$sd[i]<-0
        Resultados_Res$ari[i]<-0
        next
    }else{
        # Calculate ARI
        original_clusters <- data$Harmony_Log_res.0.7
        sub_clusters <- Idents(data)
        ari_value <- adjustedRandIndex(original_clusters, sub_clusters)

        # Calculate silhouette
        pc_data <- Embeddings(data, "pca")[, 1:dim.final]
        distance_matrix <- dist(pc_data)
        silhouette_res <- silhouette(x=as.integer(Idents(data)), dist=distance_matrix)
        Resultados_Res$mean[i]<-mean(silhouette_res[, "sil_width"])
        Resultados_Res$sd[i] <- sd((silhouette_res[, "sil_width"])) 
        Resultados_Res$ari[i] <- ari_value
    }
    
}

saveRDS(Resultados_Res,paste(path.guardar,"Resultados_Res_NotSumbSample.rds",sep="/"))
writexl::write_xlsx(Resultados_Res,paste(path.guardar,"Resultados_Res_NotSumbSample.xlsx",sep="/"))

max_mean_resolution <- Resultados_Res$resolution[which.max(Resultados_Res$mean)]
max_mean_value<-Resultados_Res$mean[which.max(Resultados_Res$mean)]

Signature_Total<-ggplot(Resultados_Res, aes(x = resolution, y = mean)) +
geom_point() +
geom_line()+
scale_x_continuous(breaks = sort(unique(Resultados_Res$resolution)))+
geom_segment(aes(x = max_mean_resolution, y = max_mean_value, yend = max_mean_resolution,xend=max_mean_resolution), 
               color = "red", size = 1)+
theme_bw() +  # Optional: a cleaner theme for the plot
labs(title = "Silhouete Score by Resolution",
    x = "Resolution",
    y = "Silhouete Mean Score") +
theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels if needed
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.background = element_blank())+
annotate("text", x = max_mean_resolution, y = max_mean_value, label = sprintf("Max Mean: %.2f at Res: %.2f", max_mean_value, max_mean_resolution), 
           hjust = -5, vjust = -0.5, color = "black", size = 4, angle = 0)
ggsave(plot=Signature_Total,paste(path.guardar,"ClusteringEvaluation_Final_NotSumbSample.png",sep="/"),height = 6,width = 15)

ARI_Plot <- ggplot(Resultados_Res, aes(x = resolution, y = ari)) +
    geom_point(col="#5F90BB") +
    geom_line(col="#3F6E9A") +
    #geom_errorbar(aes(ymin = ari_mean - ari_sd, ymax = ari_mean + ari_sd), width = 0.02) +
    theme_bw() +
    labs(title = "Adjusted Rand Index (ARI) by Resolution",
         x = "Resolution",
         y = "ARI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
    geom_vline(xintercept = 0.7, linetype = "dashed", size = 1,col="#E84746")

ggsave(plot = ARI_Plot, paste0(path.guardar, "/ARI_Evaluation_NotSumbSample.png"), height = 5, width = 10)
