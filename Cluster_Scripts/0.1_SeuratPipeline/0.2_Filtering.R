################################################################################

# SCRIPT: Filtering Seurat Object 
# AUTHOR: SOFIA PETRISSANS MOLL
# DATE: 06-06-2025

################################################################################

# LIBRARIES ....................................................................

source('/ijc/USERS/spetrissans/ANNA_MARIA/3.MERGE/0.1_Scripts/0.1_SeuratPipeline/0.0_Paths.R')
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(writexl)


# CONFIGURACIÓN ................................................................

# JOB ARRAY
args <- commandArgs(trailingOnly = TRUE)
# Definir los thresholds que quieres probar
thresholds <- c(5, 10, 15, 20, 25, 30, 35, 40)

if (length(args) == 0) {
  stop("Debes proporcionar un índice de SLURM array como argumento.")
}
array_index <- as.integer(args[1])

if (array_index < 1 || array_index > length(thresholds)) {
  stop("Índice fuera del rango de thresholds definidos.")
}

mito_thresh <- thresholds[array_index]
cat(paste0("Ejecutando filtrado con threshold: ", mito_thresh, "%\n"))




# Crear carpeta específica para cada threshold
path.guardar_original <- paths_AnnaMaria$Path_guardar
folder_name <- paste0("mt", mito_thresh)
path.guardar <- file.path(path.guardar_original, "2.Filtering", folder_name)
dir.create(path.guardar, showWarnings = FALSE, recursive = TRUE)
# Crear una carpeta para los plots de QC
path.guardar.plots <- file.path(path.guardar, "plots")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)
dir.create(path.guardar.plots, recursive = TRUE, showWarnings = FALSE)

################################################################################

# ..............................................................................
# CARGAR OBJETO SEURAT
# ..............................................................................
main_path_objeto <- paths_AnnaMaria$Path_seurat_object
path_objeto <- file.path(main_path_objeto, "1.QC")
seurat_object <- readRDS(file.path(path_objeto, "RawObject_Total_QC.rds"))

getPalette <- colorRampPalette(brewer.pal(9, "Set3"))
colors <- getPalette(length(unique(seurat_object$Sample)))
# ..............................................................................
# FILTRADO 
# ..............................................................................


cat(paste0("Filtrando con threshold percent.mt < ", mito_thresh, "\n"))

seurat_filtered <- subset(seurat_object,
                          subset = (nFeature_RNA >= 250) &
                                   (percent.mt < mito_thresh) &
                                   (log10GenesPerUMI > 0.8))


# ..............................................................................
# GUARDAR OBJETO FILTRADO
# ..............................................................................

saveRDS(seurat_filtered, file.path(path.guardar, "QC_Filtered.rds"))


# ..............................................................................
# SUMMARY
# ..............................................................................


results_filtered <- data.frame(
  Parameters_Filtering = paste0("(nFeature_RNA>=250) & (percent.mt<", mito_thresh, ") & (log10GenesPerUMI>0.8)"),
  Previous_Number_Cells = ncol(seurat_object),
  Previous_Number_Genes = nrow(seurat_object),
  NumberCells_AfterFiltering = ncol(seurat_filtered),
  NumberGenes_AfterFiltering = nrow(seurat_filtered),
  CellsDiscarded = ncol(seurat_object) - ncol(seurat_filtered),
  MeanGenes_AfterFiltering = mean(seurat_filtered$nFeature_RNA)
)
write_xlsx(results_filtered, file.path(path.guardar, "Results_Filtered.xlsx"))


# Media de genes por sample
media_genes_por_paciente <- seurat_filtered@meta.data %>%
  group_by(Sample) %>%
  summarise(MediaGenes = mean(nFeature_RNA))
write_xlsx(media_genes_por_paciente, file.path(path.guardar, "MeanGenesPerPatient.xlsx"))

cat("Filtrado finalizado. Objeto guardado en:", file.path(path.guardar, "QC_Filtered.rds"), "\n")


# ..............................................................................
# PLTOS QUALITY CONTROL 
# ..............................................................................

# Library Size Histogram .......................................................
print("Library size histograma")
library_size_plot <- seurat_filtered@meta.data %>%
  ggplot(aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  labs(x = "Total UMI per Cell", y = "Number of Cells") +
  theme_light(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Library_size_hist.png",sep="/"), plot = library_size_plot, width = 25, height = 10)


# Library Size by Sample .......................................................
print("Library size by Sample")
library_size_sample <- seurat_filtered@meta.data %>%
  ggplot(aes(x = Sample, y = nCount_RNA, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.4, color = NA) +
  geom_jitter(aes(color = Sample), width = 0.25, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_continuous(trans = "log2") +
  labs(x = "", y = "Total UMI per Cell") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Library_size_by_Sample.png",sep="/"),plot = library_size_sample,width = 25,height = 10)


# Number of Genes ..............................................................
print("Number of genes")
number_genes <- seurat_filtered@meta.data %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  labs(x = "Number of Detected Genes per Cell", y = "Number of Cells") +
  theme_light(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Number_of_genes.png",sep="/"),plot = number_genes,width = 25,height = 10)


# Number of Genes by Sample ...................................................
print("Number of genes by sample")
number_genes_by_sample <- seurat_filtered@meta.data %>%
  ggplot(aes(x = Sample, y = nFeature_RNA, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.4, color = NA) +
  geom_jitter(aes(color = Sample), width = 0.25, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_continuous(trans = "log2") +
  labs(x = "", y = "Number of Genes per Cell") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Number_genes_by_Sample.png",sep="/"),plot = number_genes_by_sample,width = 25,height = 10)
 

# Complexity ....................................................................
print("Complexity")
Complexity <- seurat_filtered@meta.data %>% ggplot(aes(x=log10GenesPerUMI, color = Phenotype, fill=Phenotype)) +
  geom_density(alpha = 0.2) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Complexity.png",sep="/"),plot = Complexity,width = 25,height = 10)


# Mitocondrial VlnPlot .........................................................
print("Mitocondrial Violin Plot")
mitocondrial_VlnPlot <- seurat_filtered@meta.data %>%
  ggplot(aes(x = Sample, y = percent.mt, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.4, color = NA) +
  geom_jitter(aes(color = Sample), width = 0.25, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_continuous(trans = "log2", labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "", y = "Mitochondrial % per Cell") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Mitocondrial_ViolinPlot.png",sep="/"),plot = mitocondrial_VlnPlot,width = 25,height = 10)


#  Mitocondrial Density ......................................................
print("Mitocondrial density")
mitocondrial_density <- seurat_filtered@meta.data %>%
  ggplot(aes(color = Sample, x = percent.mt, fill = Sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = scales::label_number(accuracy = 0.1)) +
  geom_vline(xintercept = 20, linetype = "dashed") +
  labs(x = "Mitochondrial Content (%) [log10]", y = "Density") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Mitocondrial_density.png",sep="/"),plot = mitocondrial_density,width = 25,height = 10)


# Ribosomal VlnPlot ...........................................................
print("Ribosomal Violin Plot")
ribosomal_VlnPlot <- seurat_filtered@meta.data %>%
  ggplot(aes(x = Sample, y = percent.rb, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.4, color = NA) +
  geom_jitter(aes(color = Sample), width = 0.25, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_continuous(trans = "log2", labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "", y = "Ribosomal % per Cell") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = paste(path.guardar.plots,"Ribosomal_VlnPlot.png",sep="/"),plot = ribosomal_VlnPlot,width = 25,height = 10)


# Cell counts .................................................................
print("Cell counts per Sample")
cell_counts <- seurat_filtered@meta.data %>%
  group_by(Sample) %>%
  summarise(CellCount = n()) %>%
  ggplot(aes(x = Sample, y = CellCount, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = CellCount), vjust = -0.5, size = 4) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Number of Cells") +
  ggtitle("Cell Counts per Sample") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = file.path(path.guardar.plots, "Cell_counts_sample.png"), plot = cell_counts, width = 12, height = 6)


# UMICount ............................................................
print("UMI counts")
UMI_count <- seurat_filtered@meta.data %>%
  ggplot(aes(x = nCount_RNA, fill = Sample, color = Sample)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  geom_vline(xintercept = 500, linetype = "dashed", color = "black", linewidth = 0.7) +
  labs(
    x = "UMIs per cell (log10 scale)",
    y = "Cell density",
    title = "UMI Count Distribution"
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = file.path(path.guardar.plots, "UMI_count.png"), plot = UMI_count, width = 10, height = 6)


# GenesPerCell ........................................................
print("Genes per cell")
genes_per_cell <- seurat_filtered@meta.data %>%
  ggplot(aes(x = nFeature_RNA, fill = Sample, color = Sample)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  geom_vline(xintercept = 300, linetype = "dashed", color = "black", linewidth = 0.7) +
  labs(
    x = "Number of detected genes per cell (log10 scale)",
    y = "Cell density",
    title = "Gene Count Distribution"
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = file.path(path.guardar.plots, "Genes_cell.png"), plot = genes_per_cell, width = 10, height = 6)



qc_data <- data.frame(
  n_genes_by_counts = seurat_filtered@meta.data$nFeature_RNA,
  total_counts = seurat_filtered@meta.data$nCount_RNA,
  patient=seurat_filtered@meta.data$Sample
)

# Total count per cell .......................................................
print("Total counts per cell")
total_counts_per_cell <- ggplot(qc_data, aes(x = patient, y = total_counts, fill = patient)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, color = NA) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.5, color = "black") +
  scale_y_log10() +
  labs(
    title = "Distribution of UMIs per cell",
    x = "Sample",
    y = "Total UMIs (log10)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24)
  )
ggsave(filename = file.path(path.guardar.plots, "Total_count_cell.png"), plot = total_counts_per_cell, width = 15, height = 7)


# Total counts vs number og cells .......................................................
print("Total Count vs Number of Genes")
total_count_vs_number_genes <- ggplot(qc_data, aes(x = total_counts, y = n_genes_by_counts, color = patient)) +
  geom_point(alpha = 0.4, size = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "UMIs vs. Genes detected per cell",
    x = "Total UMIs (log10)",
    y = "Genes detected (log10)"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
ggsave(filename=paste(path.guardar.plots,"TotalCount_vs_NumberGenes.png",sep="/"), plot=total_count_vs_number_genes, width=10, height=10)


