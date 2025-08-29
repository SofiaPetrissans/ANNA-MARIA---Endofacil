
seurat_object <- readRDS("/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.2_Annotation_All_Populations/3.Annotation_Layer/mt15/Seurat_annotated.rds")
table(seurat_object$Annotation_Layer1)


path.guardar <- "/ijc/LABS/GRAUPERA/RAW/SPETRISSANS/ANNA_MARIA/3.MERGE/0.2_Annotation_All_Populations/3.Annotation_Layer/mt15"

# Calcular proporciones
df1 <- as.data.frame(table(seurat_object$Annotation_Layer1))
colnames(df1) <- c("CellType", "Count")
df1 <- df1 %>% mutate(Percentage = Count / sum(Count) * 100)

# Ordenar por porcentaje (opcional)
df1$CellType <- factor(df1$CellType, levels = df1$CellType[order(df1$Percentage, decreasing = TRUE)])

# Plot horizontal (stacked)
p1 <- ggplot(df1, aes(x = 1, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = Annotation_Colors_L1) +
  labs(x = NULL, y = "Percentage", fill = "Cell Type") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(file.path(path.guardar, "Barplot_Stacked_L1.png"), plot = p1, width = 8, height = 2)

# Calcular porcentajes por Phenotype
df2 <- seurat_object@meta.data %>%
  group_by(Annotation_Layer1, Phenotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Annotation_Layer1) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Ordenar factores según total de células
celltype_order <- df2 %>%
  group_by(Annotation_Layer1) %>%
  summarise(total = sum(Count)) %>%
  arrange(desc(total)) %>%
  pull(Annotation_Layer1)

df2$Annotation_Layer1 <- factor(df2$Annotation_Layer1, levels = celltype_order)

# Plot horizontal
p2 <- ggplot(df2, aes(x = Annotation_Layer1, y = Percentage, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  scale_fill_manual(values = c("vehicle" = "#D9D9D9", "CHOP" = "#8DD3C7")) +
  labs(x = "Cell Type", y = "Percentage", fill = "Phenotype") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(path.guardar, "Barplot_L1_by_Phenotype.png"), plot = p2, width = 7, height = 5)
