# cargar paquetes necesarios
library(SummarizedExperiment)
library(readr)
library(tibble)
library(psych)
library(dplyr)
library(ggplot2)
library(knitr)
library(pheatmap)

# cargar dataset metabolomica
df <- read_csv("C:/Users/soliz-j/Downloads/Master/human_cachexia.csv")

# implementar una seed para crear el objeto 
set.seed(1234)

# obtener información del dataset
assay_matrix <- as.matrix(df[, -(1:2)])  # eliminar columnas ID y grupo
rownames(assay_matrix) <- df$`Patient ID`  # cada fila = muestra
assay_matrix <- t(assay_matrix)

# crear colData del dataset
col_data <- DataFrame(
  PatientID = df$`Patient ID`,
  Group = factor(df$`Muscle loss`)  # cachexic / control
)
rownames(col_data) <- df$`Patient ID`

# crear rowData para los metabolitos
row_data <- DataFrame(
  Metabolite = rownames(assay_matrix)
)
rownames(row_data) <- row_data$Metabolite

# cargar descripción del dataset desde description.md
description_lines <- readLines("C:/Users/soliz-j/Downloads/Master/description.md")
description_text <- paste(description_lines, collapse = "\n")

# crear objeto de clase SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(metabolites = assay_matrix),
  colData = col_data,
  rowData = row_data,
  metadata = list(description = description_text)
)

# visualizar resumen del objeto
se
save(se, file = "C:/Users/soliz-j/Downloads/Master/mi_objeto.Rda")

# visualizar metadata de las muestras
colData(se)

# visualizar metabolitos 
rowData(se)

# visualizar descripciones generales del experimento
metadata(se)

# preparacion de datos para análisis exploratorio
matriz <- assay(se)
grupos <- colData(se)$Group
names(grupos) <- rownames(colData(se))

df_long <- as.data.frame(t(matriz))
df_long$Sample <- rownames(df_long)
df_long$Group <- grupos[rownames(df_long)]

df_long <- df_long |>
  pivot_longer(cols = -c(Sample, Group),
               names_to = "Metabolite",
               values_to = "Value")


# crear tabla resumen de los datos por grupo y metabolito
resultados <- lapply(unique(df_long$Metabolite), function(met) {
  datos_met <- subset(df_long, Metabolite == met)
  res <- describeBy(datos_met$Value, group = datos_met$Group, mat = TRUE)
  res$Metabolite <- met
  res$Group <- res$group1
  res$ID <- paste0(res$Metabolite, "-", res$Group)
  res
})
resultados_df <- bind_rows(resultados)
head(resultados_df)

# análisis evaluativo de los datos

# lista de metabolitos
metabs <- unique(df_long$Metabolite)

# generar histograma para cada metabolito
for (m in metabs) {
  print(
    ggplot(
      dplyr::filter(df_long, Metabolite == m),
      aes(x = Value, fill = Group)
    ) +
      geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
      theme_minimal() +
      labs(title = paste("Histograma de", m)) +
      theme(legend.position = "bottom")
  )
}

# generar bloxplot por cada metabolito

for (m in unique(df_long$Metabolite)) {
  print(
    ggplot(filter(df_long, Metabolite == m),
           aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot() +
      theme_minimal() +
      labs(title = paste("Boxplot de", m)) +
      theme(legend.position = "none")
  )
}


# PCA
# escalar y transponer matriz para generar el PCA
scaled_expr <- scale(t(matriz))
pca <- prcomp(scaled_expr)

# agrupar datos
pca_df <- as.data.frame(pca$x)
pca_df$Group <- colData(se)$Group

# generar PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA de metabolomica de control vs cachexia")

# HEATMAP
# calcular coeficiente de variación
cv <- apply(assay(se), 1, function(x) sd(x) / mean(x))
top_metabs <- names(sort(cv, decreasing = TRUE))[1:25]

# generar heatmap
pheatmap(assay(se)[top_metabs, ],
         annotation_col = as.data.frame(colData(se)),
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Top 25 metabolitos más variables",
         fontsize_row = 6,   # Tamaño de letra de las filas
         fontsize_col = 2     # Tamaño de letra de las columnas
)


