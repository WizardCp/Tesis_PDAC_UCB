# Instalación de paquetes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Paquetes de Bioconductor
BiocManager::install(c("DESeq2", "biomaRt", "pheatmap", "ComplexHeatmap", "EnhancedVolcano"))

# Paquetes de CRAN
install.packages(c("tidyverse", "readr"))

# Cargar librerías
library(tidyverse)
library(readr)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)

# Definir el directorio de trabajo
setwd("C:/Users/davidcp/OneDrive - UCB/Tesis/Datos/2da Pre defensa/R")
# Cargar el archivo CSV en un DataFrame
conteos <- read.csv("Conteos_preparados.csv")

# Verificar la carga del archivo
head(conteos)

########## Agregamos los símbolos a los genes.
conteos$GENE_ID <- sub("\\.\\d+$", "", conteos$GENE_ID)  # Remover sufijos numéricos
conteos$GENE_ID <- gsub("_.*", "", conteos$GENE_ID)  # Remover sufijos no estándar

# Usa la columna GENE_ID ya limpia para crear
ensg_ids_clean <- conteos$GENE_ID

# Verificar si todavía hay identificadores con sufijos largos
genes_with_suffixes <- ensg_ids_clean[grepl("_", ensg_ids_clean)]
if (length(genes_with_suffixes) > 0) {
  cat("Hay genes con sufijos restantes:\n")
  print(genes_with_suffixes)
} else {
  cat("No hay genes con sufijos restantes.\n")
}

# Contar el número de duplicados en la columna GENE_ID
num_duplicados <- sum(duplicated(conteos$GENE_ID))

# Imprimir el número de duplicados
cat("Número de identificadores GENE_ID duplicados:", num_duplicados, "\n")


# Contar el número de genes únicos en el DataFrame
num_genes_unicos <- length(unique(conteos$GENE_ID))

# Función para conectar a Ensembl probando diferentes espejos
connect_to_ensembl <- function(mirrors) {
  for (mirror in mirrors) {
    message(paste("Intentando conectar al espejo:", mirror))
    ensembl <- tryCatch(
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = mirror),
      error = function(e) {
        message(paste("Error al conectar con el espejo:", mirror, "- Probando otro..."))
        return(NULL)
      }
    )
    if (!is.null(ensembl)) {
      message(paste("Conectado exitosamente al espejo:", mirror))
      return(ensembl)
    }
  }
  stop("No se pudo conectar a ningún servidor de Ensembl. Intenta nuevamente más tarde.")
}

# Espejos a probar
mirrors <- c("www", "useast", "asia")

# Conectar a Ensembl
ensembl <- connect_to_ensembl(mirrors)

# Realizar el mapeo de genes
gene_ids <- conteos$GENE_ID  # Usar la columna GENE_ID del DataFrame "conteos"

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_ids,  # Los identificadores ya deben estar limpios
  mart = ensembl
)

# Imprimir cuántos genes fueron mapeados
cat("Número de genes mapeados:", nrow(gene_mapping), "\n")

# Identificar cuántos no fueron mapeados
genes_no_mapeados <- setdiff(gene_ids, gene_mapping$ensembl_gene_id)
cat("Número de genes no mapeados:", length(genes_no_mapeados), "\n")

# Mostrar los genes no mapeados
cat("Lista de genes no mapeados:\n")
print(genes_no_mapeados)

# Unir el resultado de mapeo con el DataFrame original para agregar la columna de nombres de genes Hugo
conteos <- merge(conteos, gene_mapping, by.x = "GENE_ID", by.y = "ensembl_gene_id", all.x = TRUE)

# Mostrar un resumen del DataFrame actualizado
head(conteos)

# Verificar si hay duplicados en la columna de nombres de genes mapeados (hgnc_symbol)
genes_mapeados_duplicados <- conteos$hgnc_symbol[duplicated(conteos$hgnc_symbol) & !is.na(conteos$hgnc_symbol)]

# Imprimir los genes mapeados que están duplicados
if (length(genes_mapeados_duplicados) > 0) {
  cat("Genes mapeados duplicados:\n")
  print(unique(genes_mapeados_duplicados))
} else {
  cat("No hay genes mapeados duplicados.\n")
}

# Eliminar las filas donde el nombre del gen mapeado (hgnc_symbol) es NA
conteos <- conteos[!is.na(conteos$hgnc_symbol), ]

# Verificar el número de genes restantes después de eliminar los no mapeados
cat("Número de genes después de eliminar los no mapeados:", nrow(conteos), "\n")


# Crear un nuevo DataFrame con "hgnc_symbol" como la primera columna y sin la columna "GENE_ID"
conteos2 <- conteos[, !names(conteos) %in% "GENE_ID"]  # Eliminar la columna "GENE_ID"
conteos2 <- conteos2[, c("hgnc_symbol", setdiff(names(conteos2), "hgnc_symbol"))]  # Reordenar para que "hgnc_symbol" sea la primera columna

# Verificar que el nuevo DataFrame esté en el formato correcto
head(conteos2)

# Cargar la librería DESeq2
library(DESeq2)

# Crear una matriz de conteos sin la columna de nombres de genes (hgnc_symbol)
conteos_matriz <- as.matrix(conteos2[, -which(names(conteos2) == "hgnc_symbol")])
rownames(conteos_matriz) <- conteos2$hgnc_symbol  # Usar los nombres de genes como nombres de filas


# Crear un objeto de diseño para DESeq2 indicando los grupos (Sensible o Resistente)
# Ajusta este vector según tus datos (por ejemplo, las primeras 8 columnas como "Sensible" y el resto como "Resistente")
grupo <- c(rep("Sensible", 8), rep("Resistente", ncol(conteos_matriz) - 8))

is.factor(grupo) #Verificamos en que formato se encuentra

grupo <- as.factor(grupo) #Lo volvemos un factor para definir el grupo de referencia

grupo <- factor(grupo, levels = c("Sensible", "Resistente"))

levels(grupo) #El que aparezca primero es el grupo de referencia. 

# Crear un objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = conteos_matriz, colData = data.frame(grupo), design = ~ grupo)

# Realizar la normalización y análisis de DESeq2
dds <- DESeq(dds)

# Obtener los resultados normalizados
conteos_normalizados <- counts(dds, normalized = TRUE)

# Ver un resumen de los datos normalizados
head(conteos_normalizados)

# Convertir la matriz de conteos normalizados a un DataFrame
conteos_normalizados_df <- as.data.frame(conteos_normalizados)

# Ver un resumen del DataFrame normalizado
head(conteos_normalizados_df)

# Asegurarnos de que tenemos el objeto DESeq2 y los resultados
res <- results(dds)

# Convertir los resultados a un DataFrame si no lo has hecho aún
res_df <- as.data.frame(res)

# Asegurarte de que tienes una columna con los nombres de los genes
res_df$gene <- rownames(res_df)

# Filtrar los genes que tienen nombres identificados y eliminar aquellos que no fueron mapeados correctamente
res_df_filtrado <- res_df %>%
  filter(!grepl("^\\.", gene) & gene != "")  # Elimina genes que comienzan con un punto o que están vacíos

# Verificar el nuevo DataFrame filtrado
head(res_df_filtrado)

# Cargar la librería ComplexHeatmap
library(ComplexHeatmap)

# Filtrar los 50 genes más significativamente expresados usando pvalue < 0.05 y ordenando por pvalue
genes_significativos <- res_df_filtrado %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>%
  slice_head(n = 50)  # Seleccionar los primeros 30 genes

# Extraer los nombres de los genes significativos
genes_interes <- genes_significativos$gene  # Asegúrate de que esta columna contiene los nombres de genes

# Filtrar la matriz de datos normalizados para los genes significativos
conteos_heatmap <- conteos_normalizados_df[rownames(conteos_normalizados_df) %in% genes_interes, ]

# Normalizar los datos con escala Z si es necesario
conteos_heatmap_scaled <- t(scale(t(conteos_heatmap)))

# Crear un vector que identifique los grupos de las muestras
grupo_colores <- c(rep("Sensible", 8), rep("Resistente", 15))  # Ajusta estos números según tus datos

# Asignar colores a los grupos
colores_grupo <- setNames(c("orange", "yellow"), c("Resistente", "Sensible"))

# Crear una anotación para los grupos con una cuadrícula
anotacion_colores <- HeatmapAnnotation(
  Grupo = grupo_colores,
  col = list(Grupo = colores_grupo),
  show_legend = TRUE,  # Mostrar la leyenda de la anotación
  annotation_legend_param = list(title = "Grupo"),  # Agregar el título de la leyenda
  show_annotation_name = FALSE,  # Ocultar el nombre de la fila "Grupo"
  gp = gpar(col = "black", lwd = 0.5)  # Agregar la cuadrícula a la anotación
)

# Crear un esquema de colores más suaves para el heatmap
color_palette <- colorRampPalette(c("navy", "white", "firebrick"))(100)

# Crear el heatmap con ComplexHeatmap sin el cladograma del eje Y
Heatmap(
  conteos_heatmap_scaled,
  name = "Escala Z",  # Título de la escala de colores
  col = color_palette,  # Aplicar la paleta de colores suaves
  top_annotation = anotacion_colores,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,  # Desactivar agrupación de las filas para eliminar el cladograma
  cluster_columns = TRUE,
  row_names_side = "right",  # Mostrar los nombres de los genes en la derecha
  column_names_side = "bottom",
  rect_gp = gpar(col = "black", lwd = 0.5),  # Agregar cuadriculado para divisiones de filas y columnas
  heatmap_legend_param = list(title = "Escala Z")  # Agregar título a la leyenda de la escala de colores
)

library(ggplot2)
library(dplyr)
library(ggrepel)

# Crear columna de expresiondiferencial si no existe
res_df_filtrado$expresiondiferencial <- "No"
res_df_filtrado$expresiondiferencial[res_df_filtrado$log2FoldChange > 1 & res_df_filtrado$pvalue < 0.05] <- "UP"
res_df_filtrado$expresiondiferencial[res_df_filtrado$log2FoldChange < -1 & res_df_filtrado$pvalue < 0.05] <- "DOWN"


genes_interes <- c("SOX30", "CENPA", "PLD4", "PRSS12", "AXIN2", "NKX3-1", "CYSLTR1")

# Filtrar el dataframe res_df para mantener solo los genes de interés
genes_seleccionados <- res_df_filtrado %>% filter(gene %in% genes_interes)

# Crear el dotplot detallado sin título de la gráfica ni título de la leyenda
ggplot(data = genes_seleccionados, aes(x = gene, y = log2FoldChange, color = expresiondiferencial, size = -log10(pvalue))) +
  geom_point(alpha = 0.7) +  # Tamaño de los puntos y transparencia para mejor visualización
  geom_text_repel(aes(label = gene), size = 4, box.padding = 0.5, point.padding = 0.3) + # Añadir etiquetas a los genes sin que se superpongan
  scale_color_manual(values = c("UP" = "darkred", "No" = "grey", "DOWN" = "royalblue"),
                     labels = c("Sobrexpresado", "No Significativo", "Subexpresado")) +
  scale_size_continuous(range = c(3, 6), name = expression("-log"[10]*" p-value")) + # Cambiar el rango de los tamaños de los puntos para reflejar la significancia
  theme_classic(base_size = 15) +
  labs(title = NULL,  # Eliminar el título de la gráfica
       x = "Genes",
       y = expression("log"[2]*" Fold Change"),
       color = NULL) +  # Eliminar el título de la leyenda de expresión diferencial
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") + # Mover la leyenda a la derecha
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Línea horizontal en y = 0 para representar la no expresión diferencial
  geom_vline(xintercept = 1:length(genes_interes), linetype = "dotted", color = "grey") # Líneas verticales para delimitar cada gene


# Crear el dotplot detallado
ggplot(data = genes_seleccionados, aes(x = gene, y = log2FoldChange, color = expresiondiferencial, size = -log10(pvalue))) +
  geom_point(alpha = 0.7) +  # Tamaño de los puntos y transparencia para mejor visualización
  geom_text_repel(aes(label = gene), size = 4, box.padding = 0.5, point.padding = 0.3) + # Añadir etiquetas a los genes sin que se superpongan
  scale_color_manual(values = c("UP" = "darkred", "No" = "grey", "DOWN" = "royalblue"),
                     labels = c("Sobrexpresado", "No Significativo", "Subexpresado")) +
  scale_size_continuous(range = c(3, 6), name = expression("-log"[10]*" p-value")) + # Cambiar el rango de los tamaños de los puntos para reflejar la significancia
  theme_classic(base_size = 15) +
  labs(title = "Dotplot Detallado de Genes Seleccionados",
       x = "Genes",
       y = expression("log"[2]*" Fold Change"),
       color = "Expresión Diferencial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") + # Mover la leyenda a la derecha
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Línea horizontal en y = 0 para representar la no expresión diferencial
  geom_vline(xintercept = 1:length(genes_interes), linetype = "dotted", color = "grey") # Líneas verticales para delimitar cada gene
