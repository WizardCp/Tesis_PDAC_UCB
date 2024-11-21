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

#----------------------------------------------------
# Ver cuántos genes  tienen un padj < 0.05
genes_pvalue_significativos <- sum(res_df_filtrado$pvalue < 0.05, na.rm = TRUE)
genes_padj_significativos <- sum(res_df_filtrado$padj < 0.05, na.rm = TRUE)

# Imprimir el número de genes significativamente expresados
cat("Número de genes con padj < 0.05:", genes_padj_significativos, "\n")

# Calcular la proporción de genes significativamente expresados
total_genes <- nrow(res_df_filtrado)
cat("Proporción de genes con padj < 0.05:", round(genes_padj_significativos / total_genes * 100, 2), "%\n")

# Imprimir el número de genes significativamente expresados según el p-value no ajustado
cat("Número de genes con p-value < 0.05:", genes_pvalue_significativos, "\n")

# Calcular la proporción de genes con p-value < 0.05
proporcion_pvalue <- round(genes_pvalue_significativos / total_genes * 100, 2)
cat("Proporción de genes con p-value < 0.05:", proporcion_pvalue, "%\n")

#-------------------------------------------------

# Filtrar los genes con pvalue < 0.05
res_df_filtrado_significativos <- res_df_filtrado %>%
  filter(pvalue < 0.05)

# Identificar los 150 genes más sobreexpresados en los pacientes sensibles (log2FoldChange > 0)
genes_sensibles <- res_df_filtrado_significativos %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(150)

# Identificar los 150 genes más sobreexpresados en los pacientes resistentes (log2FoldChange < 0)
genes_resistentes <- res_df_filtrado_significativos %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  head(150)

# Mostrar un resumen de los nuevos data frames
cat("Genes sobreexpresados en pacientes sensibles:\n")
print(genes_sensibles)

cat("\nGenes sobreexpresados en pacientes resistentes:\n")
print(genes_resistentes)


# Cargar la librería biomaRt
library(biomaRt)

# Conectarte a Ensembl usando biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extraer las listas de genes de los data frames genes_sensibles y genes_resistentes
genes_sensibles_hugo <- unique(genes_sensibles$gene)
genes_resistentes_hugo <- unique(genes_resistentes$gene)

# Obtener información sobre procesos biológicos relacionados a los genes sobreexpresados en pacientes sensibles
anotaciones_sensibles <- getBM(
  attributes = c("hgnc_symbol", "go_id", "name_1006"),
  filters = "hgnc_symbol",
  values = genes_sensibles_hugo,
  mart = ensembl
)

# Obtener información sobre procesos biológicos relacionados a los genes sobreexpresados en pacientes resistentes
anotaciones_resistentes <- getBM(
  attributes = c("hgnc_symbol", "go_id", "name_1006"),
  filters = "hgnc_symbol",
  values = genes_resistentes_hugo,
  mart = ensembl
)

# Ver las anotaciones de procesos biológicos para genes sobreexpresados en pacientes sensibles
cat("Procesos biológicos para genes sobreexpresados en pacientes sensibles:\n")
print(anotaciones_sensibles)

# Ver las anotaciones de procesos biológicos para genes sobreexpresados en pacientes resistentes
cat("\nProcesos biológicos para genes sobreexpresados en pacientes resistentes:\n")
print(anotaciones_resistentes)

