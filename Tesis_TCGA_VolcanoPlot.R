


library(tidyverse) #incluye ggplot2 para visualización de datos, y dplyr para manipulación de los datos.
library (RColorBrewer) #Para ampliar nuestra gama de colores
library(ggrepel) #Para lindas anotaciones

res_df <- read.csv("C:/Users/david/OneDrive - UCB/Tesis/Datos/2da Pre defensa/R/volcanoplot/sensibles/TCGA_V3_res_df.csv")

#Primero ponemos de menor a mayor en función al p-value
res_df <- res_df[order(res_df$pvalue), ]

#Tema
theme_set(theme_classic(base_size = 20)+
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), colour = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), colour = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

#Graficar el esqueleto principal
ggplot(data = res_df,aes(x=log2FoldChange, y= -log10(pvalue))) +
  geom_vline(xintercept = c(-0.6,0.6), col='gray', linetype='dashed') + #para crear treshold lines que me separe cuando considero up and down regulated
  geom_hline(yintercept = 1.3, col='gray', linetype='dashed') +
  geom_point()

#Para poder diferenciar por color los Up y Down regulated primero debemos clasificarlos creando una nueva columna
res_df$expresiondiferencial <- 'No' 
res_df$expresiondiferencial[res_df$log2FoldChange > 1 & res_df$pvalue < 0.05] <- 'UP'
res_df$expresiondiferencial[res_df$log2FoldChange < -1 & res_df$pvalue < 0.05] <- 'DOWN'

#Graficar el esqueleto principal costumizado
ggplot(data = res_df,aes(x=log2FoldChange, y= -log10(pvalue), colour = expresiondiferencial)) +
  geom_vline(xintercept = c(-1,1), col='gray', linetype='dashed') + #para crear treshold lines que me separe cuando considero up and down regulated
  geom_hline(yintercept = 1.3, col='gray', linetype='dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("royalblue","grey","darkred"),
                     labels = c("Genes Subexpresados", "No Significativo", "Genes Sobrexpresados")) +
  coord_cartesian(ylim = c(0,7), xlim = c(-5,5)) + #Ahora los ejes y los limites
  scale_x_continuous(breaks = seq(-10,10,2)) +
  labs(color='Expresión Diferencial', #Titulo leyenda
       x= expression("log"[2]*"FC"), y= expression("-log"[10]*"p-value")) + #Ponemos el nombre de los ejes
  ggtitle('Expresión Diferencial entre Pacientes Respondedores y No Respondedores')  #Ponemos el titulo

#Añadir etiquetas para los genes significativos usando 'delabel'
#Seleccionamos los primeros 50 genes según el p-value 
top50seg <- head(res_df[order(res_df$pvalue, na.last = NA), "gene"], 50)

res_df$delabel <- ifelse(res_df$gene %in% top50seg, res_df$gene, NA)

# Modificar el volcánoplot para eliminar el título de la gráfica y el título de la leyenda

volcanoplot <- ggplot(data = res_df, aes(x=log2FoldChange, y= -log10(pvalue), colour = expresiondiferencial, label= delabel)) +
  geom_vline(xintercept = c(-1,1), col='gray', linetype='dashed') + #para crear treshold lines que me separe cuando considero up and down regulated
  geom_hline(yintercept = 1.3, col = 'gray', linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("royalblue","grey","darkred"),
                     labels = c("Genes Subexpresados", "Genes No Significativos", "Genes Sobrexpresados")) +
  coord_cartesian(ylim = c(0,15), xlim = c(-10,10)) + #Ahora los ejes y los limites
  scale_x_continuous(breaks = seq(-10,10,2)) +
  labs(color=NULL,  # Eliminar el título de la leyenda
       x= expression("log"[2]*"FC"), y= expression("-log"[10]*"p-value")) +  # Mantener nombres de los ejes
  ggtitle(NULL)  # Eliminar el título principal

# Agregar etiquetas usando geom_text_repel como antes
volcano_plot <- volcanoplot +
  geom_text_repel(
    data = subset(res_df, !is.na(delabel)),
    aes(label = delabel),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,  # O un valor mayor si es necesario
    force = 1,           # Ajusta este valor para cambiar la repulsión
    box.padding = 0.5,   # Espacio alrededor de las etiquetas
    point.padding = 0.3, # Espacio entre etiquetas y puntos
    min.segment.length = 0.1  # Controla la longitud mínima de las líneas de segmentación
  )

dev.new()
print(volcano_plot)



library(dplyr)

# Filtrar los genes que tienen un pvalue menor a 0.05
genes_significativos <- res_df_filtrado %>% filter(pvalue < 0.05)

# Contar el número de genes que cumplen con la condición
n_genes_significativos <- nrow(genes_significativos)


# Determinar cuántos de estos genes son sobreexpresados (UP) y cuántos son subexpresados (DOWN)
genes_sobreexpresados_total <- genes_significativos %>% filter(expresiondiferencial == "UP")
genes_subexpresados_total <- genes_significativos %>% filter(expresiondiferencial == "DOWN")

# Contar el número de genes que cumplen con cada condición
n_sobreexpresados <- nrow(genes_sobreexpresados_total)
n_subexpresados <- nrow(genes_subexpresados_total)

# Mostrar los resultados
cat("Número de genes con p-value < 0.05 que son sobreexpresados (UP):", n_sobreexpresados, "\n")
cat("Número de genes con p-value < 0.05 que son subexpresados (DOWN):", n_subexpresados, "\n")


# Seleccionar los primeros 50 genes según el p-value
top50_genes <- head(res_df[order(res_df$pvalue, na.last = NA), ], 50)

# Determinar cuántos y cuáles son los genes sobreexpresados (UP) y subexpresados (DOWN)
genes_sobreexpresados <- top50_genes %>% filter(expresiondiferencial == "UP")
genes_subexpresados <- top50_genes %>% filter(expresiondiferencial == "DOWN")

# Contar cuántos genes son sobreexpresados y cuántos son subexpresados
n_sobreexpresados <- nrow(genes_sobreexpresados)
n_subexpresados <- nrow(genes_subexpresados)

# Mostrar el resultado
cat("Número de genes sobreexpresados (UP):", n_sobreexpresados, "\n")
cat("Genes sobreexpresados:\n")
print(genes_sobreexpresados$gene)

cat("\nNúmero de genes subexpresados (DOWN):", n_subexpresados, "\n")
cat("Genes subexpresados:\n")
print(genes_subexpresados$gene)

