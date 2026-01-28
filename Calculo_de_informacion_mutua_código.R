# Mutual infromation calculations models Ranier

require(vroom)
library(tidyverse)
library(MASS)         # para usar binomiales negativas si se requiere
library(ggplot2)      # para gráficas lindas
library(infotheo)     # para los cálculos de información mutua
library(igraph)       # para manejar las redes
library(plot.matrix)  # para visualizar matrices
library(dplyr)        # Para manejar con más facilidad queries
library(maigesPack)   # Para calcular significancia redes por Boostrap

setwd("Archivos_CSV")

# Leemos la tabla de datos 

nac_all=read.table(file="/Archivos_CSV/Rwd_all.csv", header = FALSE, sep=",")

# Discretizamos en bines para el cálculo de información mutua

disc_nac <-  apply(X = nac_all, MARGIN = 1, FUN = infotheo::discretize)

# Calculamos la matríz de valores de información mutua
mi <- 
  sapply(disc_nac, FUN = function(i){
    sapply(disc_nac, FUN = function(j){
      infotheo::mutinformation(i,j)
    })
  })

# Analizamos la distribución de valores de información mutua

#density(mi)

denmi <-density(mi)

plot(denmi)

max(mi)

min(mi)

hist(mi, breaks=1000)

heatmap(mi,col=rainbow(7), main = "Mutual information distribution Rwd_all")

heatmap(mi,Colv = NA, Rowv = NA,col=rainbow(7), main = "Mutual information distribution Rwd_all no dendrogram")

# De este histograma notamos muchos valores de MI, relativamente bajos


# Esta es la matríz completa de valores de infomación mutua.
# Son 1527 x 1527 = 2,331,729. 

# Quizá valga la pena convertir la información a bits que 
# es una unidad más intuitiva

emei= natstobits(mi)

min(emei)

max(emei)

hist(emei, breaks=1000,main="Mutual Information Probability Distribution Rwd_all", xlab="Mutual Information (bits)")

denemei <-density(emei)

plot(denemei,main="Mutual Information Probability Density Distribution Rwd_all", xlab="Mutual Information (bits)")

write.csv(emei, row.names=F,file="MI_Full_Rwd_all.csv")


# Si vemos el histograma, hay un pico "extraño" con MI como de 3 bits
# Es que en el cálculo también tenemos a cada neurona consigo misma

# Para construir una matríz de adyacencia a partir de esta 
# matríz de información mutua es necesario establecer un umbral
# Exploremos esto. Por ejemplo si nos quedamos con el 0.1% más
# alto de los valores (nos quedarán 2,332 links).

th1m <- quantile(emei,0.999)

th1m


# Extraigo una matríz de adyacencia con los 2,332 enlaces mayores de MI

netmi <- ifelse(emei < quantile(emei, 0.999), 0, 1)

dim(netmi)

# emei[5,5]
# 
# emei[1,1]
# 
# emei[1276,1187]
# 
# netmi[4,4]
# 
# netmi[1276,1]

dim(netmi)

hist(netmi)

sum(netmi==1)

sum(netmi==0)

density(netmi)

heatmap(netmi,Colv = NA, Rowv = NA,col=c("white","black"), main = "Adjacency matrix")


# Como vemos la mayoría de las entradas de la matríz de adyacencia son cero 
# (de hecho, el 99.9 % de ellos).

denmitt <-density(netmi)

plot(denmitt)

g <- igraph::graph_from_adjacency_matrix(netmi)

edge_density(g, loops=F)

ecount(g)
vcount(g)

#plot(g)


write_graph(g, "Rwd_all_network_top_999.txt", "edgelist")

## Cálculo de significancia estadistica por bootstrap
## Es un cálculo muy largo que vale la pena correr en un cluster

# bootstrapMI(nac_all, y=NULL, bRep=100, ret="p-value")

## Para extraer redes con diferentes niveles de significancia

emeix=read.csv(file="MI_Full_Rwd_All.csv", header = TRUE)

##  emeix (em) y emei son exactamente la misma matríz, esto lo hice
## para trabajar en los umbrales varios días sin mantener en la memoria
## RAM de mi laptop una matríz tan grande

em=as.matrix.data.frame(emeix)

dim(em)

em[1,1]

## para generar las matrices con los diferentes umbrales variamos 
# los cuantiles 

# Para el top 99%

netmix <- ifelse(em < quantile(emei, 0.99), 0, 1)

gx <- igraph::graph_from_adjacency_matrix(netmix)

ecount(gx)
vcount(gx)

write_graph(gx, "Rwd_all_network_top_99.txt", "edgelist")

# Para el top 95%

netmiy <- ifelse(em < quantile(emei, 0.95), 0, 1)

gy <- igraph::graph_from_adjacency_matrix(netmiy)

ecount(gy)
vcount(gy)

write_graph(gy, "Rwd_all_network_top_95.txt", "edgelist")

# Para el top 90%

netmiz <- ifelse(em < quantile(emei, 0.9), 0, 1)

gz <- igraph::graph_from_adjacency_matrix(netmiz)

ecount(gz)
vcount(gz)

write_graph(gz, "Rwd_all_network_top_90.txt", "edgelist")

# Para el top 80%

netmiw <- ifelse(em < quantile(emei, 0.8), 0, 1)

gw <- igraph::graph_from_adjacency_matrix(netmiw)

ecount(gw)
vcount(gw)

write_graph(gw, "Rwd_all_network_top_80.txt", "edgelist")

# Para el top 75%

netmiv <- ifelse(em < quantile(emei, 0.75), 0, 1)

gv <- igraph::graph_from_adjacency_matrix(netmiv)

ecount(gv)
vcount(gv)

write_graph(gv, "Rwd_all_network_top_75.txt", "edgelist")

# Valores de corte de información mutua para los diferentes cuantiles (bits)

quantile(emei,0.99)

quantile(emei,0.95)

quantile(emei,0.90)

quantile(emei,0.80)

quantile(emei,0.75)
