# TFM_Claudia_Arce
Este repositorio contiene el código empleado en el TFM "Análisis integrativo de datos ómicos para el desarrollo de metabotipos" realizado por Claudia Arce Cuesta


## 0. Obtención de los datos
Antes de importar los datos se vio que había algunos pacientes de los que faltaba algún tipo de dato o variable, por lo que dichos pacientes han sido descartados del dataset previamente antes de importarse a R.

Para poder trabajar con el argoritmo SNF es importante que todas las variables sean cuantitativas por lo que los datos de variables cualitativas de los que se disponía en el dataset también han sido descartados antes de empezar a trabajar con ellos.

Instalación de los paquetes que se van a a utilizar:

```{r paquetes, echo=TRUE}
library(SNFtool)
```

```{r datos}
#importamos otus de la microbiota intestinal
otus<-data.frame(read.table("C:/Users/ddcla/Documents/Master/TFM/datos/bact_1.txt",  sep='\t', header=T, row.names = 1))

#Se importan los metabolitos de los pacientes
metab<-data.frame(read.table("C:/Users/ddcla/Documents/Master/TFM/datos/metabs_1.txt",  sep='\t', dec=',',header=T, row.names = 1))

#Se importan los datos de los pacientes
samples<-data.frame(read.table("C:/Users/ddcla/Documents/Master/TFM/datos/samples_1.txt", sep = "\t", dec=',', header=T, row.names = 1))
```

## 1. Preparación de los datos

Se escalan las variables antes de trabajar con ellas. Las variables deben ser las columnas y los pacientes las filas.

```{r snf}
set.seed(1234)
metab_sc <- scale(metab)
samples_sc<-scale(samples)
otus_sc<-scale(otus)
```

## 2. Distancias y affinitymatrix

Se calculan distancias en los tres grupos ded atos. La distancia Euclidea es la mas habitual pero pueden usarse otras distancias. Después se transforman las matrices de distancia en matrices de afinidad o similaridad. k es numero de vecinos y suele estar entre 10 y 30 Sigma es el kernel de similaridad, se utilizará en un valor entre 0.3 y 0.8

```{r dist}
#Se crean las matrices de distancias
dist_metab <- as.matrix(dist(metab_sc))
dist_samples <- as.matrix(dist(samples_sc))
dist_otus <- as.matrix(dist(otus_sc))

#Creamos las affinityMatrix
k<-20
sigma<-0.8
Wmetab <- affinityMatrix(dist_metab, K = k, sigma = sigma)
Wsamples <- affinityMatrix(dist_samples, K = k, sigma = sigma)
Wotus <- affinityMatrix(dist_otus, K = k, sigma = sigma)
```

## 3. Agrupamientos en los datos

Antes de fusionar las redes veremos si hay algún tipo de agrupamiento en estas y lo mostramos mediante heatmaps.

```{r clustering y heatmaps1}
#Se estima el numero de clusters óptimo
estimateNumberOfClustersGivenGraph(Wmetab) # el numero de clusters para los metabolitos será 2 
estimateNumberOfClustersGivenGraph(Wsamples) # el numero de clusters para samples será 3
estimateNumberOfClustersGivenGraph(Wotus) # el numero de clusters para la microbiota será 2

#Se crean los clusters
clus_metab <-spectralClustering(Wmetab, K = 2)
clus_samples <-spectralClustering(Wsamples, K = 3)
clus_otus <-spectralClustering(Wotus, K = 2)

#Se muestran los heatmaps de los clusters
displayClustersWithHeatmap(Wmetab,clus_metab )
displayClustersWithHeatmap(Wsamples, clus_samples)
displayClustersWithHeatmap(Wotus,clus_otus)
```
![image](https://user-images.githubusercontent.com/104385965/212497443-eb7a40ca-ac3f-4e41-8188-7742492a9bff.png)

![image](https://user-images.githubusercontent.com/104385965/212497442-254ef5a0-755b-407f-824d-b260c2449024.png)

![image](https://user-images.githubusercontent.com/104385965/212497438-e9eefb38-cb3a-447b-bd61-b07c76888b12.png)

## 4. Fusión de matrices de similitud

Juntamos las tres matrices de similitud en una sola mas compleja mediante Similarity Network Fusion. Se realiza clustering y se muestra mediante un heatmap.

```{r network}
#Se unen las tres matrices previamente
matrices<- list(Wmetab,Wsamples,Wotus)
network = SNF(matrices, K= 20, t=20)

#Realizamos el clustering de la red
estimateNumberOfClustersGivenGraph(network) 
#el número de clusters optimo será 3

#Se realiza clustering y se muestra heatmap
clusteringSNF<-spectralClustering(network, K = 3)
displayClustersWithHeatmap(network,clusteringSNF)

```
![image](https://user-images.githubusercontent.com/104385965/212497428-0a195c7d-0c5f-4237-9307-70d37d62251f.png)

## 5. Comprobar concordancia

Con la función concordanceNetworkkNMI podemos ver la concordancia entre la red creada y las matrices iniciales.

```{r concordance}
ConcordanceMatrix = concordanceNetworkNMI(list(network,Wmetab,Wsamples,Wotus),3)
ConcordanceMatrix

```
![Captura de pantalla (70)](https://user-images.githubusercontent.com/104385965/212497905-15c09b06-8a8e-4b82-91a4-8837b9fa260c.png)
