# Clase Introducción a R

- Instructora: Dámaris Esquén
- Fecha: 12/02/2025

## Preliminar: Instlación de RStudio

### Descargar el instalador de Rstudio:

```BASH
cd Bioprograms
wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2024.12.0-467-amd64.deb
```

### Instalar RStudio mediante archivo .deb

```BASH
sudo dpkg -i rstudio-2024.12.0-467-amd64.deb
```

Si ocurriera un problema como este: 

```
dpkg: dependency problems prevent configuration of rstudio:
 rstudio depends on libssl-dev; however:
  Package libssl-dev is not installed.
 rstudio depends on libclang-dev; however:
  Package libclang-dev is not installed. <-------------------ESTA ES LA DEPENDENCIA QUE FALTA.

dpkg: error processing package rstudio (--install):
 dependency problems - leaving unconfigured
Processing triggers for mailcap (3.70+nmu1ubuntu1) ...
Processing triggers for gnome-menus (3.36.0-1ubuntu3) ...
Processing triggers for desktop-file-utils (0.26-1ubuntu3) ...
Processing triggers for hicolor-icon-theme (0.17-2) ...
Processing triggers for shared-mime-info (2.1-2) ...
Errors were encountered while processing:
 rstudio
```

Para lo cual tenemos que instalar la dependencia:

```BASH
sudo apt install libclang-dev
```

Y obtendremos:

```
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
You might want to run 'apt --fix-broken install' to correct these.
The following packages have unmet dependencies:
 libclang-dev : Depends: libclang-14-dev (>= 14~) but it is not going to be installed
 rstudio : Depends: libssl-dev but it is not going to be installed
E: Unmet dependencies. Try 'apt --fix-broken install' with no packages (or specify a solution).
```

Para lo que tenemos que ejecutar:

```BASH
sudo apt --fix-broken install
```

Esperamos a que termine la instalación anterior para ejecutar nuevamente:

```BASH
sudo dpkg -i rstudio-2024.12.0-467-amd64.deb
```

Además vamos a instalar libcurl en el sistema para facilitar la descargar de paquetes
y verificar pkg-config esta instalado
```BASH
sudo apt update
sudo apt install libcurl4-openssl-dev
sudo apt install pkg-config
```

De esa forma la instalación del Rstudio estará completa. Podemos hacer la búsqueda del aplicativo en los `puntos` en la parte inferior de barra de trabajo.

## 1. INTRODUCCION A R /2. INSTALACION ----

En mi primer script realizare la suma de dos numeros

```R
1 + 2
```

### 1.1 Directorio de trabajo

```R
# Consultamos el directorio de trabajo
getwd()
# Establecemos el directorio de trabajo
setwd("/home/ins_user/cursoR")
```

#### Creamos y borramos carpetas

```R
#  Creamos Carpetas
dir.create("/home/ins_user/cursoR/AQUA")
# Borramos carpetas
unlink("/home/ins_user/cursoR/AQUA", recursive = TRUE)
# ver que carpetas tenemos 
dir("/home/ins_user/cursoR/")
```

**Ejercicio**  Crea las carpetas data, scripts, results y verifica con dir()

```R
dir.create("/home/ins_user/cursoR/data")
dir.create("/home/ins_user/cursoR/scripts")
dir.create("/home/ins_user/cursoR/results")
```

#### Creamos un objeto

```R
# Creamos los objetos a, b, x
a = 2
b <- 4
x <- "mouse"
a
```
**Ejercicio** Crea un objeto con un nombre (de tu elección) que comience con un número. ¿Qué sucede?

### Función, Argumentos y opciones

```R
# Función raíz cuadrada de 
sqrt(16)       
# Usamos una función y sus opciones
round(3.141516) # Redondea a 0 decimales
?round          # Explica como usar round
round(3.141516, digits = 5) # el resultado de tener 5 decimales
```

### Crear un R Script organizado

```R
# Análisis de datos de iris
# Autor: Damaris Esquen
# Fecha: 21-01-2025
# Este script ejecuta xxx utilizando como input yyy y
# genera un output zzz
# Input dir
data_dir <-  "/home/ins_user/cursoR/data"
# Output dir
results_dir = "/home/ins_user/cursoR/results/"
# cargar librerías
library(datasets)
# Cargar input
datos_iris = iris
```

## 3 INSTALACION DE LIBRERIAS ----

### Instalación desde CRAN

```R
# Desde source
install.packages("cowsay")
# Llamado del paquete
library(cowsay)
say( what = "Hola soy michimichi", by = "cat")
```
**Ejercicio** Que argumentos y opciones tiene la funcion say del paquete cowsay

```R
# usando pacman
install.packages("pacman")
pacman::p_load(cowsay, dplyr)
# desvincular paquetes
detach("package:dplyr", unload = TRUE) #unload = TRUE, además de desvincular el paquete, se eliminarán todas sus funciones de la memoria.
```
### Instalación desde Bioconductor 
```R
# Instalamos el administrador de Bioconductor: BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

```

Si ocurriera un problema como este: 
```R
ERROR: dependencies ‘GenomicRanges’, ‘SummarizedExperiment’, ‘genefilter’, ‘geneplotter’ are not available for package ‘DESeq2’
* removing ‘/home/ins_user/R/x86_64-pc-linux-gnu-library/4.1/DESeq2’
The downloaded source packages are in
	‘/tmp/RtmpuYGIRo/downloaded_packages’
Installation paths not writeable, unable to update packages
  path: /usr/lib/R/library
  packages:
    boot, class, cluster, codetools, foreign, KernSmooth, lattice, mgcv, nlme, nnet,
    rpart, spatial, survival
There were 13 warnings (use warnings() to see them)
```

Ejecutar este código 
```R
BiocManager::install("DESeq2", lib = "~/R/x86_64-pc-linux-gnu-library/4.1")
install.packages("curl")
install.packages("RCurl")
install.packages("httr")
```


## 4. IMPORTACIÓN/GUARDADO DE DATOS EN R----

### Input y Output dir

```R
# Imput dir
data_dir <-  "/home/ins_user/cursoR/data/"
# Output dir
results_dir <-  "/home/ins_user/cursoR/results/"
```

### ruta del archivo con función paste0

```R
# Trabajemos con la data iris
datos_iris <- iris
# Forma corta
write.csv(datos_iris, file = paste0(data_dir, "datos_iris.csv"))
# Forma larga
write.csv(datos_iris, file = "/home/ins_user/cursoR/data/datos_iris_2.csv")
```

### Guardar archivos, objetos y sesión  

Guardar un archivo csv
```R
# Guardar un archivo csv - reescribir
write.csv(datos_iris, file = paste0(data_dir, "datos_iris.csv"))
```
Guardar un objeto y la sesión R
```R
# Guardar un objeto
save(a, file = paste0(results_dir,"result_a.RData"))
# Obtener información de las versiones de los programas
sessionInfo()
# Guardar la sesión en R
save.image(file <- paste0(results_dir, "Mi_sesion.RData" ))
```
Guardar en formato de excel
```R
# Instalamos la librería openxlsx
install.packages("openxlsx")
# Llamamos a la librería
library(openxlsx)
# Guardar excel
write.xlsx(datos_iris, "data/datos.xlsx")
write.xlsx(datos_iris, "data/datos_table.xlsx", asTable = TRUE)
```

### Cargar archivos 
Archivos R.Data
```R
#  cargar archivos tipo R.Data
load(file <-  paste0(results_dir, "result_a.RData"))
load(file <-  paste0(results_dir, "Mi_sesion.RData"))

Archivos separados por un valor
```R
# Cargar archivos separados por un valor
read.table("data/datos_iris.csv", header = TRUE, sep = ",")
read.table("data/Araport11_gene", header = TRUE, sep = "\t")
read.csv(paste0(data_dir, "datos_iris_2.csv"))
data_new <-  read.csv(paste0(data_dir, "datos_iris_2.csv"))
```
Archivos excel
```R
# Importar datos desde Excel con código R
# cargar la librería readxl 
install.packages("readxl")
library(readxl)
# como mirar cuantas hojas tiene tu Excel
excel_sheets("data/Concentracion_data.xlsx")
# Importar excel con código de R: caso 1
caso_1 <- read_excel("data/Concentracion_data.xlsx")
# Importar excel con codigo de R: caso 2
caso_2 <- read_excel("data/Concentracion_data.xlsx", sheet = "Hoja 2")
# Importar excel con codigo de R: caso 3
caso_3 <- read_excel("data/Concentracion_data.xlsx", sheet = "Hoja 3", range = "B3:Q12")
```

### Importar archivos fasta

```R
BiocManager::install("Biostrings")
library(Biostrings)
gene_YY1 <- readDNAStringSet("data/gene_YY1.fasta")
gene_YY1
View(gene_YY1)
```

## 5. ESTRUCTURA DE DATOS 
Define qué tipo de valores existen en R

### 5.1 Tipos de datos 

Los tipos de datos principales son
	a) numérico (numeric) que incluye **double** e **integer**
 	b) caracter (character)
  	c) lógicos (logical)

```R
# Tipo numerico
x <- 1.1
x
typeof(x)
y <- 3
y
typeof(y)
y <- 3L
y
typeof(y)

# Tipo caracter
z <- "Estas mejorando en esto"
z
typeof(z)

# Tipo lógico
p <- TRUE
p
typeof(p)
```
### 5.2 Clases de datos
Define cómo R Interpreta el Valor

```R
class(x)
class(y)
class(z)
class(p)
```

#### 5.2.1 Factores
```R
sex_f <- factor(c("Femenino", "Masculino", "Femenino", "Femenino"))
# ver el tipo 
typeof(sex_f)
# ver la clase
class(sex_f)
# ver la estructura
str(sex_f)
sex_f
```
### 5.3 Variables
Tenemos variables continuas, discretas, categóricas y lógicas

### 5.4 Estructura de datos

### 5.4.1 Escalar
```R
x
y
z
p

```

### 5.4.2. Vector 
Un vector es simplemente una combinación de varios escalares del mismo tipo almacenados como un único objeto.
Veamos como crear vectores

#### Crear vectores

```R
# Crear vectores con la funcion c()
edad <- c(16, 18, 24, 29) # vector numerico
edad
alumnos <- c("Ines", "Ana", "Luis", "Teresa") # vector con caracteres
alumnos
honor <- c(FALSE, TRUE, FALSE, FALSE) # vector lógico
honor

# Crear vectores con a:b
Notas1 <- 16:20
Notas1

# Vectores con seq(from, to, by y/o length.out)
Notas2 <- seq(11, 20, by = 3)
Notas2
Notas3 <- seq(11, 20, length.out = 3)
Notas3

# vectores con rep(x, times, each, length.out)
Notas4 <- rep(x = 1:2, times = 3, each = 4)
Notas4
Notas5 <- rep(c(1,3,5), each = 2)
Notas5
```

#### Operadores

```R
5 == 4
5 <= 6
?Syntax
x <- 1
y <- c(1:8)
x %in% y
```

## 5.4.3 Matrices ----
Las matrices no son más que un conjunto de vectores del mismo tipo apilados, las matrices tienen **dos dimensiones**: filas y columnas

```R
# Crear matrices con la función matrix( data, nrow, ncol)
matrix1 <- matrix(data = 1:10, nrow = 2, ncol = 5)
matrix2 <- matrix(data = 1:10, nrow = 5, ncol = 2)

# Crear matrices a partir de vectores con rbind() y cbind()
nickname <- c("pud", "gab", "Lu")
animal <- c("perro", "raton", "gato")
# rbind() combina vectores como filas.
matrix3 <- rbind(nickname, animal)
# cbind() combina vectores como columnas.
matrix4 <- cbind(nickname, animal)
matrix4

# Transponemos matrices con la funcion t()
matrix5 <- t(matrix4)
matrix5
dim(matrix4)
dim(matrix5)
```
## 5.4.4 array ----
Un array es una estructura que puede contener vectores, matrices, todos del mismo tipo en un **número arbitrario de dimensiones**
Creo un array con la función array(data, dim)


## 5.4.5 Dataframes ----
Es una estructura de datos de **dos dimensiones** que puede almacenar vectores que pueden ser de diferente tipo: **heterogénea**

Analicemos el dataframe iris
```R
# carguemos nuestro dataframe iris
datos_iris <- iris
# Ver el objeto en una ventana separada
View(datos_iris)
# Ver las primeras 6 filas del dataframe dentro de R Studio
head(datos_iris)
# puedes especificar el numero de filas a mostrar
head(datos_iris, 3)
# ver la estructura del dataframe
str(datos_iris)
colnames(datos_iris)
row.names(datos_iris)
nrow(datos_iris)
ncol(datos_iris)
```

#### Crear un dataframe
Pueden ser del mismo tipo o diferente tipo de datos, creamos con la función data.frame()

```R
# Mismo tipo de datos
df_num <- data.frame(a = 1:3, b = 4:6) 
row.names(df_num)
colnames(df_num)
str(df_num)

# Diferentes tipos de datos
df <- data.frame(genes = paste0("Gen", c(1:8)), 
                 expresion = c(5.8, 5.5, 6.3, 6.1, 7.8, rep(8,3)), 
                 condicion_exp = c(rep("Control", 4), rep("Tratramiento",4)))
df
```

**Ejercicio** cual es el nombre de las filas y columnas de df, y cuantas filas tiene

```R
# Puedes usar
row.names(df)
colnames(df)
str(df)
nrow(df)
```

Agregar información de edad al df, con las siguientes opciones:

```R
#  Opción A
# Declarar un vector
log2FoldChange <- c(rep("0.000",4), 0.396, 0.433, 0.433, 0.433)
log2FoldChange
# agregarlo al df
datos <- data.frame(df, log2FoldChange)
datos

# Opción B
# Con el símbolo $
datos$funcion_gen <- rep(c("Metabolismo", "Señalización", "Apoptosis", "Inflamación"), 2)
datos
```

## 5.4.5.1 Tibble 

Recordemos que estructura tiene la data iris
```R
df_iris <- iris
df_iris
```

Para manipular tibble necesitamos el paquete tibble

```R
install.packages("tibble")
library(tibble)
# Obtener un tibble
Podemos coercionar con la función as_tibble() o crear con la función tibble()
# Obtenemos un tibble mediante coerción
tibble_iris <- as_tibble(iris)
tibble_iris
```

## 5.4.6 Listas ----
Son estructuras de datos que pueden almacenar, vectores, arrays, dataframes: **heterogéneas**, 
pero son **unidimensionales** ya que los elementos se almacenan de forma ordenada y no tiene filas ni columnas

Creamos una lista

```R
mi_lista <- list(sex_f, matrix1, df)
mi_lista
```

Analizo una lista

```R
# largo de una lista : número de elementos que contiene
length(mi_lista)
# Dimensión de una lista
dim(mi_lista)
# clase de una lista
class(mi_lista)
```

## 6. Manipulación de datos con R base ----

Recordemos como obtenemos información de una estructura de datos

```R
Creamos nuestro dataframe df
df <- data.frame(genes = paste0("Gen", c(1:8)), 
                 expresion = c(5.8, 5.5, 6.3, 6.1, 7.8, rep(8,3)), 
                 condicion_exp = c(rep("Control", 4), rep("Tratramiento",4)))
df

dim(df) 	# dimensiones [fila, columna] 
length(df) 	# largo, número de columnas 
row.names(df)	# nombre de las filas
colnames(df)	# nombre de las columnas
ncol(df)	# número de columnas
nrow(df) 	# número de filas 
names(df)	# nombre de las columnas 
str(df)	        # Estructura
head(df)	# muestra las 6 primeras filas
head(df, 2)	# se puede seleccionar numero de filas en n
tail(df)	# muestra las 6 ultimas filas
View(df) 	# Mostrar objetos en una ventana separada dentro de RStudio
```

### 6.1 Indexación de objetos en R
Acceder o filtrar los datos según ciertos criterios.

#### 6.1.1 Indexación en vectores
Para indexar un elemento dentro de un vector usando [ ], necesitamos escribir el número de posición del elemento dentro de los corchetes [ ]

```R
# Creamos el vector v
v <- c(10, 20, 30, 40, 50)
# Para seleccionar un solo elemento: vector[posición_del_elemento]
v[2]
# Para seleccionar múltiples elementos
v[c(1, 3, 5)]
# Para excluir elementos: vector[-posición_del_elemento]
v[-2]
# Para indexar con condiciones: vector[elemento y condicion]
v[v > 25]
```

#### 6.1.2 Indexación de matrices y dataframes
Para indexar un data.frame, debes especificar la posición de los valores dentro de sus dos dimensiones [filas , columnas].

```R
# Indexamos usando [numero o nombre de fila, numero o nombre de columna]
df
df[1,]
df[ , 1]
df[ , 1, drop = FALSE]
df[1,1]
df[-8,]
df[,-3]
df[,c("genes", "condicion_exp")]
```

#### 6.1.3 Indexación usando $
El operador $ se puede utilizar para seleccionar, asignar nuevos valores o crear nuevas variables (columnas) en data.frames, así como elementos en listas.

```R
# $ nos sirve para seleccionar = indexing
colnames(df)
(df$expresion)
# $ nos sirve para asignar nuevos valores
df$condicion_exp <- c(rep("T_1",4), rep("T_2",4))
df
# $ nos sirve para crear nuevas variables
df$numero_exones <- c(4, 6, 2, 4, 7, 3, 5, 6)
df
# asigna un mismo valor a toda una nueva columna
df$funcion_gen <- "metabolismo"
df
# cambiamos del nombre del dataframe
df_mod <- df
```

#### 6.1.4 Indexamos según condiciones
Podemos usar operadores de comparación para obtener valores TRUE o FALSE.
Podemos hacer subsetting en función de condiciones usando operadores lógicos

```R
# Creamos el vector test
test <- seq(from = 1, to = 9, by = 2)
test
# evaluo valores mayores que 3 en el vector test
test > 3
# evaluo valores menores iguales que 5 en el vector test
test <= 5
# evaluo valores exactamente iguales a 7
test == 7
# evaluo valores no iguales a 9
test != 9

# Puedo usar indexing para extraer elementos de un vector que cumplen ciertas condiciones lógicas
# evaluo x O y utilizando pipe |
test[test < 3 | test > 5]
# evaluo x Y y utilizando &
test[test >= 3 & test < 5]
# evaluo x %in% y (match) 
test[test %in% c(3,5)]
```

Podemos evaluar si los valores dentro de un vector son numericos o caracteres o que estructura es.
```R
# usando la funcion is.clase() consulto la clase cierto objeto y obtengo resultados TRUE o FALSE
is.numeric(frutas)
is.numeric(test)
is.character(frutas)
is.character(test)
# Podemos evaluar si es vector u otro estructura
is.vector(frutas)
is.list(frutas)
```
**ejercicios**
E. 1. Crea un dataframe llamado test2, donde la primera columna sea el vector test, la segunda columna llamada multip, sea la multiplicación de los valores de test por 5. Finalmente evalúa si los valores obtenidos de multip son mayores de 30

```R
# Ejercicio 1
test2 <- data.frame(test, 
                    multip = test * 5 ) 
test2
test2[, 2] > 30
```
E.2. Convierte la columna tratamiento  a factor. Renombra las filas con los nombres de los genes y elimina la columna 1. Visualizar los primeros 3 datos de la columna 1 manteniendo la estructura de dataframe. Guarda el dataframe en un archivo csv.

```R
# Ejercio 2
# Convierte la columna condicion_exp a factor
str(df_mod)
df_mod$condicion_exp <- as.factor(df$condicion_exp)
str(df_mod)
# Renombra las filas con los nombres de los genes
row.names(df_mod)
df_mod
row.names(df_mod) <- c(genes = paste0("Gen", c(1:8)))
df_mod 
# elimina columna 1
df_mod <- df_mod[,-1]
head(df_mod)
# Visualizar los primeros 3 datos de la columna 1 manteniendo la estructura de dataframe   
df_mod[1:3,1, drop = FALSE] # drop=FALSE nos muestra en formato dataframe

# Guardamos df_mod
write.csv(df_mod, file = "/home/ins_user/cursoR/results/df_mod.csv")
```
#### 6.1.5 Datos faltantes ----

```R
concentracion <- c(43, 47, NA, 48, 53, NA)
concentracion
# IGNORAR
# Primero verifico si hay NA en el vector
is.na(concentracion)
# media de concentración eliminando NA
mean(concentracion)
mean(concentracion, na.rm = TRUE)
# ELIMINAR
# extraer todos los elementos que no son NA
is.na(concentracion)
concentracion[!is.na(concentracion)]
# NULL es la ausencia total de datos 
jovenes <- data.frame(nombre = c("Ana", "Juan", "Luis"), edad = c(25, NA, 30)) 
jovenes
jovenes$altura 

```
### 6.2 Funciones de Rbase

```R
# Aprendamos funciones con el dataframe df_mod
df_mod
```

#### 6.2.1 length(): Devuelve el número total de elementos en un vector o lista.

Sintaxis: length(x)
```R
# length(): número total de elementos 
# si colocamos el dataframe nos da el numero de columnas
length(df_mod)
# si indicamos la columna nos devuelve el numero de elementos en la columna
length(df_mod$condicion_exp)
```

#### 6.2.2 unique() : Devuelve los valores únicos de un vector o columna.

Sintaxis: unique(x)
```R
# unique() : valores únicos 
unique(df_mod$condicion_exp)
```

#### 6.2.3 sort() : Ordena los valores de un vector.

Sintaxis: sort(x, decreasing = TRUE, na.last = TRUE)
```R
# sort() : Ordena y muestra elementos ordenados
sort(df_mod$numero_exones)
# decreasing = TRUE ordena en forma descendente y na.last = TRUE coloca NA al final
sort(df_mod$numero_exones, decreasing = TRUE)
```

#### 6.2.4 order() : Ordenar un dataframe por una columna específica.

Sintaxis: order(x, decreasing = FALSE)
```R
# order() : Ordena y muestra indices de posición de elementos ordenados
order(df_mod$numero_exones)
# Podemos indexar estos elementos, es bastante util en dataframes
df_mod[order(df_mod$numero_exones),]
```

#### 6.2.5 table() : Crear una tabla de frecuencias de un vector.

Sintaxis: table(x, exclude = NULL, useNA = c("no", "ifany", "always"))
```R
# table() : genera tabla de frecuencias
table(df_mod$numero_exones)
# El argumento useNA,  "no" (No se incluirán NA), "ifany" (se incluirán NA solo si existen en el vector)
# y "always" (Siempre incluir los valores NA, aunque no haya ninguno)
table(concentracion, useNA = "ifany")
```

#### 6.2.6 subset() seleccionamos una parte de un conjunto de datos (segun condiciones lógicas) o seleccionamos columnas especificas.

Sintaxis: subset(x, subset, select, drop = FALSE, ...)
```R
# con el argumento subset selecciona un dataframe según una condición lógica
subset(df_mod, subset = expresion < 6)
# con el argumento select selecciono columnas especificas
subset(df_mod, select = c(expresion, numero_exones))
# Puedo usar ambos argumentos
subset(df_mod, subset = expresion < 6, select = c(expresion, numero_exones))
# Puedo eliminar columnas si le coloco - adelante de la columna a eliminar
subset(df_mod, select = -numero_exones)
```

#### 6.2.7 sample() : Muestra aleatoria de un vector.

Sintaxis: sample(x, size, replace = FALSE, prob=NULL)
```R
# size es el número de elementos que deseas seleccionar de x.
sample(df_mod$expresion, 2)
```

#### 6.2.8 summary() : Resumen estadístico de un vector o dataframe.

Sintaxis: summary(x)
```R
# Para vectores numéricos: Muestra el mínimo, primer cuartil, mediana, media, tercer cuartil y máximo.
summary(df_mod$expresion)
# Para factor muestra frecuencias absolutas para cada categoría 
summary(df_mod$condicion_exp)
# Para dataframes: Muestra un resumen para cada columna según su tipo (numérica o categórica).
summary(df_mod$funcion_gen)
```

#### 6.2.9 duplicated() Verificar si hay elementos duplicados en un vector o dataframe.

Sintaxis: duplicated(x)
```R
# Identifica elementos duplicados en un vector, lista o dataframe. Devuelve un vector LÓGICO 
duplicated(df_mod$expresion)
# Pero podemos extraer la información haciendo indexing
df_mod[duplicated(df_mod$expresion), ]
```

#### 6.2.10 which(): Devuelve los índices de las posiciones donde una condición lógica es verdadera

Sintaxis: which(x, arr.ind = FALSE)
```R
# x es un vector lógico o una condición que evaluará si cada elemento cumple con la condición.
which(df_mod$expresion < 7)
# arr.ind devuelve una matriz con las posiciones por fila y columna
which(df == "Gen4", , arr.ind = TRUE) 
```

**Ejercicio** ¿Que gen contiene una expresión de 6.1? 

```R
# Opcion A
which(df_mod == 6.1, arr.ind = TRUE)
# Opcion B
df_mod[which(df_mod == 6.1),]
```

## Instalemos paquetes para la siguiente clase

## Instalar estos paquetes 

```BASH
sudo apt update
sudo apt ungrade
sudo apt-get install libudunits2-dev libgdal-dev libgeos-dev libproj-dev
```


```R
pacman::p_load(dplyr, ggplot2, RColorBrewer, paletteer, cowplot)
pacman::p_load(rnaturalearth, rnaturalearthdata, sf, viridis)
```




