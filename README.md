# Clase Introducción a R"

- author: Damaris Esquen
- Fecha: 12/02/2025

## 1. INTRODUCCION A R /2. INSTALACION ----

```R
getwd()
setwd("D:/cursoR/")
```

### En mi primer script realizare la suma de dos numeros

```R
1 + 2
```

### Consultamos el directorio de trabajo

```R
getwd()
```
### Establecemos el directorio de trabajo

```R
setwd("D:/cursoR")
```

### Crear carpetas

```R
dir.create("D:/cursoR/AQUA")
```

### Borrar carpetas

```R
unlink("D:/cursoR/AQUA", recursive = TRUE)
dir.create("D:/cursoR/data")
dir.create("D:/cursoR/results")
# ver que carpetas tenemos 
dir("/cursoR/")
```

### Creamos un objeto

```R
a = 2
b <- 4
x <- "mouse"
a
```
### Usamos una función

```R
sqrt(16)        # devuelve la raíz cuadrada de 16
```

### Usamos una función y sus opciones

```R
round(3.141516) # Redondea a 0 decimales
args(round)     # Muestra los argumentos de round()
round(3.141516, digits = 5) # el resultado de tener 5 decimales
?round          # Explica como usar round
```

### Crear un R Script organizado

```R
# Análisis de datos de iris
# Autor: Damaris Esquen
# Fecha: 21-01-2025
# Este script ejecuta xxx utilizando como input yyy y
# genera un output zzz
# Input dir
data_dir <-  "D:/cursoR/data"
# Output dir
results_dir = "D:/cursoR/results/"
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
library(cowsay)
say( what = "Hola soy michimichi", by = "cat")
?say
```

```R
# usando pacman
install.packages("pacman")
pacman::p_load(cowsay, dplyr)
```

```R
# Instalación desde Bioconductor 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
```

## 4. IMPORTACIÓN/GUARDADO DE DATOS EN R----

### Input dir

```R
data_dir <-  "/cursoR/data/"
```

### Output dir

```R
results_dir <-  "/cursoR/results/"
```

### ruta del archivo con ejemplo de guardado csv

```R
datos_iris <- iris
```

### Forma corta

```R
write.csv(datos_iris, file = paste0(data_dir, "datos_iris.csv"))
```

### Forma larga

```
write.csv(datos_iris, file = "D:/courses_bioinfo/Curso INS/cursoR/data/datos_iris_2.csv")
```

### Guardar un archivo csv

```R
write.csv(datos_iris, file = paste0(data_dir, "datos_iris.csv"))
```

### Guardar un objeto

```R
save(a, file = paste0(results_dir,"result_a.RData"))
```

### Guardar la sesión en R

```R
save.image(file <- paste0(results_dir, "Mi_sesion.RData" ))
```

### Obtener información de las versiones de los programas

```R
sessionInfo()
```

### Guardar en formato de excel

```R
library(openxlsx)
write.xlsx(datos_iris, "data/datos.xlsx")
write.xlsx(datos_iris, "data/datos_table.xlsx", asTable = TRUE)
```

### Cargar archivos tipo R.Data

```R
load(file <-  paste0(results_dir, "result_a.RData"))
load(file <-  paste0(results_dir, "Mi_sesion.RData"))
```

### Cargar archivos separados por un valor

```R
read.table("data/datos_iris.csv", header = TRUE, sep = ",")
read.table("data/Araport11_gene", header = TRUE, sep = "\t")
read.csv(paste0(data_dir, "datos_iris_2.csv"))
data_new <-  read.csv(paste0(data_dir, "datos_iris_2.csv"))
```

### Importar datos desde Excel con código R

```R
# cargar la librería readxl 
install.packages("readxl")
library(readxl)
```

```R
# como mirar cuantas hojas tiene tu Excel
excel_sheets("data/Concentracion_data.xlsx")
```

```R
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

## 5. ESTRUCTURA DE DATOS ----

### Tipos de datos ----

#### Tipo numérico

```R
x <- 1.1
x
typeof(x)
y <- 3
y
typeof(y)
y <- 3L
y
typeof(y)
```

#### Tipo caracter

```R
z <- "Estas mejorando en esto"
z
typeof(z)
```

#### Tipo Logico

```R
p <- TRUE
p
typeof(p)
class(p)
```

### 5.1. Vector ----

#### vectores con c

```R
edad <- c(16, 18, 24, 29) # vector numerico
edad
alumnos <- c("Ines", "Ana", "Luis", "Teresa") # vector con caracteres
alumnos
honor <- c(FALSE, TRUE, FALSE, FALSE) # vector lógico
honor
```

#### vectores con a:b

```R
Notas1 <- 16:20
Notas1
```

#### Vectores con seq(from, to, by y/o length.out)

```R
Notas2 <- seq(11, 20, by = 3)
Notas2
Notas3 <- seq(11, 20, length.out = 3)
Notas3
```

#### vectores con rep(x, times, each, length.out)

```R
Notas4 <- rep(x = 1:2, times = 3, each = 4)
Notas4
Notas5 <- rep(c(1,3,5), each = 2)
Notas5
```

#### Tipo factor : factor especial

```R
sex <- c("Femenino", "Masculino", "Femenino", "Femenino")
sex_f <- factor(c("Femenino", "Masculino", "Femenino", "Femenino"))
# ver el tipo del dato
typeof(sex_f)
# ver la clase del objeto
class(sex_f)
# ver la estructura del objeto
str(sex_f)
sex_f
# si queremos reordenar los niveles
sex_f  <- factor(sex_f, levels = c("Masculino", "Femenino"))
sex_f
```

#### Ver la lista de operadores y su significado

```R
5 == 4
5 <= 6
?Syntax
x <- 1
y <- c(1:8)
x %in% y
```

## 5.2 Matrices ----

### matrix() crea matrices desde cero con datos, filas y columnas.

```R
matrix1 <- matrix(data = 1:10, nrow = 2, ncol = 5)
matrix2 <- matrix(data = 1:10, nrow = 5, ncol = 2)
```

### Crear matrices con rbind() y cbind()

```R
nickname <- c("pud", "gab", "Lu")
animal <- c("perro", "raton", "gato")
```

### rbind() combina vectores como filas.

```R
matrix3 <- rbind(nickname, animal)
```

### cbind() combina vectores como columnas.

```R
matrix4 <- cbind(nickname, animal)
matrix4
```

### tansponer matrix

```R
matrix5 <- t(matrix4)
matrix5
dim(matrix5)
```
## 5.3 array ----

```R
mi_array <- array(data = 1:16, dim = c(2, 2, 2, 2))
mi_array
dim(mi_array)
```

## 5.4 Dataframes ----

```R
datos_iris <- iris
```

### ver las primeras 6 filas del dataframe

```R
head(datos_iris)
```

### puedes especificar el numero de filas a mostrar

```R
head(datos_iris, 3)
```

### ver la estructura del dataframe

```R
str(datos_iris)
```

### Crear un dataframe

#### Mismo tipo de datos

```R
df_num <- data.frame(a = 1:3, b = 4:6) 
row.names(df_num)
colnames(df_num)
str(df_num)
```

#### Diferentes tipos de datos

```R
df <- data.frame(genes = paste0("Gen", c(1:8)), 
                 expresion = c(5.8, 5.5, 6.3, 6.1, 7.8, rep(8,3)), 
                 condicion_exp = c(rep("Control", 4), rep("Tratramiento",4)))
df
```

### Conozco mi estructura de datos

```R
colnames(df)
row.names(df)
str(df)
nrow(df)
ncol(df)
```

### Agregar información de edad al df

#### Opción A

```R
# declarar un vector
log2FoldChange <- c(rep("0.000",4), 0.396, 0.433, 0.433, 0.433)
log2FoldChange
```

```R
# agregarlo al df
datos <- data.frame(df, log2FoldChange)
datos
```

#### Opción B

```R
datos$funcion_gen <- rep(c("Metabolismo", "Señalización", "Apoptosis", "Inflamación"), 2)
datos
```

## 5.5 Tibble ----

### Recordemos que estructura tiene la data iris

```R
df_iris <- iris
df_iris
str(iris)
```

### Para manipular tibble necesitamos el paquete tidyverse

```R
install.packages("tidyverse")
library(tidyverse)

```

### Obtener un tibble

#### Opción A

```R
tibble_iris <- as_tibble(iris)
tibble_iris
```

#### Opción B 

```R
tibble(x = 1:5, y = 1, z = x^2 + y)
```

## 5.6 Listas ----

### Creamos una lista

```R
mi_lista <- list(sex, matrix1, df)
mi_lista
```

### largo de una lista : número de elementos que contiene

```R
length(mi_lista)
```

### Dimensión de una lista

```R
dim(mi_lista)
```

### clase de una lista

```R
class(mi_lista)
```

## 6. Manipulación de datos con R base ----

### Recordemos como obtenemos información de la estructura

```R
dim(df) 		# dimensiones [fila, columna] 
length(df) 	# largo, número de columnas 
row.names(df)	# nombre de las filas
colnames(df)	# nombre de las columnas
ncol(df)	 	# número de columnas
nrow(df) 	# número de filas 
names(df)	# nombre de las columnas 
str(df)		 # Estructura
head(df)	# muestra las 6 primeras filas
head(df, 2)	# se puede seleccionar numero de filas en n
tail(df)		# muestra las 6 ultimas filas
```

### Se puede adornar la salida colocando notas

```R
cat("Las dimensiones son:", dim(df), "\n")
cat("El número de columnas es:", ncol(df), "\n")
```

### Indexación de objetos en R

#### En vectores

```R
v <- c(10, 20, 30, 40, 50)
```

#### Para seleccionar un solo elemento

```R
v[2]
```

#### Para seleccionar múltiples elementos

```R
v[c(1, 3, 5)]
``` 

#### Para excluir elementos

```R
v[-2]
``` 

#### Para indexar con condiciones

```R
v[v > 25]
```

#### Indexación de matrices y dataframes

```R
df <- data.frame(genes = paste0("Gen", c(1:8)), 
                 expresion = c(5.8, 5.5, 6.3, 6.1, 7.8, rep(8,3)), 
                 condicion_exp = c(rep("Control", 4), rep("Tratramiento",4)))
```

```R
df
df[1,]
df[ , 1]
df[ , 1, drop = FALSE]
df[1,1]
df[-8,]
df[,-3]
df[,c("genes", "condicion_exp")]
```

### Alternativamente podemos usar $

#### $ nos sirve para seleccionar = indexing

```R
colnames(df)
(df$expresion)
```
#### $ nos sirve para asignar nuevos valores

```R
df$condicion_exp <- c(rep("T_1",4), rep("T_2",4))
df
```

#### $ nos sirve para crear nuevas variables

```R
df$numero_exones <- c(4, 6, 2, 4, 7, 3, 5, 6)
df
```

#### asigna un mismo valor a toda una nueva columna

```R
df$funcion_gen <- "metabolismo"
df
```

#### guardamos el dataframe

```R
df_mod <- df
```

#### index de una lista: sublistas

```R
mi_lista[1] # sublista
mi_lista[[1]]        # vector
mi_lista[[1]][[1]]   # elementos del vector
mi_lista[[2]]        # matrix
mi_lista[[2]][[1,3]] # en elemento de matrix
mi_lista[[3]]        # dataframe
mi_lista[[3]][[1,3]] # en elemento de dataframe
```

#### Podemos acceder a los elementos de la lista utilizando $

```R
mi_lista$vector
mi_lista[["vector"]]
```

#### mi_lista[[c(sublista, elemento)]] 

```
mi_lista[[c(3,3)]]
```

#### Podemos asignar nombres a cada sublista

```R
names(mi_lista) <- c("vector", "matrices", "dataframe")
str(mi_lista)
```

#### Pruebas de declaración utilizando operadores lógicos

```R
test <- seq(from = 1, to = 9, by = 2)
test
```

#### evaluo valores mayores que 3 en el vector test

```R
test > 3
```

#### evaluo valores menores iguales que 5 en el vector test

```R
test <= 5
```

#### evaluo valores exactamente iguales a 7

```R
test == 7
```

#### evaluo valores no iguales a 9

```R
test != 9
```

#### evaluo x O y utilizando pipe |

```R
test[test < 3 | test > 5]
```

#### evaluo x Y y utilizando &

```R
test[test >= 3 & test < 5]
```

#### evaluo x %in% y (match) 

```R
test[test %in% c(3,5)]
```

#### Se pueden evaluar numeros y tambien caracteres 

```R
frutas <- c("uva", "manzana", "fresa", "sandía")
frutas
frutas[frutas == "fresa"]
```

#### Podemos evaluar si los valores dentro de un vector son numericos o caracteres

```R
is.numeric(frutas)
is.numeric(test)
is.character(frutas)
is.character(test)
```

#### Podemos evaluar si es vector u otro estructura

```R
is.vector(frutas)
is.list(frutas)
# La diferencia entre 
frutas == "uva"
frutas[frutas == "uva"]
```

```R
# Ejercicio 1
test2 <- data.frame(test, 
                    multip = test * 5 ) 
test2
test2[, 2] > 30
```

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
```

```R
# Datos faltantes ----
concentracion <- c(43, 47, NA, 48, 53, NA)
concentracion
# verificar si hay NA en el vector
is.na(concentracion)
# media de concentración eliminando NA
mean(concentracion)
mean(concentracion, na.rm = TRUE)
# Opción A
# extraer todos los elementos que no son NA
is.na(concentracion)
concentracion[!is.na(concentracion)]
# Opción B
# Elimina NA y guarda la posición de los valores eliminados
na.omit(concentracion)
#Opción C:
# genera un vector lógico (TRUE para valores completos y FALSE para los NA).
complete.cases(concentracion)
# Filtra los valores completos directamente
concentracion_cleaned <- concentracion[complete.cases(concentracion)]
concentracion_cleaned
```

### Condicionales (if, else)

```R
temperatura <- c(15, 30, 32, 28, 37, 42)
if (temperatura[1] < 20) {
  print("Es un clima fresco")
} else if (temperatura[1] >= 20 & temperatura[1] <= 30) {
  print("El clima es cálido")
} else {
  print("El clima es caluroso")
}

# condicional for, while, repeat ----
# ejemplo for
for (i in 1:5) {
  print(i)
}
# ejemplo while
i <- 1
while (i <= 5) { # esta es la condición
  print(i);
  i <- i + 1
}
```

### Condicional for, hagamos un loop para analizar todos los valores del vector

```R
for (i in 1:length(temperatura)) {
  if (temperatura[i] < 20) {
    print(paste("Temperatura de", temperatura[i],": es un clima fresco"))
  } else if (temperatura[i] >= 20 & temperatura[i] <= 30) {
    print(paste("Temperatura de", temperatura[i],": es un clima cálido"))
  } else {
    print(paste("Temperatura de", temperatura[i],": es un clima caluroso"))
  }
}
```

### Funciones de Rbase

```R
df_mod
```

#### length(): Devuelve el número total de elementos en un vector o lista.

```R
length(df_mod$condicion_exp)
length(df_mod) # si colocamos el dataframe nos da el numero de columnas
```

#### unique() : Devuelve los valores únicos de un vector o columna.

```R
unique(df_mod$condicion_exp)
```

#### sort() : Ordena los valores de un vector en orden ascendente (por defecto) o descendente.

```R
sort(df_mod$numero_exones)
sort(df_mod$numero_exones, decreasing = TRUE)
```

#### order() : Ordenar un dataframe por una columna específica.

```R
order(df_mod$numero_exones)
df_mod[order(df_mod$numero_exones),]
```

#### table() : Crear una tabla de frecuencias de un vector.

```R
table(df_mod$numero_exones)
table(concentracion, useNA = "ifany") #
```

#### subset() Filtrar un dataframe según una condición lógica y seleccionar columnas.

```R
subset(df_mod, subset = expresion < 6)
subset(df_mod, select = c(expresion, numero_exones))
subset(df_mod, subset = expresion < 6, select = c(expresion, numero_exones))
subset(df_mod, select = -numero_exones)
```

#### sample() : Muestra aleatoria de un vector.

```R
sample(df_mod$expresion, 2, replace = FALSE, prob = c(0.6, 0.4))
```

#### summary() : Resumen estadístico de un vector o dataframe.

```R
summary(df_mod$expresion)
summary(df_mod$condicion_exp)
summary(df_mod$funcion_gen)
```

#### duplicated() Verificar si hay elementos duplicados en un vector o dataframe.

```R
duplicated(df_mod$expresion)
df_mod[duplicated(df_mod$expresion), ]
```

#### which(): Devuelve los índices de las posiciones donde una condición lógica es verdadera

```R
which(df_mod$expresion < 7)
which(honor)
which.max(df_mod$expresion) 
```

#### Ejercicio 3

```R
which(df_mod == 6.1, arr.ind = TRUE)
df_mod[which(df_mod == 6.1),]
```
