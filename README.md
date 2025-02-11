# Clase_Introduccion_a_R
## author: "Damaris Esquen"
## Fecha: "12/02/2025"

# 1. INTRODUCCION A R /2. INSTALACION ----
getwd()
setwd("D:/cursoR/")

# En mi primer script realizare la suma de dos numeros
1 + 2

# Consultamos el directorio de trabajo
getwd()

# Establecemos el directorio de trabajo
setwd("D:/cursoR")

# crear carpetas
dir.create("D:/cursoR/AQUA")
# borrar carpetas
unlink("D:/cursoR/AQUA", recursive = TRUE)
dir.create("D:/cursoR/data")
dir.create("D:/cursoR/results")
# ver que carpetas tenemos 
dir("/cursoR/")

# Creamos un objeto
a = 2
b <- 4
x <- "mouse"
a

# Usamos una función
sqrt(16)        # devuelve la raíz cuadrada de 16
# Usamos una función y sus opciones
round(3.141516) # Redondea a 0 decimales
args(round)     # Muestra los argumentos de round()
round(3.141516, digits = 5) # el resultado de tener 5 decimales
?round          # Explica como usar round


# Crear un R Script organizado
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

# 3 INSTALACION DE LIBRERIAS ----
## Instalación desde CRAN
# Desde source
install.packages("cowsay")
library(cowsay)
say( what = "Hola soy michimichi", by = "cat")

# usando pacman
install.packages("pacman")
pacman::p_load(cowsay, dplyr)

## Instalación desde Bioconductor 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# 4. IMPORTACIÓN DE DATOS EN R----
# Input dir
data_dir <-  "/cursoR/data/"
# Output dir
results_dir <-  "/cursoR/results/"

# ruta del archivo con ejemplo de guardado csv
datos_iris <- iris
# Forma corta
write.csv(datos_iris, file = paste0(data_dir, "datos_iris.csv"))
#Forma larga
write.csv(datos_iris, file = "D:/courses_bioinfo/Curso INS/cursoR/data/datos_iris_2.csv")

# Guardar un archivo csv
write.csv(datos_iris, file = paste0(data_dir, "datos_iris.csv"))
# Guardar un objeto
save(a, file = paste0(results_dir,"result_a.RData"))
#Guardar la sesión en R
save.image(file <- paste0(results_dir, "Mi_sesion.RData" ))
# Obtener información de las versiones de los programas
sessionInfo()

# Guardar en formato de excel
library(openxlsx)
write.xlsx(datos_iris, "data/datos.xlsx")
write.xlsx(datos_iris, "data/datos_table.xlsx", asTable = TRUE)

# Cargar archivos tipo R.Data
load(file <-  paste0(results_dir, "result_a.RData"))
load(file <-  paste0(results_dir, "Mi_sesion.RData"))

# Cargar archivos separados por un valor
read.table("data/datos_iris.csv", header = TRUE, sep = ",")
read.table("data/Araport11_gene", header = TRUE, sep = "\t")
read.csv(paste0(data_dir, "datos_iris_2.csv"))
data_new <-  read.csv(paste0(data_dir, "datos_iris_2.csv"))

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

# Importar archivos fasta
BiocManager::install("Biostrings")
library(Biostrings)
gene_YY1 <- readDNAStringSet("data/gene_YY1.fasta")
gene_YY1
View(gene_YY1)


# 5. ESTRUCTURA DE DATOS ----
# Tipos de datos ----
### Tipo numérico
x <- 1.1
x
typeof(x)
y <- 3
y
typeof(y)
y <- 3L
y
typeof(y)
### Tipo caracter
z <- "Estas mejorando en esto"
z
typeof(z)
### Tipo Logico
p <- TRUE
p
typeof(p)
class(p)

## 5.1. Vector ----
# vectores con c
edad <- c(16, 18, 24, 29) # vector numerico
edad
alumnos <- c("Ines", "Ana", "Luis", "Teresa") # vector con caracteres
alumnos
honor <- c(FALSE, TRUE, FALSE, FALSE) # vector lógico
honor
# vectores con a:b
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

### Tipo factor : factor especial
sex <- c("Femenino", "Masculino", "Femenino", "Masculino")
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

# Ver la lista de operadores y su significado
5 == 4
5 <= 6
?Syntax
x <- 1
y <- c(1:8)
x %in% y

## 5.2 Matrices ----
# matrix() crea matrices desde cero con datos, filas y columnas.
matrix1 <- matrix(data = 1:10, nrow = 2, ncol = 5)
matrix2 <- matrix(data = 1:10, nrow = 5, ncol = 2)
# Crear matrices con rbind() y cbind()
nickname <- c("pud", "gab", "Lu")
animal <- c("perro", "raton", "gato")
# rbind() combina vectores como filas.
matrix3 <- rbind(nickname, animal)
# cbind() combina vectores como columnas.
matrix4 <- cbind(nickname, animal)
matrix4
# tansponer matrix
matrix5 <- t(matrix4)
matrix5
dim(matrix5)

## 5.3 array ----
mi_array <- array(data = 1:16, dim = c(2, 2, 2, 2))
mi_array
dim(mi_array)

## 5.4 Dataframes ----
datos_iris <- iris
# ver las primeras 6 filas del dataframe
head(datos_iris)
# puedes especificar el numero de filas a mostrar
head(datos_iris, 3)
# ver la estructura del dataframe
str(datos_iris)

# Crear un dataframe
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
# Conozco mi estructura de datos
colnames(df)
row.names(df)
str(df)
nrow(df)
ncol(df)

# Agregar información de edad al df
## Opción A
# declarar un vector
log2FoldChange <- c(rep("0.000",4), 0.396, 0.433, 0.433, 0.433)
log2FoldChange
# agregarlo al df
datos <- data.frame(df, log2FoldChange)
datos
## Opción B
datos$funcion_gen <- rep(c("Metabolismo", "Señalización", "Apoptosis", "Inflamación"), 2)
datos


## 5.5 Tibble ----
# Recordemos que estructura tiene la data iris
df_iris <- iris
df_iris
str(iris)
# Para manipular tibble necesitamos el paquete tidyverse
install.packages("tidyverse")
library(tidyverse)
# Obtener un tibble
# Opción A
tibble_iris <- as_tibble(iris)
tibble_iris
# Opción B 
tibble(x = 1:5, y = 1, z = x^2 + y)

## 5.6 Listas ----
# Creamos una lista
mi_lista <- list(sex, matrix1, df)
mi_lista
# largo de una lista : número de elementos que contiene
length(mi_lista)
# Dimensión de una lista
dim(mi_lista)
# clase de una lista
class(mi_lista)
# elementos de esta lista: sublistas
mi_lista[1] # sublista
mi_lista[[1]]        # vector
mi_lista[[1]][[1]]   # elementos del vector
mi_lista[[2]]        # matrix
mi_lista[[2]][[1,3]] # en elemento de matrix
mi_lista[[3]]        # dataframe
mi_lista[[3]][[1,3]] # en elemento de dataframe
# Podemos acceder a los elementos de la lista utilizando $
mi_lista$vector
mi_lista[["vector"]]
# mi_lista[[c(sublista, elemento)]] 
mi_lista[[c(3,3)]]
# Podemos asignar nombres a cada sublista
names(mi_lista) <- c("vector", "matrices", "dataframe")
str(mi_lista)

# 6. Manipulación de datos con R base ----
# Recordemos como obtenemos información de la estructura
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
# Se puede adornar la salida colocando notas
cat("Las dimensiones son:", dim(df), "\n")
cat("El número de columnas es:", ncol(df), "\n")


# Indexación de objetos en R
# En vectores
v <- c(10, 20, 30, 40, 50)
# Para seleccionar un solo elemento
v[2]
# Para seleccionar múltiples elementos
v[c(1, 3, 5)] 
# Para excluir elementos
v[-2] 
# Para indexar con condiciones
v[v > 25]
# Indexación de matrices y dataframes
df <- data.frame(genes = paste0("Gen", c(1:8)), 
                 expresion = c(5.8, 5.5, 6.3, 6.1, 7.8, rep(8,3)), 
                 condicion_exp = c(rep("Control", 4), rep("Tratramiento",4)))
df
df[1,]
df[ , 1]
df[ , 1, drop = FALSE]
df[1,1]
df[-8,]
df[,-3]
df[,c("genes", "condicion_exp")]

# Alternativamente podemos usar $
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
# guardamos el dataframe
df_mod <- df

# Pruebas de declaración utilizando operadores lógicos
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
# evaluo x O y utilizando pipe |
test[test < 3 | test > 5]
# evaluo x Y y utilizando &
test[test >= 3 & test < 5]
# evaluo x %in% y (match) 
test[test %in% c(3,5)]
# Se pueden evaluar numeros y tambien caracteres 
frutas <- c("uva", "manzana", "fresa", "sandía")
frutas
frutas[frutas == "fresa"]
# Podemos evaluar si los valores dentro de un vector son numericos o caracteres
is.numeric(frutas)
is.numeric(test)
is.character(frutas)
is.character(test)
# Podemos evaluar si es vector u otro estructura
is.vector(frutas)
is.list(frutas)
# La diferencia entre 
frutas == "uva"
frutas[frutas == "uva"]

# Ejercicio 1
test2 <- data.frame(test, 
                    multip = test * 5 ) 
test2
test2[, 2] > 30

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

# Condicionales (if, else)
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

# Condicional for, hagamos un loop para analizar todos los valores del vector
for (i in 1:length(temperatura)) {
  if (temperatura[i] < 20) {
    print(paste("Temperatura de", temperatura[i],": es un clima fresco"))
  } else if (temperatura[i] >= 20 & temperatura[i] <= 30) {
    print(paste("Temperatura de", temperatura[i],": es un clima cálido"))
  } else {
    print(paste("Temperatura de", temperatura[i],": es un clima caluroso"))
  }
}

# Funciones de Rbase
df_mod
# length(): Devuelve el número total de elementos en un vector o lista.
length(df_mod$condicion_exp)
length(df_mod) # si colocamos el dataframe nos da el numero de columnas
# unique() : Devuelve los valores únicos de un vector o columna.
unique(df_mod$condicion_exp)
# sort() : Ordena los valores de un vector en orden ascendente (por defecto) o descendente.
sort(df_mod$numero_exones)
sort(df_mod$numero_exones, decreasing = TRUE)
# order() : Ordenar un dataframe por una columna específica.
order(df_mod$numero_exones)
df_mod[order(df_mod$numero_exones),]
# table() : Crear una tabla de frecuencias de un vector.
table(df_mod$numero_exones)
table(concentracion, useNA = "ifany") #
# subset() Filtrar un dataframe según una condición lógica y seleccionar columnas.
subset(df_mod, subset = expresion < 6)
subset(df_mod, select = c(expresion, numero_exones))
subset(df_mod, subset = expresion < 6, select = c(expresion, numero_exones))
subset(df_mod, select = -numero_exones)
# sample() : Muestra aleatoria de un vector.
sample(df_mod$expresion, 2, replace = FALSE, prob = c(0.6, 0.4))
# summary() : Resumen estadístico de un vector o dataframe.
summary(df_mod$expresion)
summary(df_mod$condicion_exp)
summary(df_mod$funcion_gen)
# duplicated() Verificar si hay elementos duplicados en un vector o dataframe.
duplicated(df_mod$expresion)
df_mod[duplicated(df_mod$expresion), ]
# which(): Devuelve los índices de las posiciones donde una condición lógica es verdadera
which(df_mod$expresion < 7)
which(honor)
which.max(df_mod$expresion) 
# Ejercicio 3
which(df_mod == 6.1, arr.ind = TRUE)
df_mod[which(df_mod == 6.1),]

# Funciones de dyplr ----
library(tidyverse)
library(datasets)
iris = iris
View(iris)
str(iris)

## select() Seleccionar columnas ----
# seleccionar una columna
iris %>% select(Sepal.Length) %>% head(4)
# seleccionar dos columnas
iris %>% select(Sepal.Length,Species) %>% head(4)
# seleccionar columnas desde : hasta
iris %>% select(Sepal.Length:Petal.Width) %>% head(4)
# seleccionar columnas por índice
iris %>% select(1,3,5) %>% head(4)
# eliminar columnas 
iris %>% select(-Sepal.Length) %>% head(4)
# eliminar varias columnas
iris %>% select(-c(Sepal.Length, Species)) %>% head(4)
iris %>% select(-c(1,5)) %>% head(4)
# Selecciona todo excepto Species
iris %>% select(!Species) %>% head(4)
# select(contain(")): selecciona columnas que contienen un texto dado
iris %>% select(contains("Sepal")) %>% head(4)
# select(starts_with("")) : columnas que comienzan con el prefijo dado
iris %>% select(starts_with("Petal")) %>% head(4)
# select(ends_with("")) : columnas que terminan con el prefijo dado
iris %>% select(ends_with("Length")) %>% head(4)
# select(where()): columnas que cumplen con una condición dada
iris %>% select(where(is.factor)) %>% head(4)
# select(last_col()): Selecciona la última columna
iris %>% select(last_col()) %>% head(4)

# filter() : filtrar filas de un data frame en base a una o varias condiciones.
# filtrar filas en base a una condición con operador de comparación
iris %>% filter(Sepal.Length > 7.5) %>% head()
iris %>% filter(Species == "versicolor") %>% head()
# filtrar filas en base a una función
iris %>% filter(Petal.Length <= mean(Petal.Length)) %>% head()
# filtrar filas en base a una condición con operador lógico
iris %>% filter(Sepal.Length %in% 5) %>% head()
iris %>% filter(Sepal.Length %in% c(4.5, 5))
# filtrar filas que sean el opuesto a la condición : NO : !
iris %>% filter(!(Sepal.Length %in% c(4.5, 5)))  %>% head()
# filtrar filas que tengan un texto específico con grepl
iris %>% filter(grepl("1.4", Petal.Length)) %>% head()
# filtrar con varias condiciones usan & , |
iris %>% filter(Sepal.Length > 5 & Petal.Length < 1.5)
iris %>% filter(Sepal.Length > 5 | Petal.Length < 1.5)

# slice() : filtra filas basándose en su posición/índice.
# filtra segun posicion de filas
iris %>% slice(1:2)
# filtra la primera fila
iris %>% slice_head()
# filtra la ultima fila 
iris %>% slice_tail()

# arrange(): reordenar filas
# ordena filas de una columna en orden ascendente
iris %>% arrange(Sepal.Length) %>% head()
# ordena filas de una columna en orden descendennte
iris %>% arrange(desc(Sepal.Length)) %>% head()
iris %>% arrange(Species) %>% head()
# ordena filas de varias columnas
iris %>% arrange(Sepal.Length) %>% head()
iris %>% arrange(Sepal.Length, Sepal.Width) %>% head()
# ordena filas por grupos
iris %>% arrange(Species, .by_group = TRUE)

# mutate(): crear nuevas columnas, modificar existentes
# crear una nueva columna
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width) %>% head()
# crear varias columnas
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width,
                 Area_Petal = Petal.Length * Petal.Width) %>% head()
# Modificar columnas existentes
iris %>% mutate(Sepal.Length = Sepal.Length/5) %>% head()
# Modificar columnas, especificando que columnas usando across
iris %>% mutate(across(c(Petal.Length, Petal.Width), log2)) %>% head()
iris %>% mutate(across(!Species, log2)) %>% head()
# Pocisión de las nuevas columnas usando .before y .after en mutate
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width, .before = Species) %>% head()
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width, .after = Sepal.Width) %>% head()
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width, .before = 1) %>% head()
# Mantener o eliminar columnas con .keep en mutate 
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width, .keep = "used") %>% head()
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width, .keep = "unused") %>% head()
iris %>% mutate( Area_Sepal = Sepal.Length * Sepal.Width, .keep = "none") %>% head()

# group_by() para agrupar 
iris %>% group_by(Species)
# Puedo hacer calculos de los datos agrupados?
iris %>% group_by(Species) %>% mean()
iris %>% group_by(Species) %>% mean(Petal.Length)

# summarize(): resumir datos
# resumir datos (media) de una columna
iris %>% group_by(Species) %>% summarise(promedio = mean(Petal.Length))
# resumir datos (media) de varias columnas
iris %>%
  group_by(Species) %>%
  summarize(
    mean_PL = mean(Petal.Length),
    mean_SL = mean(Sepal.Length),
    mean_PW = mean(Petal.Width),
    mean_SW = mean(Sepal.Width)
  )
# si quisieramos resumir de todas las columnas con valores numéricos
iris %>% group_by(Species) %>% summarise_if(is.numeric, mean)
# resumir datos de todas las columnas con summarise_all
iris %>% group_by(Species) %>% summarise_all(mean)
# resumir columnas especificas colocando en un vector summarise_at
iris %>% group_by(Species) %>% summarise_at(c("Sepal.Length", "Sepal.Width"), mean) 

# count()
# Contar cuántas filas hay por cada columna
iris %>% count(Species)
iris %>% count(Petal.Length)
# Puede renombrar la columna n por otro nombre
iris %>% count(Species, name = "número_por_grupo")

# rename() : 
# Cambiar el nombre de las columnas
iris %>% rename(Longitud.Sepalo = Sepal.Length) %>% head()
iris %>% rename(Longitud.Sepalo = 1) %>% head()
# renombrar varias columnas
iris %>% rename(
    Longitud.Sepalo = Sepal.Length,
    Anchura.Sepalo = Sepal.Width,
    Largo.Petalo = Petal.Length,
    Ancho.Petalo = Petal.Width) %>% head()
# renombrar columnas en base a una función con rename_with
iris %>% rename_with(toupper) %>% head()
# rennombrar columnas en base a una funcion para columnas con un texto especifico
iris %>% rename_with(toupper, .cols = contains("Sep")) %>% head()
# Agregar un prefijo a las columnas
iris %>% rename_with(~paste0("NEW_", .)) %>% head()
iris = iris
library(dplyr)
# distinct(): Eliminar filas duplicadas
iris %>% distinct(Species)
iris %>% distinct(Species, .keep_all = TRUE)
# comparemos con unique, duplicated de R base
unique(iris$Species)
duplicated(iris$Species)
sum(duplicated(iris$Species))

# join() : combinar dataframes
d1 <- data.frame(ID = 1:2, X1 = c("a1", "a2"))
d1
d2 <- data.frame(ID = 2:3,X2 = c("b1", "b2"))
d2
# left_join mantiene todas las filas de d1 y agrega la info coincidentente de d2
left_join(d1, d2, by = "ID")
# right_join mantiene todas las filas de d2 y agrega la info coincidente de d1
right_join(d1, d2, by = "ID")
# inner_join unimos filas reteniendo solo filas que estan en ambas df
inner_join(d1, d2, by = "ID")
# full_join unimos filas reteniendo todos los valores, todas las filas
full_join(d1, d2, by = "ID")

# Utilizar select_if() con tipos de clase
iris %>% select_if(is.numeric) %>% head()
# usando condiciones lógicas
iris %>%  select_if(~ is.numeric(.) && all(. < 7))













