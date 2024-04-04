#Creador: MARIANO GUERRA
#Llamado de bases de datos

#Llamar la base de datos desde la computadora

install.packages("readr")
library(readr)

library(readr)
datos_titulacion_2_ <- read_csv("Downloads/datos_titulacion(2).csv")
View(datos_titulacion_2_)

#Llamar desde un repositorio (internet)

repositorio <- read_csv(file = "https://raw.githubusercontent.com/ManuelLaraMVZ/titulacion_amino_acidos/main/datos_titulacion%20(2).csv")
head(repositorio) #Para ver encabezado
View(repositorio)


#Gráfica

install.packages("ggplot2")
library(ggplot2)


grafica <- ggplot(repositorio, aes(x=Volumen, y=pH))+ geom_line()+ 
  labs(title = "Titulación cisteína",
       x= "Volumen de ácido (uL)", 
       y="Valor de pH")+
  theme_dark()

grafica
