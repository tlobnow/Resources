library(data.table)
library(stackoverflow)
library(dplyr)
library(magrittr)
library(ggplot2)



Trap <- read.csv("HZ19_Trap_eco.csv", header = T)
Dissection <- read.csv("HZ19_Dissections.csv", header = T)
Extractions <- read.csv("DNA extractions - results.csv", header = T)
head(Extractions)

    Sample_ID     <- Dissection$Mouse_ID
    Transect      <- Dissection$Transect
    State         <- Dissection$State
    Address       <- Dissection$Address
    Lat           <- Dissection$Latitude
    Long          <- Dissection$Longitude
    Trap_Date     <- Dissection$Trap_date
    D_Date        <- Dissection$Dissection_date
    Sex           <- Dissection$Sex
    Status        <- Dissection$Status
    Species       <- Dissection$Species
    year          <- Dissection$Year

  myData <- data.frame(Sample = Sample_ID,
                         Transect = Area,
                         Address = DataPoint,
                         Lat = Latitude,
                         Long = Longitude
                      )
    
    
ggplot(myData, aes(x = Long, y = Lat, size = 190)) + 
  geom_point(alpha = 1)
