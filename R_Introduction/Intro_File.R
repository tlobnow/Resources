# HOW TO GET AND SET YOUR WORKING DIRECTORY (especially when working locally and you want to save files later)
getwd()
setwd("/Users/finnlo/Documents/Github/Resources/R_Introduction/")


# HOW TO INSTALL PACKAGES AND LOAD LIBRARIES -----------------------------------
library(tidyverse)  #install.packages('tidyverse')
library(data.table) #install.packages('data.table')
library(tibble)     #install.packages('tibble')
library(readr)
library(readxl)



# THERE ARE DIFFERENT DATA CLASSES ---------------------------------------------
165     # numr = numeric
'Hello' # chr  = character
T       # logical = True/False, you can use short or long version
FALSE   


# DATA STRUCTURES --------------------------------------------------------------
  # Vector     = single row of data, same class
              Vector  <- c(1:100)
              Vector2 <- Vector*2
              # subsetting vectors:
              Days_Vector <- c("Mon", "Tue", "Wed", "Thurs", "Fri")
              Days_Vector[1]   # just Monday
              Days_Vector[2:5] # Tue - Friday
              Days_Vector[-5]  # all but Friday
              Days_Vector[10]  # not available
  # List       = single row of data, but DIFFERENT data types possible
              List <- c(1, 'apple', T, 'BANANA')
  # matrix     = 2D version of a vector
  # data frame = many vectors pasted together as columns

rm(Days_Vector, List, Vector, Vector2) # remove one or more things (lists, vectors, data frames) from your environment

# YOU CAN USE DIFFERENT BASE FUNCTIONS -----------------------------------------
Values <- 1:100

mean(Values)
median(Values)
min(Values)
max(Values)
sum(Values)
sd(Values)
class(Values)
length(Values)
log(Values) # natural log by default
log10(Values)
mysqrt <- sqrt(Values)
rm(mysqrt, Values)


# HOW TO READ A FILE (preferably you make your life easier and have csv files)
# Read the same table with different file types and make sure you get the same result for all of them:
library(readr)
library(readxl)

dirty_csv <- read.csv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/dirty_files/Dirty_Data.csv', skip = 2)
dirty_tsv <- read_tsv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/dirty_files/Dirty_Data.tsv', skip = 2)
dirty_xlsx <- read.xlsx('dirty_files/Dirty_Data.xlsx', sheetIndex = 1, startRow = 3)
rm(dirty_tsv, dirty_xlsx) # I just want you to see that you CAN open different file types,
                          # but we will remove the other files and work with the data frame we made from the csv file


# CORRECT THE COLNAMES (so they are all written in the same way Xxxx) ----------
# Correct all names at once:
library(data.table)
# EXAMPLE: setnames(df, old = c("old1", "old2", "old3"), new = c("new1", "new2", "new3"), skip_absent = T)
setnames(dirty_csv, old = c("gender", "AGE", "fruit", "sport"), new = c("Gender", "Age", "Fav_Fruit", "Fav_Sport"), skip_absent = T)

# correct colnames individually:
# EXAMPLE: colnames(df)[colnames(df)%in%"Colname1_old"] <- "Colname1_new"
colnames(dirty_csv)[colnames(dirty_csv) %in% 'gender'] <- 'Gender'
colnames(dirty_csv)[colnames(dirty_csv) %in% 'AGE']    <- 'Age'
colnames(dirty_csv)[colnames(dirty_csv) %in% 'fruit']  <- 'Fav_Fruit'
colnames(dirty_csv)[colnames(dirty_csv) %in% 'sport']  <- 'Fav_Sport'


# 


















Pokemon <- read.csv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/Pokemon.csv') # online path
Pokemon <- read.csv('~/Documents/Github/Resources/R_Introduction/Pokemon.csv') # my local path



















# Understand your data 
dim(Pokemon)      # 251 rows, 15 columns
head(Pokemon)     # look at first 6 rows by default
head(Pokemon, 10) # look at first 10 rows
tail(Pokemon)     # look at last 2 rows by default
str(Pokemon)      # shows structure of a data frame = 
                    # variables (names of the columns)
                    # data classes (int, factor, numeric)
summary(Pokemon)  # dumps a load of data on you:
                    # for numerical values:       Min, 1st Qu, Median, Mean, 3rd Qu, Max 
                    # for categorical variables:  length, class, mode

glimpse(Pokemon) # need package 'tibble' for that <- I like this best!



SOTA    <- read.csv('https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv')
Jarda   <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Mouse_data/HZ19_GenotypingJarda.csv", na.strings = c('', ' ', 'NA'))

# find columns in your data using data$...
SOTA$Mouse_ID

# select columns for specific purposes by putting them as a vector
basics           <- c("Mouse_ID", "Address", "Sex", "Longitude", "Latitude", "Year", "HI", "HI_NLoci")

gen.loci        <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
                     "Es1", "Gpd1", "Idh1", "Mpi", "Np", "Sod1", "Es1C", "Gpd1C",
                     "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci", "HI", "Zfy2", "Y", "Region")

Crypto_qPCR.cols <- c("ILWE_Crypto_Ct", "Oocyst_Predict_Crypto", "ILWE_DNA_Content_ng.microliter")



# CORRECT COLNAMES -------------------------------------------------------------

# multiple names like this: setnames(df, old = c("Colname1_old", "Colname2_old", "..."), new = c("Colname1_new", "Colname2_new", "..."), skip_absent = T)
setnames(Jarda, old = c("PIN", "X_Longit", "Y_Latit"), new = c("Mouse_ID", "Longitude", "Latitude"), skip_absent = T)

# single colnames
colnames(df)[colnames(df)%in%"Colname1_old"] <- "Colname1_new"



# CORRECT INPUT OF SPECIFIC ROWS -----------------------------------------------
SOTA$Year[ SOTA$Mouse_ID %in% c("SK_2903", "SK_2904")] <- 2014
    

# CORRECT ROW NAME INPUT BY PATTERN SEARCH -------------------------------------
df$colname <- gsub(pattern = "pattern", replacement = "xxx", x = df$colname)
df$colname[grep("pattern*.", df$colname)] <- "F" # replace the entire input that follows a specific pattern

Jarda_21$HI_NLoci <- gsub(pattern = "HI ", replacement = "", x = Jarda_21$HI_NLoci) # replace part of the colname
    


# Subset data
subSOTA <- SOTA[colnames(SOTA) %in% c(basics, Crypto_qPCR.cols, 'delta_ct_cewe_MminusE')]










