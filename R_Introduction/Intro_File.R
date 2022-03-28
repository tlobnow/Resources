# DATA CLASSES -----------------------------------------------------------------
165     # numr = numeric
'Hello' # chr  = character
T       # logical = True/False, you can use short or long version
FALSE   


# DATA STRUCTURES --------------------------------------------------------------
  # Vector     = single row of data, same class

  # List       = single row of data, but DIFFERENT data types possible
  # matrix     = 2D version of a vector
  # data frame = many vectors pasted together as columns

Vector  <- c(1:100)
Vector2 <- Vector*2
  # subsetting vectors:
    Days_Vector <- c("Mon", "Tue", "Wed", "Thurs", "Fri")
    days[1]   # just Monday
    days[2:5] # Tue - Friday
    days[-5]  # all but Friday
    days[10]  # not available
List <- c(1, 'apple', T, 'BANANA')
rm(Days_Vector, List, Vector, Vector2, Vector_subset)

# USING FUNCTIONS --------------------------------------------------------------
Values <- 1:100

mean(Values)
median(Values)
min(Values)
max(Values)
sum(Values)
sd(Values)
class(Values)
length(Values)
log(Values)
log10(Values)
mysqrt <- sqrt(Values)




# HOW TO INSTALL PACKAGES AND LOAD LIBRARIES -----------------------------------

library(ggplot2)
#install.packages('ggplot2')
library(dplyr)
library(mosaic)







