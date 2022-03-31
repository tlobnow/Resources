# HOW TO GET AND SET YOUR WORKING DIRECTORY ------------------------------------

#This is mostly important if you want to read local files (the computer will in your current directory to find stuff)
#If the computer cannot find the files there, it will of course throw you an error message

getwd()
setwd("/Users/finnlo/Documents/Github/Resources/R_Introduction/") # my personal path on my computer



# HOW TO INSTALL PACKAGES AND LOAD LIBRARIES -----------------------------------

# Every package can be installed using install.packages('xxx').
# Every package can be loaded by using library(xxx).
# Here we will check if you have all the packages installed for the tutorial:
  
library(tidyverse)  #install.packages('tidyverse')
library(data.table) #install.packages('data.table') # for the setnames() function
library(tibble)     #install.packages('tibble')     # to use tibbles
library(readr)      #install.packages('readr')      # to read tsv files
library(xlsx)       #install.packages('readxl')     # to read xlsx files
library(knitr)      #install.packages('knitr')      # for making tables in markdown
library(visdat)     #install.packages('visdat')     # for looking at missing data
library(scales)     #install.packages('scales')     # for the percent() function



# DATA CLASSES -----------------------------------------------------------------

#In R you will encounter data of different classes:
  
# Logical  (e.g. TRUE/T or FALSE/F)
#Numeric   (num, e.g. 1, 2.0, 3, 44.5, 100)
#double    (dbl)
#integer   (int)
#Character (chr, e.g. 'Hello', 'a', '13.4', 'True')

class(165)     
class('Hello')
class(T)     
class(FALSE)




# DATA STRUCTURES  -------------------------------------------------------------

#Vector = single row of data, same class
Vector  <- c(1:100)
Vector2 <- Vector*2
Vector + Vector2

# you can subset your vectors as well
Days_Vector <- c("Mon", "Tue", "Wed", "Thurs", "Fri")
Days_Vector[1]   # just Monday
Days_Vector[2:5] # Tue - Friday
Days_Vector[-5]  # all but Friday
Days_Vector[10]  # not available

#Lists = single row of data, but DIFFERENT data types possible
List <- c(1, 'apple', T, 'BANANA')
List
#matrix     = 2D version of a vector
#data frame = many vectors pasted together as columns


# SOME BASE FUNCTIONS ----------------------------------------------------------
# There are tons of functions in Base R (meaning that they are not in extra packages you need to load)
# Here are just a few that might come in handy!
  
# Don't forget, if you ever need help with a function and how to use it: type ?function and the help page will pop up

Values <- 1:100

?mean
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
# remove one or more things (lists, vectors, data frames) from your environment
rm(Days_Vector, List, Vector, Vector2, mysqrt, Values) 

# HOW TO READ A FILE -----------------------------------------------------------

#Read the same table with different file types and make sure you get the same result for all of them:
#(preferably you make your life easier and have csv files)

#1. We want to skip unnecessary rows (in our data the first two rows don't belong to actually usable data)
#2. We want to take care of how NAs are read in the data. Some people will:
  
#leave a row input blank ()
#type a space ( )
#or actually put 'NA' in the data

#To make sure that the table we add is correctly addressing the NA-situation, we add the na.strings = c(...) or equivalents to each of the read functions (beware, there is no such thing for reading xlsx, so you need to take care of this later)

library(readr)
library(xlsx)

dirty_csv <- read.csv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/dirty_files/Dirty_Data.csv', skip = 2, na.strings = c('', ' ', 'NA'))
dirty_tsv <- read_tsv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/dirty_files/Dirty_Data.tsv', skip = 2, na = c('', ' ', 'NA'))

# I will load the .xlsx file locally, so be aware of your workding directory:
dirty_xlsx <- read.xlsx('dirty_files/Dirty_Data.xlsx', sheetIndex = 1, startRow = 3)
rm(dirty_tsv, dirty_xlsx) 

# I just want you to see that you CAN open different file types,
# but we will remove the other files and work with the data frame we made from the csv file




# HOW TO CORRECT DATA INPUT ----------------------------------------------------

## Correct the colnames
#Student --> Student_ID
#gender  --> Gender
#AGE     --> Age
#fruit   --> Fav_Fruit
#sport   --> Fav_Sport

#We want all column names to have an uppercase first letter, rest can be lower case (Xxxx).
#You can do this by either correcting all in one or change them individually.
#For the 'all in one' version you have to load the 'data.table' package and follow this pattern:

    # EXAMPLE: setnames(df, old = c("old1", "old2", "old3"), new = c("new1", "new2", "new3"), skip_absent = T)

#Alternatively, you can change the colnames individually, for example using this pattern:
colnames(df)[colnames(df)%in%"Colname1_old"] <- "Colname1_new"


library(data.table)

# all at once
setnames(dirty_csv, old = c("Student", "gender", "AGE", "fruit", "sport"), new = c("Student_ID", "Gender", "Age", "Fav_Fruit", "Fav_Sport"), skip_absent = T)

# individually
colnames(dirty_csv)[colnames(dirty_csv) %in% 'Student'] <- 'Student_ID'
colnames(dirty_csv)[colnames(dirty_csv) %in% 'gender']  <- 'Gender'
colnames(dirty_csv)[colnames(dirty_csv) %in% 'AGE']     <- 'Age'
colnames(dirty_csv)[colnames(dirty_csv) %in% 'fruit']   <- 'Fav_Fruit'
colnames(dirty_csv)[colnames(dirty_csv) %in% 'sport']   <- 'Fav_Sport'

dirty_csv

## Correct the row inputs per column -------------------------------------------
### Student_ID #################################################################

#Character instead of number
#wrong input _0010 instead of _010

#Let's look at two ways to do this: you can either select the column and change input for several rows at once:

    # EXAMPLE: df$colname[df$colname %in% c("old1", "old2")] <- c("new1", "new2")

#Or you can change the row input using the gsub function, which is generally better when you also look for a pattern,
#but it can be an individual row input like in our case (generally the first option is better, but both get the job done!):

    # EXAMPLE: df$colname <- gsub(pattern = "old1", replacement = "new1", x = df$colname)


# all in one
dirty_csv$Student_ID[dirty_csv$Student_ID %in% c("Student_00eight", "Student_0010")] <- c("Student_008", "Student_010")
dirty_csv$Student_ID[ dirty_csv$Name %in% "Thomas"] <- 'Student_011'

# or individually
dirty_csv$Student_ID <- gsub(pattern = "Student_00eight", replacement = "Student_008", x = dirty_csv$Student_ID)
dirty_csv$Student_ID <- gsub(pattern = "Student_0010",    replacement = "Student_010", x = dirty_csv$Student_ID)
dirty_csv$Student_ID[ dirty_csv$Name %in% "Thomas"] <- 'Student_011'

dirty_csv


### Name #######################################################################

#exclude NA students (no name, we exclude the input)
#unified name pattern (correct names to follow Xxxx)


#We will now talk about piping: 
#You can use the pipe symbol **%>%** to start with data, narrow it down, change it, add columns, etc.
#For that, you have to specify the data you start with and tell R where you want to save the piped things (assign object).
#This object can be the same data frame you start with, but it in that case you will *override* the data from before. 

# df   <- df %>% bla bla bla

#If you want to keep the data frame as it was before, you must assign to a new object (make it descriptive, but short)

# df_1 <- df %>% bla bla bla

#The way that you can find NAs is by using is.na()
#If you want to find all row inputs of a column that are NOT NA, then you can use !is.na()
#The '!' is useful in many cases to do the opposite of said call, e.g.:
  #!= is not equal to something
  #!is.na()
  #!duplicated()
  #!(df$colname %in% colset)


# exclude Student_IDs that have no assigned name, not useful!
dirty_csv <- dirty_csv %>% filter(!is.na(Student_ID))


# change the student names as we have learned in previous steps:
dirty_csv$Name[dirty_csv$Name %in% c("laura", "selma", "ELISA", "max")] <- c("Laura", "Selma", "Elisa", "Max")

dirty_csv



### Gender #####################################################################

#correct name pattern (F and M for all inputs)
#correct Gender of Priscilla (f)
#wrong input (tennis)


#To change the gender accordingly, we will use the case_when() function:

#This is a function that allows you to adjust/make columns according to different conditions.

    # EXAMPLE: df <- df %>% mutate(colname = casewhen(colname == 'old_pattern1' ~ 'new_pattern1',
    # EXAMPLE:                                        colname == 'old_pattern2' ~ 'new_pattern2'))
                                           
dirty_csv <- dirty_csv %>% mutate(Gender = case_when(Gender == 'f' ~ 'f',
                                                     Gender == 'FEMALE' ~ 'f',
                                                     Gender == 'fem' ~ 'f',
                                                     is.na(Gender) ~ 'f',
                                                     Gender == 'm' ~ 'm',
                                                     Gender == 'mal' ~ 'm',
                                                     Gender == 'male' ~ 'm',
                                                     Gender == 'tennis' ~ 'm'))
# don't forget that the last input (tennis) was mistakenly in this column and should be added to Fav_Sport later!
  
dirty_csv$Gender[ dirty_csv$Name %in% 'Priscilla'] <- 'f'

dirty_csv


### Age ########################################################################

#numbers out of range (1000)
#chr instead of number (ten)
#change class to numeric

dirty_csv$Age[dirty_csv$Age %in% 'ten'] <- 10
dirty_csv$Age <- as.numeric(dirty_csv$Age)
dirty_csv <- dirty_csv %>% mutate(Age = ifelse(Age > 20,  # logical expression
                                               Age / 100, # what happens if the logical condition is true
                                               Age))      # what happens if the logical condition is false

dirty_csv



### Fav_Fruit ##################################################################

#wrong input (soccer)
dirty_csv$Fav_Fruit[dirty_csv$Fav_Fruit %in% "soccer"] <- "strawberry" # don't forget to change in the Fav_Sport column
dirty_csv



### Fav_Sport ##################################################################

#wrong input (MALE, strawberry)
dirty_csv$Fav_Sport[dirty_csv$Fav_Sport %in% c("MALE", "strawberry")] <- c("tennis", "soccer")
dirty_csv



### Do_you_like_Math ###########################################################

#Here we want to change unspecific answers --> make them definitive yes OR no
#We will look at the grep() function, which allows you to 'grab' a pattern and apply changes to all rows that follow it. 
#This is more generalized and is better for columns that have few options already (male/female, yes/no). 
#- casewhen() gives you more freedom if the row input is super diverse and you have to correct a lot. 
#- grep() pattern search is best combined with '*' as a 'wild card' ( = find stuff related to my pattern)

# EXAMPLE: df $ colname [ grep ("old_pattern*.", df$colname) ] <- "new_pattern"

dirty_csv$Do_you_like_Math[grep("yes*", dirty_csv$Do_you_like_Math)] <- "yes"
dirty_csv$Do_you_like_Math[grep("no*",  dirty_csv$Do_you_like_Math)] <- "no"

dirty_csv


## JOINING TABLES --------------------------------------------------------------

#There are different joins you can perform. 
#First of all, let's rename our prior data frame and load the other table we want to join with!
## let's rename our first data frame --> Student_Data1

Student_Data1 <- dirty_csv
# let's add the second data set
Student_Data2 <- read.csv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/dirty_files/Dirty_Data2.csv', na.strings = c('', ' ', 'NA'))
Student_Data2


#For simplicity, this table is clean already, we just want to look what happens when you join with different functions.
    #1. left_join() = join matching rows from b to a (includes all rows in a)
    #2. right_join() = join matching rows from a to b (includes all rows in b)
    #3. inner_join() = join data and retain only rows in both sets (includes all rows in a AND b)
    #4. full_join() = join data, retain all values, all rows (includes all rows in a OR b)


library(knitr)

# left join
Student_Data <- left_join(Student_Data1, Student_Data2) %>% arrange(Student_ID)
Student_Data
vis_miss(Student_Data)
# right join
Student_Data <- right_join(Student_Data1, Student_Data2) %>% arrange(Student_ID)
Student_Data
vis_miss(Student_Data)
# inner join
Student_Data <- inner_join(Student_Data1, Student_Data2) %>% arrange(Student_ID)
Student_Data
vis_miss(Student_Data)
# full join
Student_Data <- full_join(Student_Data1, Student_Data2) %>% arrange(Student_ID)
Student_Data
vis_miss(Student_Data)




## ADD THE MISSING DATA --------------------------------------------------------

#There are two csv files that contain the remaining, missing data:
#  - Student_Data3
#  - Student_Data4

Student_Data3 <- read.csv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/dirty_files/Dirty_Data3.csv', na.strings = c('', ' ', 'NA'))
Student_Data4 <- read.csv('https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/dirty_files/Dirty_Data4.csv', na.strings = c('', ' ', 'NA'))


# let's use full join
Student_Data <- full_join(Student_Data, Student_Data3) %>% arrange(Student_ID)
Student_Data <- full_join(Student_Data, Student_Data4) %>% arrange(Student_ID)
rm(Student_Data1, Student_Data2, Student_Data3, Student_Data4)

Student_Data
vis_miss(Student_Data)



#Let's see how we can get rid of the duplicate Student_IDs and fill the missing data
  #1. arrange by Student_ID with arrange()
  #2. group by Student_ID with group_by()
  #3. fill the empty rows with data from above/below with fill(c(everything()), .direction = 'downup')
  #4. ungroup with ungroup()
  #5. retain the distinct rows with distinct() and we want to keep all rows

Student_Data <- Student_Data %>% 
  arrange(Student_ID) %>% 
  group_by(Student_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>%  
  distinct(Student_ID, .keep_all = T)
 


#Let's look at the the table.
#Do we still have missing data?

vis_miss(Student_Data)
 

## CLEAN UP THE DUPLICATE COLUMNS ----------------------------------------------

#Fav_Animal = Fav_ANIMAL


#Fav_Animal is the correct name of the column, so we want to end up with that column name and ALL the data in that column. 
#Currently, we have two partially filled columns - how can we solve that?
  
#1. Let's start with our Student_Data
#2. Select the column we use for identification (Student_ID) and the columns we want to merge
#3. we use the pivot_longer() function to create a looooong table:
    #one column is Student_ID
    #the second column is 'Temp' and contains the previous names of the columns we want data from ('Fav_Animal' and 'Fav_ANIMAL')
    #the third column contains the actual data and that's the column we refer our values to (values_to = 'Fav_Animal')

#4. we now do the same thing as before, we look at duplicate Student_ID rows, fill the blanks with data from above/below, and remove complete duplicates
#5. We now have ONE column with ALL our data
#6. Now we only merge with the Student_Data (deselect the unneeded column Fav_ANIMAL) and have everything in our data frame

Temp_df <- Student_Data %>% select(Student_ID, Fav_Animal, Fav_ANIMAL) %>% 
  pivot_longer(names_to = "Temp", values_to = "Fav_Animal", cols = c(Fav_Animal, Fav_ANIMAL)) %>% 
  arrange(Student_ID) %>% 
  group_by(Student_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  select(Student_ID, Fav_Animal) %>% 
  distinct(Student_ID, .keep_all = T)

Temp_df

Student_Data <- full_join(Student_Data, Temp_df) %>% 
  select(-c(Fav_ANIMAL)) %>% 
  arrange(Student_ID) %>% 
  group_by(Student_ID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  distinct(Student_ID, .keep_all = T)
rm(Temp_df)

#Any blanks left?
vis_miss(Student_Data)




## CALCULATE VALUES FROM EXISTING COLUMNS AND MUTATE NEW COLUMNS ---------------

#calculate the mean of the students grades (Performance_Mean)

    # Btw, you can combine the case_when() function with different conditionals:
    # <   smaller than  
    # >   bigger than
    # <=  smaller or equal to
    # >=  bigger or equal to
    # ==  equal to
    # !=  not equal to
    # &   logical 'AND' condition
    # |   logical 'OR'  condition

Student_Data <- Student_Data %>% mutate(
  Performance_Mean = (Grade_English + Grade_Math + Grade_Biology)/3,
  Student_Eval = case_when(Performance_Mean <  1.5 ~ 'excellent',
                           Performance_Mean >= 1.5 & Performance_Mean < 2  ~ 'very good',
                           Performance_Mean >= 2   & Performance_Mean < 3  ~ 'good',
                           Performance_Mean >= 3   & Performance_Mean < 4  ~ 'average',
                           Performance_Mean >= 4   ~ 'bad'))


# DATA VISUALIZATION WITH GGPLOT2 ----------------------------------------------

#Data Visualization will be done using two main packages 'ggplot2' and 'leaflet'.
#ggplot2 is used for general visualization like dot plots, bar plots, histograms, linear regression, ..

#1. Start with ggplot and give the data frame you want to visualize from
#2. define the aesthetics using aes()
    # aes() wants to know what you want on the x-axis and on the y-axis
    # here you can also define shape, color, alpha (color transparency)
#3. after that you can add different geoms (plot types) using + geom_xxx(), e.g:
    #geom_point()    = dot plot, compare two variables x and y
    #geom_bar()      = bar plot, takes only x or y
    #geom_smooth()   = regression line (method = 'lm' will give you the linear regression fit, a straight line, without 'lm' it'll go through all points + sd ribbon)
    #geom_abline()   = line with x = y
    #geom_boxplot()  = boxplot
    #geom_violin()   = sometimes better than boxplot (shows median, ranges, variability effectively, allows comparison of groups of different sizes) 
    # ...
  
#More info on violin plots you can find [here](https://blog.bioturing.com/2018/05/16/5-reasons-you-should-use-a-violin-graph/)
#A more in-depth introduction to data visualization with ggplot is [here](https://r4ds.had.co.nz/data-visualisation.html)


## IRIS DATA -------------------------------------------------------------------

#We will start by looking at the 'iris' dataset provided by R
#It has 5 columns and 150 rows
    #Sepal.Length
    #Sepal.Width
    #Petal.Length
    #Petal.Width
    #Species

glimpse(iris)

# Dot Plot #####################################################################

# compare Sepal.Length and Sepal.Width, colored by Species
# the regression lines are per species, no standard deviation ribbon
iris %>%
  ggplot(aes(Sepal.Length, Sepal.Width, col = Species)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

plot(iris$Sepal.Length, iris$Sepal.Width) # base functions are always an alternative!


# violin plot & boxplot together ###############################################
iris %>%
  ggplot(aes(Species, Petal.Width)) +
  geom_violin(aes(col = Species, fill = Species, alpha = 0.3)) +
  geom_boxplot(aes(alpha = 0)) +
  ggtitle('Iris Boxplot + Violin plot') +
  xlab('Species') +
  ylab('Petal Width')



# Histogram ####################################################################
iris %>%
  ggplot(aes(Sepal.Width)) +
  geom_histogram(bins = 12)

hist(iris$Sepal.Width) # base functions are always an alternative!



## ANSCOMBE AND CORRELATION ----------------------------------------------------

#Anscombe (1973) has a nice example where he uses a constructed dataset to emphasize 
#the importance of using graphs in statistical analysis. 
#There are 8 variables, representing four pairings of an outcome and a predictor. 
#All sets have 11 observations, the same mean of x (9) and y (7.5), the same fitted 
#regression line (y = 3 + 0.5 x), the same regression and residual sum of squares and 
#therefore the same multiple R-squared of 0.67.
#If you want to learn more about Correlation and Regression, there is a great website 
#with different datasets and problem sets based on different regressions [here](https://data.princeton.edu/wws509/datasets).
#If you are interested in the Datacamp course on Correlation and Regression, 
#there is a free course available  [here](https://stat-ata-asu.github.io/CorrelationRegression/)


### Look at the different sets #################################################

data <- datasets::anscombe # already a base R data set!

# we have to prepare this data set a tiny bit to separate the data into sets
Anscombe <- data.frame(
  set = rep(1:4, each = 11),
  x = unlist(data[ ,c(1:4)]),
  y = unlist(data[ ,c(5:8)])
)
rownames(Anscombe) <- NULL
head(Anscombe)

Anscombe %>%
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  facet_wrap(~set) +
  theme_gray() # you can adjust themes as well


### Calculate the number of observations, means of x/y, Standard deviations of x/y
Anscombe %>%
  group_by(set) %>%
  summarize(N = n(),
            mean_x = mean(x),
            sd_x = sd(x),
            mean_y = mean(y),
            sd_y = sd(y),
            cor_between_x_and_y = cor(x,y))

# Same same but different!
# SO: always remember that visual analysis is important!



## SOTA VISUALIZATION - HMHZ Wild Samples --------------------------------------

#Let's do some data visualization using the actual SOTA Data product now.
#The data frame is pretty big, it currently has 117 columns and 1921 rows.
#We want to look at visualization of some data (I will take advantage of this and show you Crypto-related stuff)

### Bar Plot (geom_col()) of Crypto Prevalence since 2016

SOTA <- read.csv('https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv')
SOTA <- SOTA[SOTA$Mouse_ID %like% "AA_", ] # we want Brandenburg HMHZ samples only

# we want to compare the catching rate overall to the mice that were Crypto-positive
All_Samples <- SOTA %>% group_by(Year) %>% count()
Pos_Samples <- SOTA %>% filter(ILWE_Crypto_Ct > 0, Year >= 2016) %>% group_by(Year) %>% count()

# let's combine both and calculate the prevalence
Samples_Yr <- full_join(Pos_Samples, All_Samples, by = "Year") %>% 
  mutate(Prevalence = (n.x / n.y)) %>% 
  filter(!is.na(Prevalence))

# let's visualize that with a bar plot (geom_col is like geom_bar, but can take both x AND y)
Samples_Yr %>%
  ggplot(aes(x = Year, label = Prevalence)) +
  geom_col(aes(y = n.y, fill = "blue")) +
  geom_col(aes(y = n.x, fill = "red")) +
  geom_text(aes(label = percent(Prevalence),
                y = (n.x / n.y)), nudge_y = 9) +
  labs(y = "Samples [n]") +
  theme(legend.position = "none") +
  ggtitle("Cryptosporidium spp. prevalence in the HMHZ since 2016")                                       
                                       
                                
### Cryptosporidium Oocyst Prediction ##########################################

#To predict oocysts (or better oocyst equivalents, as we technically quantify oocyst DNA), we need a couple of things:
#(This is predicted not on SOTA, but on my Crypto Data product and joined with SOTA later)

#1. A Standard Curve with determined quantities of oocyst DNA
#2. A linear regression model
#3. A Prediction function
#3. A correction step for Ct = 0 (ultra high value, actually 0, just to differentiate nothing measured (0) from NAs)


# 1 Load the Standard_Curve
Standard_Curve <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Cryptosporidium/Crypto_Standard_Curve.csv")
colnames(Standard_Curve)[colnames(Standard_Curve)%in%"Ct_mean"] <- "ILWE_Crypto_Ct" # rename Ct_mean
Standard_Curve     <-  filter(Standard_Curve, ILWE_Crypto_Ct > 0) # only values above > 0

Crypto_Detection <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv") %>% select(-X)

# 2 Design of the linear Model 
# I had a 2-fold dilution series, therefore I set it up like that, usually it's a 10-fold dilution
Linear_Model     <- lm(log2(Amount_Oocysts) ~ ILWE_Crypto_Ct, data = Standard_Curve)

# 3 Prediction, depending on your dilution steps
Oocyst_Predict_Crypto    <- 2^predict(Linear_Model, newdata = Crypto_Detection)
Crypto_Detection <- data.frame(Crypto_Detection, Oocyst_Predict_Crypto) # add the Oocyst_Predict_Crypto vector to our existing data frame
# replace the value predicted for Ct = 0
Crypto_Detection <- Crypto_Detection %>% mutate(Oocyst_Predict_Crypto.1 = replace(Oocyst_Predict_Crypto.1, Oocyst_Predict_Crypto.1 == "4292821751815.77", "0"))
# make sure that the prediction is of the correct class (we want integers, not doubles)
Crypto_Detection$Oocyst_Predict_Crypto.1 <- as.integer(Crypto_Detection$Oocyst_Predict_Crypto.1)

Crypto_Detection %>% arrange(desc(Oocyst_Predict_Crypto.1)) %>% select(Mouse_ID, Year, HI, Oocyst_Predict_Crypto.1) %>% head(10)



# DATA VISUALIZATION WITH LEAFLET - MAPS ---------------------------------------

## COLOR PALETTES ##############################################################

#There are great color palettes already available in R, 
#[here](https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf) 
#is a summary of the most useful ones.

#I tend to use the RColorBrewer package a lot, as it offers different kinds of palettes:
#   - Sequential (e.g. 'YlOrRd', 'Oranges', 'Greens', 'Blues', ...)
#   - Qualitative (e.g. 'Set2', 'Paired', 'Dark2', 'Accent', ...)
#   - Diverging (e.g. 'Spectral', 'RdYlBl', ...)

#To check out the different options of the package, type display.brewer.all() after loading library(RColorBrewer)
#To use a palette, you have to add the brewer palette to your graphic (e.g. in a ggplot), as follows:

library(RColorBrewer) #display.brewer.all()

iris %>%
  ggplot(aes(Sepal.Length, Sepal.Width, col = Species)) +
  geom_point() +
  scale_color_brewer(palette = 'Dark2',
                     direction = 1) # direction = -1 would be inverted color order



#Just for the sake of showing it to you:
#You can also set up your own palettes 
#I did this to get the HI gradient we have used for quite some time in publications.

#1. You can find the RGB code for each color on different sites, e.g. [here](https://rgbcolorcode.com/)
#2. Separate each color into the according vector r,g and b (as many colors as you wish, but 3-8 should usually do the trick)
#3. follow the setup code
#4. Btw, you can use HEX codes (e.g. #3399FF = 'dodger blue') instead of color names for almost everything! So you're not stuck with 'red', 'green', 'yellow', ..
                               
                               
# Load Palette 
r <- c(0, 64, 128, 179, 217, 255)
g <- c(0, 12, 25, 25, 12,  0)
b <- c(255, 249, 243, 191,  95,   0)

myPal <- function (n, name = c("myPal.colors")) 
{
  myPal.colors = rgb(r,g,b,maxColorValue = 255)
  name = match.arg(name)
  orig = eval(parse(text = name))
  rgb = t(col2rgb(orig))
  temp = matrix(NA, ncol = 3, nrow = n)
  x = seq(0, 1, , length(orig))
  xg = seq(0, 1, , n)
  for (k in 1:3) {
    hold = spline(x, rgb[, k], n = n)$y
    hold[hold < 0] = 0
    hold[hold > 255] = 255
    temp[, k] = round(hold)
  }
  palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
  palette
}
                               
                                                                    
                                                                    
                                                                    
## HYBRID INDEX MAP ------------------------------------------------------------
                                                                    
#Let's look at mapping with leaflet now.
#This is a fantastic tool to get a geographic overview of your data and include specific information using color codes and legends.
#Leaflet maps are interactive and you can load the SOTA coordinates of our mice quite easily.
#We need a couple of packages to make this work, but again - just an overview of possibilities for you!
#To check out other maps from third-party providers, type 'addProviderTiles(providers$)' and push tab after the dollar sign.
#You can add different markers, I used *circle markers* that allow popup information, but you may want to use other ones.



library(RColorBrewer) #display.brewer.all()
library(leaflet)
library(htmltools)

SOTA <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv")  %>% filter(!is.na(Latitude), !is.na(Longitude), !is.na(HI), !Mouse_ID %in% c("SK_2697", "SK_3173", "SK_905"))
Crypto_Detection <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv") 
Crypto_Detection_tested <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv")  %>% filter(ILWE_Crypto_Ct >= 0, !is.na(Latitude), !is.na(Longitude))
Crypto_Positive <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv")  %>% filter(ILWE_Crypto_Ct > 0, !is.na(Latitude), !is.na(Longitude))
Crypto_Detection_21 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv")  %>% filter(Year == 2021, ILWE_Crypto_Ct > 0, !is.na(Latitude), !is.na(Longitude))

High_Infection_Samples <- SOTA %>% 
  select(Mouse_ID, ILWE_Crypto_Ct, Oocyst_Predict_Crypto, Year, Latitude, Longitude, HI, Sex) %>% 
  filter(ILWE_Crypto_Ct > 0,
         Year >= 2016) %>% 
  arrange(ILWE_Crypto_Ct) %>% 
  head(20)

#- EIM_INFECTED         ... This location had samples with Eimeria Infection, radius == Infection intensity (delta_ct_cewe_MminusE)
SOTA <- SOTA %>% 
  mutate(Eim_Species = ifelse(eimeriaSpecies == "E_falciformis", "E_falciformis",
                              ifelse(eimeriaSpecies == "E_ferrisi", "E_ferrisi",
                                     ifelse(eimeriaSpecies == "Eimeria_alorani", "Eimeria_alorani",
                                            ifelse(eimeriaSpecies == "Eimeria_apionodes", "Eimeria_apionodes",
                                                   ifelse(eimeriaSpecies == "Eimeria_falciformis", "Eimeria_falciformis",
                                                          ifelse(eimeriaSpecies == "Eimeria_sp_Apodemus", "Eimeria_sp_Apodemus",
                                                                 ifelse(eimeriaSpecies == "Eimeria_vermiformis", "Eimeria_vermiformis",
                                                                        ifelse(eimeriaSpecies == "Negative", "Negative",
                                                                               ifelse(NA))))))))),
         Eimeria_Positive = case_when(Eim_Species != "Negative" | Ct.Eimeria > 0 | OPG > 0 ~ T,
                                      Eim_Species == "Negative" | Ct.Eimeria == 0 | OPG == 0 ~ F))

# We want to look at Eim_infected samples in our map, so we extract a dataframe with samples that fulfill this condition
EIM_INFECTED <- SOTA %>% filter(Eimeria_Positive == T)

# Cutting the data into slices helps to make nice legends!
SOTA$HI <- as.numeric(SOTA$HI)
SOTA$HI_Level <-  cut(SOTA$HI, c(0, 0.001, 0.250, 0.500, 0.750, 0.999, 1), include.lowest = T , 
                      labels = c('HI = 0', 'HI < 0.25', 'HI < 0.5', 'HI < 0.75', 'HI < 1', 'HI = 1'))

# We also want to look at samples that sent into sequencing (was interesting for my Bachelor thesis)
Sequenced <- Crypto_Detection %>%
  mutate(seq = Mouse_ID %in% c("AA_0144", "AA_0325", "AA_0689", "AA_0209", "AA_0282", "AA_0793", "AA_0667", "AA_0805", "AA_0900",
                               "AA_0523", "AA_0534", "AA_0537", "AA_0545", "AA_0553", "AA_0554", "AA_0578", "AA_0580", "AA_0585",
                               "AA_0546", "AA_0589", "AA_0571", "AA_0667")) %>%
  filter(seq == T) # only retain samples that fulfill this condition


data_col_HI       = colorFactor(myPal(6), SOTA$HI) # for our samples, uncut
data_col_HI_Level = colorFactor(myPal(6), SOTA$HI_Level) # for our legend palette, cut

# We must provide a map to which we anchor the data.
# Check out other third party tile providers using: addProviderTiles(providers$)
# CartoDB worked best for the colors I used
map <- Crypto_Detection_21 %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%  
  setView(lat = 52.520007, lng =13.404954, zoom = 7)

map %>%
  addPolylines(lat = c(55.0000, 53.6000, 53.51885, 52.8875  , 52.6053, 51.8978, 45.0000), 
               lng = c(10.0000, 11.4563, 12.4464,13.8119 , 13.8756, 13.8103, 13.5000), 
               color = "purple", 
               weight = 55, 
               opacity = 0.1) %>%
  addCircleMarkers(data = EIM_INFECTED,
                   col = "#72cac3",
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 0.3,
                   radius = sqrt((EIM_INFECTED$delta_ct_cewe_MminusE)^2),
                   group = "Eim_Infected") %>%
  addCircleMarkers(data = SOTA,
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = "Samples (total)") %>%
  addCircleMarkers(data = Crypto_Detection_tested,
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = "Samples (Crypto-tested)") %>%
  addCircleMarkers(data = Crypto_Positive,
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = "Samples (Crypto-positive)") %>%
  addCircleMarkers(data = Sequenced, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 5,
                   radius = 3,
                   group = "Sequenced") %>%
  addCircleMarkers(data = High_Infection_Samples, 
                   color = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 3,
                   radius = 3,
                   group = "High_Infection_Samples") %>%
  addLegend("bottomleft", 
            pal = data_col_HI_Level, 
            title = "HI",
            values = SOTA$HI_Level, 
            group = c('= 0', '< 0.25', '< 0.5', '< 0.75', '< 1', '= 1'),
            opacity = 1) %>%
  addLayersControl(baseGroups = c("Samples (total)", "Samples (Crypto-tested)", "Samples (Crypto-positive)", "High_Infection_Samples", "Sequenced"), 
                   overlayGroups = c("Eim_Infected"),
                   options = layersControlOptions(collapsed = F))


# SAVING PICTURES AND WRITING FILES WITH R

#Let's say we want to:
                                                                      
#1. write a csv file with all the high infection samples
#2. save the picture of the Standard Curve Plot

#There are different functions for saving different types - as usual, I would recommend saving data frames as csv files
#For saving pictures, there are probably even more possibilities:
#  - My personal favorite for high resolution pictures: 
#  - Load library(cowplot), assign your plot to an object, e.g. my_Plot <- ...
#follow the saving path: 
  
  
save_plot(my_Plot, filename = 'my_Plot.jpg', base_height = 5, base_width = 7.5)

#change the base_height and _width according to your needs (plots with facet_wrap() tend to get quite unreadable if exported too small), I would go up to base_height 20 and base_weight 30 if you want really high resolution
#beware, that the font size WILL SHRINK dramatically with increasing base_height/weight conditions!
#export the picture through the *viewer* window (export let's you either export as Image, PDF or copy to clipboard)
#save the picture with the jpeg() function (or similar), which allows you to set the quality of the image as a percentage (I think this works best/better in normal R scripts, not in Markdown files :/ )

# JPEG device
jpeg('Crypto_SC.jpeg', quality = 75)
plot(Crypto_SC)
dev.off()
    

library(cowplot)

write.csv(High_Infection_Samples, 'High_Infection_Samples.csv') # don't forget to add the type of file, in this case '.csv' at the end (otherwise it will be a binary, non-readable file for R)
                                                                                                                        
Crypto_SC <- Standard_Curve %>%
  ggplot(aes(log2(Amount_Oocysts), ILWE_Crypto_Ct)) +
  geom_smooth(method = 'lm') +
  geom_point() +
  xlab('Amount Oocysts (log2)') +
  ylab('Ct') +
  ggtitle('Cryptosporidium Standard Curve')
Crypto_SC

save_plot(Crypto_SC, filename = 'Crypto_SC.jpg', base_height = 5, base_width = 7.5) 



# PICTURES IN MARKDOWN

#To display a picture in markdown, you should follow this pattern:

#![Caption](https://raw.githubusercontent.com/tlobnow/Resources/main/R_Introduction/Crypto_SC.jpg)

                                                                                                                        




