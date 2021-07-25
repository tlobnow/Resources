library(ggplot2)
library(dplyr)
library(openintro)
library(usdata)
library(XLConnect)
library(assertive)
library(stringr)
library(lubridate)
library(visdat)
library(stringdist)
library(fuzzyjoin)
library(reclin)

setwd("/Users/FinnLo/Documents/Programming/R/DataCamp/DataCamp")


# INTRODUCTION TO IMPORTING DATA IN R ----
# READING .CSV FILES INTO R
# Import swimming_pools.csv: pools
pools <- read.csv("swimming_pools.csv")

# Print the structure of pools
str(pools)

# STRINGS AS FACTORS
# Import swimming_pools.csv correctly: pools
pools <- read.csv("swimming_pools.csv", stringsAsFactors = FALSE)

# Check the structure of pools
str(pools)

# READ .DELIM OR .TABLE FILES
# READ .TXT FILES
# Import hotdogs.txt: hotdogs
hotdogs <- read.delim("hotdogs.txt",
                      header = FALSE)

# Summarize hotdogs
summary(hotdogs)

#Path to the hotdogs.txt file: path
path <- file.path("data", "hotdogs.txt")

# Import the hotdogs.txt file: hotdogs
hotdogs <- read.table(path,
                      sep = "\t", 
                      header = FALSE,
                      col.names = c("type", "calories", "sodium"))
# Call head() on hotdogs
head(hotdogs)

# IMPORT AND SELECT
# Finish the read.delim() call
hotdogs <- read.delim("hotdogs.txt", header = FALSE, col.names = c("type", "calories", "sodium"))

# Select the hot dog with the least calories: lily
lily <- hotdogs[which.min(hotdogs$calories), ]

# Select the observation with the most sodium: tom
tom <- hotdogs[which.max(hotdogs$sodium), ]


# Print lily and tom
print(lily)
print(tom)

# COL.CLASSES() --> SPECIFY COLUMN TYPES == COLUMN CLASSES
# Previous call to import hotdogs.txt
hotdogs <- read.delim("hotdogs.txt", header = FALSE, col.names = c("type", "calories", "sodium"))

# Display structure of hotdogs
str(hotdogs)
# Edit the colClasses argument to import the data correctly: hotdogs2
hotdogs2 <- read.delim("hotdogs.txt", 
                       header = FALSE, 
                       col.names = c("type", "calories", "sodium"),
                       colClasses = c("factor", "NULL", "numeric"))
# Display structure of hotdogs2
str(hotdogs2)

# DIFFERENT READ FUNCTIONS FOR REGIONAL DIFFERENCES
  # read.table = main function
  # read.csv / .delim / .csv2 / .delim2 are wrappers around main function
  # .csv2 / .delim2 are important for regional differences == when is . / , / ; as separator used? eg. in numbers can throw table off!

# USE PACKAGES TO IMPORT DATA 
#   READR, READ_CSV:

# Load the readr package
library(readr)

# Import potatoes.csv with read_csv(): potatoes
potatoes <- read_csv("potatoes.csv")

# READ_TSV
# Column names
properties <- c("area", "temp", "size", "storage", "method","texture", "flavor", "moistness")
# Import potatoes.txt: potatoes
potatoes <- read_tsv("potatoes.txt", col_names = properties)
# Call head() on potatoes
head(potatoes)

# READ_DELIM
# Column names
properties <- c("area", "temp", "size", "storage", "method", "texture", "flavor", "moistness")
# Import potatoes.txt using read_delim(): potatoes
potatoes <- read_delim("potatoes.txt",
                       delim = "\t",
                       col_names = FALSE)
# Print out potatoes
print(potatoes)

# skip specifies the number of lines you're ignoring in the flat file before actually starting to import data.
#n_max specifies the number of lines you're actually importing.
# Say for example you have a CSV file with 20 lines, and set skip = 2 and n_max = 3, you're only reading in lines 3, 4 and 5 of the file.
# Watch out: Once you skip some lines, you also skip the first line that can contain column names!
# Finish the first read_tsv() call to import observations 7, 8, 9, 10 and 11 from potatoes.txt: 
# readr is already loaded
# Column names
properties <- c("area", "temp", "size", "storage", "method", "texture", "flavor", "moistness")
# Import 5 observations from potatoes.txt: potatoes_fragment ()
potatoes_fragment <- read_tsv("potatoes.txt", skip = 6, n_max = 5, col_names = properties)

# COL TYPES
# You can also specify which types the columns in your imported data frame should have. 
# You can do this with col_types. If set to NULL, the default, functions from the readr package will try to find the correct types themselves. 
# You can manually set the types with a string, where each character denotes the class of the column: character, double, integer and logical. _ skips the column as a whole.
# readr is already loaded

# Column names
properties <- c("area", "temp", "size", "storage", "method", "texture", "flavor", "moistness")

# Import all data, but force all columns to be character: potatoes_char
potatoes_char <- read_tsv("potatoes.txt", col_types = "cccccccc", col_names = properties)

# Print out structure of potatoes_char
str(potatoes_char)

# COL_TYPES WITH COLLECTORS
# readr is already loaded

# Import without col_types
hotdogs <- read_tsv("hotdogs.txt", col_names = c("type", "calories", "sodium"))

# Display the summary of hotdogs
summary(hotdogs)

# The collectors you will need to import the data
fac <- col_factor(levels = c("Beef", "Meat", "Poultry"))
int <- col_integer()

# Edit the col_types argument to import the data correctly: hotdogs_factor
hotdogs_factor <- read_tsv("hotdogs.txt", col_names = c("type", "calories", "sodium"),
                           col_types = list(fac, int, int))

# Display the summary of hotdogs_factor
summary(hotdogs_factor)


# DATA.TABLE PACKAGE --> FREAD:
# infer column types & separators
# simply works
# extremely fast (good for large data)
# it's possible to manually specify lots of stuff
# load the data.table package using library()
library(data.table)

# Import potatoes.csv with fread(): potatoes
potatoes <- fread("potatoes.csv")

# Print out potatoes
print(potatoes)

# ADVANCED USE OF FREAD
# different ways to import different variables:
#   fread("path/to/file.txt", drop = 2:4)
#   fread("path/to/file.txt", select = c(1, 5))
#   fread("path/to/file.txt", drop = c("b", "c", "d"))
#   fread("path/to/file.txt", select = c("a", "e"))

# Import columns 6 and 8 of potatoes.csv: potatoes
potatoes <- fread("potatoes.csv", select = c(6,8))

# Plot texture (x) and moistness (y) of potatoes
plot(x = potatoes$texture, y = potatoes$moistness)

# READ EXCEL FILES --> READXL
# # Load the readxl package
library(readxl)

# type dir() to find your working directory
# load the file and print the names of all worksheets
excel_sheets("urbanpop.xlsx")


# IMPORT AN EXCEL SHEET
# The readxl package is already loaded

# Read the sheets, one by one
pop_1 <- read_excel("urbanpop.xlsx", sheet = 1)
pop_2 <- read_excel("urbanpop.xlsx", sheet = 2)
pop_3 <- read_excel("urbanpop.xlsx", sheet = 3)

# Put pop_1, pop_2 and pop_3 in a list: pop_list
pop_list <- list(pop_1, pop_2, pop_3)

# Display the structure of pop_list
str(pop_list)

# READING A WORKBOOK:
# In the previous exercise you generated a list of three Excel sheets that you imported. However, loading in every sheet manually and then merging them in a list can be quite tedious. Luckily, you can automate this with lapply().
# example:
#   my_workbook <- lapply(excel_sheets("data.xlsx"),
#                   read_excel,
#                   path = "data.xlsx")

# The readxl package is already loaded

# Read all Excel sheets with lapply(): pop_list
pop_list <- lapply(excel_sheets("urbanpop.xlsx"),
                   read_excel, 
                   path = "urbanpop.xlsx")
# Display the structure of pop_list
str(pop_list)

# THE COL_NAMES ARGUMENT
# Apart from path and sheet, there are several other arguments you can specify in read_excel(). One of these arguments is called col_names.
# By default it is TRUE, denoting whether the first row in the Excel sheets contains the column names. If this is not the case, you can set col_names to FALSE. In this case, R will choose column names for you. You can also choose to set col_names to a character vector with names for each column. It works exactly the same as in the readr package.
# You'll be working with the urbanpop_nonames.xlsx file. It contains the same data as urbanpop.xlsx but has no column names in the first row of the excel sheets.

# Import the first Excel sheet of urbanpop_nonames.xlsx (R gives names): pop_a
pop_a <- read_excel("urbanpop_nonames.xlsx",
                    sheet = 1,
                    col_names = FALSE)

# Import the first Excel sheet of urbanpop_nonames.xlsx (specify col_names): pop_b
pop_b <- read_excel("urbanpop_nonames.xlsx",
                    sheet = 1,
                    col_names = cols)
cols <- c("country", paste0("year_", 1960:1966))


# Print the summary of pop_a
summary(pop_a)

# Print the summary of pop_b
summary(pop_b)

# THE SKIP ARGUMENT
# Another argument that can be very useful when reading in Excel files that are less tidy, is skip. With skip, you can tell R to ignore a specified number of rows inside the Excel sheets you're trying to pull data from. Have a look at this example:
#   read_excel("data.xlsx", skip = 15)
# In this case, the first 15 rows in the first sheet of "data.xlsx" are ignored.
# If the first row of this sheet contained the column names, this information will also be ignored by readxl. Make sure to set col_names to FALSE or manually specify column names in this case!
# The file urbanpop.xlsx is available in your directory; it has column names in the first rows.

# Import the second sheet of urbanpop.xlsx, skipping the first 21 rows: urbanpop_sel
urbanpop_sel <- read_excel("urbanpop.xlsx", skip = 21,
                           sheet = 2,
                           col_names = FALSE)
# Print out the first observation from urbanpop_sel
head(urbanpop_sel, n = 1)

# GDATA FOR READING EXCEL
# Load the gdata package
library(gdata)
# Import the second sheet of urbanpop.xls: urban_pop
urban_pop <- read.xls("urbanpop.xls",
                      sheet = 2)
# Print the first 11 observations using head()
head(urban_pop, 11)

# READ.XLS() WRAPS AROUND READ.TABLE()
# Column names for urban_pop
columns <- c("country", paste0("year_", 1967:1974))
# Finish the read.xls call
urban_pop <- read.xls("urbanpop.xls", 
                      sheet = 2,
                      skip = 50, 
                      header = FALSE, 
                      stringsAsFactors = FALSE,
                      col.names = columns)
# Print first 10 observation of urban_pop
head(urban_pop, 10)

# CLEANING UP EXCEL DATA, USE CBIND(), REMOVE NA ROWS
# Add code to import data from all three sheets in urbanpop.xls
path <- "urbanpop.xls"
urban_sheet1 <- read.xls(path, sheet = 1, stringsAsFactors = FALSE)
urban_sheet2 <- read.xls(path, sheet = 2, stringsAsFactors = FALSE)
urban_sheet3 <- read.xls(path, sheet = 3, stringsAsFactors = FALSE)


# Extend the cbind() call to include urban_sheet3: urban
urban <- cbind(urban_sheet1, urban_sheet2[-1], urban_sheet3[-1])
# Remove all rows with NAs from urban: urban_clean
urban_clean <- na.omit(urban)
# Print out a summary of urban_clean
summary(urban_clean)

# READING DATA 
# XLConnect --> work with Excel through R


library("XLConnect")
library(rJava)

# load workbook <-loadWorkbook("file.xlsx")
# getSheets(book)
# readWorksheet(book, sheet = "...")

# urbanpop.xlsx is available in your working directory
# Load the XLConnect package
library("XLConnect")
# Build connection to urbanpop.xlsx: my_book
my_book <- loadWorkbook("urbanpop.xlsx")
# Print out the class of my_book
class(my_book)

# LIST AND READ EXCEL SHEETS
# XLConnect is already available
# Build connection to urbanpop.xlsx
my_book <- loadWorkbook("urbanpop.xlsx")
# List the sheets in my_book
getSheets(my_book)
# Import the second sheet in my_book
readWorksheet(my_book, sheet = 2)

# USE CBIND()
# XLConnect is already available
# Build connection to urbanpop.xlsx
my_book <- loadWorkbook("urbanpop.xlsx")
# Import columns 3, 4, and 5 from second sheet in my_book: urbanpop_sel
urbanpop_sel <- readWorksheet(my_book, sheet = 2, startCol = 3, endCol = 5)
# Import first column from second sheet in my_book: countries
countries <- readWorksheet(my_book, sheet = 2, startCol = 1, endCol = 1, header = TRUE)
# cbind() urbanpop_sel and countries together: selection
selection <- cbind(countries, urbanpop_sel)

# ADD A WORKSHEET TO AN EXISTING WORKBOOK:
# Build connection to urbanpop.xlsx
my_book <- loadWorkbook("urbanpop.xlsx")
# Add a worksheet to my_book, named "data_summary"
createSheet(my_book, "data_summary")
# Use getSheets() on my_book
getSheets(my_book)


# POPULATE A WORKSHEET (AKA ADD DATA FRAMES AND STORE THEM IN YOUR WORKBOOK)
# Build connection to urbanpop.xlsx
my_book <- loadWorkbook("urbanpop.xlsx")
# Add a worksheet to my_book, named "data_summary"
createSheet(my_book, "data_summary")
# Create data frame: summ
sheets <- getSheets(my_book)[1:3]
dims <- sapply(sheets, function(x) dim(readWorksheet(my_book, sheet = x)), USE.NAMES = FALSE)
summ <- data.frame(sheets = sheets,
                   nrows = dims[1, ],
                   ncols = dims[2, ])
# Add data in summ to "data_summary" sheet
writeWorksheet(my_book, summ, "data_summary")
# Save workbook as summary.xlsx
saveWorkbook(my_book, "summary.xlsx")

# RENAMING SHEETS
# Build connection to urbanpop.xlsx: my_book
my_book <- loadWorkbook("urbanpop.xlsx")
# Rename "data_summary" sheet to "summary"
renameSheet(my_book, "data_summary", "summary")
# Print out sheets of my_book
getSheets(my_book)
# Save workbook to "renamed.xlsx"
saveWorkbook(my_book, "renamed.xlsx")

# REMOVING SHEETS
# Load the XLConnect package
library(XLConnect)
# Build connection to renamed.xlsx: my_book
my_book <- loadWorkbook("renamed.xlsx")
# Remove the fourth sheet
removeSheet(my_book, sheet = 4)
# Save workbook to "clean.xlsx"
saveWorkbook(my_book, "clean.xlsx")




# INTRODUCTION TO DATA IN R ----

# Load data
data(email50)

# View the structure of the data
str(email50)

# Glimpse email50
glimpse(email50)

# Subset of emails with big numbers: email50_big
email50_big <- email50 %>%
  filter(number == "big")

# Glimpse the subset
glimpse(email50_big)

# Subset of emails with big numbers: email50_big
email50_big <- email50 %>%
  filter(number == "big")

# Table of the number variable
table(email50_big$number)

# Drop levels
email50_big$number_dropped <- droplevels(email50_big$number)

# Table of the number_dropped variable
table(email50_big$number_dropped)

# Calculate median number of characters: med_num_char
med_num_char <- median(email50$num_char)

# Create num_char_cat variable in email50
email50_fortified <- email50 %>%
  mutate(num_char_cat = ifelse(num_char < med_num_char, "below median", "at or above median"))

# Count emails in each category
email50_fortified %>%
  count(num_char_cat)

# Create number_yn column in email50
email50_fortified <- email50 %>%
  mutate(
    number_yn = case_when(
      # if number is "none", make number_yn "no"
      number == "none" ~ "no", 
      # if number is not "none", make number_yn "yes"
      number != "none" ~ "yes"  
    )
  )

# Visualize the distribution of number_yn
ggplot(email50_fortified, aes(x = number_yn)) +
  geom_bar()

# Load ggplot2
library(ggplot2)

# Scatterplot of exclaim_mess vs. num_char
ggplot(email50, aes(x = num_char, y = exclaim_mess, color = factor(spam))) +
  geom_point()

# Load data
library(gapminder)

# Glimpse data
glimpse(gapminder)

# Identify type of study: observational or experimental
type_of_study <- "observational"

# Load packages
library(dplyr)

# Count number of male and female applicants admitted
ucb_admit %>%
  count(Gender, Admit)

ucb_admission_counts %>%
  # Group by gender
  group_by(Gender) %>%
  # Create new variable
  mutate(prop = n / sum(n)) %>%
  # Filter for admitted
  filter(Admit == "Admitted")

glimpse(ucb_admission_counts)
ucb_admission_counts <- ucb_admit %>%
  # Counts by department, then gender, then admission status
  count(Dept, Gender, Admit)

# See the result
ucb_admission_counts

# Simple random sample: states_srs
states_srs <- us_regions %>%
  sample_n(size = 8)

# Count states by region
states_srs %>%
  count(region)

# Stratified sample
states_str <- us_regions %>%
  group_by(region) %>%
  sample_n(size = 2)

# Count states by region
states_str %>%
  count(region)

# BEAUTY IN THE CLASS ROOM
library(dplyr)
# Inspect evals

glimpse(evals)

# Inspect variable types
glimpse(evals)


# Remove non-factor variables from the vector below
cat_vars <- c("rank", "ethnicity", "gender", "language",
              "cls_level", "cls_profs", "cls_credits",
              "pic_outfit", "pic_color")

# Recode cls_students as cls_type
evals_fortified <- evals %>%
  mutate(
    cls_type = case_when(
      cls_students <= 18 ~ "small",
      cls_students >=19 & cls_students <= 59 ~ "midsize",
      cls_students >= 60 ~ "large"
    )
  )

# Scatterplot of score vs. bty_avg
ggplot(evals, aes(bty_avg, score)) +
  geom_point()

# Scatterplot of score vs. bty_avg colored by cls_type
ggplot(evals, aes(bty_avg, score, color = cls_type)) +
  geom_point()





# EXPLORATORY DATA ANALYSIS IN R ----

# CONTINGENCY TABLE REVIEW
# Print the first rows of the data
comics

# Check levels of align
levels(comics$align)

# Check the levels of gender
levels(comics$gender)

# Create a 2-way contingency table
table(comics$gender, comics$align)

# DROPPING LEVELS
# Load dplyr
library(dplyr)

# Print tab
print(tab)

# Remove align level
comics_filtered <- comics %>%
  filter(align != "Reformed Criminals") %>%
  droplevels()

# See the result
comics_filtered

# SIDE BY SIDE BARCHARTS
# Load ggplot2
library(ggplot2)

# Create side-by-side barchart of gender by alignment
ggplot(comics, aes(x = align, fill = gender)) + 
  geom_bar(position = "dodge")

# Create side-by-side barchart of alignment by gender
ggplot(comics, aes(x = gender, fill = align)) + 
  geom_bar(position = "dodge") +
  theme(axis.text.x = element_text(angle = 90))

# COUNTS VS. PROPORTIONS
# Plot of gender by align
ggplot(comics, aes(x = align, fill = gender)) +
  geom_bar()

# Plot proportion of gender, conditional on align
ggplot(comics, aes(x = align, fill = gender)) + 
  geom_bar(position = "fill") +
  ylab("proportion")

# MARGINAL BARCHART
# Change the order of the levels in align
comics$align <- factor(comics$align, 
                       levels = c("Bad", "Neutral", "Good"))

# Create plot of align
ggplot(comics, aes(x = align)) + 
  geom_bar()

# CONDITIONAL BARCHART
# Plot of alignment broken down by gender
ggplot(comics, aes(x = align)) + 
  geom_bar() +
  facet_wrap(~ gender)

# IMPROVE A PIECHART 
# Put levels of flavor in descending order
lev <- c("apple", "key lime", "boston creme", "blueberry", "cherry", "pumpkin", "strawberry")
pies$flavor <- factor(pies$flavor, levels = lev)

# Create barchart of flavor
ggplot(pies, aes(x = flavor)) + 
  geom_bar(fill = "chartreuse") + 
  theme(axis.text.x = element_text(angle = 90))

# FACETED HISTOGRAM
# Load package
library(ggplot2)

# Learn data structure
str(cars)

# Create faceted histogram
ggplot(cars, aes(x = city_mpg)) +
  geom_histogram() +
  facet_wrap(~ suv)

# BOXPLOTS AND DENSITY PLOTS
# Filter cars with 4, 6, 8 cylinders
common_cyl <- filter(cars, ncyl %in% c(4,6,8))

# Create box plots of city mpg by ncyl
ggplot(common_cyl, aes(x = as.factor(ncyl), y = city_mpg)) +
  geom_boxplot()

# Create overlaid density plots for same data
ggplot(common_cyl, aes(x = city_mpg, fill = as.factor(ncyl))) +
  geom_density(alpha = .3)

# MARGINAL AND CONDITIONAL HISTOGRAMS
# Create hist of horsepwr
cars %>%
  ggplot(aes(x = horsepwr)) +
  geom_histogram() +
  ggtitle("Distr of HPW for all cars")

# Create hist of horsepwr for affordable cars
cars %>% 
  filter(msrp < 25000) %>%
  ggplot(aes(x = horsepwr)) +
  geom_histogram() +
  xlim(c(90, 550)) +
  ggtitle("Distr affordable cars")

# THREE BINWIDTHS
# Create hist of horsepwr with binwidth of 3
cars %>%
  ggplot(aes(x = horsepwr)) +
  geom_histogram(binwidth = 3) +
  ggtitle("Histo Horsepower, bin 3")

# Create hist of horsepwr with binwidth of 30
cars %>%
  ggplot(aes(x = horsepwr)) +
  geom_histogram(binwidth = 30) +
  ggtitle("Histo Horsepower, bin 30")

# Create hist of horsepwr with binwidth of 60
cars %>%
  ggplot(aes(x = horsepwr)) +
  geom_histogram(binwidth = 60) +
  ggtitle("Histo Horsepower, bin 60")

# BOXPLOTS FOR OUTLIERS
# Construct box plot of msrp
cars %>%
  ggplot(aes(x = 1, y = msrp)) +
  geom_boxplot()

# Exclude outliers from data
cars_no_out <- cars %>%
  filter(msrp < 100000)

# Construct box plot of msrp using the reduced dataset
cars_no_out %>%
  ggplot(aes(x = 1, y = msrp)) +
  geom_boxplot()

# PLOT SELECTION -- BOXPLOT OR DENSITY PLOT?
# Create plot of city_mpg
cars %>%
  ggplot(aes(x = 1, y = city_mpg)) +
  geom_boxplot()

# Create plot of width
cars %>% 
  ggplot(aes(width)) +
  geom_density()

# 3 VARIABLE PLOT
# Facet hists using hwy mileage and ncyl
common_cyl %>%
  ggplot(aes(x = hwy_mpg)) +
  geom_histogram() +
  facet_grid(ncyl ~ suv, labeller = label_both) +
  ggtitle("hwy_mpg - ncyl vs. suv")

# CALCULATE CENTER MEASURES
glimpse(gapminder)
# Create dataset of 2007 data
gap2007 <- filter(gapminder, year %in% c(2007))

# Compute groupwise mean and median lifeExp
gap2007 %>%
  group_by(continent) %>%
  summarize(mean(lifeExp),
            median(lifeExp))

# Generate box plots of lifeExp for each continent
gap2007 %>%
  ggplot(aes(x = continent, y = lifeExp)) +
  geom_boxplot()

# CALCULATE SPREAD MEASURES
glimpse(gap2007)
# Compute groupwise measures of spread
gap2007 %>%
  group_by(continent) %>%
  summarize(sd(lifeExp),
            IQR(lifeExp),
            n())

# Generate overlaid density plots
gap2007 %>%
  ggplot(aes(x = lifeExp, fill = continent)) +
  geom_density(alpha = 0.3)

# CHOOSE CORRECT MEASURES FOR CENTER AND SPREAD
glimpse(gap2007)
# Compute stats for lifeExp in Americas
gap2007 %>%
  filter(continent == "Americas") %>%
  summarize(mean(lifeExp),
            sd(lifeExp))

# Compute stats for population
gap2007 %>%
  summarize(median(pop),
            IQR(pop))

# TRANSFORMATIONS
# Create density plot of old variable
gap2007 %>%
  ggplot(aes(x = pop)) +
  geom_density(alpha = .3)

# Transform the skewed pop variable
gap2007 <- gap2007 %>%
  mutate(log_pop = log(pop))

# Create density plot of new variable
gap2007 %>%
  ggplot(aes(x = log_pop)) +
  geom_density(alpha = .3)

# IDENTIFY OUTLIERS
# Filter for Asia, add column indicating outliers
gap_asia <- gap2007 %>%
  filter(continent == "Asia") %>%
  mutate(is_outlier = lifeExp < 50)

# Remove outliers, create box plot of lifeExp
gap_asia %>%
  filter(!is_outlier) %>%
  ggplot(aes(x = 1, y = lifeExp)) +
  geom_boxplot()

# INTRODUCING DATA
# SPAM AND NUM_CHAR
# Load packages
library(ggplot2)
library(dplyr)
library(openintro)

# Compute summary statistics
email %>%
  group_by(spam) %>%
  summarize(median(num_char),
            IQR(num_char))

# Create plot
email %>%
  mutate(log_num_char = log(num_char)) %>%
  ggplot(aes(x = spam, y = log_num_char)) +
  geom_boxplot()

# SPAM AND !!!!!
# Compute center and spread for exclaim_mess by spam

email %>%
  group_by(spam) %>%
  summarize(median(exclaim_mess),
            IQR(exclaim_mess))

# Create plot for spam and exclaim_mess
email %>%
  mutate(log_exclaim_mess = log(exclaim_mess)) %>%
  ggplot(aes(x = log_exclaim_mess)) +
  geom_histogram() +
  facet_wrap(~ spam)

# COLLAPSING LEVELS
# Create plot of proportion of spam by image

email %>%
  mutate(has_image = image > 0) %>%
  ggplot(aes(x = has_image, fill = spam)) +
  geom_bar(position = "fill")

# DATA INTEGRITY
# Test if images count as attachments
sum(email$num_char > 0)

sum(email$image >= 1)
sum(email$attach >= 1)
sum(email$image == email$attach)
sum(email$image > email$attach)
sum(email$image < email$attach)

# ANSWERING QUESTIONS USING CHAINS
# 1. For emails containing the word "dollar", does the typical spam email contain a greater number of occurrences of the word than the typical non-spam email? Create a summary statistic that answers this question.
# Question 1
email %>%
  filter(dollar > 0) %>%
  group_by(spam) %>%
  summarize(median(dollar))

# 2. If you encounter an email with greater than 10 occurrences of the word "dollar", is it more likely to be spam or not-spam? Create a barchart that answers this question.
# Question 2
email %>%
  filter(dollar > 10) %>%
  ggplot(aes(x = spam)) +
  geom_bar()

# WHAT'S IN A NUMBER --> REORDER LEVELS AND CONSTRUCT FACETED BAR CHART
# Reorder levels
email$number_reordered <- factor(email$number, levels = c("none", "small", "big"))

# Construct plot of number_reordered
ggplot(email, aes(number_reordered)) +
  geom_bar() +
  facet_wrap( ~ spam)









# DATA MANIPULATION WITH DPLYR ----
# Select the columns 
county_complete %>%
  select(state, county, population, poverty)
# ARRANGE OBSERVATIONS
counties_selected <- counties %>%
  select(state, county, population, private_work, public_work, self_employed)

# Add a verb to sort in descending order of public_work
counties_selected %>%
  arrange(desc(public_work))

# FILTERING FOR CONDITIONS
counties_selected <- counties %>%
  select(state, county, population)

# Filter for counties with a population above 1000000
counties_selected %>%
  filter(population > 1000000)

# Filter for counties in the state of California that have a population above 1000000
counties_selected %>%
  filter(state == "California", population > 1000000)

# FILTERING AND ARRANGING
counties_selected <- counties %>%
  select(state, county, population, private_work, public_work, self_employed)

# Filter for Texas and more than 10000 people; sort in descending order of private_work
counties_selected %>%
  filter(state == "Texas", population > 10000) %>%
  arrange(desc(private_work))

# CALCULATE THE NUMBER OF GOVERNMENT EMPLOYEES
counties_selected <- counties %>%
  select(state, county, population, public_work)

# Add a new column public_workers with the number of people employed in public work
counties_selected %>%
  mutate(public_workers = population * public_work / 100)


# CALCULATE THE PERCENTAGE OF WOMEN IN A COUNTY
# Select the columns state, county, population, men, and women
counties_selected <- counties %>%
  select(state, county, population, men, women)


# Calculate proportion_women as the fraction of the population made up of women
counties_selected %>%
  mutate(proportion_women = women / population)

# SELECT - MUTATE - FILTER - ARRANGE
counties %>%
  # Select the five columns 
  select(state, county, population, men, women) %>%
  
  # Add the proportion_men variable
  mutate(proportion_men = men / population) %>%
  
  # Filter for population of at least 10,000
  filter(population > 10000)   %>%
  
  # Arrange proportion of men in descending order 
  arrange(desc(proportion_men))

# COUNT BY REGION
# Use count to find the number of counties in each region
counties_selected %>%
  count(region, sort = TRUE)

# COUNTING CITIZENS BY STATE
# Find number of counties per state, weighted by citizens
counties_selected %>%
  count(state, sort = TRUE, wt = citizens)

# MUTATING AND COUNTING
counties_selected %>%
  # Add population_walk containing the total number of people who walk to work 
  mutate(population_walk = population * walk / 100) %>%
  
  # Count weighted by the new column
  count(state, wt = population_walk, sort = TRUE)

# SUMMARIZING
# Summarize to find minimum population, maximum unemployment, and average income
counties_selected %>%
  summarize(  min_population = min(population),
              max_unemployment = max(unemployment),
              average_income = mean(income))

# SUMMARIZE BY STATE
# Group by state and find the total area and population
counties_selected %>%
  group_by(state) %>%
  summarize(  total_area = sum(land_area),
              total_population = sum(population))

# Summarize to find the total population
counties_selected %>%
  group_by(region, state) %>%
  summarize(total_pop = sum(population))

# TOP_N() --> SELECTING A COUNTY FROM EACH REGION
# Group by region and find the greatest number of citizens who walk to work
counties_selected %>%
  group_by(region) %>%
  top_n(1, walk)

# FINDING THE HIGHEST INCOME STATE IN EACH REGION
counties_selected %>%
  group_by(region, state) %>%
  # Calculate average income
  summarize(average_income = mean(income)) %>%
  
  # Find the highest income state in each region
  top_n(1, average_income)

# USING SUMMARIZE - TOP_N() AND COUNT TOGETHER
# Find the total population for each combination of state and metro
counties_selected %>%
  group_by(state, metro) %>%
  summarize(total_pop = sum(population))

# SELECTING COLUMNS
# Glimpse the counties table
glimpse(counties)

counties %>%
  # Select state, county, population, and industry-related columns
  select(state, county, population, professional, service, office, construction, production) %>%
  
  # Arrange service in descending order 
  arrange(desc(service))

# SELECT HELPERS --> STARTS_WITH - ENDS_WITH - CONTAINS
counties %>%
  # Select the state, county, population, and those ending with "work"
  select(state, county, population, ends_with("work")) %>%
  filter(public_work > 50)

# RENAME --> RENAMING A COLUMNS AFTER COUNT
# Rename the n column to num_counties
counties %>%
  count(state) %>%
  rename(num_counties = n)

# RENAMING A COLUMN AS PART OF A SELECT
# Select state, county, and poverty as poverty_rate
counties %>%
  select(state, county, poverty_rate = poverty)

# TRANSMUTE()
counties %>%
  # Keep the state, county, and populations columns, and add a density column
  transmute(state, county, population, density = population / land_area) %>%
  
  # Filter for counties with a population greater than one million 
  filter(population > 1000000) %>%
  
  # Sort density in ascending order 
  arrange(density)

# CHOOSING CORRECT VERBS
# Change the name of the unemployment column
counties %>%
  rename(unemployment_rate = unemployment)

# Keep the state and county columns, and the columns containing poverty
counties %>%
  select(state, county, contains("poverty"))

# Calculate the fraction_women column without dropping the other columns
counties %>%
  mutate(fraction_women = women / population)

# Keep only the state, county, and employment_rate columns
counties %>%
  transmute(state, county, employment_rate = employed / population)


# THE BABYNAMES DATASET - FILTERING AND ARRANGING FOR ONE YEAR
glimpse(babynames)
babynames %>%
  # Filter for the year 1990
  filter(year == 1990) %>%
  
  # Sort the number column in descending order 
  arrange(desc(n))

glimpse(babynames)
# Find the most common name in each year
Baby_gr_b_name <- babynames %>%
  group_by(year, name) %>%
  top_n(10, n)

#VISUALIZE NAMES WITH GGPLOT2
# Filter for the names Steven, Thomas, and Matthew 
selected_names <- babynames %>%
  filter(name %in% c("Steven", "Thomas", "Matthew"))

# GROUPED MUTATE
# Calculate the fraction of people born each year with the same name
babynames %>%
  group_by(year) %>%
  mutate(year_total = sum(n)) %>%
  ungroup() %>%
  mutate(fraction = n / year_total) %>%
  # Find the year each name is most common
  group_by(name) %>%
  top_n(1, fraction)
  

# ADDING THE TOTAL AND MAXIMUM FOR EACH NAME
# Add columns name_total and name_max for each name
babynames %>%
  group_by(name) %>%
  mutate(   name_total = sum(n),
            name_max = max(n))
names_normalized <- babynames %>%
  group_by(name) %>%
  mutate(name_total = sum(n),
         name_max = max(n)) %>%
  # Ungroup the table 
  ungroup() %>%
  # Add the fraction_max column containing the number by the name maximum 
  mutate(fraction_max = n / name_max)

# VISUALIZE THE NORMALIZED CHANGE IN POPULARITY
# Filter for the names Steven, Thomas, and Matthew
names_filtered <- names_normalized %>%
  filter(name %in% c("Steven", "Thomas", "Matthew"))

# Visualize these names over time
ggplot(names_filtered, aes(x = year, y = fraction_max, color = name)) +
  geom_point()

# USING WINDOW FUNCTIONS - LAG()
# USING RATIO TO DESCRIBE THE FREQUENCY OF A NAME
babynames_fraction %>%
  # Arrange the data in order of name, then year 
  arrange(name, year) %>%
  # Group the data by name
  group_by(name) %>%
  # Add a ratio column that contains the ratio between each year 
  mutate(ratio = fraction / lag(fraction))

# BIGGEST JUMPS IN A NAME
babynames_ratios_filtered %>%
  # Extract the largest ratio from each name 
  top_n(1, ratio) %>%
  # Sort the ratio column in descending order 
  arrange(desc(ratio)) %>%
  # Filter for fractions greater than or equal to 0.001
  filter(fraction >= 0.001)



# CLEANING DATA IN R ----

# CHANGING THE TYPES OF DATA AND ASSERTING THEM (GIVES TRUE / FALSE)
  # Glimpse at bike_share_rides
  glimpse(bike_share_rides)
  # Summary of user_birth_year
  summary(bike_share_rides$user_birth_year)
  # Convert user_birth_year to factor: user_birth_year_fct
  bike_share_rides <- bike_share_rides %>%
    mutate(user_birth_year_fct = as.factor(user_birth_year))
  # Assert user_birth_year_fct is a factor
  assert_is_factor(bike_share_rides$user_birth_year_fct)
  # Summary of user_birth_year_fct
  summary(bike_share_rides$user_birth_year_fct)

# TRIMMING STRINGS
  glimpse(bike_share_rides)
  bike_share_rides <- bike_share_rides %>%
    # Remove 'minutes' from duration: duration_trimmed
    mutate(duration_trimmed = str_remove(duration, "minutes"),
           # Convert duration_trimmed to numeric: duration_mins
           duration_mins = as.numeric(duration_trimmed))
  # Glimpse at bike_share_rides
  glimpse(bike_share_rides)
  # Assert duration_mins is numeric
  assert_is_numeric(bike_share_rides$duration_mins)
  # Calculate mean duration
  mean(bike_share_rides$duration_mins)

# RANGE CONSTRAINTS
  # Create breaks
  breaks <- c(min(bike_share_rides$duration_min), 0, 1440, max(bike_share_rides$duration_min))
  # Create a histogram of duration_min
  ggplot(bike_share_rides, aes(duration_min)) +
    geom_histogram(breaks = breaks)
  # duration_min_const: replace vals of duration_min > 1440 with 1440
  bike_share_rides <- bike_share_rides %>%
    mutate(duration_min_const = replace(duration_min, duration_min > 1440, 1440))
  # Make sure all values of duration_min_const are between 0 and 1440
  assert_all_are_in_closed_range(bike_share_rides$duration_min_const, lower = 0, upper = 1440)

# EXCLUDE IMPOSSIBLE DATA (eg. FUTURE DATES, etc.)
  library(lubridate)
  # Convert date to Date type
  bike_share_rides <- bike_share_rides %>%
    mutate(date = as.Date(date))
  
  # Make sure all dates are in the past
  assert_all_are_in_past(bike_share_rides$date)
  
  # Filter for rides that occurred before or on today's date
  bike_share_rides_past <- bike_share_rides %>%
    filter(date <= today())
  
  # Make sure all dates from bike_share_rides_past are in the past
  assert_all_are_in_past(bike_share_rides_past$date)
  
# UNIQUENESS CONSTRAINTS == FILTER DUPLICATES --> DUPLICATED()
  # full vs. partial duplicates
  # duplicated()
  # distinct()
  
# REMOVING FULL DUPLICATES
  # Count the number of full duplicates
  sum(duplicated(bike_share_rides))
  # or: 
  filter(bike_share_rides, duplicated(bike_share_rides))
  # Remove duplicates
  bike_share_rides_unique <- distinct(bike_share_rides)
  # Count the full duplicates in bike_share_rides_unique
  sum(duplicated(bike_share_rides_unique))
  
# REMOVING PARTIAL DUPLICATES
  # Find duplicated ride_ids
  bike_share_rides %>% 
    count(ride_id) %>% 
    filter(n > 1)
  # Remove full and partial duplicates
  bike_share_rides_unique <- bike_share_rides %>%
    # Only based on ride_id instead of all cols
    distinct(ride_id, .keep_all = TRUE)
  # Find duplicated ride_ids in bike_share_rides_unique
  bike_share_rides_unique %>%
    # Count the number of occurrences of each ride_id
    count(ride_id) %>%
    # Filter for rows with a count > 1
    filter(n > 1)
  
# AGGREGATING PARTIAL DUPLICATES
  # compute a summary statistic of differing values (e.g. mean / median / max / min)
  bike_share_rides %>%
    # Group by ride_id and date
    group_by(ride_id, date) %>%
    # Add duration_min_avg column
    mutate(duration_min_avg = mean(duration_min) ) %>%
    # Remove duplicates based on ride_id and date, keep all cols
    distinct(ride_id, date, .keep_all = TRUE) %>%
    # Remove duration_min column
    select(-duration_min)
  
# FIND NON-MEMBERS OF A DATASET
  # Anti-joins can help you identify the rows that are causing issues
  # Semi-joins can remove the issue-causing rows
  
  # Find bad dest_size rows
  sfo_survey %>% 
    # Join with dest_sizes data frame to get bad dest_size rows
    anti_join(dest_sizes, by = "dest_size") %>%
    # Select id, airline, destination, and dest_size cols
    select(id, airline, destination, dest_size)
  # Remove bad dest_size rows
  sfo_survey %>% 
    # Join with dest_sizes
    semi_join(dest_sizes, by = "dest_size") %>%
    # Count the number of each dest_size
    count(dest_size)

# CATEGORICAL DATA PROBLEMS
  # SOLVE PROBLEMS WITH UPPER / LOWER CASE DIFFERENCES == CONVERT / CONVERSION
  # SOLVE PROBLEMS WITH WHITESPACE == TRIM
  # Add new columns to sfo_survey
  sfo_survey <- sfo_survey %>%
    # dest_size_trimmed: dest_size without whitespace and
    # cleanliness_lower: cleanliness converted to lowercase
    mutate(dest_size_trimmed = str_trim(dest_size),
           cleanliness_lower = str_to_lower(cleanliness))
  # Count values of dest_size_trimmed and
  # Count values of cleanliness_lower
  sfo_survey %>%
    count(dest_size_trimmed)
  sfo_survey %>%
    count(cleanliness_lower)

# COLLAPSING CATEGORIES
  # e.g. INSTEAD OF LABRADOR / BOXER / BEAGLE --> JUST "DOGS"
  # Count categories of dest_region
  sfo_survey %>%
    count(dest_region)
  # Categories to map to Europe
  europe_categories <- c("Europ", "EU", "eur")
  # Add a new col dest_region_collapsed
  sfo_survey %>%
    # Map all categories in europe_categories to Europe
    # and Count categories of dest_region_collapsed
    mutate(dest_region_collapsed = fct_collapse(dest_region, 
                                                Europe = europe_categories)) %>%
    count(dest_region_collapsed)
  
# CLEANING TEXT DATA
  # detect() e.g. detect wrong input like credit card data
  # str_replace_all() - e.g. replace hyphens
  # str_remove_all() - e.g. remove all hyphens and spaces
  # str_length() - e.g. find wrong length numbers (like less/more than 16 digits)
  # filter()
  # fixed() - e.g. for characters that are treated differently in a regular expression (), [], $, . , +, *, ...
  #   search for these characters requires using the following function:
  #       str_detect(column, fixed("$"))
  
  # FILTERING DIRTY DATA
  # Filter for rows with "-" in the phone column
  sfo_survey %>%
    filter(str_detect(phone, "-"))
  # Filter for rows with "(" or ")" in the phone column
  sfo_survey %>%
    filter(str_detect(phone, fixed("(")) | str_detect(phone, fixed(")")))
  
  # REPLACE AND REMOVE
  # Remove parentheses from phone column
  phone_no_parens <- sfo_survey$phone %>%
    # Remove "("s
    str_remove_all(fixed("(")) %>%
    # Remove ")"s
    str_remove_all(fixed(")"))
  # Add phone_no_parens as column
  sfo_survey %>%
    mutate(phone_no_parens = phone_no_parens,
           # Replace all hyphens in phone_no_parens with spaces
           phone_clean = str_replace_all(phone_no_parens, "-", " "))
  # FILTER AND REMOVE
  # Check out the invalid numbers
  sfo_survey %>%
    filter(str_length(phone) != 12)
  # Remove rows with invalid numbers
  sfo_survey %>%
    filter(str_length(phone) == 12)

# UNIFORMITY
  # FILTER VALUES THAT DON'T FIT THE DATA DUE TO DIFFERENT UNITS
  #   e.g. C/F, g/lb/kg, dates...
  #   there are many date formats --> get overview: ?strptime
  # parse_date_time()
  #     e.g. parse_date_time("Monday, January 3",
  #                         orders = c("%Y-%m-%d", "%m-%d-%y", "%B-%d-%Y"))
  
# DATE UNIFORMITY
  # Check out the accounts data frame
  head(accounts)
  # Define the date formats
  formats <- c("%Y-%m-%d", "%B %d, %Y")
  # Convert dates to the same format
  accounts %>%
    mutate(date_opened_clean = parse_date_time(date_opened, formats))

# CURRENCY UNIFORMITY
  # Scatter plot of opening date and total amount
  accounts %>%
    ggplot(aes(x = date_opened, y = total)) +
    geom_point()
  
  # Left join accounts to account_offices by id
  accounts %>%
    left_join(account_offices, by = "id") %>%
    # Convert totals from the Tokyo office to USD
    mutate(total_usd = ifelse(office == "Tokyo", total / 104, total)) %>%
    # Scatter plot of opening date vs total_usd
    ggplot(aes(x = date_opened, y = total_usd)) +
    geom_point()
  
# VALIDATING THE DATA 
# VALIDATING TOTALS
  # DO NUMBERS IN THE TOTAL ADD UP TO THE THEORETICAL TOTAL VALUES?
  # Find invalid totals
  accounts %>%x
  # theoretical_total: sum of the three funds
  mutate(theoretical_total =  fund_A + fund_B + fund_C)%>%
    # Find accounts where total doesn't match theoretical_total
    filter(theoretical_total != total)
  
# VALIDATING AGE
  # Find invalid acct_age
  accounts %>%
    # theoretical_age: age of acct based on date_opened
    mutate(theoretical_age = floor(as.numeric(date_opened %--% today(), "years"))) %>%
    # Filter for rows where acct_age is different from theoretical_age
    filter(theoretical_age != acct_age)
  
# COMPLETENESS -- MISSING DATA (e.g. NA, nan, 0, 99, . , ...)
  # is.na(variable) --> returns TRUE / FALSE for each observation
  # sum(is.na(variable)) --> gives general amount of NA's in all variable observations
  # visualizing:visdat package
  #   vis_miss() --> missing = black
  
# INVESTIGATE MISSINGNESS
  # Visualize the missing values by column
  vis_miss(accounts)
  
  accounts %>%
    # missing_inv: Is inv_amount missing?
    mutate(missing_inv = is.na(inv_amount)) %>%
    # Group by missing_inv
    group_by(missing_inv) %>%
    # Calculate mean age for each missing_inv group
    summarize(avg_age = mean(age))
  # Sort by age and visualize missing vals
  accounts %>%
    arrange(age) %>%
    vis_miss()
  
# TREATING MISSING DATA
  # here missing data can be replaced with 5 * inv_amount, but obviously depends highly on data
  # Create accounts_clean
  accounts_clean <- accounts %>%
    # Filter to remove rows with missing cust_id
    filter(!is.na(cust_id)) %>%
    # Add new col acct_amount_filled with replaced NAs
    mutate(acct_amount_filled = ifelse(is.na(acct_amount), inv_amount * 5, acct_amount))
  # Assert that cust_id has no missing vals
  assert_all_are_not_na(accounts_clean$cust_id)
  # Assert that acct_amount_filled has no missing vals
  assert_all_are_not_na(accounts_clean$acct_amount_filled)
  
  
# COMPARING STRINGS - HOW SIMILAR ARE WORDS?
  # TYPES OF EDIT DISTANCE
  # different algorithms avaliable --> see which one works best for you
  # here "dl" algorithm used (Damerau-Levenshtein distance)
  library(stringdist)
  library(fuzzyjoin)
  # Calculate Damerau-Levenshtein distance
  stringdist("las angelos", "los angeles", method = "dl")
  # Calculate LCS distance
  stringdist("las angelos", "los angeles", "lcs")
  # Calculate Jaccard distance
  stringdist("las angelos", "los angeles", "jaccard")
  
  
# FIXING TYPOS WITH STRING DISTANCE
  # Count the number of each city variation
  zagat %>%
    count(city)
  # Join zagat and cities and look at results
  zagat %>%
    # Left join based on stringdist using city and city_actual cols
    stringdist_left_join(cities, by = c("city" = "city_actual")) %>%
    # Select the name, city, and city_actual cols
    select(name, city, city_actual)
  
# COMPARE AND GENERATE PAIRS:
  # 1. Start with 2 data sets --> clean the datasets
  # 2. generate pairs of records
  # 3. compare pairs (compare sparate columns of each pair)
  # 4. score pairs using summing or probability
  # 5. select pairs that are matches based on their score
  # 6. link datasets together
  
# PAIR BLOCKING:
  # Load reclin
  library(reclin)
  # Generate all possible pairs
  pair_blocking(zagat, fodors)
  # Generate pairs with same city
  pair_blocking(zagat, fodors, blocking_var = "city")
  
# COMPARING PAIRS:
  # Generate pairs
  pair_blocking(zagat, fodors, blocking_var = "city") %>%
    # Compare pairs by name using lcs()
    compare_pairs(by = "name",
                  default_comparator = lcs())
  
  # Generate pairs
  pair_blocking(zagat, fodors, blocking_var = "city") %>%
    # Compare pairs by name, phone, addr
    compare_pairs(by = c("name", "phone", "addr"),
                  default_comparator = jaro_winkler())
  
# SCORING AND LINKING
  # Create pairs
  pair_blocking(zagat, fodors, blocking_var = "city") %>%
    # Compare pairs
    compare_pairs(by = "name", default_comparator = jaro_winkler()) %>%
    # Score pairs
    score_problink() %>%
    # Select pairs
    select_n_to_m() %>%
    # Link data 
    link()
  
  
  
  
  
  
  
  
# INTERMEDIATE R ----

# RELATIONAL OPERATORS
# Equality
  # Comparison of logicals
  TRUE == FALSE
  # Comparison of numerics
  -6*14 != 17 -101
  # Comparison of character strings
  "useR" == "user"
  # Compare a logical with a numeric
  TRUE == 1
  #Awesome! Since TRUE coerces to 1 under the hood, TRUE == 1 evaluates to TRUE. 
  #Make sure not to mix up == (comparison) and = (assignment), == is what need to check the equality of R objects.

# Greater than vs. less than
  # Comparison of numerics
  -6 * 5 + 2 >= -10 + 1
  # Comparison of character strings
  "raining" <= "raining dogs"
  # Comparison of logicals
  TRUE > FALSE
  
# Compare Vectors
  #On which days did the number of LinkedIn profile views exceed 15?
  # The linkedin and facebook vectors have already been created for you
  linkedin <- c(16, 9, 13, 5, 2, 17, 14)
  facebook <- c(17, 7, 5, 16, 8, 13, 14)
  # Popular days
  linkedin >15
  #When was your LinkedIn profile viewed only 5 times or fewer?
  # Quiet days
  linkedin <= 5
  #When was your LinkedIn profile visited more often than your Facebook profile?
  # LinkedIn more popular than Facebook
  linkedin > facebook
  
# Compare Matrices
  #When were the views exactly equal to 13? Use the views matrix to return a logical matrix.
  # The social data has been created for you
  linkedin <- c(16, 9, 13, 5, 2, 17, 14)
  facebook <- c(17, 7, 5, 16, 8, 13, 14)
  views <- matrix(c(linkedin, facebook), nrow = 2, byrow = TRUE)
  
  views == 13
  #For which days were the number of views less than or equal to 14? Again, have R return a logical matrix.
  views <= 14

# LOGICAL OPERATORS
  # & vs. | (AND vs. OR)
  linkedin <- c(16, 9, 13, 5, 2, 17, 14)
  last <- tail(linkedin, 1)

  # Is last under 5 or above 10?
  last <5 | last > 10
  # Is last between 15 (exclusive) and 20 (inclusive)?
  last >20 | last >= 20
  #When did LinkedIn views exceed 10 and did Facebook views fail to reach 10 for a particular day? Use the linkedin and facebook vectors.
  linkedin > 10 & facebook < 10
  #When were one or both of your LinkedIn and Facebook profiles visited at least 12 times?
  linkedin >=12 | facebook >=12
  #When is the views matrix equal to a number between 11 and 14, excluding 11 and including 14?
  views >11 & views <= 14
  
# BLENDING OPERATORS TOGETHER
  #Select the entire second column, named day2, from the li_df data frame as 
  #a vector and assign it to second.
  second <- li_df$day2
  #Use second to create a logical vector, that contains TRUE if the 
  #corresponding number of views is strictly greater than 25 or strictly lower 
  #than 5 and FALSE otherwise. Store this logical vector as extremes.
  extremes <- second > 25 | second < 5
  #Use sum() on the extremes vector to calculate the number of TRUEs in extremes
  #(i.e. to calculate the number of employees that are either very popular or 
  #very low-profile). Simply print this number to the console.
  sum(extremes)

# CONDITIONAL STATEMENTS
  if (condition1) {
    expr1
  } else if (condition2) {
    expr2
  } else if (condition3) {
    expr3
  } else {
    expr4
  }
  # use IF statement
  # Variables related to your last day of recordings
  medium <- "LinkedIn"
  num_views <- 14
  # Examine the if statement for medium
  if (medium == "LinkedIn") {
    print("Showing LinkedIn information")
  }
  # Write the if statement for num_views
  if (num_views > 15) {
    print("You are popular!")
  }
  
  # add an ELSE statement
  # Control structure for medium
  if (medium == "LinkedIn") {
    print("Showing LinkedIn information")
  } else {
    print("Unknown medium")
  }
  # Control structure for num_views
  if (num_views > 15) {
    print("You're popular!")
  } else {
    print("Try to be more visible!")
  }
  
  # ELSE IF statement
  # Variables related to your last day of recordings
  medium <- "LinkedIn"
  num_views <- 14

  # Control structure for medium
  if (medium == "LinkedIn") {
    print("Showing LinkedIn information")
  } else if (medium == "Facebook") {
    # Add code to print correct string when condition is TRUE
    print("Showing Facebook information")
  } else {
    print("Unknown medium")
  }
  # Control structure for num_views
  if (num_views > 15) {
    print("You're popular!")
  } else if (num_views <= 15 & num_views > 10) {
    # Add code to print correct string when condition is TRUE
    print("Your number of views is average")
  } else {
    print("Try to be more visible!")
  }
  
  # Combining all
  li <- 15
  fb <- 9
  
  # Code the control-flow construct
  if (li >=15 & fb >=15) {
    sms <- 2 * (li + fb)
  } else if (li <10 & fb < 10) {
    sms <- 0.5 * (li + fb)
  } else {
    sms <- li + fb
  } 
  # Print the resulting sms to the console
  print(sms)
  
# LOOPS
  # WHILE loop
  while (condition) {
    expr
  }
  #The condition of the while loop should check if speed is higher than 30.
  speed <- 64
  #Inside the body of the while loop, print out "Slow down!".
  #Inside the body of the while loop, decrease the speed by 7 units and assign 
  #this new value to speed again. This step is crucial; 
  #otherwise your while loop will never stop and your session will expire.
  while (speed > 30) {
    print("Slow down!")
    speed <- speed - 7
  }
  speed
  
  speed <- 64
  # Extend/adapt the while loop
  #If the speed is greater than 48, have R print out "Slow down big time!", 
  #and decrease the speed by 11.
  #Otherwise, have R simply print out "Slow down!", and decrease the speed by 6.
  while (speed > 30) {
    print(paste("Your speed is",speed))
    if (speed >48 ) {
      print("Slow down big time!")
      speed <- speed - 11
    } else {
      print("Slow down!")
      speed <- speed - 6
    }
  }

  # Stopping WHILE loops --> break
  speed <- 88
  while (speed > 30) {
    print(paste("Your speed is", speed))
    # Break the while loop when speed exceeds 80
    if (speed > 80) {
      break
    }
    if (speed > 48) {
      print("Slow down big time!")
      speed <- speed - 11
    } else {
      print("Slow down!")
      speed <- speed - 6
    }
  }
  
  # Initialize i as 1 
  i <- 1
  #Finish the while loop so that it:
  #prints out the triple of i, so 3 * i, at each run.
  #is abandoned with a break if the triple of i is divisible by 8, 
  #but still prints out this triple before breaking.
  while (i <= 10) {
    print(3*i)
    if ((3 * i) %% 8 == 0) {
      print(3*i)
      break
    }
    i <- i + 1
  }
  
  # FOR loops
  primes <- c(2, 3, 5, 7, 11, 13)
  # loop version 1
  for (p in primes) {
    print(p)
  }
  # loop version 2
  for (i in 1:length(primes)) {
    print(primes[i])
  }
  
  #Write a for loop that iterates over all the elements of linkedin and prints 
  #out every element separately. Do this in two ways: using the loop version 1 
  #and the loop version 2 in the example code above.
  linkedin <- c(16, 9, 13, 5, 2, 17, 14)
  # Loop version 1
  for(i in linkedin) {
    print(i)
  }
  # Loop version 2
  for(i in 1:length(linkedin)) {
    print(linkedin[i])
  }
  
  # Loop over a list
  #As in the previous exercise, loop over the nyc list in two different ways to print its elements:
  nyc <- list(pop = 8405837, 
              boroughs = c("Manhattan", "Bronx", "Brooklyn", "Queens", "Staten Island"), 
              capital = FALSE)
  #Loop directly over the nyc list (loop version 1).
  for(i in nyc) {
    print(i)
  }
  #Define a looping index and do subsetting using double brackets (loop version 2).
  for(i in 1:length(nyc)) {
    print(nyc[[i]])
  }
  
  # Loop over a matrix
  for (var1 in seq1) {
    for (var2 in seq2) {
      expr
    }
  }
  
  #The outer loop should loop over the rows, with loop index i (use 1:nrow(ttt)).
  #The inner loop should loop over the columns, with loop index j (use 1:ncol(ttt)).
  #Inside the inner loop, make use of print() and paste() to print out 
  #information in the following format: "On row i and column j the board contains x",
  #where x is the value on that position.
  # The tic-tac-toe matrix ttt has already been defined for you
  # define the double for loop
  for (i in 1:nrow(ttt)) {
    for (j in 1:ncol(ttt)) {
      print(paste("On row", i, "and column", j, "the board contains", ttt[i,j]))
    }
  }
  
  #Add code to the for loop that loops over the elements of the linkedin vector:
  #If the vector element's value exceeds 10, print out "You're popular!".
  #If the vector element's value does not exceed 10, print out "Be more visible!"
  linkedin <- c(16, 9, 13, 5, 2, 17, 14)
  for (li in linkedin) {
    if (li > 10) {
      print("You're popular!")
    } else {
      print("Be more visible!")
    }
    print(li)
  }
  
  
  #Extend the for loop with two new, separate if tests as follows:
  #If the vector element's value exceeds 16, print out "This is ridiculous, 
  #I'm outta here!" and have R abandon the for loop (break).
  #If the value is lower than 5, print out "This is too embarrassing!" 
  #and fast-forward to the next iteration (next).
  for (li in linkedin) {
    if (li > 10) {
      print("You're popular!")
    } else {
      print("Be more visible!")
    }
    # Add if statement with break
    if(li > 16) {
      print("This is ridiculous, I'm outta here!")
      break
    }
    # Add if statement with next
    if (li < 5) {
      print("This is too embarrassing!")
      next
    }
    print(li)
  }
  
  #This exercise will not introduce any new concepts on for loops.
  #We already went ahead and defined a variable rquote. This variable has been 
  #split up into a vector that contains separate letters and has been stored in 
  #a vector chars with the strsplit() function.
  #Can you write code that counts the number of r's that come before the first u in rquote?
  #Initialize the variable rcount, as 0.
  #if char equals "r", increase the value of rcount by 1.
  #if char equals "u", leave the for loop entirely with a break.
  #Finally, print out the variable rcount to the console to see if your code is correct.
  # Pre-defined variables
  rquote <- "r's internals are irrefutably intriguing"
  chars <- strsplit(rquote, split = "")[[1]]
  
  # Initialize rcount
  rcount <- 0
  
  for (char in chars) {
    if(char == "r") {
      rcount = rcount + 1
    }
    if(char == "u") {
      break
    }
  }
  print(rcount)
  
  
  
  
  
  
  
  
  
# LEAFLET #####################################################################
  
library(leaflet)
library(leaflet.extras)
library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)
library(htmltools)


# providers list 
providers
names(providers)
str_detect(names(providers), "CartoDB")
names(providers)[str_detect(names(providers), "CartoDB")]


# add custom Map Tile
leaflet() %>% 
  addProviderTiles(provider = "CartoDB")
  # Esri
  # CartoDB.PositronNoLabels

# Setting the Default Map View
## geocode()
## single point --> setView() 
## select rectangle --> fitBounds()
## stay focused --> turn pan off: dragging = F 
## e.g.:    350 5th Ave, Floor 77, New York, NY 10118
leaflet()  %>% 
  addProviderTiles("CartoDB")  %>% 
  setView(lng = -73.98575, lat = 40.74856, zoom = 6)
## e.g.:    get coordinates from file
leaflet() %>% 
  addProviderTiles("CartoDB.PositronNoLabels") %>% 
  setView(lng = dc_hq$lon[2], lat = dc_hq$lat[2], zoom = 4)

## Create a map with a narrower view
leaflet(options = leafletOptions(
  # Set minZoom and dragging 
  minZoom = 12, dragging = T))  %>% 
  addProviderTiles("CartoDB")  %>% 
  # Set default zoom level 
  setView(lng = dc_hq$lon[2], lat = dc_hq$lat[2], zoom = 14) %>% 
  # Set max bounds of map 
  setMaxBounds(lng1 = dc_hq$lon[2] + .05, 
               lat1 = dc_hq$lat[2] + .05, 
               lng2 = dc_hq$lon[2] - .05, 
               lat2 = dc_hq$lat[2] - .05) 

## Plotting Points on the map
# Plot DataCamp's NYC HQ
leaflet() %>% 
  addProviderTiles("CartoDB") %>% 
  addMarkers(lng = dc_hq$lon[1], lat = dc_hq$lat[1])
## or
# Plot DataCamp's NYC HQ with zoom of 12    
leaflet() %>% 
  addProviderTiles("CartoDB") %>% 
  addMarkers(lng = -73.98575, lat = 40.74856)  %>% 
  setView(lng = -73.98575, lat = 40.74856, zoom = 12)    
## or
# Plot both DataCamp's NYC and Belgium locations
leaflet() %>% 
  addProviderTiles("CartoDB") %>% 
  addMarkers(lng = dc_hq$lon, lat = dc_hq$lat)

## Store Map Data
# Store leaflet hq map in an object called map
map <- leaflet() %>%
  addProviderTiles("CartoDB") %>%
  # Use dc_hq to add the hq column as popups
  addMarkers(lng = dc_hq$lon, lat = dc_hq$lat,
             popup = dc_hq$hq)
# Center the view of map on the Belgium HQ with a zoom of 5 
map_zoom <- map %>%
  setView(lat = 50.881363, lng = 4.717863,
          zoom = 5)
map_zoom

# Chapter 2: 

## Cleaning up the base map
# clearMarkers() == remove 1+ features
# clearBounds() == clear bounds and automatically determine bounds based on map elements
# Remove markers, reset bounds, and store the updated map in the m object
map_clear <- map %>%
  clearMarkers() %>% 
  clearBounds()
map_clear

## IPEDS == Integrated Postsecondary Education Data System
# Remove colleges with missing sector information
ipeds <- 
  ipeds_missing %>% 
  drop_na()
# Create a list of US States in descending order by the number of colleges in each state
ipeds %>%
  group_by(state) %>%
  count() %>%
  arrange(desc(n))

## Mapping California Colleges
# Create a dataframe called `ca` with data on only colleges in California
ca <- ipeds %>%
  filter(state == "CA")

# Use `addMarkers` to plot all of the colleges in `ca` on the `m` leaflet map
map %>%
  addMarkers(lng = ca$lng, lat = ca$lat)

# Center the map on LA 
map %>% 
  addMarkers(data = ca) %>% 
  setView(lat = la_coords$lat, lng = la_coords$lon, zoom = 12)
# Set the zoom level to 8 and store in the m object
map_zoom <- 
  map %>%
  addMarkers(data = ca) %>%
  setView(lat = la_coords$lat, lng = la_coords$lon, zoom = 8)
map_zoom

# Clear the markers from the map 
map2 <- map %>% 
  clearMarkers
# Use addCircleMarkers() to plot each college as a circle
map2 %>%
  addCircleMarkers(lng = ca$lng, lat = ca$lat)
# Change the radius of each circle to be 2 pixels and the color to red
map2 %>% 
  addCircleMarkers(lng = ca$lng, lat = ca$lat,
                   radius = 2, color = "red")

# Add circle markers with popups for college names
map %>%
  addCircleMarkers(data = ca, radius = 2, popup = ~name)
# Change circle color to #2cb42c and store map in map_color object
map_color <- map %>% 
  addCircleMarkers(data = ca, radius = 2, color = "#2cb42c", popup = ~name)
map_color



## COSTUMIZING POP-UPs
# Add circle markers with popups that display both the institution name and sector
map2 %>% 
  addCircleMarkers(data = ca, radius = 2, 
                   popup = ~paste0(name, "<br/>", sector_label))
# Make the institution name in each popup bold
map2 %>% 
  addCircleMarkers(data = ca, radius = 2, 
                   popup = ~paste0("<b>", name, "</b>", "<br/>", sector_label))

# Use paste0 to add sector information to the label inside parentheses 
map %>% 
  addCircleMarkers(data = ca, radius = 2, label = ~paste0(name, " (", sector_label, ")"))


# COLOR CODING
# Make a color palette called pal for the values of `sector_label` using `colorFactor()`  
# Colors should be: "red", "blue", and "#9b4a11" for "Public", "Private", and "For-Profit" colleges, respectively
pal <- colorFactor(palette = c("red", "blue", "#9b4a11"), 
                   levels = c("Public", "Private", "For-Profit"))

# Add circle markers that color colleges using pal() and the values of sector_label
map2 <- 
  map %>% 
  addCircleMarkers(data = ca, radius = 2, 
                   color = ~pal(sector_label), 
                   label = ~paste0(name, " (", sector_label, ")"))
map2

# Customize the legend
m %>% 
  addLegend(pal = pal, 
            values = c("Public", "Private", "For-Profit"),
# opacity of .5, title of Sector, and position of topright
            opacity = 0.5, title = "Sector", 
            position = "topright")




# Chapter 3 --> leaflet.extras package
library(leaflet.extras)

leaflet() %>%
  addTiles() %>% 
  addSearchOSM() %>% 
  addReverseSearchOSM() 

## Building a Base
m2 <- 
  ipeds %>% 
  leaflet() %>% 
  # use the CartoDB provider tile
  addProviderTiles("CartoDB") %>% 
  # center on the middle of the US with zoom of 3
  setView(lat = 39.8282, lng = -98.5795, zoom = 3)
# Map all American colleges 
m2 %>% 
  addCircleMarkers()


## Overlay Groups
# plot as separate layers --> addLayersControl()
library(htmltools)
# Create data frame called public with only public colleges
public <- filter(ipeds, sector_label == "Public")  
# Create a leaflet map of public colleges called m3 
m3 <- leaflet() %>% 
  addProviderTiles("CartoDB") %>% 
  addCircleMarkers(data = public, radius = 2, label = ~htmlEscape(name),
                   color = ~pal(sector_label), group = "Public")
m3

# Create data frame called private with only private colleges
private <- filter(ipeds, sector_label == "Private")  
# Add private colleges to `m3` as a new layer
m3 <- m3 %>% 
  addCircleMarkers(data = private, 
                   radius = 2, 
                   label = ~htmlEscape(name),
                   color = ~pal(sector_label), 
                   group = "Private") %>% 
  addLayersControl(overlayGroups = c("Public", "Private"))
m3

# Create data frame called profit with only For-Profit colleges
profit <- filter(ipeds, sector_label == "For-Profit")  
# Add For-Profit colleges to `m3` as a new layer
m3 <- m3 %>% 
  addCircleMarkers(data = profit, radius = 2, label = ~htmlEscape(name),
                   color = ~pal(sector_label),   group = "For-Profit")  %>% 
  addLayersControl(overlayGroups = c("Public", "Private", "For-Profit"))  
# Center the map on the middle of the US with a zoom of 4
m4 <- m3 %>%
  setView(lat = 39.8282, lng = -98.5795, zoom = 4) 
m4


## BASE MAP GROUPS
leaflet() %>% 
  # Add the OSM, CartoDB and Esri tiles
  addTiles(group = "OSM") %>% 
  addProviderTiles("CartoDB", group = "CartoDB") %>% 
  addProviderTiles("Esri", group = "Esri") %>% 
  # Use addLayersControl to allow users to toggle between basemaps
  addLayersControl(baseGroups = c("OSM", "CartoDB", "Esri"))


## EVERYTHING TOGETHER
m4 <- leaflet() %>% 
  addTiles(group = "OSM") %>% 
  addProviderTiles("CartoDB", group = "Carto") %>% 
  addProviderTiles("Esri", group = "Esri") %>% 
  addCircleMarkers(data = public, radius = 2, label = ~htmlEscape(name),
                   color = ~pal(sector_label),  group = "Public") %>% 
  addCircleMarkers(data = private, radius = 2, label = ~htmlEscape(name),
                   color = ~pal(sector_label), group = "Private")  %>% 
  addCircleMarkers(data = profit, radius = 2, label = ~htmlEscape(name),
                   color = ~pal(sector_label), group = "For-Profit")  %>% 
  addLayersControl(baseGroups = c("OSM", "Carto", "Esri"), 
                   overlayGroups = c("Public", "Private", "For-Profit")) %>% 
  setView(lat = 39.8282, lng = -98.5795, zoom = 4) 
m4

## CLUSTERING
# Make each sector of colleges searchable 
m4_search <- m4  %>% 
  addSearchFeatures(
    targetGroups = c("Public", "Private", "For-Profit"), 
    # Set the search zoom level to 18
    options = searchFeaturesOptions(zoom = 18)) 
m4_search

# library(ipeds)
# library(htmltools)
# color palette pal
ipeds %>% 
  leaflet() %>% 
  addTiles() %>% 
  # Sanitize any html in our labels
  addCircleMarkers(radius = 2, 
                   label = ~htmlEscape(name),
                   color = ~pal(sector_label),
                   # Cluster all colleges using `clusterOptions`
                   clusterOptions = markerClusterOptions()) 


# SPATIAL DATA --> using Polygons

summary(shp)
class(shp)
slotNames(shp)
glimpse(shp@data)
class(shp@data)
shp@data$GEOID10 # Print GEOID10

## Joining Spatial Data
glimpse(nc_income)
summary(nc_income)

# Left join nc_income onto shp@data and store in shp_nc_income
shp_nc_income <- shp@data %>% 
  left_join(nc_income, by = c("GEOID10" = "zipcode"))
# Print the number of missing values of each variable in shp_nc_income
shp_nc_income  %>%
  summarize_all(funs(sum(is.na(.))))

## Mapping Polygons
shp %>% 
  leaflet() %>% 
  addTiles() %>% 
  addPolygons()
# which zips were not in the income data?
shp_na <- shp[is.na(shp$mean_income),]
# map the polygons in shp_na
shp_na %>% 
  leaflet() %>% 
  addTiles() %>% 
  addPolygons()

## Visualize quartiles
# summarize the mean income variable
summary(shp$mean_income)
# subset shp to include only zip codes in the top quartile of mean income
high_inc <- shp[!is.na(shp$mean_income) & shp$mean_income > 55917,]
# map the boundaries of the zip codes in the top quartile of mean income
high_inc %>%
  leaflet() %>%
  addTiles() %>%
  addPolygons()

# create color palette with colorNumeric()
nc_pal <- colorNumeric("YlGn", domain = high_inc@data$mean_income)
high_inc %>%
  leaflet() %>%
  addTiles() %>%
  addPolygons(weight = 1, 
              color = ~nc_pal(mean_income),   # set boundary thickness to 1 and color polygons
              label = ~paste0("Mean Income: ", dollar(mean_income)), # add labels that display mean income
              highlight = highlightOptions(weight = 5, # highlight polygons on hover
                                           color = "white", 
                                           bringToFront = TRUE))
# Use the log function to create a new version of nc_pal
nc_pal <- colorNumeric("YlGn", domain = log(high_inc@data$mean_income))
# comment out the map tile
high_inc %>%
  leaflet() %>%
  #addProviderTiles("CartoDB") %>%
  # apply the new nc_pal to the map
  addPolygons(weight = 1, color = ~nc_pal(log(mean_income)), fillOpacity = 1,
              label = ~paste0("Mean Income: ", dollar(mean_income)),
              highlightOptions = highlightOptions(weight = 5, color = "white", bringToFront = TRUE))


slotNames(wealthy_zips)
summary(wealthy_zips@data$mean_income)
# plot zip codes with mean incomes >= $200k
wealthy_zips %>% 
  leaflet() %>% 
  addProviderTiles("CartoDB") %>% 
  # set color to green and create Wealth Zipcodes group
  addPolygons(weight = 1, fillOpacity = .7, color = "green",  group = "Wealthy Zipcodes", 
              label = ~paste0("Mean Income: ", dollar(mean_income)),
              highlightOptions = highlightOptions(weight = 5, color = "white", bringToFront = TRUE))

## FINAL MAP
# Add polygons using wealthy_zips
final_map <- m4 %>% 
  addPolygons(data = wealthy_zips, 
              weight = 1, 
              fillOpacity = .5, 
              color = "Grey",  
              group = "Wealthy Zip Codes", 
              label = ~paste0("Mean Income: ", dollar(mean_income)),
              highlight = highlightOptions(weight = 5, 
                                           color = "white", 
                                           bringToFront = TRUE)) %>% 
  # Update layer controls including "Wealthy Zip Codes"
  addLayersControl(baseGroups = c("OSM", "Carto", "Esri"), 
                   overlayGroups = c("Public", "Private", "For-Profit", "Wealthy Zip Codes"))
final_map




  
  
  
  
  
  