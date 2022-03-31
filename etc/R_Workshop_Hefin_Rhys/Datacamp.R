# Data Camp -- Intro Course ----

library(dplyr)
library(ggplot2)

table(iris$Species)

# Drop unused Levels
iris$Species <- droplevels(iris$Species) 
iris$Species

# use of mutate & ifelse function
  # ifelse([logical condition], [do this if TRUE], [do this if FALSE])
  # mutate(read_cat = ifelse(
  #        read < avg_read,
  #        "below average",
  #        "at or above average"
  #      )
  # )


# Sampling Data in R ----

library(dplyr)
library(openintro)
# load county data
data(county)
# remove DC
county_noDC <- county %>%
  filter(state != "District of Columbia") %>%
  droplevels()

# Simple Random Sample ----
    county_srs <- county_noDC %>%
      sample_n(size = 150)
    glimpse(county_srs)

  # State distribution of SRS counties
    county_srs %>%
      group_by(state) %>%
      count()

# Stratified Sample of 150 counties, each state is a stratum ----
  county_str <- county_noDC %>%
      group_by(state) %>%
      sample_n(size = 3)
  # State distribution of stratified sample counties
    glimpse(county_str)

    
    
# Plot of alignment broken down by gender
    ggplot(iris, aes(x = iris$Sepal.Length)) + 
      geom_bar() +
      facet_wrap(~ Species)
  
# PIE CHART ----    
# Put levels of flavor in descending order
    lev <- c("apple", "key lime", "boston creme", "blueberry", "cherry", "pumpkin", "strawberry")
    pies$flavor <- factor(pies$flavor, levels = lev)
    
    # Create barchart of flavor
    ggplot(pies, aes(x = flavor)) + 
      geom_bar(fill = "chartreuse") + 
      theme(axis.text.x = element_text(angle = 90))
    
# CARS DATASET ----
    data(cars)
    str(cars)
# Filter cars with 4, 6, 8 cylinders
    common_cyl <- filter(cars, ncyl %in% c(4,6,8))
    
# Create box plots of city mpg by ncyl
    ggplot(common_cyl, aes(x = as.factor(ncyl), y = city_mpg)) +
      geom_boxplot()
    
# Create overlaid density plots for same data
    ggplot(common_cyl, aes(x = city_mpg, fill = as.factor(ncyl))) +
      geom_density(alpha = .3)

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
 
    
    
library(dplyr) 
library(ggplot2)

# Compute groupwise measures of spread
    county %>%
      group_by(state %in% c("Virginia")) %>%
      summarize(sd(pop2017),
                IQR(pop2017),
                n())

# Generate overlaid density plots
    county %>%
      group_by(state %in% c("Virginia")) %>%
      ggplot(aes(x = pop2017, y = median_edu)) +
      geom_boxplot()
    
    
    
    # Facet hists using hwy mileage and ncyl
    common_cyl %>%
      ggplot(aes(x = hwy_mpg)) +
      geom_histogram() +
      facet_grid(ncyl ~ suv, labeller = label_both) +
      ggtitle("hwy_mpg - ncyl vs. suv")
    
    library(gapminder)
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
    
    # Filter for Asia, add column indicating outliers
    gap_asia <- gap2007 %>%
      filter(continent == "Asia") %>%
      mutate(is_outlier = lifeExp < 50)
    
    # Remove outliers, create box plot of lifeExp
    gap_asia %>%
      filter(!is_outlier) %>%
      ggplot(aes(x = 1, y = lifeExp)) +
      geom_boxplot()
    
library(usdata)
library(ggplot2)
library(dplyr)   

    glimpse(county)
    
    mut_county <-county %>%
      filter(state == "Virginia") %>%
      mutate(is_outlier = pop2010 > 500000)
    
    mut_county %>%
      filter(!is_outlier) %>%
      ggplot(aes(x = per_capita_income , col = name)) +
      geom_density(alpha = .3)
    
    
    
    