library(dplyr)
library(ggplot2)
library(usa)
library(UScensus2000cdp)
library(UScensus2010)
library(UScensus2000tract)

count(counties)
glimpse(city.name)
glimpse(counties)
glimpse(county.name)
glimpse(people)
data(virginia.tract)

people %>%
  ggplot(aes(edu, fill = gender)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ vote)

table(people$gender, people$race)

people %>%
  ggplot(aes(race, fill = gender)) +
  geom_bar(alpha = 0.3) +
  facet_wrap(~ gender)

virginia.tract %>%
  ggplot(aes(pop2000)) +
  geom_bar(alpha = 0.3) +
  facet_wrap(~ white)
