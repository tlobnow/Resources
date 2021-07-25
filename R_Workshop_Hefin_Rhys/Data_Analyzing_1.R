library(ggplot2)
library(dplyr)
library(openintro)
library(usdata)
Pokemon <- read.csv("Pokemon.csv")
Pokemon
glimpse(county)
str(county)

plot(unemployment_rate, poverty)
us_model <- lm(county$unemployment_rate ~ county$poverty, data = county) 

county %>%
ggplot(aes(log(unemployment_rate), fill = state)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~state)
  
countyWA <- county %>%
  filter(state == "Washington") %>%
  write.csv("countyWA.csv")

plot(countyWA)

countyVA <- county %>%
  filter(state == "Virginia") %>%
  write.csv("countyVA.csv")
plot(countyVA)

countyVA %>%
  ggplot(aes(log(x = per_capita_income))) +
  geom_density(alpha = 0.3)

