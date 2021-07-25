library(ggplot2)
library(dplyr)
library(assertive)


data("storms")

glimpse(storms)


storms %>%
  ggplot(aes(year, fill = status)) +
  geom_density(alpha = 0.5)
