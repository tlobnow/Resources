library(ggplot2)
library(dplyr)
library(babynames)
library(gapminder)

# babynames long ----
properties <- c(col.names = "State", "Gender", "Year", "Name", "Number")
AK <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/AK.TXT", sep = ",", col.names = properties)
AL <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/AL.TXT", sep = ",", col.names = properties)
AR <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/AR.TXT", sep = ",", col.names = properties)
AZ <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/AZ.TXT", sep = ",", col.names = properties)
CA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/CA.TXT", sep = ",", col.names = properties)
CO <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/CO.TXT", sep = ",", col.names = properties)
CT <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/CT.TXT", sep = ",", col.names = properties)
DC <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/DC.TXT", sep = ",", col.names = properties)
DE <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/DE.TXT", sep = ",", col.names = properties)
FL <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/FL.TXT", sep = ",", col.names = properties)
GA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/GA.TXT", sep = ",", col.names = properties)
HI <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/HI.TXT", sep = ",", col.names = properties)
IA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/IA.TXT", sep = ",", col.names = properties)
ID <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/ID.TXT", sep = ",", col.names = properties)
IL <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/IL.TXT", sep = ",", col.names = properties)
IN <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/IN.TXT", sep = ",", col.names = properties)
KS <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/KS.TXT", sep = ",", col.names = properties)
KY <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/KY.TXT", sep = ",", col.names = properties)
LA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/LA.TXT", sep = ",", col.names = properties)
MA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/MA.TXT", sep = ",", col.names = properties)
MD <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/MD.TXT", sep = ",", col.names = properties)
ME <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/ME.TXT", sep = ",", col.names = properties)
MI <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/MI.TXT", sep = ",", col.names = properties)
MN <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/MN.TXT", sep = ",", col.names = properties)
MO <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/MO.TXT", sep = ",", col.names = properties)
MS <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/MS.TXT", sep = ",", col.names = properties)
MT <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/MT.TXT", sep = ",", col.names = properties)
NC <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/NC.TXT", sep = ",", col.names = properties)
ND <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/ND.TXT", sep = ",", col.names = properties)
NE <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/NE.TXT", sep = ",", col.names = properties)
NH <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/NH.TXT", sep = ",", col.names = properties)
NJ <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/NJ.TXT", sep = ",", col.names = properties)
NM <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/NM.TXT", sep = ",", col.names = properties)
NV <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/NV.TXT", sep = ",", col.names = properties)
NY <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/NY.TXT", sep = ",", col.names = properties)
OH <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/OH.TXT", sep = ",", col.names = properties)
OK <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/OK.TXT", sep = ",", col.names = properties)
OR <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/OR.TXT", sep = ",", col.names = properties)
PA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/PA.TXT", sep = ",", col.names = properties)
RI <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/RI.TXT", sep = ",", col.names = properties)
SC <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/SC.TXT", sep = ",", col.names = properties)
SD <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/SD.TXT", sep = ",", col.names = properties)
TN <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/TN.TXT", sep = ",", col.names = properties)
TX <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/TX.TXT", sep = ",", col.names = properties)
UT <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/UT.TXT", sep = ",", col.names = properties)
VA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/VA.TXT", sep = ",", col.names = properties)
VT <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/VT.TXT", sep = ",", col.names = properties)
WA <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/WA.TXT", sep = ",", col.names = properties)
WI <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/WI.TXT", sep = ",", col.names = properties)
WV <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/WV.TXT", sep = ",", col.names = properties)
WY <- read.table("~/Documents/Programming/R/DataCamp/namesbystate/WY.TXT", sep = ",", col.names = properties)


all_names <- bind_rows(list(AK, AL, AR, AZ, CA, CO, CT, DC, DE, FL, GA, HI, IA, ID, IL, IN, 
                            KS, KY, LA, MA, MD, ME, MI, MN, MO, MS, MT, NC, ND, NE, NJ, NM, 
                            NV, NY, OH, OK, OR, PA, RI, SC, SD, TN, TX, UT, VA, VT, WA, WI, 
                            WV, WY, id = NULL))

# babynames from package ----

glimpse(babynames)

babynames %>%
  filter(name %in% c("Christin", "Sarah", "Jenny", "Minnie")) %>%
  ggplot(aes(year, fill = name)) +
  geom_density(alpha = 0.3)


# GAPMINDER DATASET ----
library(gapminder)
glimpse(gapminder)

# filter by country
gapminder %>%
  filter(country %in% c("Germany", "United States", "United Kingdom", "Brazil", "Mexico", "Nigeria")) %>%
  ggplot(aes(year, lifeExp, col = country)) +
  geom_line()

glimpse(gapminder)
summary(gapminder$continent)
summary(gapminder$country)
summary(gapminder$pop)
summary(gapminder$year)
summary(gapminder$gdpPercap)
summary(gapminder$lifeExp)

Africa <- gapminder %>%
  filter(continent == "Africa") %>%
  group_by(country) %>%
  top_n(5, lifeExp)
  ggplot(Africa, aes(year, lifeExp, col = country)) +
  geom_line()

Europe <- gapminder %>%
  filter(continent == "Europe") %>%
  group_by(country)
  ggplot(Europe, aes(year, lifeExp, col = country)) +
  geom_line()
  
Asia <- gapminder %>%
  filter(continent == "Asia") %>%
  group_by(country)
  ggplot(Asia, aes(year, lifeExp, col = country)) +
  geom_line()

Americas <- gapminder %>%
  filter(continent == "Americas") %>%
  group_by(country)
ggplot(Americas, aes(year, lifeExp, col = country)) +
  geom_line()

Oceania <- gapminder %>%
  filter(continent == "Oceania") %>%
  group_by(country)
ggplot(Oceania, aes(year, lifeExp, col = country)) +
  geom_line()
  
plot(gapminder)
plot(gapminder$year, gapminder$pop, col = gapminder$conti)

pop_ratio <- gapminder %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(ratio = )
  

gapminder_fraction <- gapminder %>%
  group_by(country) %>%
  mutate(country_total = sum(pop),
         country_max = max(pop),
         fraction = pop / country_total) %>%
  ungroup() %>%
  mutate(fraction_max = pop / country_total)

gapminder_ratio <- gapminder_fraction %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(ratio = fraction / lag(fraction))

gapminder_ratio %>%
  group_by(continent, country) %>%
  filter(country %in% c("Germany", "United States", "China", "Vietnam")) %>%
  ggplot(aes(year, pop, col = country)) +
  geom_point(alpha = 2, shape = continent)

# filter by continent
gapminder %>%
  group_by(continent, country) %>%
  filter(continent == "Asia") %>%
  ggplot(aes(year, lifeExp, col = country)) +
  geom_line()

countries_count <- gapminder %>%
  count(country)
