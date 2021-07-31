#### COVID19_DATA ##############################################################
# Load the readr, ggplot2, and dplyr packages
library(readr)
library(ggplot2)
library(dplyr)

# Read datasets/confirmed_cases_worldwide.csv into confirmed_cases_worldwide
confirmed_cases_worldwide <- read_csv("datasets/confirmed_cases_worldwide.csv")

# See the result
confirmed_cases_worldwide

library(testthat) 
library(IRkernel.testthat)

soln_confirmed_cases_worldwide <- read_csv("datasets/confirmed_cases_worldwide.csv")

run_tests({
  test_that("readr is loaded", {
    expect_true(
      "readr" %in% .packages(), 
      info = "Did you load the `readr` package?"
    )
  })
  test_that("ggplot2 is loaded", {
    expect_true(
      "ggplot2" %in% .packages(), 
      info = "Did you load the `ggplot2` package?"
    )
  })
  test_that("dplyr is loaded", {
    expect_true(
      "dplyr" %in% .packages(), 
      info = "Did you load the `dplyr` package?"
    )
  })
  
  test_that("confirmed_cases_worldwide is a data.frame", {
    expect_s3_class(
      confirmed_cases_worldwide,
      "data.frame",
    )
  })
  test_that("confirmed_cases_worldwide has the correct column", {
    expect_identical(
      colnames(confirmed_cases_worldwide),
      colnames(soln_confirmed_cases_worldwide), 
      info = "The column names of the `confirmed_cases_worldwide` data frame do not correspond with the ones in the CSV file: `\"datasets/confirmed_cases_worldwide.csv\"`."
    ) 
  })
  test_that("has the correct data", {
    expect_equal(
      confirmed_cases_worldwide,
      soln_confirmed_cases_worldwide, 
      info = "The data of the `confirmed_cases_worldwide` data frame do not correspond with data in the CSV file: \"datasets/confirmed_cases_worldwide.csv\"."
    )
  })
})

# Draw a line plot of cumulative cases vs. date
# Label the y-axis
ggplot(confirmed_cases_worldwide, aes(date, cum_cases)) +
  geom_line() +
  ylab("Cumulative confirmed cases")

run_tests({
  plot <- last_plot()
  test_that("the plot is created", {
    expect_false(
      is.null(plot),
      info = "Could not find a plot created with `ggplot()`."
    )
  })
  test_that("the plot uses the correct data", {
    expect_equal(
      plot$data,
      confirmed_cases_worldwide,
      info = "The dataset used in the last plot is not `confirmed_cases_worldwide`."
    )
  })
  test_that("the plot uses the correct x aesthetic", {
    expect_equal(
      quo_name(plot$mapping$x),
      "date",
      info = "The x aesthetic used in the last plot is not `date`."
    )
  })
  test_that("the plot uses the correct y aesthetic", {
    expect_equal(
      quo_name(plot$mapping$y),
      "cum_cases",
      info = "The y aesthetic used in the last plot is not `cum_cases`."
    )
  })
  test_that("the plot uses the correct geom", {
    expect_true(
      "GeomLine" %in% class(plot$layers[[1]]$geom),
      info = "The geom used in the last plot is not `geom_line()`."
    )
  })
  test_that("the plot uses the correct y label", {
    expect_true(
      grepl("[Cc]umulative\\s+[Cc]onfirmed\\s+[Cc]ases", plot$labels$y),
      info = "The y label used in the last plot is not `\"Cumulative confirmed cases\"`."
    )
  })
})

# Read in datasets/confirmed_cases_china_vs_world.csv
confirmed_cases_china_vs_world <- read_csv("datasets/confirmed_cases_china_vs_world.csv")

# See the result
confirmed_cases_china_vs_world

# Draw a line plot of cumulative cases vs. date, colored by is_china
# Define aesthetics within the line geom
plt_cum_confirmed_cases_china_vs_world <- ggplot(confirmed_cases_china_vs_world) +
  geom_line(aes(date, cum_cases, color = is_china)) +
  ylab("Cumulative confirmed cases")

# See the plot
plt_cum_confirmed_cases_china_vs_world

soln_confirmed_cases_china_vs_world <- read_csv("datasets/confirmed_cases_china_vs_world.csv")

run_tests({
  test_that("confirmed_cases_china_vs_world is a data.frame", {
    expect_s3_class(
      confirmed_cases_china_vs_world,
      "data.frame"
    )
  })
  test_that("confirmed_cases_china_vs_world has the correct column names", {
    expect_identical(
      colnames(confirmed_cases_china_vs_world),
      colnames(soln_confirmed_cases_china_vs_world), 
      info = "The column names of the `confirmed_cases_china_vs_world` data frame do not correspond with the ones in the CSV file: `\"datasets/confirmed_cases_china_vs_world.csv\"`."
    ) 
  })
  test_that("confirmed_cases_china_vs_world has the correct data", {
    expect_equal(
      confirmed_cases_china_vs_world,
      soln_confirmed_cases_china_vs_world, 
      info = "The data of the `confirmed_cases_china_vs_world` data frame do not correspond with data in the CSV file: \"datasets/confirmed_cases_china_vs_world.csv\"."
    )
  })
  # NOTE: glimpse is not tested. Can this be done?
  test_that("plt_cum_confirmed_cases_china_vs_world is not NULL", {
    expect_false(
      is.null(plt_cum_confirmed_cases_china_vs_world),
      info = "`plt_cum_confirmed_cases_china_vs_world` is NULL."
    )
  })
  test_that("plt_cum_confirmed_cases_china_vs_world is a plot", {
    expect_true(
      "ggplot" %in% class(plt_cum_confirmed_cases_china_vs_world),
      info = "`plt_cum_confirmed_cases_china_vs_world` is not a `ggplot()` object."
    )
  })
  test_that("plt_cum_confirmed_cases_china_vs_world uses the correct data", {
    expect_equal(
      plt_cum_confirmed_cases_china_vs_world$data,
      confirmed_cases_china_vs_world,
      info = "The dataset used in `plt_cum_confirmed_cases_china_vs_world` is not `confirmed_cases_china_vs_world`."
    )
  })
  layer <- plt_cum_confirmed_cases_china_vs_world$layers[[1]]
  test_that("plt_cum_confirmed_cases_china_vs_world uses uses the correct geom", {
    expect_false(
      is.null(layer),
      info = "The geom used in `plt_cum_confirmed_cases_china_vs_world` is not `geom_line()`."
    )
  })
  test_that("plt_cum_confirmed_cases_china_vs_world uses uses the correct geom", {
    expect_true(
      "GeomLine" %in% class(layer$geom),
      info = "The geom used in `plt_cum_confirmed_cases_china_vs_world` is not `geom_line()`."
    )
  })
  test_that("plt_cum_confirmed_cases_china_vs_world uses uses the correct x aesthetic", {
    expect_equal(
      quo_name(layer$mapping$x),
      "date",
      info = "The x aesthetic used in `plt_cum_confirmed_cases_china_vs_world` is not `date`."
    )
  })
  test_that("plt_cum_confirmed_cases_china_vs_world uses uses the correct y aesthetic", {
    expect_equal(
      quo_name(layer$mapping$y),
      "cum_cases",
      info = "The y aesthetic used in `plt_cum_confirmed_cases_china_vs_world` is not `cum_cases`."
    )
  })
  test_that("plt_cum_confirmed_cases_china_vs_world uses uses the correct color aesthetic", {
    expect_equal(
      quo_name(layer$mapping$colour),
      "is_china",
      info = "The color aesthetic used in `plt_cum_confirmed_cases_china_vs_world` is not `is_china`."
    )
  })
})

who_events <- tribble(
  ~ date, ~ event,
  "2020-01-30", "Global health\nemergency declared",
  "2020-03-11", "Pandemic\ndeclared",
  "2020-02-13", "China reporting\nchange"
) %>%
  mutate(date = as.Date(date))

# Using who_events, add vertical dashed lines with an xintercept at date
# and text at date, labeled by event, and at 100000 on the y-axis
plt_cum_confirmed_cases_china_vs_world +
  geom_vline(data = who_events, aes(xintercept = date), linetype = "dashed") +
  geom_text(data = who_events, 
            aes(x = date, y = 100000, label = event)) 

run_tests({
  plot <- last_plot()
  test_that("the plot got created", {
    expect_false(
      is.null(plot),
      info = "Could not find a plot created with `ggplot()`."
    )
  })
  layer1 <- plot$layers[[2]]
  layer2 <- plot$layers[[3]]
  test_that("the plot has both geoms", {
    expect_false(
      is.null(layer1) || is.null(layer2),
      info = "Could not fin `geom_vline()` and `geom_text()` in your last plot."
    )
  })
  test_that("the plot has both geoms", {
    expect_true(
      "GeomVline" %in% class(layer1$geom) && "GeomText" %in% class(layer2$geom) ||
        "GeomText" %in% class(layer1$geom) && "GeomVline" %in% class(layer2$geom),
      info = "Could not fin `geom_vline()` and `geom_text()` in your last plot."
    )
  })
  if ("GeomVline" %in% class(layer1$geom)) {
    vline <- layer1
    text <- layer2
  } else {
    vline <- layer2
    text <- layer1
  }
  test_that("the plot uses the correct data", {
    expect_equal(
      vline$data,
      who_events,
      info = "The dataset used in the `geom_vline()` is not `who_events`."
    )
  })
  test_that("the geom uses the correct xintercept aesthetic", {
    expect_equal(
      quo_name(vline$mapping$xintercept),
      "date",
      info = "The xintercept aesthetic used in the `geom_vline()` is not `date`."
    )
  })
  test_that("the geom uses the correct lintype parameter", {
    expect_equal(
      vline$aes_params$linetype,
      "dashed",
      info = "The linetype parameter used in the `geom_vline()` is not `\"dashed\"`."
    )
  })
  test_that("the geom uses the correct data", {
    expect_equal(
      text$data,
      who_events,
      info = "The dataset used in the `geom_text()` is not `who_events`."
    )
  })
  test_that("the geom uses the correct x aesthetic", {
    expect_equal(
      quo_name(text$mapping$x),
      "date",
      info = "The x aesthetic used in the `geom_text()` is not `date`."
    )
  })
  test_that("the geom uses the correct label aesthetic", {
    expect_equal(
      quo_name(text$mapping$label),
      "event",
      info = "The label aesthetic used in the `geom_text()` is not `event`."
    )
  })
  if(!is.null(text$aes_params$y)) {
    test_that("the geom uses the correct y parameter", {
      expect_equal(
        text$aes_params$y,
        100000
      )
    })
  } else if (!is.null(quo_name(text$mapping$y))) {
    test_that("the geom uses the correct y parameter", {
      expect_equal(
        quo_name(text$mapping$y),
        '1e+05'
      )
    })
  }
})

# Filter for China, from Feb 15
china_after_feb15 <- confirmed_cases_china_vs_world %>%
  filter(is_china == "China", date >= "2020-02-15")

# Using china_after_feb15, draw a line plot cum_cases vs. date
# Add a smooth trend line using linear regression, no error bars
ggplot(china_after_feb15, aes(date, cum_cases)) +
  geom_line() +
  geom_smooth(method = "lm", se = F) +
  ylab("Cumulative confirmed cases")

run_tests({
  test_that("the data is filtered correctly", {
    soln_china_after_feb15 <- confirmed_cases_china_vs_world %>%
      filter(is_china == "China", date >= "2020-02-15")
    expect_equivalent(
      soln_china_after_feb15,
      china_after_feb15,
      info = "`china_after_feb15` has not been filtered correctly."
    )
  })
  plot <- last_plot()
  test_that("the plot is created", {
    expect_false(
      is.null(plot),
      info = "Could not find a plot created with `ggplot()`."
    )
  })
  test_that("the plot uses the correct data", {
    expect_equal(
      plot$data,
      china_after_feb15,
      info = "The dataset used in the last plot is not `soln_china_after_feb15`."
    )
  })
  test_that("the plot uses the correct x aesthetic", {
    expect_equal(
      quo_name(plot$mapping$x),
      "date",
      info = "The x aesthetic used in the last plot is not `date`."
    )
  })
  test_that("the plot uses the correct y aesthetic", {
    expect_equal(
      quo_name(plot$mapping$y),
      "cum_cases",
      info = "The y aesthetic used in the last plot is not `cum_cases`."
    )
  })
  layer1 <- plot$layers[[1]]
  layer2 <- plot$layers[[2]]
  test_that("the plot has the correct geoms", {
    expect_false(
      is.null(layer1) || is.null(layer2),
      info = "Could not fin `geom_line()` and `geom_smooth()` in your last plot."
    )
  })
  test_that("the plot has the correct geoms", {
    expect_true(
      "GeomLine" %in% class(layer1$geom) && "GeomSmooth" %in% class(layer2$geom) ||
        "GeomSmooth" %in% class(layer1$geom) && "GeomLine" %in% class(layer2$geom),
      info = "Could not fin `geom_line()` and `geom_smooth()` in your last plot."
    )
  })
  if ("GeomLine" %in% class(layer1$geom)) {
    line <- layer1
    smooth <- layer2
  } else {
    line <- layer2
    smooth <- layer1
  }
  test_that("the geom has the correct method parameter", {
    expect_equal(
      smooth$stat_params$method,
      "lm",
      info = "The method parameter used in the `geom_smooth()` is not `\"lm\"`."
      
    )
  })
  test_that("the geom has the correct se parameter", {
    expect_equal(
      smooth$stat_params$se,
      FALSE,
      info = "The se parameter used in the `geom_smooth()` is not `\"FALSE\"`."
    )
  })
})

# Filter confirmed_cases_china_vs_world for not China
not_china <- confirmed_cases_china_vs_world %>% filter(is_china == "Not China")

# Using not_china, draw a line plot cum_cases vs. date
# Add a smooth trend line using linear regression, no error bars
plt_not_china_trend_lin <- ggplot(not_china, aes(date, cum_cases)) +
  geom_line() +
  geom_smooth(method = "lm", se = F) +
  ylab("Cumulative confirmed cases")

# See the result
plt_not_china_trend_lin 

run_tests({
  test_that("the data is filtered correctly", {
    soln_not_china <- confirmed_cases_china_vs_world %>%
      filter(is_china == "Not China")
    expect_equal(
      soln_not_china,
      not_china,
      info = "`not_china` has not been filtered correctly."
    )
  })
  plot <- last_plot()
  test_that("the plot is created", {
    expect_false(
      is.null(plot),
      info = "Could not find a plot created with `ggplot()`."
    )
  })
  test_that("the plot uses the correct data", {
    expect_equal(
      plot$data,
      not_china,
      info = "The dataset used in the last plot is not `not_china`."
    )
  })
  test_that("the plot uses the correct x aesthetic", {
    expect_equal(
      quo_name(plot$mapping$x),
      "date",
      info = "The x aesthetic used in the last plot is not `date`."
    )
  })
  test_that("the plot uses the correct y aesthetic", {
    expect_equal(
      quo_name(plot$mapping$y),
      "cum_cases",
      info = "The y aesthetic used in the last plot is not `cum_cases`."
    )
  })
  layer1 <- plot$layers[[1]]
  layer2 <- plot$layers[[2]]
  test_that("the plot uses the correct geoms", {
    expect_false(
      is.null(layer1) || is.null(layer2),
      info = "Could not fin `geom_line()` and `geom_smooth()` in your last plot."
    )
  })
  test_that("the plot uses the correct geoms", {
    expect_true(
      "GeomLine" %in% class(layer1$geom) && "GeomSmooth" %in% class(layer2$geom) ||
        "GeomSmooth" %in% class(layer1$geom) && "GeomLine" %in% class(layer2$geom),
      info = "Could not fin `geom_line()` and `geom_smooth()` in your last plot."
    )
  })
  if ("GeomLine" %in% class(layer1$geom)) {
    line <- layer1
    smooth <- layer2
  } else {
    line <- layer2
    smooth <- layer1
  }
  test_that("the geom uses the correct method parameter", {
    expect_equal(
      smooth$stat_params$method,
      "lm",
      info = "The method parameter used in the `geom_smooth()` is not `\"lm\"`."
    )
  })
  test_that("the geom uses the correct se parameter", {
    expect_equal(
      smooth$stat_params$se,
      FALSE,
      info = "The se parameter used in the `geom_smooth()` is not `\"FALSE\"`."
    )
  })
})

# Modify the plot to use a logarithmic scale on the y-axis
plt_not_china_trend_lin + 
  scale_y_log10(date)

run_tests({
  plot <- last_plot()
  test_that("the plot is created", {
    expect_false(
      is.null(plot),
      info = "Could not find a plot created with `ggplot()`."
    )
  })
  scale <- plot$scales$get_scales(aes("y"))
  test_that("the plot has a scale", {
    expect_false(
      is.null(scale),
      info = "Could not find a scale in your last plot."
    )
  })
  test_that("the plot uses the correct scale", {
    expect_equal(
      scale$trans$name,
      "log-10",
      info = "Could not find a logarithmic y scale: `scale_y_log10()`."
    )
  })
})

# Run this to get the data for each country
confirmed_cases_by_country <- read_csv("datasets/confirmed_cases_by_country.csv")
glimpse(confirmed_cases_by_country)

# Group by country, summarize to calculate total cases, find the top 7
top_countries_by_total_cases <- confirmed_cases_by_country %>%
  group_by(country) %>%
  summarize(total_cases = max(cum_cases)) %>%
  top_n(7, total_cases)

# See the result
top_countries_by_total_cases

run_tests({
  test_that("the data is manipulated correctly", {
    soln_top_countries_by_total_cases <- confirmed_cases_by_country %>%
      group_by(country) %>%
      summarize(total_cases = max(cum_cases)) %>%
      top_n(7, total_cases)
    expect_equivalent(
      soln_top_countries_by_total_cases,
      top_countries_by_total_cases,
      info = "`top_countries_by_total_cases` has not been filtered correctly."
    )
  })
})

# Read in the dataset from datasets/confirmed_cases_top7_outside_china.csv
confirmed_cases_top7_outside_china <- read_csv("datasets/confirmed_cases_top7_outside_china.csv")

# Glimpse at the contents of confirmed_cases_top7_outside_china
glimpse(confirmed_cases_top7_outside_china)

# Using confirmed_cases_top7_outside_china, draw a line plot of
# cum_cases vs. date, colored by country
ggplot(confirmed_cases_top7_outside_china, aes(date, cum_cases, col = country)) +
  geom_line()

soln_confirmed_cases_top7_outside_china <- read_csv("datasets/confirmed_cases_top7_outside_china.csv")

run_tests({
  test_that('confirmed_cases_top7_outside_china is a data.frame', {
    expect_s3_class(
      confirmed_cases_top7_outside_china,
      'data.frame'
    )
  })
  test_that('confirmed_cases_top7_outside_china had the correct column names', {
    expect_identical(
      colnames(confirmed_cases_top7_outside_china),
      colnames(soln_confirmed_cases_top7_outside_china), 
      info = "The column names of the `confirmed_cases_top7_outside_china` data frame do not correspond with the ones in the CSV file: `\"datasets/confirmed_cases_top7_outside_china.csv\"`."
    ) 
  })
  test_that('confirmed_cases_top7_outside_china had the correct data', {
    expect_equal(
      confirmed_cases_top7_outside_china,
      soln_confirmed_cases_top7_outside_china,
      info = "The data of the `confirmed_cases_top7_outside_china` data frame do not correspond with data in the CSV file: \"datasets/confirmed_cases_top7_outside_china.csv\"."
    )
  })
  # NOTE: glimpse is not tested. Can this be done?
  plot <- last_plot()
  test_that('the plot is created', {
    expect_false(
      is.null(plot),
      info = "Could not find a plot created with `ggplot()`."
    )
  })
  test_that('the plot uses the correct data', {
    expect_equal(
      plot$data,
      confirmed_cases_top7_outside_china,
      info = "The dataset used in the last plot is not `not_china`."
    )
  })
  line <- plot$layers[[1]]
  test_that('the plot uses the correct geom', {
    expect_false(
      is.null(line),
      info = "Could not fin `geom_line()` in your last plot."
    )
  })
  test_that('the plot uses the correct geom', {
    expect_true(
      'GeomLine' %in% class(line$geom),
      info = "Could not fin `geom_line()` in your last plot."
    )
  })
  mapping <- plot$mapping
  geom_mapping <- line$mapping
  test_that('the plot uses the correct x aesthetic', {
    expect_true(
      !is.null(mapping$x) && quo_name(mapping$x) == "date" ||
        !is.null(geom_mapping$x) && quo_name(geom_mapping$x) == "date",
      info = "The x aesthetic used in the last plot is not `date`."
      
    )
  })
  test_that('the plot uses the correct y aesthetic', {
    expect_true(
      !is.null(mapping$y) && quo_name(mapping$y) == "cum_cases" ||
        !is.null(geom_mapping$y) && quo_name(geom_mapping$y) == "cum_cases",
      info = "The y aesthetic used in the last plot is not `cum_cases`."
    )
  })
  test_that('the plot uses the correct color aesthetic', {
    expect_true(
      !is.null(mapping$colour) && quo_name(mapping$colour) == "country" ||
        !is.null(geom_mapping$colour) && quo_name(geom_mapping$colour) == "country",
      info = "The color aesthetic used in the last plot is not `country`."
    )
  })
})

#### SEMMELWEISS ##############################################################

# Load in the tidyverse package
library(tidyverse)

# Read datasets/yearly_deaths_by_clinic.csv into yearly
yearly <- read_csv("datasets/yearly_deaths_by_clinic.csv")

# Print out yearly
print(yearly)

library(testthat) 
library(IRkernel.testthat)
run_tests({
  test_that("Read in data correctly.", {
    expect_is(yearly, "tbl_df", 
              info = 'You should use read_csv (with an underscore) to read "datasets/yearly_deaths_by_clinic.csv" into yearly.')
  })
  
  test_that("Read in data correctly.", {
    yearly_temp <- read_csv('datasets/yearly_deaths_by_clinic.csv')
    expect_equivalent(yearly, yearly_temp, 
                      info = 'yearly should contain the data in "datasets/yearly_deaths_by_clinic.csv"')
  })
})

# Adding a new column to yearly with proportion of deaths per no. births
yearly <- yearly %>%
  mutate(proportion_deaths = deaths / births)
# Print out yearly
yearly

run_tests({
  test_that("A proportion_deaths column exists", {
    expect_true("proportion_deaths" %in% names(yearly), 
                info = 'yearly should have the new column proportion_deaths')
  })
  
  test_that("Read in data correctly.", {
    yearly_temp <- read_csv('datasets/yearly_deaths_by_clinic.csv') %>% 
      mutate(proportion_deaths = deaths / births)
    expect_equivalent(yearly, yearly_temp, 
                      info = 'proportion_deaths should be calculated as deaths / births')
  })
})

# Setting the size of plots in this notebook
options(repr.plot.width=7, repr.plot.height=4)

# Plot yearly proportion of deaths at the two clinics
ggplot(yearly, aes(proportion_deaths, year, color = clinic)) +
  geom_line()

run_tests({
  test_that("The right columns are plotted", {
    mappings <- str_replace(as.character(last_plot()$mapping), "~", "")
    expect_true(all(c("year", "proportion_deaths", "clinic") %in% mappings),
                info = "year should be on the x-axis, proportion_deaths should be on the y-axis, and clinic should be mapped to color.")
  })
})


# Read datasets/monthly_deaths.csv into monthly
monthly <- read_csv("datasets/monthly_deaths.csv")

# Adding a new column with proportion of deaths per no. births
monthly <- monthly %>% mutate(proportion_deaths = deaths / births)

# Print out the first rows in monthly
head(monthly)

run_tests({
  
  test_that("Read in data correctly.", {
    expect_is(monthly, "tbl_df", 
              info = 'You should use read_csv (with an underscore) to read "datasets/monthly_deaths.csv" into monthly.')
  })
  
  test_that("Read in monthly correctly.", {
    monthly_temp <- read_csv("datasets/monthly_deaths.csv")
    expect_true(all(names(monthly_temp) %in% names(monthly)), 
                info = 'monthly should contain the data in "datasets/monthly_deaths.csv"')
  })
  
  test_that("proportion_death is calculated correctly.", {
    monthly_temp <- read_csv("datasets/monthly_deaths.csv")
    monthly_temp <- monthly_temp %>% 
      mutate(proportion_deaths = deaths / births)
    expect_equivalent(monthly, monthly_temp, 
                      info = 'proportion_deaths should be calculated as deaths / births')
  })
})

# Plot monthly proportion of deaths
ggplot(monthly, aes(date, proportion_deaths)) +
  geom_line() +
  labs(x = "Date",
       y = "Proportion of Deaths")

run_tests({
  test_that("The right columns are plotted", {        
    mappings <- str_replace(as.character(last_plot()$mapping), "~", "")
    expect_true(all(c("date", "proportion_deaths") %in% mappings), 
                info = "date should be on the x-axis, proportion_deaths on the y-axis")
  })
})

# From this date handwashing was made mandatory
handwashing_start = as.Date('1847-06-01')

# Add a TRUE/FALSE column to monthly called handwashing_started
monthly <- monthly %>% mutate(handwashing_started = monthly$date >= '1847-06-01')

# Plot monthly proportion of deaths before and after handwashing
ggplot(monthly, aes(date, proportion_deaths, color = handwashing_started)) +
  geom_line() +
  labs(x = "Date",
       y = "Proportion of Deaths")

run_tests({
  test_that("handwashing_started has been defined", {
    expect_true("handwashing_started" %in% names(monthly),
                info = 'monthly should contain the column handwashing_started.')
  })  
  
  test_that("there are 22 rows where handwashing_started is TRUE", {
    expect_equal(22, sum(monthly$handwashing_started),
                 info = 'handwashing_started should be a TRUE/FALSE column where the rows where handwashing was enforced are set to TRUE.')
  })
  
  test_that("The right columns are plotted", {        
    mappings <- str_replace(as.character(last_plot()$mapping), "~", "")
    expect_true(all(c("date", "proportion_deaths", "handwashing_started") %in% mappings), 
                info = 'date should be on the x-axis, proportion_deaths on the y-axis, and handwashing_started should be mapped to color.')
  })
})

# Calculating the mean proportion of deaths 
# before and after handwashing.

monthly_summary <- monthly %>%
  group_by(handwashing_started) %>%
  summarise(mean_proportion_deaths = mean(proportion_deaths))

# Printing out the summary.
monthly_summary

run_tests({
  test_that("mean_proportion_deaths was calculated correctly", {
    flat_summary <- as.numeric(unlist(monthly_summary))
    handwashing_start = as.Date('1847-06-01')
    monthly_temp <- read_csv("datasets/monthly_deaths.csv") %>% 
      mutate(proportion_deaths = deaths / births) %>% 
      mutate(handwashing_started = date >= handwashing_start) %>% 
      group_by(handwashing_started) %>%
      summarise(mean_proportion_deaths = mean(proportion_deaths))
    expect_true(all(monthly_temp$mean_proportion_deaths %in% flat_summary),
                info = 'monthly_summary should contain the mean monthly proportion of deaths before and after handwashing was enforced.')
  })  
})

# Calculating a 95% Confidence intrerval using t.test 
test_result <- t.test( proportion_deaths ~ handwashing_started, data = monthly)
test_result

run_tests({
  test_that("the confidence intervals match", {
    temp_test_result <- t.test( proportion_deaths ~ handwashing_started, data = monthly)
    expect_equivalent(test_result$conf.int, temp_test_result$conf.int,
                      info = 'The t-test should be calculated with proportion_deaths as a function of handwashing_started.')
  })  
})

# The data Semmelweis collected points to that:
doctors_should_wash_their_hands <- TRUE

run_tests({
  test_that("The project is finished.", {
    expect_true(doctors_should_wash_their_hands, 
                info = "Semmelweis would argue that doctors_should_wash_their_hands should be TRUE .")
  })  
})

#### PHYLLOTAXIS ###############################################################

# Set plot images to a nice size
options(repr.plot.width = 4, repr.plot.height = 4)

# Load the ggplot2 package
library(ggplot2)

library(testthat) 
library(IRkernel.testthat)

run_tests({
  test_that("Test that ggplot2 is loaded", {
    expect_true( "package:ggplot2" %in% search(), 
                 info = "The ggplot2 package should be loaded using library().")
  })
})

# Create circle data to plot
t <- seq(0, 2*pi, length.out = 50)
x <- sin(t)
y <- cos(t)
df <- data.frame(t, x, y)

# Make a scatter plot of points in a circle
p <- ggplot(df, aes(x, y))
p + geom_point()

run_tests({
  test_that("Check that a geom_point plot was plotted.", {
    expect_true( "GeomPoint" %in% class( last_plot()$layers[[1]]$geom ) , 
                 info = "Add geom_point() to produce a scatter plot.")
  })
})

# Define the number of points
points <- 500

# Define the Golden Angle
angle <- pi * (3 - sqrt(5))

t <- (1:points) * angle
x <- sin(t)
y <- cos(t)
df <- data.frame(t, x, y)

# Make a scatter plot of points in a spiral
p <- ggplot(df, aes(x*t, y*t))
p + geom_point()

run_tests({
  test_that("points are 500.", {
    expect_equal(points, 500, 
                 info = "There should be 500 points.")
  })
  
  test_that("angle is golden.", {
    expect_equal(angle, pi*(3-sqrt(5)), 
                 info = "angle should be set to the Golden Angel. Check the hint!")
  })
})

df <- data.frame(t, x, y)

# Make a scatter plot of points in a spiral and remove some plot components
p <- ggplot(df, aes(x*t, y*t))
p + geom_point() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank())

run_tests({
  test_that("Background is white.", {
    expect_equal(last_plot()$theme$panel.background$fill, "white", 
                 info = "The background should be white.")
  })
  test_that("Plot components are removed.", {
    expect_true("element_blank" %in% class(last_plot()$theme$panel.grid), 
                info = "The grid lines should be removed.")
    expect_true("element_blank" %in% class(last_plot()$theme$axis.ticks), 
                info = "The axis ticks should be removed.")
    expect_true("element_blank" %in% class(last_plot()$theme$axis.title), 
                info = "The axis titles should be removed.")
    expect_true("element_blank" %in% class(last_plot()$theme$axis.text), 
                info = "The axis text should be removed.")        
  })
})

# Change the code from Task 4 to modify the 
# size, transparency, and color of the points
p <- ggplot(df, aes(x*t, y*t))
p + geom_point(size = 8,
               alpha = 0.5,
               color = "darkgreen") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank())


run_tests({
  test_that("Point size equal to 8.", {
    expect_equal(last_plot()$layers[[1]]$aes_params$size, 8, 
                 info = "size should be set 8.")
  })
  test_that("alpha equal to 0.5.", {
    expect_equal(last_plot()$layers[[1]]$aes_params$alpha, 0.5, 
                 info = "alpha should be set 0.5.")
  })
})

# Copy the code from Task 5 and modify the 
# color, size, and shape of the points
p <- ggplot(df, aes(x*t, y*t))
p + geom_point(size = 8,
               alpha = 0.5,
               shape = 8,
               aes(size = t)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        legend.position = "none")

run_tests({
  test_that("Map size of points to t.", {
    expect_equal(last_plot()$labels$size, "t", 
                 info = "Map size of points to t. Check the hint!")
  })
  test_that("point shape is asterisk.", {
    expect_equal(last_plot()$layers[[1]]$aes_params$shape, 8, 
                 info = "Change the shape of all points to asterisks.")
  })
  test_that("Legend is removed.", {
    expect_equal(last_plot()$theme$legend.position, "none", 
                 info = "Remove the legend from the plot.")
  })
})

# Copy the code from Task 6 and modify the color and
# shape of the points, and the background color
p <- ggplot(df, aes(x*t, y*t, colour = "yellow"))
p + geom_point(aes(size=t), 
               alpha=0.5, 
               shape=17, 
               color="yellow") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "darkmagenta"),
        panel.grid = element_blank(),
        legend.position = "none")

run_tests({
  test_that("point shape is filled triangles.", {
    expect_equal(last_plot()$layers[[1]]$aes_params$shape, 17, 
                 info = "Change the shape of all points to filled triangles. Check the hint.")
  })
  test_that("The triangles are yellow", {
    expect_equal(last_plot()$layers[[1]]$aes_params$colour, "yellow", 
                 info = "The triangles are not yellow. Check the hint.")
  })
  test_that("The background is dark magenta", {
    expect_equal(last_plot()$theme$panel.background$fill, "darkmagenta", 
                 info = "The background is not dark magenta. Check the hint.")
  })
})

# Change the value of the angle
angle <- 2.0
points <- 1000

t <- (1:points)*angle
x <- sin(t)
y <- cos(t)
df <- data.frame(t, x, y)

# Copy the plotting code from Task 7
p <- ggplot(df, aes(x*t, y*t, colour = "yellow"))
p + geom_point(aes(size=t), 
               alpha=0.5, 
               shape=17, 
               color="yellow") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "darkmagenta"),
        panel.grid = element_blank(),
        legend.position = "none")

run_tests({
  test_that("angle is 2.", {
    expect_equal(angle, 2, 
                 info = "angle should be equal to 2")
  })
})

# Change the values of angle and points
angle <- 13 * pi / 180
points <- 2000

t <- (1:points)*angle
x <- sin(t)
y <- cos(t)
df <- data.frame(t, x, y)

# Adjust the plot parameters to create the magenta flower
p <- ggplot(df, aes(x*t, y*t))
p + geom_point(size = 80, 
               alpha = 0.1, 
               shape = 1, 
               color = "magenta4")+
  theme(legend.position="none",
        panel.background = element_rect(fill = "white"),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank())

run_tests({
  test_that("points is equal to 2000.", {
    expect_equal(points, 2000, 
                 info = "There should be 2000 points.")
  })
  test_that("point shape is empty circle.", {
    expect_equal(last_plot()$layers[[1]]$aes_params$shape, 1, 
                 info = "Change the shape of all points to empty circles. Check the hint!")
  })
  test_that("alpha is equal 0.1", {
    expect_equal(last_plot()$layers[[1]]$aes_params$alpha, 0.1, 
                 info = "alpha of points should be 0.1")
  })
  test_that("Background is white.", {
    expect_equal(last_plot()$theme$panel.background$fill, "white", 
                 info = "The background should be white.")
  })
  test_that("angle is 13*pi/180.", {
    expect_equal(angle, 13*pi/180, 
                 info = "angle should be set to 13*pi/180.")
  })
})

#### KAGGLE DATA SURVEY ########################################################

# Load necessary packages
library(tidyverse)

# Load the data
responses <- read_csv("datasets/kagglesurvey.csv")

# Print the first 10 rows
head(responses)

library("testthat")
library('IRkernel.testthat')

run_tests({
  test_that("Read in data correctly.", {
    expect_is(responses, "tbl_df", 
              info = 'You should use read_csv() (with an underscore) to read "datasets/kagglesurvey.csv" into responses.')
  })
  
  test_that("Read in data correctly.", {
    responses_test <- read_csv('datasets/kagglesurvey.csv')
    expect_equivalent(responses, responses_test, 
                      info = 'responses should contain the data in "datasets/kagglesurvey.csv".')
  })
  
})

# Print the first respondent's tools and languages
responses[1, 2]

# Add a new column, and unnest the new column
tools <- responses  %>% 
  mutate(work_tools = str_split(WorkToolsSelect, ","))  %>% 
  unnest(work_tools)

# View the first 6 rows of tools
head(tools, 6)

run_tests({
  test_that("Tools and Languages were Split and Unnested", {
    expect_true(nrow(tools) == 47409, 
                info = 'Make sure that you split the tools at the commas and unnested them.')
  })
  
  test_that("Tools and Languages were Unnested", {
    expect_is(tools$work_tools, "character", 
              info = 'The work_tools column should be of class "character". Make sure that you unnested the results of str_split().')
  })
  
})

# Group the data by work_tools, summarise the counts, and arrange in descending order
tool_count <- tools  %>% 
  group_by(work_tools)  %>% 
  summarise(count = n()) %>%
  arrange(desc(count))

# Print the first 6 results
head(tool_count)

run_tests({
  test_that("Tools were Grouped and Summarised", {
    expect_true(nrow(tool_count) == 50, 
                info = 'Make sure that you grouped by tools and then summarised the counts.')
  })
  
  test_that("Values were sorted correctly", {
    expect_true(tool_count[1, 2] == 6073, 
                info = 'Do not forget to sort your tool counts from largest to smallest.')
  })
  
})

# Create a bar chart of the work_tools column, most counts on the far right
ggplot(tool_count, aes(x = fct_reorder(a, b), y = b)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1))

run_tests({
  test_that("Plot is a bar chart",{
    p <- last_plot()
    q <- p$layers[[1]]
    expect_is(q$geom, "GeomBar", 
              info = "You should plot a bar chart with ggplot().")
  })
})

# Create a new column called language preference
debate_tools <- responses  %>% 
  mutate(language_preference = case_when(
    str_detect(WorkToolsSelect, "R") & ! str_detect(WorkToolsSelect, "Python") ~ "R",
    str_detect(WorkToolsSelect, "Python") & ! str_detect(WorkToolsSelect, "R")   ~ "Python",
    str_detect(WorkToolsSelect, "R") &  str_detect(WorkToolsSelect, "Python") ~ "both",
    
    TRUE ~ "neither"
  ))

# Print the first 6 rows
head(debate_tools)

debate_tools_counts <- debate_tools %>% 
  count(language_preference)

run_tests({
  test_that("New column was created", {
    expect_is(debate_tools$language_preference, "character", 
              info = 'The language_preference column should be of class "character". Make sure that you filled this new column correctly.')
  })
  test_that("Language preferences are correct", {
    expect_equal(filter(debate_tools_counts, language_preference == "both")  %>% pull(n), 3660, 
                 info = 'There is an incorrect amount of "both". Please check the case_when() statements.')
    expect_equal(filter(debate_tools_counts, language_preference == "neither")  %>% pull(n), 2860, 
                 info = 'There is an incorrect amount of "neither". Please check the case_when() statements.')
    expect_equal(filter(debate_tools_counts, language_preference == "Python")  %>% pull(n), 2413, 
                 info = 'There is an incorrect amount of "Python". Please check the case_when() statements.')
    expect_equal(filter(debate_tools_counts, language_preference == "R")  %>% pull(n), 1220, 
                 info = 'There is an incorrect amount of "R". Please check the case_when() statements.')
    
  })
  
})

# Group by language preference, calculate number of responses, and remove "neither"
debate_plot <- debate_tools  %>% 
  group_by(language_preference)  %>% 
  summarise(count = n()) %>% 
  filter(language_preference != "neither")

# Create a bar chart
ggplot(debate_plot, aes(language_preference)) +
  geom_bar(stat = "identity")

run_tests({
  test_that("Plot is a bar chart",{
    p <- last_plot()
    q <- p$layers[[1]]
    expect_is(q$geom, "GeomBar",
              info = "You should plot a bar chart with ggplot().")
  })
})

# Group by, summarise, arrange, mutate, and filter
recommendations <- debate_tools  %>% 
  group_by(language_preference, LanguageRecommendationSelect)  %>% 
  summarise(count = n())  %>% 
  arrange(language_preference, desc(count)) %>%
  mutate(row = row_number()) %>%
  filter(row <= 4)

run_tests({
  test_that("Tools have been summarised", {
    expect_true(nrow(recommendations) == 16, 
                info = 'Make sure that you are only keeping the top 4 responses for each language used.')
  })
  
})

# Create a faceted bar plot
ggplot(recommendations, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  facet_wrap(~language_preference)


run_tests({
  test_that("Plot is a bar chart",{
    p <- last_plot()
    q <- p$layers[[1]]
    expect_is(q$geom, "GeomBar",
              info = "You should plot a bar chart with ggplot().")
  })
})

# Would R users find this statement TRUE or FALSE?
R_is_number_one = TRUE

run_tests({
  test_that("The question has been answered", {
    expect_true(R_is_number_one, 
                info = 'Try again! Should R_is_number_one be set to TRUE or FALSE?')
  })
  
})

#### CANDYCRUSH DATA ###########################################################


# This sets the size of plots to a good default.
options(repr.plot.width = 5, repr.plot.height = 4)

# Loading in packages
library(readr)
library(dplyr)
library(ggplot2)

library(testthat) 
library(IRkernel.testthat)

run_tests({
  test_that("the packages are loaded", {
    expect_true( all(c("package:ggplot2", "package:readr", "package:dplyr") %in% search() ), 
                 info = "The dplyr, readr and ggplot2 packages should be loaded using library().")
  })
})

# Reading in the data
data <- read_csv("datasets/candy_crush.csv")

# Printing out the first six rows
head(data, 6)

run_tests({
  test_that("data is read in correctly", {
    correct_data <- read_csv("datasets/candy_crush.csv")
    expect_equal(correct_data, data, 
                 info = "data should countain datasets/candy_crush.csv read in using read_csv")
  })
})

# Count and display the number of unique players
print("Number of players:")
count(data, player_id, sort = T) %>% unique()

# Display the date range of the data
print("Period for which we have data:")
range(data$date)



run_tests({
  test_that("nothing", {
    expect_true(TRUE, info = "")
  })
})

# Calculating level difficulty
difficulty <- data %>%
  group_by(level) %>%
  summarise(attempts = sum(num_attempts), wins = sum(num_success)) %>%
  mutate(p_win = wins / attempts)

# Printing out the level difficulty
print(difficulty)

run_tests({
  test_that("p_win is calculated correctly", {
    correct_difficulty <- data %>%
      group_by(level) %>%
      summarise(attempts = sum(num_attempts), wins = sum(num_success)) %>%
      mutate(p_win = wins / attempts)
    expect_equal(correct_difficulty$p_win, difficulty$p_win, 
                 info = "difficulty$p_win should be estimated probability to pass each level in a single attempt")
  })
})

# Plotting the level difficulty profile
ggplot(difficulty, aes(level, p_win)) +
  geom_line() +
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(label = scales::percent)



run_tests({
  test_that("the student plotted a ggplot", {
    expect_true('ggplot' %in% class(last_plot()), 
                info = "You should plot difficulty using ggplot.")
  })
})

# Adding points and a dashed line
ggplot(difficulty, aes(level, p_win)) +
  geom_line() +
  geom_point() +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(label = scales::percent)



run_tests({
  plot_layers <- sapply(last_plot()$layers, function(layer)  class(layer$geom)[1])
  test_that("the student has plotted lines, points and a hline", {
    expect_true(all(c('GeomLine', 'GeomPoint', 'GeomHline') %in%  plot_layers), 
                info = "The plot should include lines between the datapoints, points at the datapoints and a horisontal line.")
  })
})

# Computing the standard error of p_win for each level
difficulty <- difficulty %>%
  mutate(error = sqrt(p_win * (1- p_win) / attempts))

run_tests({
  test_that("error is correct", {
    correct_difficulty <- difficulty %>%
      mutate(error = sqrt(p_win * (1 - p_win) / attempts))
    expect_equal(correct_difficulty$error, difficulty$error,
                 info = "difficulty$error should be calculated as sqrt(p_win * (1 - p_win) / attempts)")
  })
})

# Adding standard error bars
ggplot(difficulty, aes(level, p_win)) +
  geom_line() +
  geom_point() +
  geom_hline(aes(yintercept = 0.1), linetype = "dashed") +
  geom_errorbar(aes(level, ymin = p_win - error, ymax = p_win + error)) +
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(label = scales::percent)


run_tests({
  plot_layers <- sapply(last_plot()$layers, function(layer)  class(layer$geom)[1])
  test_that("the student has plotted lines, points and a hline", {
    expect_true("GeomErrorbar" %in%  plot_layers, 
                info = "The plot should include error bats using geom_errorbar.")
  })
})

# The probability of completing the episode without losing a single time
p <- prod(difficulty$p_win)

# Printing it out
p

run_tests({
  test_that("p is correct", {
    correct_p <- prod(difficulty$p_win)
    expect_equal(correct_p, p,
                 info = "p should be calculated as the product of difficulty$p_win .")
  })
})

# Should our level designer worry about that a lot of 
# players will complete the episode in one attempt?
should_the_designer_worry = FALSE # TRUE / FALSE

run_tests({
  test_that("should_the_designer_worry is FALSE", {
    expect_false(should_the_designer_worry,
                 info = "The probability is really small, so I don't think the designer should worry that much...")
  })
})

#### MUSIC DATA ################################################################

# Loading individual Tidyverse packages
library(dplyr)
library(readr)
library(ggplot2)

# Reading in the McGill Billboard chord data
bb <- read_csv("datasets/bb_chords.csv")

# Taking a look at the first rows in bb
head(bb)

# These packages need to be loaded in the first `@tests` cell. 
library(testthat) 
library(IRkernel.testthat)

run_tests({
  test_that("Read in data correctly.", {
    expect_is(bb, "tbl_df", 
              info = 'You should use read_csv (with an underscore) to read "datasets/bb_chords.csv" into bb')
  })
  
  test_that("Read in data correctly.", {
    bb_correct <- read_csv('datasets/bb_chords.csv')
    expect_equivalent(bb, bb_correct, 
                      info = 'bb should contain the data in "datasets/bb_chords.csv"')
  })
})

# Counting the most common chords
bb_count <- bb %>% count(chord, sort = T)

# Displaying the top 20 chords
bb_count %>% slice(1:20)



run_tests({
  test_that("bb_count is correct", {
    correct_bb_count <- bb %>%
      count(chord, sort = TRUE)
    expect_equivalent(bb_count, correct_bb_count, 
                      info = "bb_count should contain the count of each type of chord.")
  })
})

# Creating a bar plot from bb_count
bb_count %>%
  slice(1:20) %>%
  mutate(share = n / sum(n)*100,
         chord = reorder(chord, share))%>%
  ggplot(aes(chord, share)) +
  geom_col(aes(fill = chord)) +
  coord_flip() +
  xlab("Share of total chords") +
  ylab("Chord") +
  theme(legend.position = 'none')





run_tests({
  test_that("bb_count has some data in it", {
    expect_true(length(bb_count) > 0, 
                info = "Looks like you're missing data in `bb_count`.")
  })
})

# Wrangling and counting bigrams
bb_bigram_count <- bb %>%
  mutate(next_chord = lead(chord),
         next_title = lead(title),
         bigram = paste(chord,next_chord)) %>%
  filter(title == next_title) %>%
  count(bigram, sort = T)


# Displaying the first 20 rows of bb_bigram_count
bb_bigram_count %>% slice(1:20)





run_tests({
  test_that("bb_bigram_count is correct", {
    correct_bb_bigram_count <- bb %>%
      mutate(next_chord = lead(chord),
             next_title = lead(title),
             bigram = paste(chord, next_chord)) %>%
      filter(title == next_title) %>%
      count(bigram, sort = TRUE)
    expect_equivalent(bb_bigram_count, correct_bb_bigram_count, 
                      info = "`bb_bigram_count` should contain the count of each type of bigram. Don't forget to sort by bigram frequency!")
  })
})

# Creating a column plot from bb_bigram_count
bb_bigram_count %>%
  slice(1:20) %>%
  mutate(share = n / sum(n),
         bigram = reorder(bigram, share))%>%
  ggplot(aes(bigram, share)) +
  geom_col(aes(fill = bigram)) +
  coord_flip() +
  xlab("Change in Chords") +
  ylab("Count") +
  theme(legend.position = 'none')



run_tests({
  test_that("bb_bigram_count has some data in it", {
    expect_true(length(bb_bigram_count) > 0, 
                info = "Looks like you're missing data in `bb_bigram_count`.")
  })
})

# Finding 30 artists with the most songs in the corpus
bb_30_artists <- bb %>%
  select(artist, title) %>%
  unique() %>%
  count(artist, sort = T)

# Displaying 30 artists with the most songs in the corpus
bb_30_artists %>% slice(1:30)





run_tests({
  test_that("bb artists counted and sorted", {
    correct_bb_30_artists <- bb %>%
      select(artist, title) %>%
      unique() %>%
      count(artist, sort = TRUE)
    expect_equivalent(bb_30_artists, correct_bb_30_artists, 
                      info = "`bb_30_artists` should contain the number of soungs (not chords) by each artist in the corpus. Don't forget to sort!")
  })
})

tags <- tibble(
  artist = c('Abba', 'Billy Joel', 'Elton John', 'Stevie Wonder', 'The Rolling Stones', 'The Beatles', 'Eric Clapton'),
  instrument = c('piano', 'piano', 'piano', 'piano', 'guitar', 'guitar', 'guitar'))

# Creating a new dataframe bb_tagged that includes a new column instrument from tags
bb_tagged <- bb %>%
  inner_join(tags)

# Displaying the new data frame
bb_tagged

run_tests({
  test_that("bb artists counted and sorted", {
    correct_bb_tagged <- bb %>%
      inner_join(tags)
    expect_equivalent(bb_tagged, correct_bb_tagged, 
                      info = "`bb_tagged` should be a successful join of `bb` and `tags` that only contains records cointained in both dataframes.")
  })
})

# The top 20 most common chords
top_20 <- bb_count$chord[1:20]

# Comparing the frequency of the 20 most common chords in piano- and guitar-driven songs
bb_tagged %>%
  filter(chord %in% top_20) %>%
  count(chord, instrument, sort = T) %>%
  ggplot(aes(chord, n, fill = chord)) +
  geom_bar(stat = "identity") +
  facet_wrap(~instrument) +
  coord_flip() +
  xlab("chord") +
  ylab("count") 


run_tests({
  test_that("bb_tagged has some data in it", {
    expect_true(length(bb_tagged) > 0, 
                info = "Looks like you're missing data in `bb_tagged`.")
  })
})

# The top 20 most common bigrams
top_20_bigram <- bb_bigram_count$bigram[1:20]

# Creating a faceted plot comparing guitar- and piano-driven songs for bigram frequency
bb_tagged %>%
  mutate(next_chord = lead(chord),
         next_title = lead(title),
         bigram = paste(chord, next_chord)) %>%
  filter(title == next_title) %>%
  count(bigram, instrument, sort = T) %>%
  filter(bigram %in% top_20_bigram) %>%
  ggplot(aes(bigram, n, fill = bigram)) +
  geom_col() +
  facet_wrap(~instrument) +
  coord_flip() +
  xlab("Total bigrams") +
  ylab("Bigram") +
  theme(legend.position = 'none')



run_tests({
  test_that("bb_bigram_count has some data in it", {
    expect_true(length(bb_bigram_count) > 0, 
                info = "Looks like you're missing data in `bb_bigram_count`.")
  })
})

# Set to TRUE or FALSE to reflect your answer
hypothesis_valid <- TRUE

# Set to TRUE or FALSE to reflect your answer
more_data_needed <- TRUE

run_tests({
  test_that("hypothesis is true", {
    expect_true(hypothesis_valid, 
                info = "Are you sure the hypothesis isn't valid?!")
  })
  test_that("more_data_needed is true", {
    expect_true(more_data_needed, 
                info = "Are you sure we don't need more data?!")
  })
})


#### TAXI DATA #################################################################

# Loading the tidyverse
library(tidyverse)

# Reading in the taxi data
taxi <- read_csv("datasets/taxi.csv")

# Taking a look at the first few rows in taxi
head(taxi)

library(testthat) 
library(IRkernel.testthat)

run_tests({
  test_that("Test that tidyverse is loaded", {
    expect_true( "package:tidyverse" %in% search(), 
                 info = "The tidyverse package should be loaded using library().")
  })
  
  test_that("Read in data correctly.", {
    expect_is(taxi, "tbl_df", 
              info = 'You should use read_csv (with an underscore) to read "datasets/taxi.csv" into taxi.')
  })
  
  test_that("Read in data correctly.", {
    taxi_temp <- read_csv('datasets/taxi.csv')
    expect_equivalent(taxi, taxi_temp, 
                      info = 'taxi should contain the data in "datasets/taxi.csv".')
  })
})

glimpse(taxi)
# Renaming the location variables,
# dropping any journeys with zero fares and zero tips,
# and creating the total variable as the log sum of fare and tip
taxi <- taxi %>%
  rename(long = pickup_longitude, lat = pickup_latitude)  %>% 
  filter(fare_amount | tip_amount > 0) %>%
  mutate(total = log(fare_amount + tip_amount))




run_tests({
  test_that("rename lat", {
    expect_true(!is.null(taxi$lat), 
                info = "The taxi data frame does not contain a variable called lat. You need to rename pickup_latitude.")
  })
  test_that("rename long", {
    expect_true(!is.null(taxi$long), 
                info = "The taxi data frame does not contain a variable called long. You need to rename pickup_longitude.")
  })
  test_that("total exists", {
    expect_true(!is.null(taxi$total), 
                info = "The taxi data frame does not contain a variable called total. You need to create this as the logarithm (use the log() function) of the sum of fare_amount and tip_amount.")
  })
  test_that("Modified data correctly.", {
    taxi_temp <- read_csv('datasets/taxi.csv') %>%
      rename(long = pickup_longitude, lat = pickup_latitude)  %>% 
      filter(fare_amount > 0 | tip_amount > 0) %>%
      mutate(total = log(fare_amount + tip_amount) )
    expect_equivalent(taxi, taxi_temp, 
                      info = 'The taxi dataframe has not been modified correctly. See if you can find something is wrong with your code.')
  })
})


# Reducing the data to taxi trips starting in Manhattan
# Manhattan is bounded by the rectangle with 
# latitude from 40.70 to 40.83 and 
# longitude from -74.025 to -73.93
taxi <- taxi  %>% 
  filter(between(lat, 40.70, 40.83),
         between(long, -74.025, -73.93))

run_tests({
  test_that("The correct number of rows have been filtered away", {
    expect_equal(45766, nrow(taxi), 
                 info = "It seems you haven't filter away the taxi trips outside of Manhattan correctly.")
  })
})

# Loading in ggmap and viridis for nice colors
library(ggmap)
library(viridis)

# Retrieving a stored map object which originally was created by
# manhattan <- get_map("manhattan", zoom = 12, color = "bw")
manhattan <- readRDS("datasets/manhattan.rds")

# Drawing a density map with the number of journey start locations
ggmap(manhattan, darken = 0.5) +
  scale_fill_viridis(option = 'plasma') +
  geom_bin2d(data = taxi, aes(long, 
                              lat), 
             alpha = 0.6,
             bins = 60) +
  labs(x = "Longitude",
       y = "Latitude")


run_tests({
  
  test_that("Test that ggmap is loaded", {
    expect_true( "package:ggmap" %in% search(), 
                 info = "The ggmap package should be loaded using library().")
  })
  test_that("Test that viridis is loaded", {
    expect_true( "package:viridis" %in% search(), 
                 info = "The viridis package should be loaded using library().")
  })
  
  test_that("Check that geom_bin2d was used", {
    p <- last_plot()
    stat_classes <- as.character(sapply(p$layers, function(layer) {
      class(layer$stat)
    }))
    
    expect_true("StatBin2d" %in% stat_classes, 
                info = "You need to use geom_bin2d correctly to draw the map.")
  })
})


# Loading in the tree package
library(tree)

# Fitting a tree to lat and long
fitted_tree <- tree(total ~ lat + long, data = taxi)

# Draw a diagram of the tree structure
plot(fitted_tree) 
text(fitted_tree)

run_tests({
  test_that("Test that tree is loaded", {
    expect_true( "package:tree" %in% search(), 
                 info = "The tree package should be loaded using library().")
  })
  test_that("The tree has been fitted correctly", {
    correctly_fitted_tree <- tree(total ~ lat + long, data = taxi)
    expect_equivalent(fitted_tree, correctly_fitted_tree, 
                      info = "It seem you didn't fit the tree correctly. Check the hint, it might help!")
  })
})


# Loading in the lubridate package
library(lubridate)

# Generate the three new time variables
taxi <- taxi %>% 
  mutate(hour = hour(pickup_datetime),
         wday = wday(pickup_datetime, label = T),
         month = month(pickup_datetime, label = T))

run_tests({
  test_that("Test that lubridate is loaded", {
    expect_true( "package:lubridate" %in% search(), 
                 info = "The lubridate package should be loaded using library().")
  })
  test_that("hour is correct", {
    expect_equivalent(taxi$hour[1], 10L, 
                      info = "The `hour` column doesn't seem to be correct. Check the hint for more help.")
  })
  test_that("wday is correct", {
    expect_true(taxi$wday[1] == "Sun", 
                info = "The `wday` column doesn't seem to be correct. Check the hint for more help.")
  })
  test_that("month is correct", {
    expect_true(taxi$month[1] == "Jan", 
                info = "The `month` column doesn't seem to be correct. Check the hint for more help.")
  })
})

# Fitting a tree with total as the outcome and 
# lat, long, hour, wday, and month as predictors
fitted_tree <- tree(total ~ lat + long + hour + wday + month, data = taxi)

# draw a diagram of the tree structure
plot(fitted_tree)
text(fitted_tree)
# Summarizing the performance of the tree
summary(fitted_tree)

run_tests({
  test_that("The tree has been fitted correctly", {
    correctly_fitted_tree <- tree(total ~ lat + long + hour + wday + month, data = taxi)
    expect_equivalent(fitted_tree, correctly_fitted_tree, 
                      info = "It seem you didn't fit the tree correctly. Check the hint, it might help!")
  })
})

# Loading in the randomForest package
library(randomForest)

# Fitting a random forest
fitted_forest <- randomForest(formula = total ~ lat + long + hour + wday + month,
                              data = taxi,
                              sampsize = 10000,
                              ntree = 80)

# Printing the fitted_forest object
fitted_forest

run_tests({
  test_that("Test that randomForest is loaded", {
    expect_true( "package:randomForest" %in% search(), 
                 info = "The randomForest package should be loaded using library().")
  })
  test_that("ntree is correct.", {
    expect_true(fitted_forest$ntree == 80, 
                info = "The ntree argument to randomForest should be ntree = 80 .")
  })
  test_that("Check randomForest call was ok", {
    call_string <- paste(deparse(fitted_forest$call), collapse = " ")
    keywords <- c("total", "lat", "long", "hour", "wday", "month",
                  "ntree", "sampsize", "100")
    expect_true(all(str_detect(call_string, keywords)), 
                info = "You have not called randomForest correctly. Did you include all the predictors and the right output variable?.")
  })
})

# Extracting the prediction from fitted_forest
taxi$pred_total <- fitted_forest$predicted

# Plotting the predicted mean trip prices from according to the random forest
ggmap(manhattan, darken = 0.5) +
  scale_fill_viridis(option = 'plasma') +
  stat_summary_2d(data = taxi, aes(x = long,
                                   y = lat, 
                                   z = pred_total), 
                  alpha = 0.6,
                  bins = 60,
                  fun = mean) +
  labs(x = "Longitude",
       y = "Latitude",
       z = "Total Prediction")

run_tests({
  test_that("taxi$pred_total == fitted_forest$predicted", {
    expect_true(all(taxi$pred_total == fitted_forest$predicted), 
                info = "You should assign fitted_forest$predicted to taxi$pred_total .")
  })
  test_that("Check that stat_summary_2d was used", {
    p <- last_plot()
    stat_classes <- as.character(sapply(p$layers, function(layer) {
      class(layer$stat)
    }))
    
    expect_true("StatSummary2d" %in% stat_classes, 
                info = "You need to use geom_bin2d correctly to draw the map.")
  })
  test_that("Check that pred_total was used", {
    p <- last_plot()
    p_variables <- unlist(sapply(p$layers, function(layer) {
      as.character(layer$mapping)
    }))
    expect_true(any(str_detect(p_variables, "pred_total")), 
                info = "You need to connect pred_total to z in the aes() call correctly.")
  })
})


# Function that returns the mean *if* there are 15 or more datapoints
mean_if_enough_data <- function(x) { 
  ifelse( length(x) >= 15, mean(x), NA) 
}

# Plotting the mean trip prices from the data
ggmap(manhattan, darken = 0.5) +
  scale_fill_viridis(option = 'plasma') +
  stat_summary_2d(data = taxi, aes(x = long,
                                   y = lat, 
                                   z = total), 
                  alpha = 0.6,
                  bins = 60,
                  fun = mean_if_enough_data) +
  labs(x = "Longitude",
       y = "Latitude",
       z = "Total Prediction")

run_tests({
  test_that("Check that total was used but not pred_total", {
    p <- last_plot()
    p_variables <- unlist(sapply(p$layers, function(layer) {
      as.character(layer$mapping)
    }))
    expect_true(any(str_detect(p_variables, "total")) & 
                  !any(str_detect(p_variables, "pred_total")), 
                info = "You need to connect total to z in the aes() call correctly. Make sure you are not still using pred_total.")
  })
})

# Where are people spending the most on their taxi trips?
spends_most_on_trips <- "downtown" # "uptown" or "downtown"

run_tests({
  test_that("...", {
    expect_true(str_detect(tolower(spends_most_on_trips), "downtown"), 
                info = "Well, looking at the plot it looks like people pay more downtown.")
  })
})


#### VISUALIZING INEQUALITIES IN LIFE EXPECTANCY ###############################


# This sets plot images to a nice size
options(repr.plot.width = 6, repr.plot.height = 6)

# Loading packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Loading data
life_expectancy <- read.csv("datasets/UNdata.csv")

# Taking a look at the first few rows
life_expectancy

# These packages need to be loaded in the first `@tests` cell. 
library(testthat) 
library(IRkernel.testthat)

# Then follows one or more tests of the students code. 
# The @solution should pass the tests.
# The purpose of the tests is to try to catch common errors and to 
# give the student a hint on how to resolve these errors.

run_tests({
  test_that("Test that life_expectancy exists", {
    expect_true(exists("life_expectancy"), 
                info = "It seems that the data frame life_expectancy does not exist.")
  })
  
  test_that("Test that life_expectancy is loaded correctly", {
    expect_true(nrow(life_expectancy)==1571, 
                info = "The data frame life_expectancy is not correctly loaded.")
  })
  
  test_that("Test that life_expectancy is loaded correctly", {
    expect_true(ncol(life_expectancy)==7, 
                info = "The data frame life_expectancy is not correctly loaded.")
  })
})

# Subsetting and reshaping the life expectancy data
subdata <- life_expectancy  %>% 
  filter(Year == "2000-2005") %>%
  select(Country.or.Area, Subgroup, Value) %>%
  spread(Subgroup, Value)

# Taking a look at the first few rows
subdata

# one or more tests of the students code. 
# The @solution should pass the tests.
# The purpose of the tests is to try to catch common errors and to 
# give the student a hint on how to resolve these errors.
run_tests({
  test_that("Test that subdata exists", {
    expect_true(exists("subdata"), 
                info = "It seems that dataset subdata does not exist.")
  })
  
  test_that("Test that subdata is created correctly", {
    expect_true(nrow(subdata)==195, 
                info = "It seems that subdata is not correctly created.")
  })
  
  test_that("Test that subdata is created correctly", {
    expect_true(ncol(subdata)==3, 
                info = "It seems that subdata is not correctly created.")
  })
  
  test_that("Test that subdata is contains correct columns", {
    expect_true(sum(is.element(c("Country.or.Area", "Female", "Male"), names(subdata)))==3, 
                info = "It seems that subdata does not contain the correct columns.")
  })
})

# Plotting male and female life expectancy
ggplot(subdata, aes(Male, Female)) +
  geom_point()


run_tests({
  test_that("Check that a geom_point plot was plotted.", {
    expect_true( "GeomPoint" %in% class( last_plot()$layers[[1]]$geom ) , 
                 info = "Add geom_point() to produce a scatter plot.")
  })
  
  test_that("Check variables are correctly mapped.", {
    expect_true(deparse(last_plot()$mapping$x)=="~Male" & deparse(last_plot()$mapping$y) == "~Female",
                info = "Check that the variables are mapped to the correct axes.")
  })
  
})

# Adding an abline and changing the scale of axes of the previous plots
ggplot(subdata, aes(Male, Female)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(breaks = c(35, 45, 55, 65, 75, 85),
                     limits = c(35, 85)) +
  scale_y_continuous(breaks = c(35, 45, 55, 65, 75, 85),
                     limits = c(35, 85))





for (i in 1:length(ggplot_build(last_plot())$data)) 
{
  if ("slope" %in% colnames(ggplot_build(last_plot())$data[[i]])) i1=i
}


run_tests({
  test_that("Intercept of diagonal line is equal to 0.", {
    expect_equal(ggplot_build(last_plot())$data[[i1]]$intercept, 0, 
                 info = "Did you add the diagonal line correctly?")
  })
  test_that("Slope of diagonal line is equal to 1.", {
    expect_equal(ggplot_build(last_plot())$data[[i1]]$slope, 1, 
                 info = "Did you add the diagonal line correctly?")
  })
  test_that("Limits of x-axis.", {
    expect_equal(length(setdiff(c(39, 79), ggplot_build(last_plot())$layout$panel_scales_x[[1]]$range$range)), 0, 
                 info = "The limits of x-axis is not equal to [35, 85].")
  })
  test_that("Limits of y-axis.", {
    expect_equal(length(setdiff(c(39, 85), ggplot_build(last_plot())$layout$panel_scales_y[[1]]$range$range)), 0, 
                 info = "The limits of y-axis is not equal to [35, 85].")
  })
  
})

# Adding labels to previous plot
ggplot(subdata, aes(x=Male, y=Female))+
  geom_point(colour="white", fill="chartreuse3", shape=21, alpha=.55, size=5)+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  scale_x_continuous(limits=c(35,85))+
  scale_y_continuous(limits=c(35,85))+
  labs(title = "Life Expectancy at Birth by Country",
       subtitle = "Years. Period: 2000-2005. Average.",
       caption = "Source: United Nations Statistics Division",
       x = "Males",
       y = "Females")

run_tests({
  test_that("Title is correct.", {
    expect_equal(toupper(gsub("[[:space:]]", "", last_plot()$labels$title)), "LIFEEXPECTANCYATBIRTHBYCOUNTRY", 
                 info = "Did you add the title correctly?")
  })
  
  test_that("x-axis label is correct.", {
    expect_equal(toupper(gsub("[[:space:]]", "", last_plot()$labels$x)), "MALES", 
                 info = "Did you set the x-axis label correctly?")
  })
  
  
  test_that("y-axis label is correct.", {
    expect_equal(toupper(gsub("[[:space:]]", "", last_plot()$labels$y)), "FEMALES", 
                 info = "Did you set the y-axis label correctly?")
  })
  
  test_that("caption is correct.", {
    expect_equal(toupper(gsub("[[:space:]]", "", last_plot()$labels$caption)), "SOURCE:UNITEDNATIONSSTATISTICSDIVISION", 
                 info = "Did you set the caption correctly?")
  })
  
  
  
})

# Subseting data to obtain countries of interest
top_male <- subdata %>% arrange(Male-Female) %>% head(3)
top_female <- subdata %>% arrange(Female-Male) %>% head(3)

# Adding text to the previous plot to label countries of interest
ggplot(subdata, aes(x = Male, 
                    y = Female, 
                    label = Country.or.Area))+
  geom_point(colour="white", 
             fill="chartreuse3", 
             shape=21, alpha=.55, 
             size=5)+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  scale_x_continuous(limits=c(35,85))+
  scale_y_continuous(limits=c(35,85))+
  labs(title="Life Expectancy at Birth by Country",
       subtitle="Years. Period: 2000-2005. Average.",
       caption="Source: United Nations Statistics Division",
       x="Males",
       y="Females")+
  theme_bw() +
  geom_text(data = top_male,
            size = 3) +
  geom_text(data = top_female,
            size = 3)

texts=c()
for (i in 1:length(last_plot()$layers)) texts=c(last_plot()$layers[[i]]$data$Country.or.Area %>% as.character, texts)


run_tests({
  
  test_that("Test that countries defined by top_female and top_male are correctly labeled.", {
    expect_true(length(setdiff(texts, c("Russian Federation", "Belarus", "Estonia", "Niger", "Afghanistan", "Maldives")))==0, 
                info = "It seems that countries defined by top_female and top_male are not labeled correctly.")
  })
  
  test_that("Theme is theme_bw().", {
    expect_equal(last_plot()$theme$panel.background$fill, "white", 
                 info = "It seems that your plot does not have theme_bw().")
  })
  
  
})

# Subsetting, mutating and reshaping the life expectancy data
subdata2 <- life_expectancy %>% 
  filter(Year %in% c("1985-1990", "2000-2005")) %>% 
  mutate(Sub_Year = paste(Subgroup, Year, sep="_")) %>% 
  mutate(Sub_Year = gsub("-", "_", Sub_Year)) %>% 
  select(-Subgroup, -Year) %>% 
  spread(Sub_Year, Value) %>%
  mutate(diff_Female = Female_2000_2005 - Female_1985_1990,
         diff_Male = Male_2000_2005 - Male_1985_1990)

# Taking a look at the first few rows
subdata2


run_tests({
  test_that("Test that subdata2 is created correctly.", {
    expect_true(nrow(subdata2)==195, 
                info = "It seems that dataset subdata2 is not correctly created.")
  })
  
  test_that("Test that subdata2 is created correctly.", {
    expect_true(ncol(subdata2)==10, 
                info = "It seems that dataset subdata2 is not correctly created.")
  })
  
  test_that("Test that subdata2 is created correctly.", {
    expect_true(length(setdiff(c('diff_Female', 'diff_Male'), names(subdata2)))==0, 
                info = "It seems that subdata2 does not contain columns diff_Female or diff_Male.")
  })
  
  test_that("Test that subdata2 is created correctly.", {
    expect_true(sum(subdata2$diff_Female)==492, 
                info = "It seems that the diff_Female column is not correctly created.")
  })  
  
  test_that("Test that subdata2 is created correctly.", {
    expect_true(sum(subdata2$diff_Male)==503, 
                info = "It seems that the diff_Male column is not correctly created.")
  })
})

# Doing a nice first version of the plot with abline, scaling axis and adding labels
ggplot(subdata2, aes(x = diff_Male, 
                     y = diff_Female, 
                     label = Country.or.Area))+
  geom_point(colour="white", fill="chartreuse3", shape=21, alpha=.55, size=5)+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  scale_x_continuous(limits=c(-25,25))+
  scale_y_continuous(limits=c(-25,25))+
  labs(title="Life Expectancy at Birth by Country in Years",
       subtitle="Difference between 1985-1990 and 2000-2005. Average.",
       caption="Source: United Nations Statistics Division",
       x="Males",
       y="Females")+
  theme_bw()

run_tests({
  
  #    test_that("Check that a geom_point plot was plotted.", {
  #     expect_true( "GeomPoint" %in% class( last_plot()$layers[[1]]$geom ) , 
  #                  info = "Add geom_point() to produce a scatter plot.")
  #   })
  
  test_that("Check variables are correctly mapped.", {
    expect_true( deparse(last_plot()$mapping$x)=="~diff_Male" & deparse(last_plot()$mapping$y)=="~diff_Female", 
                 info = "Check that the variables are mapped to the correct axes.")
  })
  
  
  #     test_that("Intercept of diagonal line is equal to 0.", {
  #     expect_equal(ggplot_build(last_plot())$data[[2]]$intercept, 0, 
  #         info = "Did you add the diagonal line correctly?")
  #     })
  #     test_that("Slope of diagonal line is equal to 1.", {
  #     expect_equal(ggplot_build(last_plot())$data[[2]]$slope, 1, 
  #         info = "Did you add the diagonal line correctly?")
  #     })
  test_that("Limits of x-axis", {
    expect_equal(length(setdiff(c(-20, 15), ggplot_build(last_plot())$layout$panel_scales_x[[1]]$range$range)), 0, 
                 info = "Limits of x-axis is not equal to [-25, 25].")
  })
  test_that("Limits of y-axis", {
    expect_equal(length(setdiff(c(-24, 15), ggplot_build(last_plot())$layout$panel_scales_y[[1]]$range$range)), 0, 
                 info = "Limits of y-axis is not equal to [-25, 25]")
  })
  
  #     test_that("Intercept of diagonal line is equal to 0.", {
  #     expect_equal(toupper(gsub("[[:space:]]", "", last_plot()$labels$title)), "LIFEEXPECTANCYATBIRTHBYCOUNTRYINYEARS", 
  #         info = "Did you add the title correctly?")
  #     })
  #     test_that("Slope of diagonal line is equal to 1.", {
  #     expect_equal(last_plot()$theme$panel.background$fill, "white", 
  #         info = "It seems that your plot does not have theme_bw().")
  #     })
})

# Adding an hline and vline to previous plots
ggplot(subdata2, aes(x=diff_Male, y=diff_Female, label=Country.or.Area))+
  geom_point(colour="white", fill="chartreuse3", shape=21, alpha=.55, size=5)+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  scale_x_continuous(limits=c(-25,25))+
  scale_y_continuous(limits=c(-25,25))+
  geom_vline(xintercept = 0, linetype= "dashed") +
  geom_hline(yintercept = 0, linetype= "dashed") +
  labs(title="Life Expectancy at Birth by Country",
       subtitle="Years. Difference between 1985-1990 and 2000-2005. Average.",
       caption="Source: United Nations Statistics Division",
       x="Males",
       y="Females")+
  theme_bw()

for (i in 1:length(ggplot_build(last_plot())$data)) 
{
  if ("slope"      %in% colnames(ggplot_build(last_plot())$data[[i]])) i1=i
  if ("yintercept" %in% colnames(ggplot_build(last_plot())$data[[i]])) i2=i
  if ("xintercept" %in% colnames(ggplot_build(last_plot())$data[[i]])) i3=i
}


run_tests({
  #     test_that("Intercept of diagonal line is equal to 0.", {
  #     expect_equal(ggplot_build(last_plot())$data[[i1]]$intercept, 0, 
  #         info = "Did you add the diagonal line correctly?")
  #     })
  #     test_that("Slope of diagonal line is equal to 1.", {
  #     expect_equal(ggplot_build(last_plot())$data[[i1]]$slope, 1, 
  #         info = "Did you add the diagonal line correctly?")
  #     })
  
  test_that("Horizontal line is well defined.", {
    expect_equal(ggplot_build(last_plot())$data[[i2]]$yintercept, 0, 
                 info = "Did you add the horizontal line correctly?")
  })
  test_that("Vertical line is well defined.", {
    expect_equal(ggplot_build(last_plot())$data[[i3]]$xintercept, 0, 
                 info = "Did you add the vertical line correctly?")
  })
})

# Subseting data to obtain countries of interest
top <- subdata2 %>% arrange(diff_Male+diff_Female) %>% head(3)
bottom <- subdata2 %>% arrange(-(diff_Male+diff_Female)) %>% head(3)

# Adding text to the previous plot to label countries of interest
ggplot(subdata2, aes(x=diff_Male, y=diff_Female, label=Country.or.Area), guide=FALSE)+
  geom_point(colour="white", fill="chartreuse3", shape=21, alpha=.55, size=5)+
  geom_abline(intercept = 0, slope = 1, linetype=2)+
  scale_x_continuous(limits=c(-25,25))+
  scale_y_continuous(limits=c(-25,25))+
  geom_hline(yintercept=0, linetype=2)+
  geom_vline(xintercept=0, linetype=2)+
  labs(title="Life Expectancy at Birth by Country",
       subtitle="Years. Difference between 1985-1990 and 2000-2005. Average.",
       caption="Source: United Nations Statistics Division",
       x="Males",
       y="Females")+
  geom_text(data = top,
            size = 3) +
  geom_text(data = bottom,
            size = 3) +
  theme_bw()

texts=c()
for (i in 1:length(last_plot()$layers)) texts=c(last_plot()$layers[[i]]$data$Country.or.Area %>% as.character, texts)

run_tests({
  test_that("Test that dataset bottom exists.", {
    expect_true(exists("bottom"), 
                info = "It seems that bottom does not exist.")
  })
  
  test_that("Test that dataset bottom is correctly created.", {
    expect_true(nrow(bottom)==3, 
                info = "It seems that bottom is not correctly created.")
  })
  
  test_that("Test that countries defined by top and bottom are correctly labeled.", {
    expect_true(length(setdiff(texts, c("Timor Leste", "Bhutan", "Egypt", "Zimbabwe", "Botswana", "Swaziland")))==0, 
                info = "It seems that countries defined by top and bottom are not labeled correctly.")
  })
  
  
})

