library(utils)
browseVignettes()

library(visdat)
library(dplyr)


## VIS_DAT
    ## helps explore the data class structure and missingness
    ## provides unit testing for your data. In a similar way, 
    ## visdat provides visual tests, the idea being that first you visualise 
    ## your data (visdat), then you run tests from testdat, or a package like assertr, 
    ## to fix these errors.
    vis_dat(iris, sort_type = T)
    vis_dat(typical_data, sort_type = T)


## VIS_MISS
    ## provides a custom plot for missing data
    vis_miss(iris)
    vis_miss(typical_data, sort_miss = T, cluster = T)

## VIS_COMPARE
    ## Sometimes you want to see what has changed in your data. 
    ## vis_compare() displays the differences in two dataframes of the same size. 
    chickwts_diff <- chickwts
    chickwts_diff[sample(1:nrow(chickwts), 30),
                  sample(1:ncol(chickwts), 2)] <- NA
    vis_compare(chickwts_diff, chickwts)
  
    
## VIS_EXPECT
    ## vis_expect visualises certain conditions or values in your data. 
    ## For example, If you are not sure whether to expect values greater than 25 
    ## in your data (airquality), you could write: vis_expect(airquality, ~.x >= 25), 
    ## and you can see if there are times where the values in your data are greater 
    ## than or equal to 25.
    vis_expect(airquality, ~.x >= 25)
    
    ## This shows the proportion of times that there are values greater than 25, 
    ## as well as the missings.
    ## You could also, for example, explore a set of bad strings, or possible NA 
    ## values and visualise where they are using vis_expect(data, ~.x %in% bad_strings) 
    ## where bad_strings is a character vector containing bad strings like N A, N/A etc.
    bad_data <- data.frame(x = c(rnorm(100), rep("N/A", 10)),
                           y = c(rep("N A ", 30), rnorm(80)))
    vis_expect(bad_data, ~.x %in% c("N/A", "N A "))

## VIS_COR
    ## To make it easy to plot correlations of your data, use vis_cor:
    vis_cor(airquality)
    ## Under the hood, vis_cor is powered by the cor function in base R, 
    ## and takes a character string indicating which correlation coefficient 
    ## (or covariance) is to be computed. One of “pearson” (default), “kendall”, 
    ## or “spearman”.
    vis_cor(airquality, cor_method = "spearman")
    ## You can also specify what to do for the missing data using the na_action function, 
    ## which again borrows from the cor methods. This can be “everything”, 
    ## “all.obs”, “complete.obs”, “na.or.complete”, or “pairwise.complete.obs" (default)
    vis_cor(airquality, na_action = "complete.obs")

        
## VIS_GUESS
    ## vis_guess() takes a guess at what each cell is. It’s best illustrated 
    ## using some messy data, which we’ll make here.
    messy_vector <- c(TRUE,
                      T,
                      "TRUE",
                      "T",
                      "01/01/01",
                      "01/01/2001",
                      NA,
                      NaN,
                      "NA",
                      "Na",
                      "na",
                      "10",
                      10,
                      "10.1",
                      10.1,
                      "abc",
                      "$%TG")
    set.seed(1114)
    messy_df <- data.frame(var1 = messy_vector,
                           var2 = sample(messy_vector),
                           var3 = sample(messy_vector))
    vis_guess(messy_df)
    vis_dat(messy_df)
    ## So here we see that there are many different kinds of data in your dataframe. 
    ## As an analyst this might be a depressing finding. 
    ## We can see this comparison above.
    ## Here, you might just assume your data is weird because it’s all factors - 
    ## or worse, not notice that this is a problem.
    
    ## At the moment vis_guess is very slow. Please take this into consideration 
    ## when you are using it on data with more than 1000 rows. We’re looking into 
    ## ways of making it faster, potentially using methods from the parallel package, 
    ## or extending the c++ code from readr:::collectorGuess.
    
## COMBINE VIS_ AND GGPLOTLY
    ## You can make the plots in visdat by wrapping them in plotly::ggplotly:
    library(plotly)
    ggplotly(vis_dat(airquality))
    ggplotly(vis_miss(airquality))
    ggplotly(vis_guess(airquality))
