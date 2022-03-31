# R has three main plotting systems:
  #   1. main plotting embedded in R
  #   2. lattice plotting system
  #   3. ggplot <- good for plotting data!


  library(ggplot2)

 
# IRIS data set <- load available sets with data() ----

  data()
    data("iris")
    head(iris)
    summary(iris)
  
# PLOT IRIS ----
  
  plot(iris)
    # interesting relationship b/w petal length & sepal width
  
  plot(iris$Petal.Length, iris$Sepal.Width)

# GGPLOT explanation ----
  
      # gg stands for "grammar of graphics" = any graphical representation of data/plot/graph
      # can be produced from a series of layers
          # layer   grid lines --> plotting area
          # layer   axis ticks
          # layer   values for the axis
          # layer   labels for the axis
          # layer   plot title
          # layer   geoms --> geometric objects (can be dots, triangles, text, etc)
          # layer   additional plots (boxplots, regression lines, confidence intervals, text, ..)
          # ...
      # ggplot fct. needs:
          # 1. data frame --> here iris
          # 2. aes fct. --> aesthetic mappings = relationship b/w variable and some aspect of the plot
            # e.g. values of x-axis, so we can map petal.length values to x-axis
            # e.g. sepal.width to the y-axis
            # you can choose color, size, transparency, ...
 
  
# GGPLOT IRIS ----
  ggplot(iris, aes(x=Petal.Length, y = Sepal.Width)) +  # <-- base plot, now add on with +
    geom_point() # now we get all single dots of data 
    ?geom_point() # look at the Aesthetics options --> this can be data input
                  # are there relationships among those 3 species?
   
    
    # add color 
      # +++ good for continuous variables
  
  ggplot(iris, aes(x=Petal.Length, 
                   y = Sepal.Width, 
                   col = Species)) +
                  geom_point()    
  
    # add representation of petal.width --> add size
  
  ggplot(iris, aes(x=Petal.Length, 
                   y = Sepal.Width, 
                   col = Species, 
                   size = Petal.Width)) +
                  geom_point() 
  
    # add shape
      # +++ categorical values
      # --- continuous values don't work
  
  ggplot(iris, aes(x=Petal.Length, 
                   y = Sepal.Width, 
                   col = Species, 
                   size = Petal.Width,
                   shape = Species)) +
                  geom_point() 

    # don't overdo it.. --> 1 variable per aesthetic
  
    # transparency = alpha 
        # --- continuous variables
        # --- changes hard to see
        # +++ use color for continuous variables
  
  ggplot(iris, aes(x=Petal.Length, 
                   y = Sepal.Width, 
                   col = Species, 
                   size = Petal.Width,
                   shape = Species,
                   alpha = Sepal.Length)) +
                  geom_point() 
      
# BAR AND BOXPLOT ----

  # x-axis = species 
  # y-axis = sepal.length
  # bar at group means
  
  ggplot(iris, aes(Species)) +
    geom_bar()  # gives error, when y is already added for a bar plot, leave out for now
                # not very useful --> simply a histogram
                # scrap that, instead we will look at Sepal.Length
  
  ggplot(iris, aes(Sepal.Length)) +
    geom_bar()  # we will now change the stats of geom_bar
                # default stat = "count"
                # stat = "summary"
                # specify WHICH summary stat you want to plot --> fun = "mean" (mean function of y)
  
  ggplot(iris, aes(Species, Sepal.Length)) +
    geom_bar(stat = "summary", fun = "mean")
  
  # col refers to the border of the geom
  ggplot(iris, aes(Species, Sepal.Length, col = Species)) +
    geom_bar(stat = "summary", fun = "mean")
  
  # col of inside = fill
  ggplot(iris, aes(Species, Sepal.Length, fill = Species)) +
    geom_bar(stat = "summary", fun = "mean")
  
  # all bars blue
  ggplot(iris, aes(Species, Sepal.Length)) +
    geom_bar(stat = "summary", fun = "mean", fill = "blue")
  
  # specific color control with color codes
  ggplot(iris, aes(Species, Sepal.Length)) +
    geom_bar(stat = "summary", fun = "mean", fill = "#ff0076")
  
  # add border to your bars
  ggplot(iris, aes(Species, Sepal.Length)) +
    geom_bar(stat = "summary", fun = "mean", fill = "#ff0076", col = "black")

  # add dots for each individual data value on top of bars
  ggplot(iris, aes(Species, Sepal.Length)) +
    geom_bar(stat = "summary", fun = "mean", fill = "#ff0076", col = "black") +
    geom_point()
  
  # quite small, hard to differentiate --> position points
  ggplot(iris, aes(Species, Sepal.Length)) +
    geom_bar(stat = "summary", fun = "mean", fill = "#ff0076", col = "black") +
    geom_point(position = position_jitter())
  
  # add width to points, shapes have numbers --> 
  ggplot(iris, aes(Species, Sepal.Length)) +
    geom_bar(stat = "summary", fun = "mean", fill = "#ff0076", col = "black") +
    geom_point(position = position_jitter(0.2), size = 3, shape = 5)

  # put everything into an object
  myPlot <- ggplot(iris, aes(Species, Sepal.Length)) +
              geom_bar(stat = "summary", fun = "mean", fill = "#ff0076", col = "black") +
              geom_point(position = position_jitter(0.2), size = 3, shape = 5)
          myPlot
  
  # change stuff further --> theme options
    ?theme
  
  myPlot + theme( panel.grid = element_blank(), # removes all grid (we could also remove just x/y axis)
                  panel.background = element_rect(fill = "white"), # removes gray box
                  # axis.line.y = element_line(color = "black", size = 0.2), --> can add y-axis line
                  # axis.line.x = element_line(color = "black", size = 0.2), --> can add x-axis line
                  panel.border = element_rect(color = "black", fill = NA, size = 0.2) # marks x/y axis and gives borders all around
                )
                  
                  
    myPlot  # Starting point and add to it         
    myPlot + theme_bw()                
    myPlot + theme_classic()                  
    myPlot + theme_dark()                  
    myPlot + theme_get()                  
    myPlot + theme_gray()  
    myPlot + theme_light()
    myPlot + theme_map()  
    myPlot + theme_minimal()  
    myPlot + theme_void()
    myPlot + theme_linedraw() + theme(panel.background = element_rect(fill = "blue")) # you can add anything you want
    myPlot + theme_hefin # object --> no () needed
  
  
  # now with boxplots
    ggplot(iris, aes(Species, Sepal.Length)) +
      geom_boxplot(fill = "#ff0080", col = "black", notch = TRUE) + #notch default = FALSE, Median clearer
      geom_point()
  # order matters --> eg.:
    ggplot(iris, aes(Species, Sepal.Length)) +
      geom_point() +
      geom_boxplot(fill = "#ff0080", col = "black", notch = TRUE) #notch default = FALSE, Median clearer

    
# FINISHING TOUCHES ----
    myPlot + 
      theme_hefin +
        labs(x = "", y = "Sepal length (mm)") +  # changes axis labels --> more info, looks better
        ggtitle("Sepal Length by iris species") +
        theme(plot.title = element_text(hjust = 0.5)) # title in the middle, default = left
    
    
# SAVING OUR PLOT ----
    # use ggsave function
    # saves last produced plot with specifications
    
    
    setwd("~/Documents") # set your working directory
                # save to desktop = "~/"
                # save to documents = "~/Documents"
                
    ggsave("plotting.pdf", width = 8, height = 5) 
      # specify type with 
        # .pdf 
        # .png 
        # .jpg 
        # .svg (install package "svg.lite" for this)
        # higher png resolution for windows with --> ggsave("x", type = "cairo-png")
      # specify unit --> default = inches
    
    
# FACTORIAL DATA ----
    data(ToothGrowth) # --> this is a factorial design experiment
    head(ToothGrowth)    
    summary(ToothGrowth)    
      # R sees "dose" as a continuous variable, but it's actually a categorical/ordinal variable (can be only 0.5, 1, or 2)
    
    ggplot(ToothGrowth, aes(supp, len, fill = dose)) +
      geom_bar(stat = "summary", fun = "median", col = "black") # summary, not count
        
      # tell R that dose = factor
      # tell R that position =  "stack" (stacks each mean on top of each other)
      #                         "fill" (takes 100% of each group and plots bar width as proportion of each dose)
      #                         "dodge" (easier comparisons with separate bars)
    ggplot(ToothGrowth, aes(supp, len, fill = as.factor(dose))) +
      geom_bar(stat = "summary", fun = "median", col = "black", position = "dodge") + # default "dodge" = bars are immediately each other
      geom_point(position = position_dodge(0.9)) # must specify the dot positions

    ggplot(ToothGrowth, aes(as.factor(dose), len, group = supp, col = supp)) +
      geom_line(stat = "summary", fun = "mean") + # specify that we want summary statistics
      geom_smooth(method = "lm") # draw regression lines on your plots
    
    
    