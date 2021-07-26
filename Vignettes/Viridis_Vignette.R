library(viridis)

## For base plots, you can use viridis() to generate a function
    x <- y <- seq(-8*pi,
                  8*pi,
                  len = 40)
    r <- sqrt(outer(x^2,
                    y^2,
                    "+"))
    filled.contour(cos(r^2) * exp(-r / (2*pi)),
                   axes = F,
                   color.palette = viridis,
                   asp = 1)

## for ggplot --> use scale_color_viridis() and scale_fill_viridis()
    library(ggplot2)
    ggplot(data.frame(x = rnorm(10000),
                      y = rnorm(10000)),
                      aes(x = x, y = y)) +
      geom_hex() +
      coord_fixed() +
      scale_fill_viridis() +
      theme_bw()
             
    
## Intro
    ## color maps are designed to be:
      ## Colorful
      ## perceptually uniform
      ## robust to colour blindness

## viridisLite provides the base functions for generating the color maps in base R. 
    ## The package is meant to be as lightweight and dependency-free as possible 
    ## for maximum compatibility with all the R ecosystem. 
    ## viridis provides additional functionalities, in particular bindings for ggplot2.
    
    
## The Color Scales
    ## The package contains eight color scales: “viridis”, the primary choice, 
    ## and five alternatives with similar properties - 
    ## “magma”, “plasma”, “inferno”, “civids”, “mako”, and “rocket” -, 
    ## and a rainbow color map - “turbo”.
    
    
## USAGE
    ## Here the inferno() scale is used for a raster of U.S. max temperature:
    install.packages("rgdal")
    library(rasterVis)
    library(httr)
    
    par(mfrow=c(1,1), mar=rep(0.5, 4))
    
    temp_raster <- "http://ftp.cpc.ncep.noaa.gov/GIS/GRADS_GIS/GeoTIFF/TEMP/us_tmax/us.tmax_nohads_ll_20150219_float.tif"
    try(GET(temp_raster, write_disk("us.tmax_nohads_ll_20150219_float.tif")), silent=TRUE)
    
    us <- raster("us.tmax_nohads_ll_20150219_float.tif")
    us <- projectRaster(us, crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
    
    image(us, 
          col=inferno(256), 
          asp=1, 
          axes=FALSE, 
          xaxs="i", 
          xaxt='n', 
          yaxt='n', ann=FALSE)
    
    ## Here the “magma” scale is used for a cloropleth map of U.S. unemployment:
      unemp <- read.csv("http://datasets.flowingdata.com/unemployment09.csv", header = FALSE, stringsAsFactors = FALSE)
       names(unemp) <- c("id", 
                         "state_fips", 
                         "county_fips", 
                         "name", 
                         "year",
                         "?", 
                         "?", 
                         "?", 
                         "rate")
       p$county <- tolower(gsub(" County, [A-Z]{2}", "", unemp$name))
       unemp$county <- gsub("^(.*) parish, ..$","\\1", unemp$county)
       unemp$state <- gsub("^.*([A-Z]{2}).*$", "\\1", unemp$name)
      data(unemp)
      county_df <- map_data("county", projection = "albers", parameters = c(39, 45))
      names(county_df) <- c("long", "lat", "group", "order", "state_name", "county")
      county_df$state <- state.abb[match(county_df$state_name, tolower(state.name))]
      county_df$state_name <- NULL
      state_df <- map_data("state", projection = "albers", parameters = c(39, 45))
      choropleth <- merge(county_df, unemp, by = c("state", "county"))
      choropleth <- choropleth[order(choropleth$order), ]
      ggplot(choropleth, aes(long, 
                             lat, 
                             group = group)) +
        geom_polygon(aes(fill = rate), colour = alpha("white", 1 / 2), size = 0.2) +
        geom_polygon(data = state_df, colour = "white", fill = NA) +
        
        coord_fixed() +
        theme_minimal() +
        ggtitle("US unemployment rate by county") +
        theme(axis.line = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(), axis.title = element_blank()) +
        scale_fill_viridis(option="magma")
      
      
    ## the ggplot functions also can be used for discrete values with the discrete = T argument
    p <- ggplot(mtcars,
                aes(wt, mpg))
    p + geom_point(size = 4, aes(colour = factor(cyl))) +
      scale_color_viridis(discrete = T) +
      theme_bw()

        

    
    