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



