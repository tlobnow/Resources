library(ggplot2)
library(dplyr)

US_Covid <- read.csv("~/Documents/GitHub/covid-19-data/rolling-averages/us.csv")

glimpse(US_Covid)

US_Covid %>%
  ggplot(aes(date, cases_avg_per_100k)) +
  geom_point(alpha = 0.3)

US_Covid %>%
ggplot(aes(date, deaths_avg_per_100k)) + 
  geom_point(alpha = 0.3)

US_Cov_Counties <- read.csv("~/Documents/GitHub/covid-19-data/us-counties.csv")
glimpse(US_Cov_Counties)

VA_US_counties <-  US_Cov_Counties %>%
  filter(state %in% "Virginia")

# mutate Districts in VA_Subset ----
VA_mut <- VA_US_counties %>%
  mutate(Central_District = county %in% c("Nelson", 
                                             "Albemarle",
                                             "Charlottesville city",
                                             "Greene",
                                             "Madison",
                                             "Orange",
                                             "Louisa",
                                             "Fluvanna",
                                             "Caroline",
                                             "Hanover",
                                             "King William",
                                             "Goochland",
                                             "Powhatan",
                                             "Amelia",
                                             "Henrico",
                                             "Richmond",
                                             "New Kent",
                                             "Charles City",
                                             "Hopewell",
                                             "Prince George",
                                             "Chesterfield",
                                             "Colonial Heights",
                                             "Petersburg",
                                             "Dinwiddie") ,
          West_Central_District = county %in% c("Giles",
                                               "Pulaski",
                                               "Montgomery",
                                               "Radford city",
                                               "Floyd",
                                               "Franklin",
                                               "Roanoke",
                                               "Roanoke city",
                                               "Salem",
                                               "Craig",
                                               "Botetourt",
                                               "Bedford",
                                               "Lynchburg city",
                                               "Amherst",
                                               "Appomattox",
                                               "Campbell") ,
          Southside_District = county %in% c("Patrick",
                                            "Henry",
                                            "Martinsville city",
                                            "Pittsylvania",
                                            "Danville city",
                                            "Halifax",
                                            "Charlotte",
                                            "Prince Edward",
                                            "Buckingham",
                                            "Cumberland",
                                            "Nottoway",
                                            "Lunenburg",
                                            "Mecklenburg",
                                            "Brunswick",
                                            "Emporia city",
                                            "Greensville",
                                            "Sussex",
                                            "Southampton",
                                            "Franklin city") ,
        Hampton_Roads_District = county %in% c("Mathews",
                                                "Gloucester",
                                                "York",
                                                "James City",
                                                "Williamsburg city",
                                                "Poquoson",
                                                "Hampton",
                                                "Newport News",
                                                "Surry",
                                                "Isle of Wight",
                                                "Suffolk city",
                                                "Chesapeake city",
                                                "Portsmouth",
                                                "Norfolk city",
                                                "Virginia Beach") ,
        Eastern_District = county %in% c("Westmoreland",
                                          "Richmond",
                                          "Essex",
                                          "King and Queen",
                                          "Northumberland",
                                          "Lancaster",
                                          "Middlesex",
                                          "Accomack",
                                          "Northampton") ,
          Southwest_District = county %in% c("Lee",
                                            "Wise",
                                            "Norton city",
                                            "Scott",
                                            "Dickenson",
                                            "Buchanan",
                                            "Russell",
                                            "Washington",
                                            "Bristol city",
                                            "Smyth",
                                            "Tazewell",
                                            "Bland",
                                            "Wythe",
                                            "Grayson",
                                            "Galax city",
                                            "Carroll"),
        Northern_District = county %in% c("Clarke",
                                           "Warren",
                                           "Loudoun",
                                           "Rappahannock",
                                           "Prince William",
                                           "Fauqier",
                                           "Culpeper",
                                           "Stafford",
                                           "Fredericksburg city",
                                           "Spotsylvania",
                                           "King George",
                                           "Fairfax",
                                           "Fairfax city",
                                           "Falls Church city",
                                           "Arlington",
                                           "Alexandria city",
                                           "Manassas Park city",
                                           "Manassas city"),
              Valley_District = county %in% c("Frederick",
                                         "Winchester city",
                                         "Shenandoah",
                                         "Page",
                                         "Rockingham",
                                         "Harrisonburg city",
                                         "Augusta",
                                         "Highland",
                                         "Staunton city",
                                         "Waynesboro city",
                                         "Bath",
                                         "Rockbridge",
                                         "Buena Vista city",
                                         "Lexington city",
                                         "Alleghany",
                                         "Covington city"))
        VA_mut

        
        
# Visualize differences in Districts ----
VA_mut %>%
  filter(Central_District == T) %>%
  ggplot(aes(date, cases, col = county)) +
  geom_point(alpha = 0.3)
        
  facet_wrap(~ county)

VA_mut %>%
  filter(Eastern_District == T) %>%
  ggplot(aes(date, cases, col = county)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ county)

VA_mut %>%
  filter(Hampton_Roads_District == T) %>%
  ggplot(aes(date, cases, col = county)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ county)

VA_mut %>%
  filter(Northern_District == T) %>%
  ggplot(aes(date, cases, col = county)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ county)

VA_mut %>%
  filter(Southside_District == T) %>%
  ggplot(aes(date, cases, col = county)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ county)

# DISTRICTS ----

  # Central District ----
  Central_District <- VA_US_counties %>%
    filter(county %in% c("Nelson", 
                         "Albemarle",
                         "Charlottesville city",
                         "Greene",
                         "Madison",
                         "Orange",
                         "Louisa",
                         "Fluvanna",
                         "Caroline",
                         "Hanover",
                         "King William",
                         "Goochland",
                         "Powhatan",
                         "Amelia",
                         "Henrico",
                         "Richmond",
                         "New Kent",
                         "Charles City",
                         "Hopewell",
                         "Prince George",
                         "Chesterfield",
                         "Colonial Heights",
                         "Petersburg",
                         "Dinwiddie"
                         ))
  
      Central_District %>%
        ggplot(aes(date, cases, col = county)) +
          geom_point()
      Central_District %>%  
        ggplot(aes(date, deaths, col = county)) +
          geom_point()
  
  # West_Central District ----
      
  West_Central_District <- VA_US_counties %>%
    filter(county %in% c("Giles",
                         "Pulaski",
                         "Montgomery",
                         "Radford city",
                         "Floyd",
                         "Franklin",
                         "Roanoke",
                         "Roanoke city",
                         "Salem",
                         "Craig",
                         "Botetourt",
                         "Bedford",
                         "Lynchburg city",
                         "Amherst",
                         "Appomattox",
                         "Campbell"))
  
      West_Central %>%
        ggplot(aes(date, cases, col = county)) +
        geom_point()
      West_Central %>%  
        ggplot(aes(date, deaths, col = county)) +
        geom_point()
  
  # Southside District ----
      
  Southside_District <- VA_US_counties %>%
    filter(county %in% c("Patrick",
                         "Henry",
                         "Martinsville city",
                         "Pittsylvania",
                         "Danville city",
                         "Halifax",
                         "Charlotte",
                         "Prince Edward",
                         "Buckingham",
                         "Cumberland",
                         "Nottoway",
                         "Lunenburg",
                         "Mecklenburg",
                         "Brunswick",
                         "Emporia city",
                         "Greensville",
                         "Sussex",
                         "Southampton",
                         "Franklin city"))
  
      Southside_District %>%
        ggplot(aes(date, cases, col = county)) +
        geom_point()
      Southside_District %>%  
        ggplot(aes(date, deaths, col = county)) +
        geom_point()
  
  # Hampton Roads District ----
      
  Hampton_Roads_District <- VA_US_counties %>%
    filter(county %in% c("Mathews",
                         "Gloucester",
                         "York",
                         "James City",
                         "Williamsburg city",
                         "Poquoson",
                         "Hampton",
                         "Newport News",
                         "Surry",
                         "Isle of Wight",
                         "Suffolk city",
                         "Chesapeake city",
                         "Portsmouth",
                         "Norfolk city",
                         "Virginia Beach"
                         ))
  
        Hampton_Roads_District %>%
          ggplot(aes(date, cases, col = county)) +
          geom_point()
        Hampton_Roads_District %>%  
          ggplot(aes(date, deaths, col = county)) +
          geom_point()
    
  
  # Eastern District ----
        
  Eastern_District <- VA_US_counties %>%
    filter(county %in% c("Westmoreland",
                         "Richmond",
                         "Essex",
                         "King and Queen",
                         "Northumberland",
                         "Lancaster",
                         "Middlesex",
                         "Accomack",
                         "Northampton"))
  
      Eastern_District %>%
        ggplot(aes(date, cases, col = county)) +
        geom_point()
      Eastern_District %>%  
        ggplot(aes(date, deaths, col = county)) +
        geom_point()
  
  
  # Southwest District ----
      
  Southwest_District <- VA_US_counties %>%
      filter(county %in% c("Lee",
                           "Wise",
                           "Norton city",
                           "Scott",
                           "Dickenson",
                           "Buchanan",
                           "Russell",
                           "Washington",
                           "Bristol city",
                           "Smyth",
                           "Tazewell",
                           "Bland",
                           "Wythe",
                           "Grayson",
                           "Galax city",
                           "Carroll"))
    
    Southwest_District %>%
      ggplot(aes(date, cases, col = county)) +
      geom_point()
    Southwest_District %>%  
      ggplot(aes(date, deaths, col = county)) +
      geom_point()
  
    
    
  # Northern District ----
    
  Northern <- VA_US_counties %>%
      filter(county %in% c("Clarke",
                           "Warren",
                           "Loudoun",
                           "Rappahannock",
                           "Prince William",
                           "Fauqier",
                           "Culpeper",
                           "Stafford",
                           "Fredericksburg city",
                           "Spotsylvania",
                           "King George",
                           "Fairfax",
                           "Fairfax city",
                           "Falls Church city",
                           "Arlington",
                           "Alexandria city",
                           "Manassas Park city",
                           "Manassas city"))
    
    Northern_District %>%
      ggplot(aes(date, cases, col = county)) +
      geom_point()
    Northern_District %>%  
      ggplot(aes(date, deaths, col = county)) +
      geom_point()
  
  # Valley District ----
    
  Valley <- VA_US_counties %>%
    filter(county %in% c("Frederick",
                         "Winchester city",
                         "Shenandoah",
                         "Page",
                         "Rockingham",
                         "Harrisonburg city",
                         "Augusta",
                         "Highland",
                         "Staunton city",
                         "Waynesboro city",
                         "Bath",
                         "Rockbridge",
                         "Buena Vista city",
                         "Lexington city",
                         "Alleghany",
                         "Covington city"))
  
  Valley_District %>%
    ggplot(aes(date, cases, col = county)) +
    geom_point()
  Valley_District %>%  
    ggplot(aes(date, deaths, col = county)) +
    geom_point()
  
  sum(Central_District$deaths)
  sum(Eastern_District$deaths)
  sum(Hampton_Roads_District$deaths)
  sum(Northern_District$deaths)
  sum(Southside_District$deaths)
  sum(Southwest_District$deaths)
  sum(West_Central_District$deaths)
  sum(Valley_District$deaths)
  
  ggplot(VA_US_counties, aes(deaths))
  
  
  
  


Col_Covid <- read.csv("~/Documents/GitHub/covid-19-data/colleges/colleges.csv")
glimpse(Col_Covid)
str(Col_Covid)

Col_Covid %>%
  filter(state == "Virginia") %>%
  filter(city %in% c("Charlottesville", "Richmond", "Virginia Beach")) %>%
  group_by(college) %>%
  ggplot(aes(college, cases, col = college)) +
    geom_count() +
    theme(axis.text.x = element_text(angle = 90))
  
  
Col_Covid %>%
  filter(state == "Virginia") %>%
  ggplot(aes(cases, college, col = college)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))
  
  
  
