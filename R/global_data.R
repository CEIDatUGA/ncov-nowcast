# Code used to extract John Hopkins University data related to COVID-19 epidemic, and reconfigure data.frames to fit nowcast model. 


# Import necessary packages.

library( "dplyr" )
library( "tibbletime" )
library( "tidyverse" )

# variable to define the the beginning/ first row of the analysis 

startdate <- as.Date("2019-12-01")


#### ___________________________________________________________________________________________________________________________________________________________________________________________________


# JHU CASES  
global_jhu_cases <- readr::read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

# Reconfigure cases for nowcast 
cases.global <- select( global_jhu_cases , -c( `Province/State`, Lat , Long  ) ) %>% 
  group_by( `Country/Region` ) %>% 
  summarize_all( sum)  %>%
  rename( Country = `Country/Region` ) %>%
  tidyr::pivot_longer( cols = c( -Country ) , names_to = "Date" , values_to = "Cases" ) %>%
  filter(!Country %in% c("Diamond Princess", "Grand Princess") ) %>%
  mutate( Date = as.Date( as.character( Date ) , format = "%m/%d/%y" ) ) %>%
  group_by( Country ) %>% 
  arrange( Date ) %>%
  mutate( Daily_Cases = c( 0 , diff( Cases ) ) ) %>%
  select( -Cases) %>% 
  pivot_wider( names_from = Country , values_from = "Daily_Cases" ) %>%
  padr::pad(start_val = startdate) %>%
  replace( . , is.na(.) , 0 ) %>%
  replace( . , . < 0 , 0 ) %>%
  mutate( Global = rowSums( .[,-1] , na.rm = TRUE ) ) %>%
  tibbletime::as_tbl_time(index=Date)  


#### ___________________________________________________________________________________________________________________________________________________________________________________________________


# JHU FATALITIES
global_jhu_deaths <- readr::read_csv( "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv" )

# Reconfigure cases for nowcast
cases.deaths <- select( global_jhu_deaths , -c( `Province/State`, Lat , Long  ) ) %>% 
  group_by( `Country/Region` ) %>% 
  summarize_all( sum)  %>%
  rename( Country = `Country/Region` ) %>%
  tidyr::pivot_longer( cols = c( -Country ) , names_to = "Date" , values_to = "Cases" ) %>%
  filter(!Country %in% c("Diamond Princess", "Grand Princess") ) %>%
  mutate( Date = as.Date( as.character( Date ) , format = "%m/%d/%y" ) ) %>%
  group_by( Country ) %>% 
  arrange( Date ) %>%
  mutate( Daily_Cases = c( 0 , diff( Cases ) ) ) %>%
  select( -Cases) %>% 
  pivot_wider( names_from = Country , values_from = "Daily_Cases" ) %>%
  padr::pad(start_val = startdate ) %>%
  replace( . , is.na(.) , 0 ) %>%
  replace( . , . < 0 , 0 ) %>%
  mutate( Global = rowSums( .[,-1] , na.rm = TRUE ) ) %>%
  tibbletime::as_tbl_time( index = Date )  


#### ___________________________________________________________________________________________________________________________________________________________________________________________________


# Save data.frames

saveRDS( cases.global , "data/cases.global.rds" )
saveRDS( cases.deaths , "data/fatalities.global.rds" )







