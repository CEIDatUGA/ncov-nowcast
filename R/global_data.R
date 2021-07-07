#library(coronavirus)
library(dplyr)
library(tidyverse)
library(tibbletime)

library("tidyverse")

global_jhu_cases <- readr::read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")


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
  mutate(Global = rowSums(.[,-1],na.rm = TRUE)) %>%
  tibbletime::as_tbl_time(index=Date)  


sum(is.na(cases.global))

saveRDS( cases.global , "data/cases.global.rds" )

#### ___________________________________________________________________________________________________________________________________________________________________________________________________

global_jhu_deaths <- readr::read_csv( "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv" )

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
  mutate(Global = rowSums(.[,-1],na.rm = TRUE)) %>%
  tibbletime::as_tbl_time(index=Date)  

saveRDS( cases.global , "data/fatalities.global.rds" )
