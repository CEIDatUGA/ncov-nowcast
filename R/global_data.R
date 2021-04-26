#library(coronavirus)
library(dplyr)
library(tidyverse)
library(tibbletime)

# JHU COVID-19 DATA

# Pull data from Johns Hopkins University Center for Systems Science and Engineering (JHU CCSE) Coronavirus website

#data("coronavirus") # broken command
coronavirus <- read.csv("https://raw.githubusercontent.com/RamiKrispin/coronavirus/master/csv/coronavirus.csv",stringsAsFactors = FALSE)
colnames(coronavirus)[1] <-"Date"  
coronavirus$Date <- as.Date(coronavirus$Date)

# CONFIRMED CASES

#Extract
# Some countries are considered provinces of other countries and will require seperating them

Confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country != "Australia") %>% 
  filter(country != "Canada") %>%
  filter(country != "China") %>%
  filter(country != "Denmark") %>%
  filter(country != "France") %>%
  filter(country != "Netherlands") %>%
  filter(country != "United Kingdom") %>%
  select(Date, country, cases) %>%
  spread(key = country, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)

#  Create single country-level calculations for the removed countries 

# Australia
Australia.confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country == "Australia") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  mutate(Australia = rowSums(.[,-1],na.rm = TRUE)) %>%
  select(Date, Australia) %>%
  tibbletime::as_tbl_time(index=Date)

# Canada
Canada.confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country == "Canada") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  mutate(Canada = rowSums(.[,-1],na.rm = TRUE)) %>%
  select(Date, Canada) %>%
  tibbletime::as_tbl_time(index=Date)
  
# China
China.confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country == "China") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  mutate(China = rowSums(.[,-1],na.rm = TRUE)) %>%
  select(Date, China) %>%
  tibbletime::as_tbl_time(index=Date)


# *** DENMARK, FAROE ISLANDS, GREENLAND
# *** Need to create  indvidual objects for these countries

# Denmark & provinces
Denmark.confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country == "Denmark") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)
names(Denmark.confirmed)[2] <- "Denmark"

# France & provinces 
France.confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country == "France") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)
names(France.confirmed)[2] <- "France"

# Netherlands & provinces
Netherlands.confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country == "Netherlands") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date) 
names(Netherlands.confirmed)[2] <- "Netherlands"

# United Kingdom & provinces
UK.confirmed <- coronavirus %>%
  filter(type == "confirmed") %>%
  filter(country == "United Kingdom") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)
names(UK.confirmed)[2] <- "UK"

# COMBINE AND ALPHEBETIZE

cases.global <- merge(merge(merge(merge(merge(merge(merge(
  Confirmed,
  Australia.confirmed, by = "Date"),
  Canada.confirmed, by = "Date"),
  China.confirmed, by = "Date"),
  Denmark.confirmed, by = "Date"),
  France.confirmed, by = "Date"),
  Netherlands.confirmed, by = "Date"),
  UK.confirmed, by = "Date")
cases.global.dates <- cases.global[,1]
cases.ordered<- cases.global[,-1][,order(colnames(cases.global[,-1]))]
cases.global <- cbind(cases.global.dates, cases.ordered)
names(cases.global)[1] <- "Date"


cases.global <- cases.global %>% 
  tibbletime::as_tbl_time(index=Date) %>% 
  mutate(Global = rowSums(.[,-1],na.rm = TRUE))

saveRDS(cases.global,"data/cases.global.rds")
#--------------------------------------------------------------------------------------------

# # FATALITIES

#Extract 
Death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country != "Australia") %>% 
  filter(country != "Canada") %>%
  filter(country != "China") %>%
  filter(country != "Denmark") %>%
  filter(country != "France") %>%
  filter(country != "Netherlands") %>%
  filter(country != "United Kingdom") %>%
  select(Date, country, cases) %>%
  spread(key = country, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)

#  Create single country-level calculations for the removed countries 

# Australia
Australia.death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country == "Australia") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  mutate(Australia = rowSums(.[,-1],na.rm = TRUE)) %>%
  select(Date, Australia) %>%
  tibbletime::as_tbl_time(index=Date)

# Canada
Canada.death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country == "Canada") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  mutate(Canada = rowSums(.[,-1],na.rm = TRUE)) %>%
  select(Date, Canada) %>%
  tibbletime::as_tbl_time(index=Date)

# China
China.death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country == "China") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  mutate(China = rowSums(.[,-1],na.rm = TRUE)) %>%
  select(Date, China) %>%
  tibbletime::as_tbl_time(index=Date)


# *** DENMARK, FAROE ISLANDS, GREENLAND
# *** Need to create  indvidual objects for these countries

# Denmark & provinces
Denmark.death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country == "Denmark") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)
names(Denmark.death)[2] <- "Denmark"

# France & provinces 
France.death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country == "France") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)
names(France.death)[2] <- "France"

# Netherlands & provinces
Netherlands.death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country == "Netherlands") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date) 
names(Netherlands.death)[2] <- "Netherlands"

# United Kingdom & provinces
UK.death <- coronavirus %>%
  filter(type == "death") %>%
  filter(country == "United Kingdom") %>%
  select(Date, province, cases) %>%
  spread(key = province, value = cases) %>%
  tibbletime::as_tbl_time(index=Date)
names(UK.death)[2] <- "UK"

# COMBINE AND ALPHEBETIZE

fatalities.global <- merge(merge(merge(merge(merge(merge(merge(
  Death,
  Australia.death, by = "Date"),
  Canada.death, by = "Date"),
  China.death, by = "Date"),
  Denmark.death, by = "Date"),
  France.death, by = "Date"),
  Netherlands.death, by = "Date"),
  UK.death, by = "Date")
fatalities.global.dates <- fatalities.global[,1]
fatalities.ordered<- fatalities.global[,-1][,order(colnames(fatalities.global[,-1]))]
fatalities.global <- cbind(fatalities.global.dates, fatalities.ordered)
names(fatalities.global)[1] <- "Date"


fatalities.global <- fatalities.global %>% 
  tibbletime::as_tbl_time(index=Date) %>% 
  mutate(Global = rowSums(.[,-1],na.rm = TRUE))

saveRDS(fatalities.global,"data/fatalities.global.rds")
