library(tidyverse)
library(tibbletime)
library(padr)
library(USAboundaries)
# library(rvest)

# Get Data from GitHub repo "CEIDatUGA/COVID-19-DATA" (from Wikipedia)

# CASES

path <- "https://raw.githubusercontent.com/CEIDatUGA/COVID-19-DATA/master/US/US_wikipedia_cases_fatalities/"

states <- state_codes$state_abbr
states <- states[states != ""]
startdate <- as.Date("2019-12-01")

cases.US <- read_csv(paste0(path, "UScases_by_state_wikipedia.csv")) %>% 
  # filter(Date != c("Total") & Date != c("Date")) %>% 
  # select(-c("Conf_New", "Conf_Cml", "deaths_New", "deaths_Cml", "Rec_New", "Rec_Cml", "time_last_update")) %>% 
  select(c(Date, any_of(states))) %>% 
  filter(!is.na(Date)) %>% 
  mutate_at(vars(-Date), as.integer) %>% 
  # mutate(Date = as.Date(Date, "%Y-%m-%d")) %>% 
  padr::pad(start_val = startdate) %>%
  replace(., is.na(.), 0) %>%
  replace(., . < 0, 0) %>%
  tibbletime::as_tbl_time(index=Date) %>% 
  mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
  filter(Date < Sys.Date())

# FATALITIES

fatalities.US <- read_csv(paste0(path, "USfatalities_by_state_wikipedia.csv")) %>% 
  # filter(Date != c("Total") & Date != c("Date")) %>% 
  # select(-c("deaths_New", "deaths_Cml", "time_last_update")) %>% 
  select(c(Date, any_of(states))) %>% 
  filter(!is.na(Date)) %>% 
  mutate_at(vars(-Date), as.integer) %>% 
  # mutate(Date = as.Date(paste(Date, "2020"), "%b %d %Y")) %>% 
  padr::pad(start_val = startdate) %>%
  replace(., is.na(.), 0) %>%
  replace(., . < 0, 0) %>%
  tibbletime::as_tbl_time(index=Date) %>% 
  mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
  filter(Date < Sys.Date())

# Reconcile date range

enddate <- min(max(cases.US$Date),max(fatalities.US$Date))
cases.US <- cases.US %>% filter(Date <= enddate) 
fatalities.US <- fatalities.US %>% filter(Date <= enddate) 

# save RDS
saveRDS(cases.US,"data/cases.US.rds")
saveRDS(fatalities.US,"data/fatalities.US.rds")


# Get County level data
US.counties <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")

US.counties <- US.counties %>%
  mutate(fips = ifelse(is.na(fips),paste0(county,", ",state),fips))

cases.US.counties <- US.counties %>% select(-deaths)
fatalities.US.counties <- US.counties %>% select(-cases)

fips.10deaths <- fatalities.US.counties[fatalities.US.counties$deaths >= 10,]$fips %>% unique

fatalities.US.counties.10deaths <- fatalities.US.counties %>%
  filter(fips %in% fips.10deaths)