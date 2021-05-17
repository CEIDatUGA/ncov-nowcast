library(tidyverse)
library(tibbletime)
library(padr)
library(USAboundaries)
library(here)
# library(rvest)


# Functions -----------------------------------------------------------------------------------

get_US_nowcast_data <- function(datasource, states, startdate) {
  if(datasource == "Wikipedia") {
    # Get Data from GitHub repo "CEIDatUGA/COVID-19-DATA" (from Wikipedia)
    
    # CASES
    path <- "https://raw.githubusercontent.com/CEIDatUGA/COVID-19-DATA/master/US/US_wikipedia_cases_fatalities/"
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
  }
  
  # Get Data from JHU API
  if(datasource == "JHU") {
    
    #function to partially clean us_jhu data
    clean_us_jhu <- function(df){
      usjhu <- df %>% filter(Country_Region == "US") %>%
        dplyr::select(c(-Country_Region, -Lat, -Long_, -UID, -iso2, -iso3, -code3, -Combined_Key)) %>%
        rename(state_name = Province_State)
      return(usjhu)
    }
    
    us_jhu_cases <- readr::read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
    us_jhu_deaths <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv")
    # Clean cases
    us_jhu_cases_clean <- clean_us_jhu(us_jhu_cases) %>%
      tidyr::pivot_longer(cols = c(-state_name, -FIPS, -Admin2), names_to = "Date", values_to = "Cases") %>%
      select(-FIPS) %>% 
      left_join(USAboundaries::state_codes, by = "state_name")
    # Clean deaths
    us_jhu_deaths_clean <- clean_us_jhu(us_jhu_deaths) %>% 
      tidyr::pivot_longer(cols = c(-state_name, -FIPS, -Admin2), names_to = "Date", values_to = "Deaths") %>%
      select(-FIPS) %>% 
      left_join(USAboundaries::state_codes, by = "state_name")

    #combine cases and deaths
    us_jhu_combined <- inner_join(us_jhu_cases_clean, us_jhu_deaths_clean) %>% 
      filter(!state_name %in% c("Diamond Princess", "Grand Princess") )
    us_popsize <- readRDS(here::here("data","us_popsize.rds")) %>%
      select(state_abbr = state, pop_size = total_pop)
    us_jhu_total <- inner_join(us_jhu_combined, us_popsize ) %>%
      mutate(Date = as.Date(as.character(Date),format="%m/%d/%y")) %>%
      group_by(state_name, Admin2) %>% arrange(Date) %>%
      mutate(Daily_Cases = c(0,diff(Cases))) %>%
      mutate(Daily_Deaths = c(0,diff(Deaths))) %>% 
      ungroup() %>%
      rename(Total_Deaths = Deaths, Total_Cases = Cases, Population_Size = pop_size, county_name = Admin2)
    
    #pull state data and aggregate county values
    cases.US <- us_jhu_total %>% 
      select(Date,state_abbr,Daily_Cases) %>% 
      group_by(state_abbr, Date) %>% 
      summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
      pivot_wider(names_from = state_abbr, values_from = "Daily_Cases") %>% 
      select(c(Date, any_of(states))) %>% 
      filter(!is.na(Date)) %>% 
      mutate_at(vars(-Date), as.integer) %>% 
      padr::pad(start_val = startdate) %>%
      replace(., is.na(.), 0) %>%
      replace(., . < 0, 0) %>%
      tibbletime::as_tbl_time(index=Date) %>% 
      mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
      filter(Date < Sys.Date())
    fatalities.US <- us_jhu_total %>% 
      select(Date,state_abbr,Daily_Deaths) %>% 
      group_by(state_abbr, Date) %>% 
      summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
      pivot_wider(names_from = state_abbr, values_from = "Daily_Deaths") %>% 
      select(c(Date, any_of(states))) %>% 
      filter(!is.na(Date)) %>% 
      mutate_at(vars(-Date), as.integer) %>% 
      padr::pad(start_val = startdate) %>%
      replace(., is.na(.), 0) %>%
      replace(., . < 0, 0) %>%
      tibbletime::as_tbl_time(index=Date) %>% 
      mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
      filter(Date < Sys.Date())
  }
  
  # Reconcile date range
  
  enddate <- min(max(cases.US$Date),max(fatalities.US$Date))
  cases.US <- cases.US %>% filter(Date <= enddate) 
  fatalities.US <- fatalities.US %>% filter(Date <= enddate) 
  
  out <- list(cases = cases.US, deaths = fatalities.US)
  return(out)
}

# Get Data ------------------------------------------------------------------------------------
states <- USAboundaries::state_codes$state_abbr
states <- states[states != ""]
startdate <- as.Date("2019-12-01")

# US_nowcast_data <- get_US_nowcast_data(datasource = "Wikipedia", states = states, startdate = startdate)
US_nowcast_data <- get_US_nowcast_data(datasource = "JHU", states = states, startdate = startdate)

# save RDS
saveRDS(US_nowcast_data$cases,"data/cases.US.rds")
saveRDS(US_nowcast_data$deaths,"data/fatalities.US.rds")

# Get County level data
# US.counties <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
# 
# US.counties <- US.counties %>%
#   mutate(fips = ifelse(is.na(fips),paste0(county,", ",state),fips))
# 
# cases.US.counties <- US.counties %>% select(-deaths)
# fatalities.US.counties <- US.counties %>% select(-cases)
# 
# fips.10deaths <- fatalities.US.counties[fatalities.US.counties$deaths >= 10,]$fips %>% unique
# 
# fatalities.US.counties.10deaths <- fatalities.US.counties %>%
#   filter(fips %in% fips.10deaths)