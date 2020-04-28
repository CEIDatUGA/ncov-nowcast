library(tidyverse)
library(tibbletime)
library(padr)
# library(rvest)

# Get Data from GitHub repo "CEIDatUGA/COVID-19-DATA" (from Wikipedia)

# CASES

path <- "https://raw.githubusercontent.com/CEIDatUGA/COVID-19-DATA/master/US/US_wikipedia_cases_fatalities/"

cases.US <- read_csv(paste0(path, "UScases_by_state_wikipedia.csv")) %>% 
  # filter(Date != c("Total") & Date != c("Date")) %>% 
  select(-c("Conf_New", "Conf_Cml", "deaths_New", "deaths_Cml", "Rec_New", "Rec_Cml", "time_last_update")) %>% 
  mutate_at(vars(-Date), as.integer) %>% 
  # mutate(Date = as.Date(Date, "%Y-%m-%d")) %>% 
  padr::pad(start_val = as.Date("2019-12-01")) %>%
  replace(., is.na(.), 0) %>%
  tibbletime::as_tbl_time(index=Date) %>% 
  mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
  filter(Date < Sys.Date())

# save RDS
saveRDS(cases.US,"data/cases.US.rds")

# FATALITIES

fatalities.US <- read_csv(paste0(path, "USfatalities_by_state_wikipedia.csv")) %>% 
  # filter(Date != c("Total") & Date != c("Date")) %>% 
  select(-c("deaths_New", "deaths_Cml", "time_last_update")) %>% 
  mutate_at(vars(-Date), as.integer) %>% 
  # mutate(Date = as.Date(paste(Date, "2020"), "%b %d %Y")) %>% 
  padr::pad(start_val = as.Date("2019-12-01")) %>%
  replace(., is.na(.), 0) %>%
  tibbletime::as_tbl_time(index=Date) %>% 
  mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
  filter(Date < Sys.Date())

# save RDS
saveRDS(fatalities.US,"data/fatalities.US.rds")






# Get case reports for US states
# read data directly from google sheet.

# sheet <- "https://docs.google.com/spreadsheets/d/10NYNf4UX7R_nNdM-jkbo7A5SGkB8sUVMIiqU7e_-miU/edit#gid=0"
# 
# US <- googlesheets4::sheets_read(sheet)
# 
# US <- US %>% 
#   rename_all(list(~make.names(.))) %>% 
#   mutate_if(is.list, fixcolumns) %>% 
#   mutate(Date = as.Date(Date)) %>% 
#   padr::pad(start_val = as.Date("2019-12-01")) %>%
#   replace(., is.na(.), 0) %>%
#   tibbletime::as_tbl_time(index=Date) %>% 
#   mutate(US = rowSums(.[,-1],na.rm = TRUE))


# Scrape Wikipedia

# url <- "https://en.wikipedia.org/wiki/Template:2019%E2%80%9320_coronavirus_pandemic_data/United_States_medical_cases"
# 
# cases.US.wikipedia <- url %>%
#   html() %>%
#   html_nodes(xpath='//*[@id="mw-content-text"]/div/div[2]/table') %>%
#   html_table(fill = TRUE)
# cases.US.wikipedia <- cases.US.wikipedia[[1]]
# names(cases.US.wikipedia) <- cases.US.wikipedia[1,]
# cases.US.wikipedia <- cases.US.wikipedia[2:(nrow(cases.US.wikipedia)-5),
#                                          1:(ncol(cases.US.wikipedia)-6)]
# 
# cases.US.wikipedia$Date <- paste(cases.US.wikipedia$Date, "2020") %>% 
#   as.Date("%b %d %Y")
# 
# cases.US <- cases.US.wikipedia %>%
#   mutate_if(sapply(cases.US.wikipedia, is.character), as.numeric) %>% 
#   rename_all(list(~make.names(.))) %>% 
#   padr::pad(start_val = as.Date("2019-12-01")) %>%
#   replace(., is.na(.), 0) %>%
#   tibbletime::as_tbl_time(index=Date) %>% 
#   mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
#   filter(Date < Sys.Date())
# 
# remove(cases.US.wikipedia)
# saveRDS(cases.US,"data/cases.US.rds")

# Get fatalities reports for US states

# sheets_file <- "https://docs.google.com/spreadsheets/d/1SC28cM52m6s1gTJutpvFxadT9GGu-li90tsqkGuaM48/edit#gid=219891324"
# sheets_tab <- "Cumulative fatalities reported by country"
#   
# cum.fatalities <- googlesheets4::sheets_read(ss = sheets_file, sheet = sheets_tab) %>%
#   dplyr::select(-Sources, -Contributor) %>% 
#   rename_all(list(~make.names(.))) %>%
#   # mutate_if(is.list, fixcolumns) %>%
#   mutate(Date = as.Date(Date)) %>%
#   padr::pad(start_val = as.Date("2019-12-01")) %>%
#   replace(., is.na(.), 0) %>%
#   tibbletime::as_tbl_time(index=Date) %>%
#   mutate(World = rowSums(.[,-1],na.rm = TRUE))
# 
# cum.fatalities.diff <- cum.fatalities %>% select(-Date) %>% 
#   as.matrix() %>% diff() %>% as_tibble()
# 
# # New deaths by day for all countries
# fatalities <- bind_rows(cum.fatalities[1,-1],cum.fatalities.diff) %>% bind_cols(cum.fatalities[,1],.)
# 
# # New deaths by day for the US
# fatalities.US <- fatalities %>% dplyr::select(Date, deaths = US)


# url <- "https://en.wikipedia.org/wiki/Template:2019%E2%80%9320_coronavirus_pandemic_data/United_States_medical_cases"
# 
# fatalities.US.wikipedia <- url %>%
#   html() %>%
#   html_nodes(xpath='//*[@id="mw-content-text"]/div/div[3]/table') %>%
#   html_table(fill = TRUE)
# fatalities.US.wikipedia <- fatalities.US.wikipedia[[1]]
# names(fatalities.US.wikipedia) <- fatalities.US.wikipedia[1,]
# fatalities.US.wikipedia <- fatalities.US.wikipedia[2:(nrow(fatalities.US.wikipedia)-5),
#                                                    1:(ncol(fatalities.US.wikipedia)-2)]
# 
# fatalities.US.wikipedia$Date <- paste(fatalities.US.wikipedia$Date, "2020") %>% 
#   as.Date("%b %d %Y")
# 
# fatalities.US <- fatalities.US.wikipedia %>%
#   mutate_if(sapply(fatalities.US.wikipedia, is.character), as.numeric) %>% 
#   rename_all(list(~make.names(.))) %>% 
#   padr::pad(start_val = as.Date("2019-12-01")) %>%
#   replace(., is.na(.), 0) %>%
#   tibbletime::as_tbl_time(index=Date) %>% 
#   mutate(US = rowSums(.[,-1],na.rm = TRUE)) %>% 
#   filter(Date < Sys.Date())
# 
# remove(fatalities.US.wikipedia)
# saveRDS(fatalities.US,"data/fatalities.US.rds")

# Get County level data
US.counties <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")

US.counties <- US.counties %>%
  mutate(fips = ifelse(is.na(fips),paste0(county,", ",state),fips))

cases.US.counties <- US.counties %>% select(-deaths)
fatalities.US.counties <- US.counties %>% select(-cases)

fips.10deaths <- fatalities.US.counties[fatalities.US.counties$deaths >= 10,]$fips %>% unique

fatalities.US.counties.10deaths <- fatalities.US.counties %>%
  filter(fips %in% fips.10deaths)