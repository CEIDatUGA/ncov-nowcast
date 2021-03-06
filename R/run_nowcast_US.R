start.time <- proc.time()
library(tidyverse)
source('R/nowcast.R')
source('R/data.R')

# Cores
# Set number of cores used for nowcast to 4 less than number of cores on system, and no less than one.
assign("nowcast_cores", max(1,detectCores()-4L), envir = .GlobalEnv)


# Helper functions ----------------------------------------------------------------------------

get_nowcast <- function(admin, cases, fatalities, params, samplesize, chunksize) {
  cat(as.character(Sys.time()), "###### Starting nowcast for", admin, "\n") # for debugging
  nowcast <- cases %>% select(Date, cases = admin) %>% 
    tbl_time(index = Date) %>% 
    nowcast_from_case_reports(params=params, samplesize = samplesize)

  nowcast$deaths <- pull(fatalities, admin)
  nowcast$cum.deaths <- cumsum(nowcast$deaths)
  cat(as.character(Sys.time()), "Finished nowcast for", admin, "\n") # for debugging
  return(nowcast)
}

do_nowcast <- function(admin, cases, fatalities, params, samplesize = 1.0, chunksize = 30){
  nowcast <- get_nowcast(admin=admin, cases=cases, fatalities=fatalities, params=params, 
                         samplesize=samplesize, chunksize=chunksize)
  saveRDS(nowcast, paste0("data/output/", admin,".nowcast.from.cases.rds"))
}


# get reports -------------------------------------------------------------

cases.US <- readRDS('data/cases.US.rds')
fatalities.US <- readRDS('data/fatalities.US.rds')

# # Cutoff date # for debugging
# lastdate <- '2020-04-15'
# cases.US <- cases.US %>% filter(Date <= lastdate)
# fatalities.US <- fatalities.US %>% filter(Date <= lastdate)

# Cutoff
cases.US.all <- cases.US %>% dplyr::select(Date, cases = US)
fatalities.US.all <- fatalities.US %>% dplyr::select(Date, deaths = US)


# params ------------------------------------------------------------------

params.US <- list(
  admin = "US", # state or country
  
  # https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf
  IFR = .009, 
  
  # ascertainment (calculated and added later)
  # q = 1,
  
  # accepted value
  infectious.period = list(dist="exponential", mean=7),   
  
  # our analysis onset-to-reporting
  effective.infectious.period = list(dist="gamma", mean=5.95, shape=2.2788163), 
  
  # our analysis of US data
  incubation.period = list(dist="gamma", mean=5.89, shape=2.4265511), 
  
  # Our analysis of US data (not used)
  # report.to.death.period = list(dist="skewnormal", mean = 1.33, sd = sqrt(11.06),
  #                               location = -2.388990, # xi
  #                               scale = 4.322, # omega
  #                               shape = 1.447), # alpha
  # 
  # from China data: lognormal distribution with location 2.84, scale 0.52 
  # (Jung et al.) https://doi.org/10.3390/jcm9020523 
  # onset.to.death.period = list(dist="lognormal", mean = 19.9, sd = 11.4)
  
  # from China data: Kenji Mizumoto, Gerardo Chowell. "Estimating the risk of 2019 Novel Coronavirus death during the course of the outbreak in China, 2020" medRxiv preprint. https://doi.org/10.1101/2020.02.19.20025163
  onset.to.death.period = list(dist="gamma", mean = 16, shape = 4)
)


# ascertainment ------------------------------------------------------

# calculate and load ascertainment (time varying q with upper and lower bounds)

params.US$q <- get_ascertainment(cases.US.all, 
                                 fatalities.US.all, 
                                 params.US, 
                                 window = 7, # smoothing window for forecasting 
                                 samplesize = 0.1,  # sample size as fraction of reports
                                 chunksize = 30 # process in 30 day chunks
                                 )
saveRDS(params.US,"data/params.US.rds")

# ## alternate contant q
# 
# params.US$q <- mean(params.US$q$mean.raw)
# 
# ## alternate q = grand mean
# 
# params.US$q <- list(mean = mean(params.US$q$mean.raw), 
#                     upper = mean(params.US$q$upper.raw), 
#                     lower = mean(params.US$q$lower.raw))



# US ----------------------------------------------------------------------

start.time <- proc.time()
do_nowcast(admin = "US", cases = cases.US, fatalities = fatalities.US, params.US, 
           samplesize = 0.05, chunksize = 30)
(proc.time()-start.time)['elapsed']/60

# states ------------------------------------------------------------------

statesvector <- cases.US %>% select(-Date, -US) %>% names()

for(s in statesvector){
  do_nowcast(admin = s, cases = cases.US, fatalities = fatalities.US, params.US, 
             samplesize = 1.0, chunksize = 30)
  }

elapsed <- round(as.numeric((proc.time() - start.time)[3]))
sink("log.txt")
cat(as.character(Sys.time()), 
    "\nNowcast for US completed.\nElapsed time =",
    as.character(lubridate::seconds_to_period(elapsed)))
sink()