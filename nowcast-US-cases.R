## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------

source('R/nowcast.R')
source('R/data.R')



## ----params, include=FALSE, message = FALSE, warning = FALSE, cache = FALSE-----------------------------------------------------------------------------------------------
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



## ----data, include=FALSE, message = FALSE, warning = FALSE, cache = FALSE-------------------------------------------------------------------------------------------------

cases.US <- readRDS('data/cases.US.rds')
fatalities.US <- readRDS('data/fatalities.US.rds')

# # Cutoff date # for debugging
lastdate <- '2020-04-15'
cases.US <- cases.US %>% filter(Date <= lastdate)
fatalities.US <- fatalities.US %>% filter(Date <= lastdate)

# Cutoff
cases.US.all <- cases.US %>% dplyr::select(Date, cases = US)
fatalities.US.all <- fatalities.US %>% dplyr::select(Date, deaths = US)



## ----ascertainment, include=FALSE, message = FALSE, warning = FALSE, cache = FALSE----------------------------------------------------------------------------------------


# calculate and load ascertainment (time varying q with upper and lower bounds)

params.US$q <- get_ascertainment(cases.US.all, fatalities.US.all, params.US, window = 7)
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



## ----nowcast_US, include=FALSE, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE-------------------------------------------------------------------------------

params.US <- readRDS("data/params.US.rds")

get_nowcast <- function(admin, cases, fatalities, params) {
  nowcast <- cases %>% select(Date, cases = admin) %>% 
    tbl_time(index = Date) %>% nowcast_from_case_reports(params)
  
  nowcast$deaths <- pull(fatalities, admin)
  nowcast$cum.deaths <- cumsum(nowcast$deaths)
  return(nowcast)
}

# start.time <- proc.time()
US <- get_nowcast(admin = "US", cases = cases.US, fatalities = fatalities.US, params.US)
saveRDS(US,"data/US.nowcast.from.cases.rds")
# (proc.time()-start.time)['elapsed']/60

# Export single plot
p <- plot_nowcast_from_case_reports(US, plotcumulative = TRUE, maxy = 10^8, legend = FALSE)
htmlwidgets::saveWidget(p, "USnowcast_plot.html")


## ----US, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
US <- readRDS("data/US.nowcast.from.cases.rds")
# Plot_R_effective(US, legend = FALSE)
plot_nowcast_from_case_reports(US, maxy = 10^8, legend = FALSE)


## ----nowcast_states, include=FALSE, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE---------------------------------------------------------------------------
# start.time <- proc.time()

# cases.US <- readRDS('data/cases.US.rds')
# fatalities.US <- readRDS('data/fatalities.US.rds')
# params.US <- readRDS("data/params.US.rds")

AK <- get_nowcast("AK", cases.US, fatalities.US, params.US)
AL <- get_nowcast("AL", cases.US, fatalities.US, params.US)
AR <- get_nowcast("AR", cases.US, fatalities.US, params.US)
AZ <- get_nowcast("AZ", cases.US, fatalities.US, params.US)
CA <- get_nowcast("CA", cases.US, fatalities.US, params.US)
CO <- get_nowcast("CO", cases.US, fatalities.US, params.US)
CT <- get_nowcast("CT", cases.US, fatalities.US, params.US)
DC <- get_nowcast("DC", cases.US, fatalities.US, params.US)
DE <- get_nowcast("DE", cases.US, fatalities.US, params.US)
FL <- get_nowcast("FL", cases.US, fatalities.US, params.US)
GA <- get_nowcast("GA", cases.US, fatalities.US, params.US)
saveRDS(GA,"data/GA.nowcast.from.cases.rds")
# Export single plot
p <- plot_nowcast_from_case_reports(GA, plotcumulative = TRUE, maxy = 10^7, legend = FALSE)
htmlwidgets::saveWidget(p, "GAnowcast_plot.html")

HI <- get_nowcast("HI", cases.US, fatalities.US, params.US)
IA <- get_nowcast("IA", cases.US, fatalities.US, params.US)
ID <- get_nowcast("ID", cases.US, fatalities.US, params.US)
IL <- get_nowcast("IL", cases.US, fatalities.US, params.US)
IN <- get_nowcast("IN", cases.US, fatalities.US, params.US)
KS <- get_nowcast("KS", cases.US, fatalities.US, params.US)
KY <- get_nowcast("KY", cases.US, fatalities.US, params.US)
LA <- get_nowcast("LA", cases.US, fatalities.US, params.US)
MA <- get_nowcast("MA", cases.US, fatalities.US, params.US)
MD <- get_nowcast("MD", cases.US, fatalities.US, params.US)
ME <- get_nowcast("ME", cases.US, fatalities.US, params.US)
MI <- get_nowcast("MI", cases.US, fatalities.US, params.US)
MN <- get_nowcast("MN", cases.US, fatalities.US, params.US)
MO <- get_nowcast("MO", cases.US, fatalities.US, params.US)
MS <- get_nowcast("MS", cases.US, fatalities.US, params.US)
MT <- get_nowcast("MT", cases.US, fatalities.US, params.US)
NC <- get_nowcast("NC", cases.US, fatalities.US, params.US)
ND <- get_nowcast("ND", cases.US, fatalities.US, params.US)
NE <- get_nowcast("NE", cases.US, fatalities.US, params.US)
NH <- get_nowcast("NH", cases.US, fatalities.US, params.US)
NJ <- get_nowcast("NJ", cases.US, fatalities.US, params.US)
NM <- get_nowcast("NM", cases.US, fatalities.US, params.US)
NV <- get_nowcast("NV", cases.US, fatalities.US, params.US)
NY <- get_nowcast("NY", cases.US, fatalities.US, params.US)
OH <- get_nowcast("OH", cases.US, fatalities.US, params.US)
OK <- get_nowcast("OK", cases.US, fatalities.US, params.US)
OR <- get_nowcast("OR", cases.US, fatalities.US, params.US)
PA <- get_nowcast("PA", cases.US, fatalities.US, params.US)
RI <- get_nowcast("RI", cases.US, fatalities.US, params.US)
SC <- get_nowcast("SC", cases.US, fatalities.US, params.US)
SD <- get_nowcast("SD", cases.US, fatalities.US, params.US)
TN <- get_nowcast("TN", cases.US, fatalities.US, params.US)
TX <- get_nowcast("TX", cases.US, fatalities.US, params.US)
UT <- get_nowcast("UT", cases.US, fatalities.US, params.US)
VA <- get_nowcast("VA", cases.US, fatalities.US, params.US)
VT <- get_nowcast("VT", cases.US, fatalities.US, params.US)
WA <- get_nowcast("WA", cases.US, fatalities.US, params.US)
WI <- get_nowcast("WI", cases.US, fatalities.US, params.US)
WY <- get_nowcast("WY", cases.US, fatalities.US, params.US)
GU <- get_nowcast("GU", cases.US, fatalities.US, params.US)
MP <- get_nowcast("MP", cases.US, fatalities.US, params.US)
PR <- get_nowcast("PR", cases.US, fatalities.US, params.US)
VI <- get_nowcast("VI", cases.US, fatalities.US, params.US)
# (proc.time()-start.time)['elapsed']/60


## ----AK, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(AK,  maxy = 10^7)


## ----AL, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(AL,  maxy = 10^7)


## ----AR, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(AR,  maxy = 10^7)


## ----AZ, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(AZ,  maxy = 10^7)


## ----CA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(CA,  maxy = 10^7)


## ----CO, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(CO,  maxy = 10^7)


## ----CT, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(CT,  maxy = 10^7)


## ----DC, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(DC,  maxy = 10^7)


## ----DE, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(DE,  maxy = 10^7)


## ----FL, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(FL,  maxy = 10^7)


## ----GA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(GA,  maxy = 10^7)


## ----HI, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(HI,  maxy = 10^7)


## ----IA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(IA,  maxy = 10^7)


## ----ID, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(ID,  maxy = 10^7)


## ----IL, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(IL,  maxy = 10^7)


## ----IN, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(IN,  maxy = 10^7)


## ----KS, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(KS,  maxy = 10^7)


## ----KY, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(KY,  maxy = 10^7)


## ----LA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(LA,  maxy = 10^7)


## ----MA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MA,  maxy = 10^7)


## ----MD, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MD,  maxy = 10^7)


## ----ME, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(ME,  maxy = 10^7)


## ----MI, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MI,  maxy = 10^7)


## ----MN, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MN,  maxy = 10^7)


## ----MO, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MO,  maxy = 10^7)


## ----MS, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MS,  maxy = 10^7)


## ----MT, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MT,  maxy = 10^7)


## ----NC, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(NC,  maxy = 10^7)


## ----ND, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(ND,  maxy = 10^7)


## ----NE, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(NE,  maxy = 10^7)


## ----NH, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(NH,  maxy = 10^7)


## ----NJ, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(NJ,  maxy = 10^7)


## ----NM, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(NM,  maxy = 10^7)


## ----NV, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(NV,  maxy = 10^7)


## ----NY, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(NY,  maxy = 10^7)


## ----OH, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(OH,  maxy = 10^7)


## ----OK, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(OK,  maxy = 10^7)


## ----OR, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(OR,  maxy = 10^7)


## ----PA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(PA,  maxy = 10^7)


## ----RI, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(RI,  maxy = 10^7)


## ----SC, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(SC,  maxy = 10^7)


## ----SD, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(SD,  maxy = 10^7)


## ----TN, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(TN,  maxy = 10^7)


## ----TX, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(TX,  maxy = 10^7)


## ----UT, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(UT,  maxy = 10^7)


## ----VA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(VA,  maxy = 10^7)


## ----VT, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(VT,  maxy = 10^7)


## ----WA, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(WA,  maxy = 10^7)


## ----WI, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(WI,  maxy = 10^7)


## ----WY, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(WY,  maxy = 10^7)


## ----GU, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(GU,  maxy = 10^7)


## ----MP, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(MP)


## ----PR, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(PR,  maxy = 10^7)


## ----VI, echo=FALSE, message = FALSE, warning = FALSE, cache = FALSE, out.width = '100%'----------------------------------------------------------------------------------
plot_nowcast_from_case_reports(VI,  maxy = 10^7)

