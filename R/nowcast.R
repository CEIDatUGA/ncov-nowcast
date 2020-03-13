library(tidyverse)
library(tibbletime)
library(padr)
library(tvReg)
library(forecast)
source('R/package.R')
source('R/deconcolve.R')
source('R/simple_tdar.R')
source('R/package_tv.R')

# functions

# function to trim NA's from a vector
trim_na <- function(x) {x[!is.na(x)]}

# function to pad a vector with NA's
pad_na <- function(x,n,front=FALSE) {
  if(front==TRUE){
    x <- c(rep.int(NA,n),x)
  }else{
    length(x) <- length(x)+n
  }
  return(x)
}

# function to replace tail of vector or dataframe with NA's
tail_na <- function(data, tail, columns=NULL) {
  d <- data
  if (is.vector(d)) {
    n <- length(d)
    # replace(data, (n-tail):n, NA)
    d[(n-tail[1]+1):n] <- NA
  }else{
    if (is.data.frame(d)) {
      n <- nrow(d)
      if(!is.null(columns)){c <- columns}else{c <- 1:length(data)}
      d[(n-tail[1]+1):n,c] <- NA
    }
  }
  return(d)
}

# function to calculate confidence intervals across rows of a dataframe
rowwise_confidence_intervals <- function(df,interval) {
  out <- list()
  for (i in 1:nrow(df)) {
    v <- unlist(df[i,])
    v_sd <- sd(v)
    v_n <- length(v)
    v_mean <- mean(v)
    error <- qt((interval + 1)/2, df = v_n - 1) * v_sd / sqrt(v_n)
    out$lower[i] <- v_mean - error
    out$upper[i] <- v_mean + error
  }
  return(out)
}

# Simulate from time varying distribution

# ## Exponential distribution
# backcast_exponential <- function(x,mean,rate=NULL) {
#   
#   offsets <- lapply(x$cases, rexp, rate=1/mean)
#   
# }
# 
# backcast_exponential(cases.china, mean = gamma)
# 
# ## Gamma distribution
# backcast_gamma <- function(x,n,mean,shape) {
#   
#   
#   rgamma(n=n,shape=shape,scale=mean/shape)
# }
# get_exposure_onset(10, mean = incubation.mean, shape=incubation.shape)


## line list

backward_simulate_linelist <- function(dates, counts, interval) {
  
  input <- tibble(exit.date = dates, count = counts)

  if(is.numeric(interval)){
    interval <- interval[1]
    get_intervals <- function(interval,x) {
      rep(x,rate=1/interval)
    }
  }else{
    if(interval$dist == "exponential") {
      get_intervals <- function(x,interval) {
        rexp(x,rate=1/interval$mean)
      }
    }
    if(interval$dist == "gamma") {
      get_intervals <- function(x,interval) {
        rgamma(x,shape=interval$shape,scale=interval$mean/interval$shape)
      }
    }
  }
  
  intervals <- lapply(input$count, get_intervals, interval) %>% unlist()
  
  linelist <- input %>% 
    uncount(count) %>% 
    mutate(interval = intervals) %>% 
    mutate(onset.date = exit.date - interval + 1) %>% 
    select(onset.date, interval, exit.date)

  linelist
}

## Get onset curve for range of dates from linelist

get_onset_curve <- function(dates, linelist, interval) {
  dates <- range(dates)
  onset.curve <- linelist %>% 
    mutate(onset.date = as.Date(lubridate::floor_date(onset.date))) %>% 
    count(onset.date) %>% 
    padr::pad(interval = "day", start_val = dates[1L], end_val = dates[2L])
  
  if(is.numeric(interval)){
    interval <- list(mean=interval[1], sd=0)
  }else{
    if(interval$dist == "exponential") {
      interval$sd <- sqrt((1/interval$mean)^2) # exponential
    }
    if(interval$dist == "gamma") {
      interval$scale <- interval$mean/interval$shape # gamma
      interval$sd <- sqrt(interval$shape*interval$scale^2) # gamma
    }
  }
  
  # set last several valsues to NA
  onset.curve <- tail_na(onset.curve,
                         round(interval$mean+interval$sd*2) 
                         )
  onset.curve
}


## Get value of state at a single "date" from "linelist"
get_state <- function(date, linelist) {
  linelist[as.Date(lubridate::floor_date(linelist$onset.date)) <= date
           & date < linelist$exit.date,] %>% nrow
}

## Get state curve for range of dates from linelist
get_state_curve <- function(dates, linelist, interval){
  values <- lapply(X = dates, FUN = get_state, linelist) %>% unlist
  
  if(is.numeric(interval)){
    interval <- list(mean=interval[1], sd=0)
  }else{
    if(interval$dist == "exponential") {
      interval$sd <- sqrt((1/interval$mean)^2) # exponential
    }
    if(interval$dist == "gamma") {
      interval$scale <- interval$mean/interval$shape # gamma
      interval$sd <- sqrt(interval$shape*interval$scale^2) # gamma
    }
  }
  
  values <- tail_na(values, round(interval$mean+interval$sd*2))
  
  state.curve <- tibble(
    date = dates, 
    value = values
    )
  state.curve
}

