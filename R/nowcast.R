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

# fix list columns
fixcolumns <- function(x) {
  if(is.list(x)){
    is.na(x) <- lengths(x) == 0
    x <- gsub("[()]","",x)
    as.numeric(unlist(x))
  }
}

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

### Store interval distribution in rows of database


## Simulate a line list from dates, counts (such as case reports), and interval distribution

backward_simulate_linelist <- function(dates, counts, interval) {
  
  input <- tibble(exit.date = dates, count = counts)

  if(is.numeric(interval)){
    interval <- interval[1]
    get_intervals <- function(x,interval) {
      rep(x,interval)
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

## back propagate a line list using an interval distribution

backward_propagate_linelist <- function (linelist, interval) {
  
  if(is.numeric(interval)){
    interval <- interval[1]
    get_intervals <- function(x,interval) {
      rep(x,interval)
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
  
  intervals <- get_intervals(nrow(linelist),interval)
  
  out.linelist <- linelist %>% 
    mutate(exit.date = onset.date) %>% 
    mutate(interval = intervals) %>% 
    mutate(onset.date = exit.date - interval)
  
  out.linelist
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

# forecast to present

tvar_forecast_to_present <- function(curve,lag=1) {
  curve.trimmed <- trim_na(curve)
  trimmed.length <- length(curve.trimmed)
  forecast.length <- length(curve)-trimmed.length
  
  # log curve
  log.curve.trimmed <- log(curve.trimmed)
  
  # replace infinities
  log.curve.trimmed[is.infinite(log.curve.trimmed)] <- 0

  log.curve.model.tvar <- tvReg::tvAR(log.curve.trimmed, 
                                  p = lag,  # number of lags
                                  type = "none",  # model does not contain intercept
                                  est="ll", # "local linear" non parametric estimation method
                                  tkernel = "Gaussian")
  
  log.curve.forecast <- forecast::forecast(log.curve.model.tvar$fitted,forecast.length)
  
  # antilog fit 
  # antilog forecast
  
  forecast <- tibble(fit = exp(log.curve.model.tvar$fitted) %>% 
                       pad_na(lag, front = TRUE) %>% 
                       pad_na(forecast.length),
                     mean = exp(log.curve.forecast$mean) %>% pmax(0) %>% 
                       pad_na(trimmed.length, front = TRUE),
                     upper95 = exp(log.curve.forecast$upper[,'95%']) %>% pmax(0) %>% 
                       pad_na(trimmed.length, front = TRUE),
                     lower95 = exp(log.curve.forecast$lower[,'95%']) %>% pmax(0) %>%
                       pad_na(trimmed.length, front = TRUE)
                     )
  forecast
}

# nowcast

nowcast_from_case_reports <- function(casereports, params) {
  database <- casereports
  database$q <- US.params$q
  database$a <- US.params$a
  database$c <- US.params$c
  
  # I linelist
  I.linelist <- backward_simulate_linelist(dates = database$Date,
                                           counts = database$cases,
                                           interval = params$effective.infectious.period)
  # I onset
  I.onset <- get_onset_curve(dates = database$Date,
                             linelist = I.linelist,
                             interval = params$effective.infectious.period)
  database$I.onset <- I.onset$n

  # I (detected)
  I_d <- get_state_curve(dates = database$Date,
                         linelist = I.linelist,
                         interval = params$effective.infectious.period)
  database$I_d <- I_d$value

  # I
  database <- database %>% mutate(I = I_d/(1-a))
  
  # E linelist
  E.linelist <- backward_propagate_linelist(I.linelist,
                                            interval = params$incubation.period)
  
  # E onset
  E.onset <- get_onset_curve(dates = database$Date,
                             linelist = E.linelist,
                             interval = params$incubation.period)  
  database$E.onset <- E.onset$n
  
  # E (detected)
  E_d <- get_state_curve(dates = database$Date,
                         linelist = E.linelist,
                         interval = params$incubation.period)
  database$E_d <- E_d$value
  
  # E
  database <- database %>% mutate(E = E_d/(1-a))
  
  # forecasting I
  I.forecast <- tvar_forecast_to_present(database$I)
  database$I.fit <- I.forecast$fit
  database$I.forecast.mean <- I.forecast$mean
  database$I.forecast.lower95 <- I.forecast$lower95
  database$I.forecast.upper95 <- I.forecast$upper95
  
  # forecasting E
  E.forecast <- tvar_forecast_to_present(database$E)
  database$E.fit <- E.forecast$fit
  database$E.forecast.mean <- E.forecast$mean
  database$E.forecast.lower95 <- E.forecast$lower95
  database$E.forecast.upper95 <- E.forecast$upper95
  
  # nowcast
  
  database <- database %>% 
    mutate(nowcast.mean = rowSums(
      dplyr::select(., I, I.forecast.mean, E, E.forecast.mean), na.rm = TRUE)
    ) %>% 
    mutate(nowcast.lower95 = rowSums(
      dplyr::select(., I, I.forecast.lower95, E, E.forecast.lower95), na.rm = TRUE)
    ) %>% 
    mutate(nowcast.upper95 = rowSums(
      dplyr::select(., I, I.forecast.upper95, E, E.forecast.upper95), na.rm = TRUE)
    )
  return(database)
}

# plot nowcast

plot_nowcast_from_case_reports <- function(database) {
  
  col.cases <- 'rgba(0, 0, 0, .75)'
  col.I <- 'rgba(230, 7, 7, .75)'
  col.I.ci <- 'rgba(230, 7, 7, .25)'
  col.E <- 'rgba(7, 164, 181, 0.75)'
  col.E.ci <- 'rgba(7, 164, 181, 0.0)'
  col.nowcast <- 'rgba(7, 7, 230, 0.75)'
  col.nowcast.ci <- 'rgba(7, 7, 230, 0.25)'
  
  ci.lwd <- .5
  mean.lwd <- 1
  data.lwd <- 2
  
  p_nowcast <- plotly::plot_ly(data = database, x = ~Date , y = ~cases, type = 'scatter',
                               name = 'confirmed', mode = 'lines',
                               line = list(color = col.cases, width = data.lwd)
                               ) %>% 
    plotly::add_trace(y = ~I, 
                      name = 'I', mode = 'lines',
                      line = list(color = col.I, width = data.lwd)) %>%
    plotly::add_trace(y = ~I.forecast.mean, 
                      name = 'I (forecast mean)', mode = 'lines',
                      line = list(color = col.I, width = mean.lwd, dash = 'dot')) %>% 
    plotly::add_ribbons(ymin = ~I.forecast.lower95, ymax = ~I.forecast.upper95,
                        name = 'I (forecast 95% confidence)', mode='lines',
                        line = list(color = col.I, width = ci.lwd),
                        fillcolor = col.I.ci) %>% 
    
    plotly::add_trace(y = ~E, 
                      name = 'E', mode = 'lines',
                      line = list(color = col.E, width = data.lwd)) %>%
    plotly::add_trace(y = ~E.forecast.mean, 
                      name = 'E (forecast mean)', mode = 'lines',
                      line = list(color = col.E, width = mean.lwd, dash = 'dot')) %>% 
    plotly::add_ribbons(ymin = ~E.forecast.lower95, ymax = ~I.forecast.upper95,
                        name = 'E (forecast 95% confidence)', mode='lines',
                        line = list(color = col.E, width = ci.lwd),
                        fillcolor = col.E.ci) %>% 
    
    plotly::add_trace(y = ~nowcast.mean, 
                      name = 'nowcast (mean)', mode = 'lines',
                      line = list(color = col.nowcast, width = data.lwd, dash = 'dot')) %>% 
    plotly::add_ribbons(ymin = ~nowcast.lower95, ymax = ~nowcast.upper95,
                        name = 'nowcast (95% confidence)', mode='lines',
                        line = list(color = col.nowcast, width = ci.lwd),
                        fillcolor = col.nowcast.ci
    )
  p_nowcast_logy <- p_nowcast %>% plotly::layout(yaxis = list(type = "log", range=c(-.25,5)))
  p_nowcast_logy
}




