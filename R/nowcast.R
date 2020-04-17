library(parallel)
library(tidyverse)
library(tibbletime)
library(padr)
library(tvReg)
library(forecast)
# source('R/package.R')
# source('R/deconcolve.R')
# source('R/simple_tdar.R')
# source('R/package_tv.R')

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

# set infities and NaN's to NA

na_not_finite <- function(x){
  x[!is.finite(x)] <- NA
  return(x)
}

# get sd of distribution

get_sd <- function(...){
  interval <- list(...)[[1]]

  if(is.numeric(interval)){
    if(length(interval) == 1){
      sd <- 0
    }else{
        sd <- sd(interval)
    }
    return(sd)
  }

  if(interval$dist == "exponential") {
    if(is.null(interval$mean)) {stop("exponential distribution requires a mean")}
    sd <- sqrt((1/interval$mean)^2)
    return(sd)
  }
  
  if(interval$dist == "gamma") {
    if(is.null(interval$mean) | is.null(interval$shape)) {
      stop("gamma distribution requires mean and shape")
    }
    scale <- interval$mean/interval$shape # gamma
    sd <- sqrt(interval$shape*scale^2) # gamma
    return(sd)
  }
  
  if(interval$dist == "skewnormal") {
    if(is.null(interval$mean) | is.null(interval$shape) | is.null(interval$scale)) {
      stop("skewnormal distribution requires mean, shape and scale to find standard deviation")
    }
    delta <- interval$shape / (sqrt(1 + interval$shape^2))
    variance <- interval$scale^2 * (1 - (2 * delta)/pi )
    sd <- sqrt(variance)
    return(sd)
  }

  if(interval$dist == "lognormal") {
    if(is.null(interval$sd)) {
      stop("Standard deviation must be supplied with lognormal distribution.")
    }
    return(interval$sd)
  }

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

# Smoother

# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, window=1, centered=TRUE) {
  
  if (centered) {
    before <- floor  ((window-1)/2)
    after  <- ceiling((window-1)/2)
  } else {
    before <- window-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

# Simulate from time varying distribution

## Simulate a line list from dates, counts (such as case reports), and interval distribution

backward_simulate_linelist <- function(dates, counts, interval) {
  
  input <- tibble(exit.date = dates, count = counts)

  if(is.numeric(interval)){
    interval <- interval[1]
    get_intervals <- function(x, interval) {
      rep(x, interval)
    }
  }else{
    if(interval$dist == "exponential") {
      get_intervals <- function(x, interval) {
        rexp(x, rate = 1/interval$mean)
      }
    }
    if(interval$dist == "gamma") {
      get_intervals <- function(x,interval) {
        rgamma(x, shape = interval$shape, scale = interval$mean/interval$shape)
      }
    }
    if(interval$dist == "skewnormal") {
      get_intervals <- function(x, interval) {
        sn::rsn(x, xi = interval$location, alpha = interval$shape, omega = interval$scale)
      }
    }
    if(interval$dist == "lognormal") {
      get_intervals <- function(x, interval) {
        if(is.null(interval$location)){
          interval$location <- log(interval$mean^2 / sqrt(interval$sd^2 + interval$mean^2))
          interval$shape <- sqrt(log(1 + (interval$sd^2 / interval$mean^2)))
        }
        rlnorm(x, interval$location, interval$shape)
      }
    }
  }
  
  # serial
  # intervals <- lapply(input$count, get_intervals, interval) %>% unlist()

  # parallel
  intervals <- mclapply(input$count, get_intervals, interval, mc.cores = detectCores()-1L) %>% unlist()
  
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
    if(interval$dist == "skewnormal") {
      get_intervals <- function(x, interval) {
        sn::rsn(x, xi = interval$location, alpha = interval$shape, omega = interval$scale)
      }
    }
    if(interval$dist == "lognormal") {
      get_intervals <- function(x, interval) {
        if(is.null(interval$location)){
          interval$location <- log(interval$mean^2 / sqrt(interval$sd^2 + interval$mean^2))
          interval$shape <- sqrt(log(1 + (interval$sd^2 / interval$mean^2)))
        }
        rlnorm(x, interval$location, interval$shape)
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

get_onset_curve <- function(dates, linelist, interval, next_intervals=NULL) {
  # dates <- range(dates)
  onset.curve <- linelist %>% 
    mutate(onset.date = as.Date(lubridate::floor_date(onset.date))) %>% 
    count(onset.date) %>% 
    padr::pad(interval = "day", start_val = range(dates)[1L], end_val = range(dates)[2L])
  
  # if(is.numeric(interval)){
  #   interval <- list(mean=interval[1], sd=0)
  # }else{
  #   if(interval$dist == "exponential") {
  #     interval$sd <- sqrt((1/interval$mean)^2) # exponential
  #   }
  #   if(interval$dist == "gamma") {
  #     interval$scale <- interval$mean/interval$shape # gamma
  #     interval$sd <- sqrt(interval$shape*interval$scale^2) # gamma
  #   }
  #   if(interval$dist == "skewnormal") {
  #     # if SD provided as a parameter of interval, do nothing
  #     # if sd not provided, 
  #     if(is.null(interval$sd)){
  #       delta <- interval$shape / (sqrt(1 + interval$shape^2))
  #       variance <- interval$scale^2 * (1 - (2 * delta)/pi )
  #       interval$sd <- sqrt(variance) 
  #     }
  #   }
  #   if(interval$dist == "lognormal") {
  #     # Do nothing. sd must be provided
  #   }
  # }
  
  interval$sd <- get_sd(interval)
  natail <- interval$mean+interval$sd*2
  if(!is.null(next_intervals)){
    for(i in 1:length(next_intervals)){
      next_intervals[[i]]$sd <- get_sd(next_intervals[[i]])
      natail <- natail + next_intervals[[i]]$mean+next_intervals[[i]]$sd*2
    }
  }
  
  # NA replace
  onset.curve$n <- replace_na(onset.curve$n, 0)
  
  # set last several values to NA
  onset.curve$n <- tail_na(onset.curve$n,round(natail))
  
  onset.curve <- onset.curve %>% rename(date = onset.date, value = n)

  onset.curve
}

## Get value of state at a single "date" from "linelist"
get_state <- function(date, linelist) {
  linelist[as.Date(lubridate::floor_date(linelist$onset.date)) <= date
           & date < linelist$exit.date,] %>% nrow
}

## Get state curve for range of dates from linelist
get_state_curve <- function(dates, linelist, interval, next_intervals=NULL){
  # serial
  # values <- lapply(X = dates, FUN = get_state, linelist) %>% unlist
  
  # parallel
  values <- mclapply(X = dates, FUN = get_state, linelist, mc.cores = detectCores()-1L) %>% unlist
  
  # if(is.numeric(interval)){
  #   interval <- list(mean=interval[1], sd=0)
  # }else{
  #   if(interval$dist == "exponential") {
  #     interval$sd <- sqrt((1/interval$mean)^2) # exponential
  #   }
  #   if(interval$dist == "gamma") {
  #     interval$scale <- interval$mean/interval$shape # gamma
  #     interval$sd <- sqrt(interval$shape*interval$scale^2) # gamma
  #   }
  #   if(interval$dist == "skewnormal") {
  #     # if SD provided as a parameter of interval, do nothing
  #     # if sd not provided, 
  #     if(missing(interval$sd)){
  #       delta <- interval$shape / (sqrt(1 + interval$shape^2))
  #       variance <- interval$scale^2 * (1 - (2 * delta)/pi )
  #       interval$sd <- sqrt(variance) 
  #     }
  #   }
  #   if(interval$dist == "lognormal") {
  #     # Do nothing. sd must be provided
  #   }
  # }
  
  sd <- get_sd(interval)
  natail <- interval$mean+sd*2
  if(!is.null(next_intervals)){
    for(i in 1:length(next_intervals)){
      natail <- natail + next_intervals[[i]]$mean+get_sd(next_intervals[[i]])*2
    }
  }
  
  # NA replace
  values <- values %>% replace_na(0)
  
  # set last several values to NA
  values <- tail_na(values, round(natail))
  
  state.curve <- tibble(
    date = dates, 
    value = values
    )
  state.curve
}

# forecast to present

tvar_forecast_to_present <- function(curve,lag=1,bw=NULL) {
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
                                  tkernel = "Gaussian",
                                  bw=NULL)
  ### TEMP PLOT
  # plot(log.curve.trimmed, type="l")
  # log.curve.model.tvar$fitted %>% lines(type="l", col="red")
  # title(main=paste('bw', log.curve.model.tvar$bw, 'p', log.curve.model.tvar$p))
  ### END TEMP PLOT
  
  log.curve.forecast <- forecast::forecast(log.curve.model.tvar$fitted,forecast.length)
  
  # antilog fit 
  # antilog forecast
  
  forecast <- tibble(fit = exp(log.curve.model.tvar$fitted) %>% 
                       pad_na(lag, front = TRUE) %>% 
                       pad_na(forecast.length),
                     mean = exp(log.curve.forecast$mean) %>% pmax(0) %>% 
                       pad_na(trimmed.length, front = TRUE),
                     upper80 = exp(log.curve.forecast$upper[,'80%']) %>% pmax(0) %>% 
                       pad_na(trimmed.length, front = TRUE),
                     lower80 = exp(log.curve.forecast$lower[,'80%']) %>% pmax(0) %>%
                       pad_na(trimmed.length, front = TRUE)
                     )
  forecast
}

# nowcast from case reports

nowcast_from_case_reports <- function(casereports, params, tvar.bandwidth=NULL, minimal=FALSE) {
  database <- casereports
  
  # Ascertainment
  if(is.null(params$q)){
    params$q <- 1
  }
  
  if(is.data.frame(params$q)){
    if(nrow(database) != nrow(params$q)) {
      stop("nowcast_from_case_reports: q must be a constant, a list, or a dataframe with the same number of rows as casereports.")
    }
    database$q <- params$q$mean
    database$q.upper <- params$q$upper
    database$q.lower <- params$q$lower
  }else{
    if(is.list(params$q)){
      database$q <- params$q$mean
      database$q.upper <- params$q$upper
      database$q.lower <- params$q$lower
    }else{
      database$q <- params$q[1]
      database$q.upper <- params$q[1]
      database$q.lower <- params$q[1]
      }
  }
  
  database$q[is.na(database$q)] <- 1
  database$q.upper[is.na(database$q.upper)] <- 1
  database$q.lower[is.na(database$q.lower)] <- 1

  # database$a <- params$a # no longer used
  # database$c <- params$c # no longer used
  
  
  # Compensate for ascertainment (q)
  database$cases_over_q <- database$cases / database$q
  
  # Not using upper and lower estimates of q
  # if(is.list(params$q)){
  #   database$cases_over_q.upper <- database$cases / database$q.lower
  #   database$cases_over_q.lower <- database$cases / database$q.upper
  # }
  # 
  
  # I linelist
  I.linelist <- backward_simulate_linelist(dates = database$Date,
                                           counts = database$cases_over_q,
                                           interval = params$effective.infectious.period)
  
  # Not using upper and lower estimates of q
  # if(is.list(params$q)) {
  #   # I linelist.upper
  #   I.linelist.upper <- backward_simulate_linelist(dates = database$Date,
  #                                            counts = database$cases_over_q.upper,
  #                                            interval = params$effective.infectious.period)
  #   
  #   # I linelist.lower
  #   I.linelist.lower <- backward_simulate_linelist(dates = database$Date,
  #                                                  counts = database$cases_over_q.lower,
  #                                                  interval = params$effective.infectious.period)
  # }
 
  
  # # I onset
  # I.onset <- get_onset_curve(dates = database$Date,
  #                            linelist = I.linelist,
  #                            interval = params$effective.infectious.period)
  # database$I.onset <- I.onset$n

  # # I (detected)
  # I_d <- get_state_curve(dates = database$Date,
  #                        linelist = I.linelist,
  #                        interval = params$effective.infectious.period)
  # database$I_d <- I_d$value
  # 
  # # I
  # database <- database %>% mutate(I = I_d/q,
  #                                 I.upper = I_d/q.lower,
  #                                 I.lower = I_d/q.upper)
  
  # I
  I <- get_state_curve(dates = database$Date,
                         linelist = I.linelist,
                         interval = params$effective.infectious.period)
  database$I <- I$value

  # Not using upper and lower estimates of q
  # # I.upper
  # if(is.list(params$q)){
  #   I.upper <- get_state_curve(dates = database$Date,
  #                              linelist = I.linelist.upper,
  #                              interval = params$effective.infectious.period)
  #   database$I.upper <- I.upper$value
  # }else{
  #   database$I.upper <-  database$I
  # }
  # 
  # # I.lower
  # if(is.list(params$q)){
  #   I.lower <- get_state_curve(dates = database$Date,
  #                              linelist = I.linelist.lower,
  #                              interval = params$effective.infectious.period)
  #   database$I.lower <- I.lower$value
  # }else{
  #   database$I.lower <-  database$I
  # }
  
  # E linelist
  E.linelist <- backward_propagate_linelist(I.linelist,
                                            interval = params$incubation.period)
  
  # Not using upper and lower estimates of q
  # if(is.list(params$q)){
  #   # E linelist.upper
  #   E.linelist.upper <- backward_propagate_linelist(I.linelist.upper,
  #                                             interval = params$incubation.period)
  # 
  #   # E linelist.lower
  #   E.linelist.lower <- backward_propagate_linelist(I.linelist.lower,
  #                                             interval = params$incubation.period)
  # }
  
  # E onset
  if(minimal == FALSE){
    E.onset <- get_onset_curve(dates = database$Date,
                               linelist = E.linelist,
                               interval = params$incubation.period,
                               next_intervals = list(params$effective.infectious.period))
    database$E.onset <- E.onset$value
  }
  
  # # E (detected)
  # E_d <- get_state_curve(dates = database$Date,
  #                        linelist = E.linelist,
  #                        interval = params$incubation.period)
  # database$E_d <- E_d$value
  # 
  # # E
  # database <- database %>% mutate(E = E_d/q,
  #                                 E.upper = E_d/q.lower,
  #                                 E.lower = E_d/q.upper)
  
  # E
  E <- get_state_curve(dates = database$Date,
                       linelist = E.linelist,
                       interval = params$incubation.period,
                       next_intervals = list(params$effective.infectious.period))
  database$E <- E$value
  
  # Not using upper and lower estimates of q
  # # E.upper
  # if(is.list(params$q)){
  #   E.upper <- get_state_curve(dates = database$Date,
  #                      linelist = E.linelist.upper,
  #                      interval = params$incubation.period,
  #                      next_intervals = list(params$effective.infectious.period))
  #   database$E.upper <- E.upper$value
  # }else{
  #   database$E.upper <- database$E
  # }
  # 
  # # E.lower
  # if(is.list(params$q)){
  #   E.lower <- get_state_curve(dates = database$Date,
  #                      linelist = E.linelist.lower,
  #                      interval = params$incubation.period,
  #                      next_intervals = list(params$effective.infectious.period))
  #   database$E.lower <- E.lower$value
  # }else{
  #   database$E.lower <- database$E
  # }
  
  # forecasting I
  I.forecast <- tvar_forecast_to_present(database$I, bw=tvar.bandwidth)
  database$I.fit <- I.forecast$fit
  database$I.forecast.mean <- I.forecast$mean
  database$I.forecast.upper80 <- I.forecast$upper80
  database$I.forecast.lower80 <- I.forecast$lower80
  
  # Not using upper and lower estimates of q
  # if(is.list(params$q)){
  #   I.upper.forecast <- tvar_forecast_to_present(database$I.upper)
  #   database$I.upper.fit <- I.upper.forecast$fit
  #   database$I.upper.forecast.mean <- I.upper.forecast$mean
  #   database$I.upper.forecast.lower80 <- I.upper.forecast$lower80
  #   database$I.upper.forecast.upper80 <- I.upper.forecast$upper80
  #   
  #   I.lower.forecast <- tvar_forecast_to_present(database$I.lower)
  #   database$I.lower.fit <- I.lower.forecast$fit
  #   database$I.lower.forecast.mean <- I.lower.forecast$mean
  #   database$I.lower.forecast.lower80 <- I.lower.forecast$lower80
  #   database$I.lower.forecast.upper80 <- I.lower.forecast$upper80
  # }else{
  #   database$I.upper.forecast.mean <- database$I.forecast.mean
  #   database$I.upper.forecast.lower80 <- database$I.forecast.lower80 
  #   database$I.upper.forecast.upper80 <- database$I.forecast.upper80 
  #   database$I.lower.forecast.mean <- database$I.forecast.mean
  #   database$I.lower.forecast.lower80 <- database$I.forecast.lower80 
  #   database$I.lower.forecast.upper80 <- database$I.forecast.upper80 
  # }
  
  
  # forecasting E
  E.forecast <- tvar_forecast_to_present(database$E, bw=tvar.bandwidth)
  database$E.fit <- E.forecast$fit
  database$E.forecast.mean <- E.forecast$mean
  database$E.forecast.upper80 <- E.forecast$upper80
  database$E.forecast.lower80 <- E.forecast$lower80
  
  # Not using upper and lower estimates of q
  # if(is.list(params$q)){
  #   E.upper.forecast <- tvar_forecast_to_present(database$E.upper)
  #   database$E.upper.fit <- E.upper.forecast$fit
  #   database$E.upper.forecast.mean <- E.upper.forecast$mean
  #   database$E.upper.forecast.lower80 <- E.upper.forecast$lower80
  #   database$E.upper.forecast.upper80 <- E.upper.forecast$upper80
  #   
  #   E.lower.forecast <- tvar_forecast_to_present(database$E.lower)
  #   database$E.lower.fit <- E.lower.forecast$fit
  #   database$E.lower.forecast.mean <- E.lower.forecast$mean
  #   database$E.lower.forecast.lower80 <- E.lower.forecast$lower80
  #   database$E.lower.forecast.upper80 <- E.lower.forecast$upper80
  # }else{
  #   database$E.upper.forecast.mean <- database$E.forecast.mean
  #   database$E.upper.forecast.lower80 <- database$E.forecast.lower80 
  #   database$E.upper.forecast.upper80 <- database$E.forecast.upper80 
  #   database$E.lower.forecast.mean <- database$E.forecast.mean
  #   database$E.lower.forecast.lower80 <- database$E.forecast.lower80 
  #   database$E.lower.forecast.upper80 <- database$E.forecast.upper80 
  # }

  # forecasting E.onset
  if(minimal == FALSE){
    E.onset.forecast <- tvar_forecast_to_present(database$E.onset, bw=tvar.bandwidth)
    database$E.onset.fit <- E.onset.forecast$fit
    database$E.onset.forecast.mean <- E.onset.forecast$mean
    database$E.onset.forecast.upper80 <- E.onset.forecast$upper80
    database$E.onset.forecast.lower80 <- E.onset.forecast$lower80
  }
  
  # nowcast
  
  # Not using upper and lower estimates of q
  # ## I composite bounds
  # 
  # database <- database %>% 
  #   mutate(I.nowcast.upper = rowSums(dplyr::select(., I.upper, I.upper.forecast.upper80), na.rm = TRUE),
  #          I.nowcast.lower = rowSums(dplyr::select(., I.lower, I.lower.forecast.lower80), na.rm = TRUE)
  #   )
  # 
  # ## E composite bounds
  # database <- database %>% 
  #   mutate(E.nowcast.upper = rowSums(dplyr::select(., E.upper, E.upper.forecast.upper80), na.rm = TRUE),
  #          E.nowcast.lower = rowSums(dplyr::select(., E.lower, E.lower.forecast.lower80), na.rm = TRUE)
  #   )  
  # 
  # ## nowcast and bounds
  # database <- database %>%
  #   mutate(nowcast.mean = rowSums(dplyr::select(., I, I.forecast.mean, E, E.forecast.mean), na.rm = TRUE),
  #          nowcast.upper = rowSums(dplyr::select(., I.upper, I.upper.forecast.upper80, E.upper, E.upper.forecast.upper80), na.rm = TRUE),
  #          nowcast.lower = rowSums(dplyr::select(., I.lower, I.lower.forecast.lower80, E.lower, E.lower.forecast.lower80), na.rm = TRUE)
  #   )

  # alternate nowcast without propagating uncertainty in q to forecasts.
  
  database <- database %>%
    mutate(nowcast.mean = rowSums(dplyr::select(., I, I.forecast.mean, E, E.forecast.mean), na.rm = TRUE),
           nowcast.upper = rowSums(dplyr::select(., I, I.forecast.upper80, E, E.forecast.upper80), na.rm = TRUE),
           nowcast.lower = rowSums(dplyr::select(., I, I.forecast.lower80, E, E.forecast.lower80), na.rm = TRUE)
    )
  
  # cumulative 
  if(minimal == FALSE){
    database <- database %>%
      mutate(cum.cases = cumsum(replace_na(cases,0)),
             cum.infections.mean = cumsum(
               rowSums(dplyr::select(., E.onset, E.onset.forecast.mean), na.rm = TRUE)
               ),
             cum.infections.upper80 = cumsum(
               rowSums(dplyr::select(., E.onset, E.onset.forecast.upper80), na.rm = TRUE)
               ),
             cum.infections.lower80 = cumsum(
               rowSums(dplyr::select(., E.onset, E.onset.forecast.lower80), na.rm = TRUE)
               )
             )
  }
  
  
  # R effective
  if(minimal == FALSE){
    
    database <- database %>%
      mutate(R_eff.mean = na_not_finite(params$effective.infectious.period$mean * 
                                        rowSums(dplyr::select(., E.onset, E.onset.forecast.mean), na.rm = TRUE) /
                                        rowSums(dplyr::select(., I, I.forecast.mean), na.rm = TRUE)
                                        ),
             R_eff.lower = na_not_finite(params$effective.infectious.period$mean *
                                         rowSums(dplyr::select(., E.onset, E.onset.forecast.lower80), na.rm = TRUE) / 
                                         rowSums(dplyr::select(., I, I.forecast.lower80), na.rm = TRUE)
                                        ),
             R_eff.upper = na_not_finite(params$effective.infectious.period$mean *
                                         rowSums(dplyr::select(., E.onset, E.onset.forecast.upper80), na.rm = TRUE) /
                                         rowSums(dplyr::select(., I, I.forecast.upper80), na.rm = TRUE)
                                        )
             )
  }
  
  return(database)
}

# nowcast from death reports using symtpom-onset-to-death

nowcast_from_deaths_with_onset_to_death <- function(deathreports, params, tvar.bandwidth=NULL) {
  database <- deathreports
  # # Ascertainment - not used for death reports
  # if(is.null(params$q)){
  #   params$q <- 1
  # }
  # if(is.data.frame(params$q)){
  #   if(nrow(database) != nrow(params$q)) {s
  #     stop("nowcast_from_case_reports: q must be a constant, a list, or a dataframe with the same number of rows as casereports.")
  #   }
  #   database$q <- params$q$mean
  #   database$q.upper <- params$q$upper
  #   database$q.lower <- params$q$lower
  # }else{
  #   if(is.list(params$q)){
  #     database$q <- params$q$mean
  #     database$q.upper <- params$q$upper
  #     database$q.lower <- params$q$lower
  #   }else{
  #     database$q <- params$q[1]
  #     database$q.upper <- params$q[1]
  #     database$q.lower <- params$q[1]
  #   }
  # }
  # database$a <- params$a # no longer used
  # database$c <- params$c # no longer used
  database$IFR <- params$IFR

  # Compensate for IFR
  database$deaths_over_IFR <- database$deaths / database$IFR
    
  # I linelist
  I_to_Death.linelist <- backward_simulate_linelist(dates = database$Date,
                                           counts = database$deaths_over_IFR,
                                           interval = params$onset.to.death.period)
  # # I onset
  # I.onset <- get_onset_curve(dates = database$Date,
  #                            linelist = I_to_Death.linelist,
  #                            interval = params$onset.to.death.period)
  # database$I.onset <- I.onset$n
  
  # I
  I <- get_state_curve(dates = database$Date,
                         linelist = I_to_Death.linelist,
                         interval = params$onset.to.death.period)
  database$I <- I$value
  
  # E linelist
  E_to_I.linelist <- backward_propagate_linelist(I_to_Death.linelist,
                                                 interval = params$incubation.period)
  
  # # E onset
  # E.onset <- get_onset_curve(dates = database$Date,
  #                            linelist = E_to_I.linelist,
  #                            interval = params$incubation.period,
  #                            next_intervals = list(params$onset.to.death.period))  
  # database$E.onset <- E.onset$n

  # E
  E <- get_state_curve(dates = database$Date,
                         linelist = E_to_I.linelist,
                         interval = params$incubation.period,
                         next_intervals = list(params$onset.to.death.period))
  database$E <- E$value
  
  # forecasting I
  I.forecast <- tvar_forecast_to_present(database$I, bw=tvar.bandwidth)
  database$I.fit <- I.forecast$fit
  database$I.forecast.mean <- I.forecast$mean
  database$I.forecast.lower80 <- I.forecast$lower80
  database$I.forecast.upper80 <- I.forecast$upper80
  
  # forecasting E
  E.forecast <- tvar_forecast_to_present(database$E, bw=tvar.bandwidth)
  database$E.fit <- E.forecast$fit
  database$E.forecast.mean <- E.forecast$mean
  database$E.forecast.lower80 <- E.forecast$lower80
  database$E.forecast.upper80 <- E.forecast$upper80
  
  # nowcast
  
  database <- database %>% 
    mutate(nowcast.mean = rowSums(dplyr::select(., I, I.forecast.mean, E, E.forecast.mean), na.rm = TRUE),
           nowcast.upper = rowSums(dplyr::select(., I, I.forecast.upper80, E, E.forecast.upper80), na.rm = TRUE),
           nowcast.lower = rowSums(dplyr::select(., I, I.forecast.lower80, E, E.forecast.lower80), na.rm = TRUE)
    ) 
  return(database)
}

# calculate ascertainment

get_ascertainment <- function(cases, deaths, params, window = 7) {
  params$q <- 1
  
  nowcast_from_cases <- nowcast_from_case_reports(cases, params, minimal=TRUE)
  nowcast_from_deaths <- nowcast_from_deaths_with_onset_to_death(deaths, params)
  
  maxspan <- max(nrow(nowcast_from_cases), nrow(nowcast_from_deaths))
  minspan <- min(nrow(nowcast_from_cases), nrow(nowcast_from_deaths))
  
  nowcast_from_cases <- nowcast_from_cases[1:minspan,]
  nowcast_from_deaths <- nowcast_from_deaths[1:minspan,]
  
  nowcasts <- tibble(Date = nowcast_from_cases$Date,
                     nowcast.from.cases = nowcast_from_cases$nowcast.mean,
                     nowcast.from.cases.upper = nowcast_from_cases$nowcast.upper,
                     nowcast.from.cases.lower = nowcast_from_cases$nowcast.lower,
                     nowcast.from.deaths = nowcast_from_deaths$nowcast.mean,
                     nowcast.from.deaths.upper = nowcast_from_deaths$nowcast.upper,
                     nowcast.from.deaths.lower = nowcast_from_deaths$nowcast.lower)
  
  ascertainment <- nowcasts %>%
    transmute(Date = Date,
              mean.raw = nowcast.from.cases / nowcast.from.deaths,
              upper.raw = nowcast.from.cases.upper / nowcast.from.deaths.upper,
              lower.raw = nowcast.from.cases.lower / nowcast.from.deaths.lower) 
  
  # centred box moving average
  ascertainment$mean <- movingAverage(ascertainment$mean.raw, window = window, centered=TRUE)
  ascertainment$mean[is.nan(ascertainment$mean)] <- NA
  
  ascertainment$upper <- movingAverage(ascertainment$upper.raw, window = window, centered=TRUE)
  ascertainment$upper[is.nan(ascertainment$upper)] <- NA
  
  ascertainment$lower <- movingAverage(ascertainment$lower.raw, window = window, centered=TRUE)
  ascertainment$lower[is.nan(ascertainment$lower.smooth)] <- NA
  
  lastvalue <- function(x) {tail(x[!is.na(x)],1)}
  
  if(nrow(ascertainment) < nrow(cases)){
    # mean.pad <- tail(ascertainment$mean,1)
    # upper.pad <- tail(ascertainment$upper,1)
    # lower.pad <- tail(ascertainment$lower,1)
    # mean.smooth.pad <- tail(ascertainment$mean.raw,1)
    # upper.smooth.pad <- tail(ascertainment$upper.raw,1)
    # lower.smooth.pad <- tail(ascertainment$lower.raw,1)
    ascertainment <- ascertainment %>%
      padr::pad(end_val = max(cases$Date)) 
  }
  
  ascertainment$mean[is.na(ascertainment$mean)] <- lastvalue(ascertainment$mean)
  ascertainment$upper[is.na(ascertainment$upper)] <- lastvalue(ascertainment$upper)
  ascertainment$lower[is.na(ascertainment$lower)] <- lastvalue(ascertainment$lower)
  
  # # smooth
  # smoothed <- ascertainment$mean %>% movingAverage(window = 7, centered=TRUE)
  # plot(ascertainment$mean, log='y', type = "l")
  # lines(smoothed, col="green", lwd=2)
  
  # replace zeros with NA
  
  ascertainment[ascertainment == 0] <- NA
  return(ascertainment)
}

## function to wrap text in an html span tag styled with a serif font


# plot nowcast from case reports

plot_nowcast_from_case_reports <- function(database, plotdeaths = TRUE, plotcumulative = TRUE, maxy = 10^7, logy = TRUE, legend = FALSE) {
  
  col.cases <- 'rgba(0, 0, 0, .35)'
  col.deaths <- 'rgba(0, 0, 0, 1.0)'
  col.I <- 'rgba(230, 7, 7, .75)'
  col.I.ci <- 'rgba(230, 7, 7, .15)'
  col.E <- 'rgba(7, 164, 181, 0.75)'
  col.E.ci <- 'rgba(7, 164, 181, 0.15)'
  col.nowcast <- 'rgba(7, 7, 230, 0.75)'
  col.nowcast.ci <- 'rgba(7, 7, 230, 0.15)'
  col.invisible <- 'rgba(7, 7, 230, .01)'
  
  ci.lwd <- .5
  mean.lwd <- 1
  data.lwd <- 2
  
  serif <- function(x) {
    htmltools::tags$span(x, style = htmltools::css(font.family = "serif"))
  }
  
  display <- function(x) {
    format(round(x), big.mark=",", trim = TRUE)
  }
  
  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  
  p_nowcast <- plotly::plot_ly(data = database, x = ~Date) %>% 
    plotly::add_trace(y = ~cases, type = 'bar',
                      name = 'New case notifications',
                      marker = list(color = col.cases),
                      legendgroup = 'group1',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(cases),serif("new case reports")) ) %>% 
    plotly::add_trace(y = ~I, 
                      name = 'Symptomatic cases', type = 'scatter', mode = 'lines',
                      line = list(color = col.I, width = data.lwd),
                      legendgroup = 'group2',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(I),serif("symptomatic")) ) %>%
    plotly::add_trace(y = ~I.forecast.mean, 
                      name = '(forecast average)', type = 'scatter', mode = 'lines',
                      line = list(color = col.I, width = mean.lwd, dash = 'dot'),
                      legendgroup = 'group2',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(I.forecast.mean),"symptomatic\n",
                                    "(range:", display(I.forecast.lower80), "to",
                                    display(I.forecast.upper80),")") ) %>% 
    plotly::add_ribbons(ymin = ~I.forecast.lower80, ymax = ~I.forecast.upper80,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.I, width = ci.lwd),
                        fillcolor = col.I.ci,
                        legendgroup = 'group2',
                        hoverinfo='none') %>% 
    
    plotly::add_trace(y = ~E, 
                      name = 'Latent cases', type = 'scatter', mode = 'lines',
                      line = list(color = col.E, width = data.lwd),
                      legendgroup = 'group3',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(E),serif("latent")) ) %>%
    plotly::add_trace(y = ~E.forecast.mean, 
                      name = '(forecast average)', type = 'scatter', mode = 'lines',
                      line = list(color = col.E, width = mean.lwd, dash = 'dot'),
                      legendgroup = 'group3',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(E.forecast.mean),"latent\n",
                                    "(range:", display(E.forecast.lower80), "to",
                                    display(E.forecast.upper80),")") ) %>% 
    plotly::add_ribbons(ymin = ~E.forecast.lower80, ymax = ~E.forecast.upper80,
                        name = '(prediction interval)', type = 'scatter', mode='lines',
                        line = list(color = col.E, width = ci.lwd),
                        fillcolor = col.E.ci,
                        legendgroup = 'group3',
                        hoverinfo='none') %>% 
    
    plotly::add_trace(y = ~nowcast.mean, 
                      name = 'Total unnotified cases', type = 'scatter', mode = 'lines',
                      line = list(color = col.nowcast, width = data.lwd, dash = 'dot'),
                      legendgroup = 'group4',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~if_else(nowcast.mean==nowcast.upper,
                                      paste(display(nowcast.mean),"total latent + symptomatic"),
                                      paste(display(nowcast.mean),"total latent + symptomatic\n",
                                            "(range:", display(nowcast.lower), "to",
                                            display(nowcast.upper),")")
                                      )
                      ) %>% 
    plotly::add_ribbons(ymin = ~nowcast.lower, ymax = ~nowcast.upper,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.nowcast, width = ci.lwd),
                        fillcolor = col.nowcast.ci,
                        legendgroup = 'group4',
                        hoverinfo='none')
  
    if(plotcumulative == TRUE){
      p_nowcast <- p_nowcast %>% 
       plotly::add_trace(y = maxy,
                      name = 'CUMULATIVE NUMBERS HELPER TRACE', type = 'scatter', mode = 'lines',
                      line = list(color = col.invisible),
                      visible = TRUE, showlegend = FALSE,
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste("Case notifications to date:", display(cum.cases),"\n",
                                    "Deaths to date:", display(cum.deaths),"\n",
                                    "Estimated infections to date:", display(cum.infections.mean),"\n",
                                    "(range:",display(cum.infections.lower80),"to",display(cum.infections.upper80),")"
                                    )
                      )
    }
  
  if(plotdeaths == TRUE) {
    p_nowcast <- p_nowcast %>% 
      plotly::add_trace(y = ~deaths, 
                        name = 'New death notifications', type = 'scatter', mode = 'lines',
                        line = list(color = col.deaths, width = data.lwd, shape = 'hvh'),
                        legendgroup = 'group1',
                        hoverinfo = "x+text",
                        hoverlabel = hoverlabel,
                        text = ~paste(display(deaths),serif("new deaths")) )
  }
  
  
  
  # cross reference state codes
  states <- read_csv("data/states.csv")
  state <- deparse(substitute(database))
  statename <- states %>% filter(code == state) %>% pull(name)

  p_nowcast <- p_nowcast %>% 
    plotly::layout(yaxis = list(type = ifelse(logy==TRUE,"log","linear"),
                                range = ifelse(logy==TRUE,c(-.25,log10(maxy)),c(0,maxy)),
                                title = "Number",
                                spikethickness = 0),
                   xaxis = list(
                     spikethickness = 1,
                     spikedash = "dot",
                     spikecolor = "black",
                     spikemode = "across+marker",
                     spikesnap = "cursor"
                   ),
                   hovermode = 'x',
                   hoverdistance = 1,
                   # legend = list()
                   showlegend = legend,
                   title = paste("Nowcast for",statename)
                   )
  if(plotcumulative == TRUE){
    p_nowcast <- p_nowcast %>% 
      plotly::layout(
        annotations = list(yref = 'paper', xref = 'paper',
                           y = 1, x = .1, align = "left",
                           text = paste("Cumulative case notifications through:",format(max(database$Date),"%B %d, %Y"),
                                        display(tail(database$cum.cases,1)),"\n",
                                        "Cumulative deaths through:",format(max(database$Date),"%B %d, %Y"),
                                        display(tail(database$cum.deaths,1)),"\n",
                                        "Estimated cumulative infections through:",format(max(database$Date),"%B %d, %Y"),
                                        display(tail(database$cum.infections.mean,1)),"\n",
                                        "(range of estimate:",
                                        display(tail(database$cum.infections.lower80,1)),"to",
                                        display(tail(database$cum.infections.upper80,1)),")"
                           ),
                           showarrow = FALSE,
                           xanchor = "left")
      )
  }else{
    p_nowcast <- p_nowcast %>% 
      plotly::layout(
        annotations = list(yref = 'paper', xref = 'paper',
                           y = 1, x = .1, align = "left",
                           text = paste("Cumulative case notifications through:",format(max(database$Date),"%B %d, %Y"),
                                        display(tail(database$cum.cases,1)),"\n",
                                        "Cumulative deaths through:",format(max(database$Date),"%B %d, %Y"),
                                        display(tail(database$cum.deaths,1))
                           ),
                           showarrow = FALSE,
                           xanchor = "left")
      )
  }
  p_nowcast
}

# plot nowcast from death reports

plot_nowcast_from_death_reports <- function(database, maxy = 10^7) {
  
  col.cases <- 'rgba(0, 0, 0, .35)'
  col.I <- 'rgba(230, 7, 7, .75)'
  col.I.ci <- 'rgba(230, 7, 7, 0.15)'
  col.E <- 'rgba(7, 164, 181, 0.75)'
  col.E.ci <- 'rgba(7, 164, 181, .15)'
  col.nowcast <- 'rgba(7, 7, 230, 0.75)'
  col.nowcast.ci <- 'rgba(7, 7, 230, 0.15)'
  
  ci.lwd <- .5
  mean.lwd <- 1
  data.lwd <- 3
  
  p_nowcast <- plotly::plot_ly(data = database, x = ~Date) %>% 
    plotly::add_trace(y = ~deaths, type = 'bar', 
                      name = 'New death notifications', 
                      marker = list(color = col.cases),
                      legendgroup = 'group1'
  ) %>% 
    plotly::add_trace(y = ~I, 
                      name = 'Symptomatic cases', mode = 'lines',
                      line = list(color = col.I, width = data.lwd),
                      legendgroup = 'group2') %>%
    plotly::add_trace(y = ~I.forecast.mean, 
                      name = '(forecast average)', mode = 'lines',
                      line = list(color = col.I, width = mean.lwd, dash = 'dot'),
                      legendgroup = 'group2') %>% 
    plotly::add_ribbons(ymin = ~I.forecast.lower80, ymax = ~I.forecast.upper80,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.I, width = ci.lwd),
                        fillcolor = col.I.ci,
                        legendgroup = 'group2') %>% 
    
    plotly::add_trace(y = ~E, 
                      name = 'Latent cases', mode = 'lines',
                      line = list(color = col.E, width = data.lwd),
                      legendgroup = 'group3') %>%
    plotly::add_trace(y = ~E.forecast.mean, 
                      name = '(forecast average)', mode = 'lines',
                      line = list(color = col.E, width = mean.lwd, dash = 'dot'),
                      legendgroup = 'group3') %>% 
    plotly::add_ribbons(ymin = ~E.forecast.lower80, ymax = ~E.forecast.upper80,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.E, width = ci.lwd),
                        fillcolor = col.E.ci,
                        legendgroup = 'group3') %>% 
    
    plotly::add_trace(y = ~nowcast.mean, 
                      name = 'Total unnotified cases', mode = 'lines',
                      line = list(color = col.nowcast, width = data.lwd, dash = 'dot'),
                      legendgroup = 'group4') %>% 
    plotly::add_ribbons(ymin = ~nowcast.lower, ymax = ~nowcast.upper,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.nowcast, width = ci.lwd),
                        fillcolor = col.nowcast.ci,
                        legendgroup = 'group4')
  
  p_nowcast_logy <- p_nowcast %>% plotly::layout(yaxis = list(type = "log", 
                                                              range=c(-.25,log10(maxy)),
                                                              title = "Number"),
                                                 # legend = list()
                                                 showlegend = FALSE
  )
  p_nowcast_logy
}

plot_ascertainment <- function(ascertainment){
  col.cases <- 'rgba(0, 0, 0, 1)'
  col.I <- 'rgba(230, 7, 7, .75)'
  col.I.ci <- 'rgba(230, 7, 7, .25)'
  col.E <- 'rgba(7, 164, 181, 0.75)'
  col.E.ci <- 'rgba(7, 164, 181, 0.0)'
  col.nowcast <- 'rgba(7, 7, 230, 0.75)'
  col.nowcast.ci <- 'rgba(7, 7, 230, 0.25)'
  col.other = 'rgb(164, 0, 181, .75)'         
  col.other.ci = 'rgb(164, 0, 181, .25)'
  
  ci.lwd <- .5
  mean.lwd <- 1
  data.lwd <- 2
  
  
  p_ascertainment <- plotly::plot_ly(data = ascertainment, x = ~Date , y = ~mean, type = 'scatter',
                                     name = 'Ascertainment %', mode = 'lines',
                                     line = list(color = col.other, width = data.lwd)
  ) %>%
    plotly::add_trace(y = ~lower,
                      name = 'upper 80% confidence', mode = 'lines',
                      line = list(color = col.other, width = ci.lwd, dash = 'dot')) %>%
    plotly::add_trace(y = ~upper,
                      name = 'lower 80% confidence', mode = 'lines',
                      line = list(color = col.other, width = ci.lwd)) %>%
    plotly::layout(yaxis = list(type = "log",tickformat = ".2%",title="Ascertainment %"))
  p_ascertainment
}


Plot_R_effective <- function(database, legend=TRUE) {
  col.cases <- 'rgba(0, 0, 0, 1)'
  col.I <- 'rgba(230, 7, 7, .75)'
  col.I.ci <- 'rgba(230, 7, 7, .25)'
  col.E <- 'rgba(7, 164, 181, 0.75)'
  col.E.ci <- 'rgba(7, 164, 181, 0.0)'
  col.nowcast <- 'rgba(7, 7, 230, 0.75)'
  col.nowcast.ci <- 'rgba(7, 7, 230, 0.25)'
  col.other = 'rgb(164, 0, 181, .75)'         
  col.other.ci = 'rgb(164, 0, 181, .25)'
  
  ci.lwd <- .5
  mean.lwd <- 1
  data.lwd <- 2
  
  serif <- function(x) {
    htmltools::tags$span(x, style = htmltools::css(font.family = "serif"))
  }
  
  display <- function(x) {
    format(round(x,2), big.mark=",", trim = TRUE)
  }
  
  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  
  p_Reff <- plotly::plot_ly(data = database, x = ~Date , y = ~R_eff.mean, type = 'scatter',
                                     name = 'R (mean)', mode = 'lines',
                                     line = list(color = col.other, width = data.lwd),
                            hoverinfo = "x+text",
                            hoverlabel = hoverlabel,
                            text = ~paste("mean R =", display(R_eff.mean),"\n",
                                          "(range:", display(R_eff.lower), "to",
                                          display(R_eff.upper),")")) %>%
    plotly::add_trace(y = ~R_eff.upper,
                      name = 'upper 80% confidence', mode = 'lines',
                      line = list(color = col.other, width = ci.lwd, dash = 'dot'),
                      hoverinfo = "none") %>%
    plotly::add_trace(y = ~R_eff.lower,
                      name = 'lower 80% confidence', mode = 'lines',
                      line = list(color = col.other, width = ci.lwd),
                      hoverinfo = "none") %>%
    plotly::layout(yaxis = list(# type = "log", 
                                # tickformat = ".2%", 
                                title="Effective Reproduction Number (R)")) %>% 
    plotly::layout(shapes = list(list(type = "line", 
                                      line = list(color = "black", width=2, dash = "dot"), 
                                      xref = "x", yref = "y",
                                      x0 = ~min(Date), x1 = ~max(Date), 
                                      y0 = 1, y1 = 1)),
                   annotations = list(yref = 'y', xref = 'paper',
                                      y = 1, x = .1, align = "left",
                                      text = "R = 1",
                                      showarrow = FALSE,
                                      xanchor = "left", yanchor = "bottom"),
                   yaxis = list(range=c(0,40), spikethickness = 0),
                   xaxis = list(spikethickness = 1,
                                spikedash = "dot",
                                spikecolor = "black",
                                spikemode = "across+marker",
                                spikesnap = "cursor"),
                   hovermode = 'x',
                   hoverdistance = 1,
                   # legend = list()
                   showlegend = legend
    )
  
  
  
  
  
  p_Reff
}
  
  


# subplots <- list(p_Reff, p_nowcast_logy)
# dash <- plotly::subplot(subplots, nrows = length(subplots), shareX = TRUE, titleY = TRUE) %>%
#   plotly::layout(
#     xaxis = list(
#       spikethickness = 1,
#       spikedash = "dot",
#       spikecolor = "black",
#       spikemode = "across+marker",
#       spikesnap = "cursor"
#     ),
#     yaxis = list(spikethickness = 0)
#   )
# dash

