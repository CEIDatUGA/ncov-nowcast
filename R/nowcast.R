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

# Helper Functions ----------------------------------------------------------------------------

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

# sum rows, keeping NA's
rowSums_na <- function(x){
  ifelse(apply(is.na(x),1,all), # if all columns of x == NA 
         NA,                    # set to NA
         rowSums(x,na.rm=TRUE)) # else sum columns, treating NA's as 0
}  

# set infities and NaN's to NA

na_not_finite <- function(x){
  x[!is.finite(x)] <- NA
  return(x)
}

# Smoother

# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, window=1, centered=TRUE) {
  
  if (centered) {
    before <- floor((window-1)/2)
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

# Get standard deviation of distribution ------------------------------------------------------

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

# Calculate confidence intervals across rows of a dataframe -----------------------------------

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




# Simulate linelists from time varying distributions ------------------------------------------

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
  ncores <- get0('nowcast_cores', envir = .GlobalEnv, mode = "any", inherits = TRUE, ifnotfound = 1)
  intervals <- mclapply(input$count, get_intervals, interval, mc.cores = ncores) %>% unlist()
  
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


# Get onset curve for range of dates from linelist --------------------------------------------

get_chunk_onset_curve <- function(dates, linelist, interval, next_intervals=NULL) {
  onset.curve <- linelist %>% 
    mutate(onset.date = lubridate::floor_date(onset.date,"days")) %>% 
    count(onset.date) %>% 
    padr::pad(interval = "day", start_val = range(dates)[1L], end_val = range(dates)[2L])
  
  # NA replace
  onset.curve$n <- replace_na(onset.curve$n, 0)
  # cat("After na replace, onset curve nrows = ", nrow(onset.curve), "\n")  # for debugging
  
  onset.curve <- onset.curve %>% rename(date = onset.date, value = n)

  return(onset.curve)
}

## Get onset curve for range of dates from linelist  # DEPRICATED

get_onset_curve <- function(dates, linelist, interval, next_intervals=NULL) {
  onset.curve <- linelist %>% 
    mutate(onset.date = lubridate::floor_date(onset.date, "days")) %>% 
    count(onset.date) %>% 
    padr::pad(interval = "day", start_val = range(dates)[1L], end_val = range(dates)[2L])
  # cat("onset curve nrows = ", nrow(onset.curve), "\n")  # for debugging
  
  interval$sd <- get_sd(interval)
  natail <- interval$mean+interval$sd*2
  if(!is.null(next_intervals)){
    for(i in 1:length(next_intervals)){
      next_intervals[[i]]$sd <- get_sd(next_intervals[[i]])
      natail <- natail + next_intervals[[i]]$mean+next_intervals[[i]]$sd*2
    }
  }
  # cat("After next intervals, onset curve nrows = ", nrow(onset.curve), "\n")  # for debugging
  
  # NA replace
  onset.curve$n <- replace_na(onset.curve$n, 0)
  # cat("After na replace, onset curve nrows = ", nrow(onset.curve), "\n")  # for debugging
  
  # set last several values to NA
  onset.curve$n <- tail_na(onset.curve$n,round(natail))
  # cat("After tail_na, onset curve nrows = ", nrow(onset.curve), "\n")  # for debugging
  
  onset.curve <- onset.curve %>% rename(date = onset.date, value = n)
  
  return(onset.curve)
}


# Get state curve for range of dates from linelist --------------------------------------------

## Get value of state at a single "date" from "linelist"
get_state <- function(date, linelist) {
  linelist[lubridate::floor_date(linelist$onset.date,"days") <= date
           & date < linelist$exit.date,] %>% nrow()
}

## Get state curve for range of dates from linelist
get_chunk_state_curve <- function(dates, linelist, interval, next_intervals=NULL){
  # serial
  # values <- lapply(X = dates, FUN = get_state, linelist) %>% unlist
  
  # parallel
  ncores <- get0('nowcast_cores', envir = .GlobalEnv, mode = "any", inherits = TRUE, ifnotfound = 1)
  values <- mclapply(X = dates, FUN = get_state, linelist, mc.cores = ncores) %>% unlist

  # NA replace
  values <- values %>% replace_na(0)
  
  chunk.state.curve <- tibble(
    date = dates, 
    value = values
  )
  chunk.state.curve
}

## Get state curve for range of dates from linelist  # DEPRICATED
get_state_curve <- function(dates, linelist, interval, next_intervals=NULL){
  # serial
  # values <- lapply(X = dates, FUN = get_state, linelist) %>% unlist
  
  # parallel
  ncores <- get0('nowcast_cores', envir = .GlobalEnv, mode = "any", inherits = TRUE, ifnotfound = 1)
  values <- mclapply(X = dates, FUN = get_state, linelist, mc.cores = ncores) %>% unlist
  # browser() # for dubugging
  # cat("n values = ", length(values), "\n")  # for debugging
  
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
  # cat("after na_tail, n values = ", length(values), "\n")  # for debugging
  
  state.curve <- tibble(
    date = dates, 
    value = values
    )
  state.curve
}



# Sum curve chunks ----------------------------------------------------------------------------

# Sum curve chunks and set tail to NA, depending on interval(s)
get_curve_from_chunks <- function(chunks, interval, next_intervals=NULL){
  
  values <- rowSums(chunks[,-1], na.rm = TRUE)
  
  sd <- get_sd(interval)
  natail <- interval$mean+sd*2
  if(!is.null(next_intervals)){
    for(i in 1:length(next_intervals)){
      natail <- natail + next_intervals[[i]]$mean+get_sd(next_intervals[[i]])*2
    }
  }
  
  # set last several values to NA
  values <- tail_na(values, round(natail))
  # cat("after na_tail, n values = ", length(values), "\n")  # for debugging
  
  curve <- tibble(
    date = chunks$date, 
    value = values
  )
  return(curve)
}


# Forecast to present -------------------------------------------------------------------------

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



# Nowcast from case reports -------------------------------------------------------------------


nowcast_from_case_reports <- function(casereports, 
                                      params, 
                                      tvar.bandwidth=NULL, 
                                      minimal=FALSE, # calculate E.onset and R effective
                                      samplesize=1.0, # proportion of cases to use in generating linelist
                                      chunksize=30, # rows per chunk.
                                      savechunks=FALSE  # save chunk dataframes for debugging
                                      ) {
  database <- casereports
  if(is.null(chunksize)) {
    chunksize <- nrow(database)
    }
  
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
  

  # Sample cases
  sample <- database$cases_over_q * samplesize
  
  # CHUNK PROCESSING
  
  nchunks <- ceiling(nrow(database)/chunksize)
  I <- E <- E.onset <- tibble(date=database$Date)  # init curves
  date.chunks <- split(database$Date, ceiling(seq_along(database$Date)/chunksize))
  sample.chunks <- split(sample, ceiling(seq_along(sample)/chunksize))
  
  for (i in 1:(nchunks)) {
    cat(paste("Processing chunk", i, "...\n")) # for debugging
    
    # offset = (i-1) * chunksize
    dates <- date.chunks[[i]]
    counts <- sample.chunks[[i]]
    
    if(sum(counts, na.rm = TRUE) > 0) {  
      # I linelist for sample of cases for chunk
      chunk.I.linelist <- backward_simulate_linelist(dates = dates,
                                               counts = counts,
                                               interval = params$effective.infectious.period)
      # E linelist for sample of cases for chunk
      chunk.E.linelist <- backward_propagate_linelist(chunk.I.linelist,
                                                      interval = params$incubation.period)
      
      # I for sample of cases for chunk
      chunk.I.date.min <- lubridate::floor_date(min(chunk.I.linelist$onset.date),"days")
      chunk.I.dates <- seq(chunk.I.date.min, max(database$Date), 1)
      chunk.I <- get_chunk_state_curve(dates = chunk.I.dates,
                                 linelist = chunk.I.linelist,
                                 interval = params$effective.infectious.period)
      names(chunk.I)[2] <- paste0("chunk",i)
      I <- full_join(I,chunk.I,by="date")
      
      # E for sample of cases for chunk
      chunk.E.date.min <- lubridate::floor_date(min(chunk.E.linelist$onset.date),"days")
      chunk.E.dates <- seq(chunk.E.date.min, max(database$Date), 1)
      chunk.E <- get_chunk_state_curve(dates = chunk.E.dates,
                           linelist = chunk.E.linelist,
                           interval = params$incubation.period,
                           next_intervals = list(params$effective.infectious.period))
      names(chunk.E)[2] <- paste0("chunk",i)
      E <- full_join(E,chunk.E,by="date")
      
      # E onset for sample of cases for chunk
      if(minimal == FALSE){
        # cat("Getting Exposure onset curve... \n")
        chunk.E.onset <- get_chunk_onset_curve(dates = chunk.E.dates,
                                   linelist = chunk.E.linelist,
                                   interval = params$incubation.period,
                                   next_intervals = list(params$effective.infectious.period))
        names(chunk.E.onset)[2] <- paste0("chunk",i)
        E.onset <- full_join(E.onset,chunk.E.onset,by="date")
      }else{
        E.onset[paste0("chunk",i)] <- NA
      }
    }else{
      I[paste0("chunk",i)] <- NA
      E[paste0("chunk",i)] <- NA
      E.onset[paste0("chunk",i)] <- NA
    }
  }  # end for loop
  
  # Sum chunks
  I$value <- get_curve_from_chunks(I, 
                                   interval = params$effective.infectious.period
                                   )$value
  E$value <- get_curve_from_chunks(E, 
                                   interval = params$incubation.period,
                                   next_intervals = list(params$effective.infectious.period)
                                   )$value
  E.onset$value <- get_curve_from_chunks(E.onset, 
                                         interval = params$incubation.period,
                                         next_intervals = list(params$effective.infectious.period)
                                         )$value

  # Rescale curves and store in database
  database$I <- I$value / samplesize
  database$E <- E$value / samplesize
  if(minimal==FALSE) {
    database$E.onset <- E.onset$value / samplesize
  }else{
    database$E.onset <- NA
  }
  
  
  # END CHUNK PROCESSING
  
  
  # # NON-CHUNK PROCESSING
  # 
  # # I linelist for sample of cases
  # I.linelist <- backward_simulate_linelist(dates = database$Date,
  #                                          counts = sample,
  #                                          interval = params$effective.infectious.period)
  # cat("Done simulating I linelist.\n")  # for debugging
  # # print(head(I.linelist)) # for debugging
  # # print(tail(I.linelist)) # for debugging
  # # Not using upper and lower estimates of q
  # # if(is.list(params$q)) {
  # #   # I linelist.upper
  # #   I.linelist.upper <- backward_simulate_linelist(dates = database$Date,
  # #                                            counts = database$cases_over_q.upper,
  # #                                            interval = params$effective.infectious.period)
  # #
  # #   # I linelist.lower
  # #   I.linelist.lower <- backward_simulate_linelist(dates = database$Date,
  # #                                                  counts = database$cases_over_q.lower,
  # #                                                  interval = params$effective.infectious.period)
  # # }
  # 
  # # # I onset
  # # I.onset <- get_onset_curve(dates = database$Date,
  # #                            linelist = I.linelist,
  # #                            interval = params$effective.infectious.period)
  # # database$I.onset <- I.onset$n
  # 
  # # # I (detected)
  # # I_d <- get_state_curve(dates = database$Date,
  # #                        linelist = I.linelist,
  # #                        interval = params$effective.infectious.period)
  # # database$I_d <- I_d$value
  # #
  # # # I
  # # database <- database %>% mutate(I = I_d/q,
  # #                                 I.upper = I_d/q.lower,
  # #                                 I.lower = I_d/q.upper)
  # 
  # # I
  # cat("Getting state curve for I...\n") # for debugging
  # I <- get_state_curve(dates = database$Date,
  #                        linelist = I.linelist,
  #                        interval = params$effective.infectious.period)
  # cat("Done getting state curve for I.\n")  # for debugging
  # # print(head(I)) # for debugging
  # # print(tail(I)) # for debugging
  # # cat("I curve date range:", range(I$date), "\n") # for debugging
  # # cat("I curve nrows:", nrow(I), "\n") # for debugging
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # database$I <- I$value / samplesize
  # 
  # # Not using upper and lower estimates of q
  # # # I.upper
  # # if(is.list(params$q)){
  # #   I.upper <- get_state_curve(dates = database$Date,
  # #                              linelist = I.linelist.upper,
  # #                              interval = params$effective.infectious.period)
  # #   database$I.upper <- I.upper$value
  # # }else{
  # #   database$I.upper <-  database$I
  # # }
  # #
  # # # I.lower
  # # if(is.list(params$q)){
  # #   I.lower <- get_state_curve(dates = database$Date,
  # #                              linelist = I.linelist.lower,
  # #                              interval = params$effective.infectious.period)
  # #   database$I.lower <- I.lower$value
  # # }else{
  # #   database$I.lower <-  database$I
  # # }
  # 
  # # E linelist for sample of cases
  # cat("Back propagating linelist to Exposures...\n") # for debugging
  # E.linelist <- backward_propagate_linelist(I.linelist,
  #                                           interval = params$incubation.period)
  # cat("Done back propagating linelist.\n")  # for debugging
  # # print(head(E.linelist)) # for debugging
  # # print(tail(E.linelist)) # for debugging
  # 
  # # Not using upper and lower estimates of q
  # # if(is.list(params$q)){
  # #   # E linelist.upper
  # #   E.linelist.upper <- backward_propagate_linelist(I.linelist.upper,
  # #                                             interval = params$incubation.period)
  # #
  # #   # E linelist.lower
  # #   E.linelist.lower <- backward_propagate_linelist(I.linelist.lower,
  # #                                             interval = params$incubation.period)
  # # }
  # 
  # # E onset
  # if(minimal == FALSE){
  #   cat("Getting Exposure onset curve... \n")
  #   E.onset <- get_onset_curve(dates = database$Date,
  #                              linelist = E.linelist,
  #                              interval = params$incubation.period,
  #                              next_intervals = list(params$effective.infectious.period))
  # 
  #   database$E.onset <- E.onset$value / samplesize
  #   cat("Done getting Exposure onset curve. \n")
  #   # print(head(E.onset)) # for debugging
  #   # print(tail(E.onset)) # for debugging
  #   # cat("E onset curve date range:", range(E.onset$date), "\n")  # for debugging
  #   # cat("E onset curve nrows:", nrow(E.onset), "\n")  # for debugging
  # }
  # 
  # # # E (detected)
  # # E_d <- get_state_curve(dates = database$Date,
  # #                        linelist = E.linelist,
  # #                        interval = params$incubation.period)
  # # database$E_d <- E_d$value
  # #
  # # # E
  # # database <- database %>% mutate(E = E_d/q,
  # #                                 E.upper = E_d/q.lower,
  # #                                 E.lower = E_d/q.upper)
  # 
  # # E
  # cat("Getting Exposure state curve... \n")
  # E <- get_state_curve(dates = database$Date,
  #                      linelist = E.linelist,
  #                      interval = params$incubation.period,
  #                      next_intervals = list(params$effective.infectious.period))
  # 
  # database$E <- E$value / samplesize
  # cat("Done getting Exposure state curve. \n")
  # # print(head(E)) # for debugging
  # # print(tail(E)) # for debugging
  # # cat("E state curve date range:", range(E$date), "\n")
  # # cat("E state curve nrows:", nrow(E), "\n")
  # 
  # # Not using upper and lower estimates of q
  # # # E.upper
  # # if(is.list(params$q)){
  # #   E.upper <- get_state_curve(dates = database$Date,
  # #                      linelist = E.linelist.upper,
  # #                      interval = params$incubation.period,
  # #                      next_intervals = list(params$effective.infectious.period))
  # #   database$E.upper <- E.upper$value
  # # }else{
  # #   database$E.upper <- database$E
  # # }
  # #
  # # # E.lower
  # # if(is.list(params$q)){
  # #   E.lower <- get_state_curve(dates = database$Date,
  # #                      linelist = E.linelist.lower,
  # #                      interval = params$incubation.period,
  # #                      next_intervals = list(params$effective.infectious.period))
  # #   database$E.lower <- E.lower$value
  # # }else{
  # #   database$E.lower <- database$E
  # # }
  # 
  # 
  # # END NON-CHUNK PROCESSING
  
  
  # forecasting I
  cat("Forecasting I state curve... \n")
  I.forecast <- tvar_forecast_to_present(database$I, bw=tvar.bandwidth)
  database$I.fit <- I.forecast$fit
  database$I.forecast.mean <- I.forecast$mean
  database$I.forecast.upper80 <- I.forecast$upper80
  database$I.forecast.lower80 <- I.forecast$lower80
  cat("Done forecasting I state curve. \n")
  
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
  cat("Forecasting E state curve... \n")
  E.forecast <- tvar_forecast_to_present(database$E, bw=tvar.bandwidth)
  database$E.fit <- E.forecast$fit
  database$E.forecast.mean <- E.forecast$mean
  database$E.forecast.upper80 <- E.forecast$upper80
  database$E.forecast.lower80 <- E.forecast$lower80
  cat("Done forecasting E state curve... \n")
  
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
    cat("Forecasting E onset curve... \n")
    E.onset.forecast <- tvar_forecast_to_present(database$E.onset, bw=tvar.bandwidth)
    database$E.onset.fit <- E.onset.forecast$fit
    database$E.onset.forecast.mean <- E.onset.forecast$mean
    database$E.onset.forecast.upper80 <- E.onset.forecast$upper80
    database$E.onset.forecast.lower80 <- E.onset.forecast$lower80
    cat("Done forecasting E onset curve... \n")
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
  
  ### START HERE
  # get R_effective
  # burnin <- get_sd(interval)*2
  # get_R_eff <- function(dates,E.onset.curve,I.state,gamma) {
  #   R_effective <- E.onset.curve/(I.state*gamma)
  # }
  
  
  # get_R_eff <- function(database, burnin = 0) {
  #   # dates <- range(dates)
  #   onset.curve <- linelist %>%
  #     mutate(onset.date = lubridate::floor_date(onset.date,"days")) %>%
  #     count(onset.date) %>%
  #     padr::pad(interval = "day", start_val = range(dates)[1L], end_val = range(dates)[2L])
  # 
  #   interval$sd <- get_sd(interval)
  #   natail <- interval$mean+interval$sd*2
  #   if(!is.null(next_intervals)){
  #     for(i in 1:length(next_intervals)){
  #       next_intervals[[i]]$sd <- get_sd(next_intervals[[i]])
  #       natail <- natail + next_intervals[[i]]$mean+next_intervals[[i]]$sd*2
  #     }
  #   }
  # 
  #   # NA replace
  #   onset.curve$n <- replace_na(onset.curve$n, 0)
  # 
  #   # set last several values to NA
  #   onset.curve$n <- tail_na(onset.curve$n,round(natail))
  # 
  #   onset.curve <- onset.curve %>% rename(date = onset.date, value = n)
  # 
  #   onset.curve
  # }
  ### END HERE
  
  # Save
  if(savechunks==TRUE){
    return(list(database = database, chunks = list(I = I, E = E, E.onset = E.onset))) # for debugging.
  }else{
    return(database)
  }  
  
}


# Nowcast from death reports using symtpom-onset-to-death -------------------------------------

nowcast_from_deaths_with_onset_to_death <- function(deathreports, 
                                                    params, 
                                                    tvar.bandwidth=NULL,
                                                    minimal=FALSE, # calculate E.onset and R effective
                                                    samplesize=1.0, # proportion of cases to use in generating linelist
                                                    chunksize=30, # rows per chunk.
                                                    savechunks=FALSE  # save chunk dataframes for debugging
                                                    ) {
  database <- deathreports
  if(is.null(chunksize)) {
    chunksize <- nrow(database)
  }
  # Infection Fatality Rate
  database$IFR <- params$IFR

  # Compensate for IFR
  database$deaths_over_IFR <- database$deaths / database$IFR
  
  # Sample deaths
  sample <- database$deaths_over_IFR * samplesize
  
  
  # CHUNK PROCESSING*******************************************
  
  nchunks <- ceiling(nrow(database)/chunksize)
  I <- E <- E.onset <- tibble(date=database$Date)  # init curves
  date.chunks <- split(database$Date, ceiling(seq_along(database$Date)/chunksize))
  sample.chunks <- split(sample, ceiling(seq_along(sample)/chunksize))
  
  for (i in 1:(nchunks)) {
    cat(paste("Processing chunk", i, "...\n")) # for debugging
    
    # offset = (i-1) * chunksize
    dates <- date.chunks[[i]]
    counts <- sample.chunks[[i]]
    
    if(sum(counts, na.rm = TRUE) > 0) {  
      # I to death linelist for sample of cases for chunk
      chunk.I.linelist <- backward_simulate_linelist(dates = dates,
                                                     counts = counts,
                                                     interval = params$onset.to.death.period)
      
      # E linelist for sample of cases for chunk
      chunk.E.linelist <- backward_propagate_linelist(chunk.I.linelist,
                                                      interval = params$incubation.period)
      
      # I for sample of deaths for chunk
      chunk.I.date.min <- lubridate::floor_date(min(chunk.I.linelist$onset.date),"days")
      chunk.I.dates <- seq(chunk.I.date.min, max(database$Date), 1)
      chunk.I <- get_chunk_state_curve(dates = chunk.I.dates,
                                       linelist = chunk.I.linelist,
                                       interval = params$onset.to.death.period)
      names(chunk.I)[2] <- paste0("chunk",i)
      I <- full_join(I,chunk.I,by="date")
      
      # E for sample of deaths for chunk
      chunk.E.date.min <- lubridate::floor_date(min(chunk.E.linelist$onset.date),"days")
      chunk.E.dates <- seq(chunk.E.date.min, max(database$Date), 1)
      chunk.E <- get_chunk_state_curve(dates = chunk.E.dates,
                                       linelist = chunk.E.linelist,
                                       interval = params$incubation.period,
                                       next_intervals = list(params$onset.to.death.period))
      names(chunk.E)[2] <- paste0("chunk",i)
      E <- full_join(E,chunk.E,by="date")
      
      # E onset for sample of cases for chunk
      if(minimal == FALSE){
        # cat("Getting Exposure onset curve... \n")
        chunk.E.onset <- get_chunk_onset_curve(dates = chunk.E.dates,
                                               linelist = chunk.E.linelist,
                                               interval = params$incubation.period,
                                               next_intervals = list(params$onset.to.death.period))
        names(chunk.E.onset)[2] <- paste0("chunk",i)
        E.onset <- full_join(E.onset,chunk.E.onset,by="date")
      }else{
        E.onset[paste0("chunk",i)] <- NA
      }
    }else{
      I[paste0("chunk",i)] <- NA
      E[paste0("chunk",i)] <- NA
      E.onset[paste0("chunk",i)] <- NA
    }
  }  # end for loop
  
  # Sum chunks
  I$value <- get_curve_from_chunks(I, 
                                   interval = params$onset.to.death.period
                                   )$value
  E$value <- get_curve_from_chunks(E, 
                                   interval = params$incubation.period,
                                   next_intervals = list(params$effective.infectious.period)
                                   )$value
  E.onset$value <- get_curve_from_chunks(E.onset, 
                                         interval = params$incubation.period,
                                         next_intervals = list(params$onset.to.death.period)
                                         )$value
  
  # Rescale curves and store in database
  database$I <- I$value / samplesize
  database$E <- E$value / samplesize
  if(minimal==FALSE) {
    database$E.onset <- E.onset$value / samplesize
  }else{
    database$E.onset <- NA
  }  
  
  
  # END CHUNK PROCESSING***************************************

  
  # # NON-CHUNK PROCESSING
  # 
  # # I linelist
  # I_to_Death.linelist <- backward_simulate_linelist(dates = database$Date,
  #                                          counts = database$deaths_over_IFR,
  #                                          interval = params$onset.to.death.period)
  # # # I onset
  # # I.onset <- get_onset_curve(dates = database$Date,
  # #                            linelist = I_to_Death.linelist,
  # #                            interval = params$onset.to.death.period)
  # # database$I.onset <- I.onset$n
  # 
  # # I
  # I <- get_state_curve(dates = database$Date,
  #                        linelist = I_to_Death.linelist,
  #                        interval = params$onset.to.death.period)
  # database$I <- I$value
  # 
  # # E linelist
  # E_to_I.linelist <- backward_propagate_linelist(I_to_Death.linelist,
  #                                                interval = params$incubation.period)
  # 
  # # # E onset
  # # E.onset <- get_onset_curve(dates = database$Date,
  # #                            linelist = E_to_I.linelist,
  # #                            interval = params$incubation.period,
  # #                            next_intervals = list(params$onset.to.death.period))  
  # # database$E.onset <- E.onset$n
  # 
  # # E
  # E <- get_state_curve(dates = database$Date,
  #                        linelist = E_to_I.linelist,
  #                        interval = params$incubation.period,
  #                        next_intervals = list(params$onset.to.death.period))
  # database$E <- E$value
  # 
  # # END NON-CHUNK PROCESSING
  
  # forecasting I
  cat("Forecasting I state curve... \n")
  I.forecast <- tvar_forecast_to_present(database$I, bw=tvar.bandwidth)
  database$I.fit <- I.forecast$fit
  database$I.forecast.mean <- I.forecast$mean
  database$I.forecast.lower80 <- I.forecast$lower80
  database$I.forecast.upper80 <- I.forecast$upper80
  cat("Done forecasting I state curve. \n")
  
  # forecasting E
  cat("Forecasting E state curve... \n")
  E.forecast <- tvar_forecast_to_present(database$E, bw=tvar.bandwidth)
  database$E.fit <- E.forecast$fit
  database$E.forecast.mean <- E.forecast$mean
  database$E.forecast.lower80 <- E.forecast$lower80
  database$E.forecast.upper80 <- E.forecast$upper80
  cat("Done forecasting E state curve. \n")
  
  # forecasting E.onset
  if(minimal == FALSE){
    cat("Forecasting E onset curve... \n")
    E.onset.forecast <- tvar_forecast_to_present(database$E.onset, bw=tvar.bandwidth)
    database$E.onset.fit <- E.onset.forecast$fit
    database$E.onset.forecast.mean <- E.onset.forecast$mean
    database$E.onset.forecast.upper80 <- E.onset.forecast$upper80
    database$E.onset.forecast.lower80 <- E.onset.forecast$lower80
    cat("Done forecasting E onset curve... \n")
  }
  
  # nowcast
  
  database <- database %>% 
    mutate(nowcast.mean = rowSums(dplyr::select(., I, I.forecast.mean, E, E.forecast.mean), na.rm = TRUE),
           nowcast.upper = rowSums(dplyr::select(., I, I.forecast.upper80, E, E.forecast.upper80), na.rm = TRUE),
           nowcast.lower = rowSums(dplyr::select(., I, I.forecast.lower80, E, E.forecast.lower80), na.rm = TRUE)
    ) 
  
  # cumulative 
  if(minimal == FALSE){
    database <- database %>%
      mutate(cum.deaths = cumsum(replace_na(deaths,0)),
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

  # Save
  if(savechunks==TRUE){
    return(list(database = database, chunks = list(I = I, E = E, E.onset = E.onset))) # for debugging.
  }else{
    return(database)
  }  
  
}


# Calculate ascertainment ---------------------------------------------------------------------

get_ascertainment <- function(cases, 
                              deaths, 
                              params, 
                              window = 7, 
                              tvar.bandwidth=NULL, 
                              samplesize = 1.0, 
                              chunksize = 30) {
  params$q <- 1
  
  nowcast_from_cases <- nowcast_from_case_reports(cases, 
                                                  params, 
                                                  minimal=TRUE, 
                                                  tvar.bandwidth=tvar.bandwidth, 
                                                  samplesize = samplesize,
                                                  chunksize = chunksize)
  nowcast_from_deaths <- nowcast_from_deaths_with_onset_to_death(deaths, 
                                                                 params, 
                                                                 minimal=TRUE, 
                                                                 tvar.bandwidth=tvar.bandwidth, 
                                                                 samplesize = samplesize,
                                                                 chunksize = chunksize)
  
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


# Plotting helper functions -------------------------------------------------------------------

## function to wrap text in an html span tag styled with a serif font

serif <- function(x) {
  htmltools::tags$span(x, style = htmltools::css(font.family = "serif"))
}

## format numbers for display

display <- function(x) {
  format(round(x), big.mark=",", trim = TRUE)
}


# Visual Theme -------------------------------------------------------------------

nowcast_theme <- function (name='nowcast_theme_default', 
                           col.cases = 'rgba(0, 0, 0, .35)',
                           col.deaths = 'rgba(64, 0, 0, .5)',
                           col.I = 'rgba(230, 7, 7, .75)',
                           col.I.ci = 'rgba(230, 7, 7, .15)',
                           col.E = 'rgba(7, 164, 181, 0.75)',
                           col.E.ci = 'rgba(7, 164, 181, 0.15)',
                           col.nowcast = 'rgba(7, 7, 230, 0.75)',
                           col.nowcast.ci = 'rgba(7, 7, 230, 0.15)',
                           col.invisible = 'rgba(7, 7, 230, .01)',
                           col.R = 'rgb(164, 0, 181, .75)',
                           col.R.ci = 'rgb(164, 0, 181, .25)',
                           col.other = 'rgb(0, 0, 0, .75)',
                           col.other.ci = 'rgb(0, 0, 0, .25)',
                           lwd.ci = .5,
                           lwd.mean = 1,
                           lwd.data = 2){
  
  # translate colors
  formatcolor <- function(x){
    if(is.numeric(x)) {
      if(length(x) != 4) stop ("Color must be specificed as a 4 element numeric vector (r,g,b,a) 
                               with max (255,255,255,1), or as a string. 
                               See https://plotly-r.com/working-with-colors.html")
      y = paste0('rgba(',paste(x,collapse = ','),')')
    }else{
      y = x
    }
    return(y)
  }
  
  theme <- list(
    col.cases = formatcolor(col.cases),
    col.deaths = formatcolor(col.deaths),
    col.I = formatcolor(col.I),
    col.I.ci = formatcolor(col.I.ci),
    col.E = formatcolor(col.E),
    col.E.ci = formatcolor(col.E.ci),
    col.nowcast = formatcolor(col.nowcast),
    col.nowcast.ci = formatcolor(col.nowcast.ci),
    col.invisible = formatcolor(col.invisible),
    col.R = formatcolor(col.R),
    col.R.ci = formatcolor(col.R.ci),
    col.other = formatcolor(col.other),
    col.other.ci = formatcolor(col.other.ci),
    lwd.ci = lwd.ci,
    lwd.mean = lwd.mean,
    lwd.data = lwd.data
  )
  
  assign(name,theme,envir = .GlobalEnv)
}

# Plot nowcast from case reports --------------------------------------------------------------

plot_nowcast_from_case_reports <- function(database, 
                                           plotcases=TRUE, 
                                           plotdeaths = TRUE, 
                                           plotcumulative = TRUE, 
                                           annotations = FALSE, 
                                           maxy = 10^7, 
                                           logy = TRUE, 
                                           legend = FALSE,
                                           theme = 'default') {

  # VISUAL VARIABLES ***************************
  if(is.character(theme)) {
    if(theme=='default') {
      nctheme <- get0('nowcast_theme_default', envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())
    }else{
      nctheme <- get0(theme, envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())    
    }
  }else{
    if(is.list(theme)) { 
      nctheme <- nowcast_theme()
      nctheme[[names(theme)[i]]] <- theme[[i]]
      }
  }
  for(i in seq_along(nctheme)) { 
    assign(names(nctheme)[i], nctheme[[i]]) 
  }

  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  # END VISUAL VARIABLES ***********************
  
  # Plot
  p_nowcast <- plotly::plot_ly(data = database, x = ~Date) %>% 
    plotly::add_trace(y = ~nowcast.mean, 
                      name = 'Total unnotified cases', type = 'scatter', mode = 'lines',
                      line = list(color = col.nowcast, width = lwd.data, dash = 'dot'),
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
                        line = list(color = col.nowcast, width = lwd.ci),
                        fillcolor = col.nowcast.ci,
                        legendgroup = 'group4',
                        hoverinfo='none') %>% 
  
    plotly::add_trace(y = ~E, 
                      name = 'Latent cases', type = 'scatter', mode = 'lines',
                      line = list(color = col.E, width = lwd.data),
                      legendgroup = 'group3',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(E),serif("latent")) ) %>%
    plotly::add_trace(y = ~E.forecast.mean, 
                      name = '(forecast average)', type = 'scatter', mode = 'lines',
                      line = list(color = col.E, width = lwd.mean, dash = 'dot'),
                      legendgroup = 'group3',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(E.forecast.mean),"latent\n",
                                    "(range:", display(E.forecast.lower80), "to",
                                    display(E.forecast.upper80),")") ) %>% 
    plotly::add_ribbons(ymin = ~E.forecast.lower80, ymax = ~E.forecast.upper80,
                        name = '(prediction interval)', type = 'scatter', mode='lines',
                        line = list(color = col.E, width = lwd.ci),
                        fillcolor = col.E.ci,
                        legendgroup = 'group3',
                        hoverinfo='none') %>% 
    
    plotly::add_trace(y = ~I, 
                      name = 'Symptomatic cases', type = 'scatter', mode = 'lines',
                      line = list(color = col.I, width = lwd.data),
                      legendgroup = 'group2',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(I),serif("symptomatic")) ) %>%
    plotly::add_trace(y = ~I.forecast.mean, 
                      name = '(forecast average)', type = 'scatter', mode = 'lines',
                      line = list(color = col.I, width = lwd.mean, dash = 'dot'),
                      legendgroup = 'group2',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(I.forecast.mean),"symptomatic\n",
                                    "(range:", display(I.forecast.lower80), "to",
                                    display(I.forecast.upper80),")") ) %>% 
    plotly::add_ribbons(ymin = ~I.forecast.lower80, ymax = ~I.forecast.upper80,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.I, width = lwd.ci),
                        fillcolor = col.I.ci,
                        legendgroup = 'group2',
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
  
  if(plotcases == TRUE) {
    p_nowcast <- p_nowcast %>% 
      plotly::add_trace(y = ~cases, type = 'bar',
                        name = 'New case notifications', 
                        marker = list(color = col.cases),
                        legendgroup = 'group1',
                        hoverinfo = "x+text",
                        hoverlabel = hoverlabel,
                        text = ~paste(display(cases),serif("new case reports")) )
  }
  
  
  if(plotdeaths == TRUE) {
    p_nowcast <- p_nowcast %>% 
      plotly::add_trace(y = ~deaths, 
                        name = 'New death notifications', type = 'scatter', mode = 'lines',
                        line = list(color = col.deaths, width = lwd.data, shape = 'hvh'),
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
                                title = "Nowcast",
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
  if(plotcumulative == TRUE & annotations == TRUE){
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
  }
  # }else{
  #   p_nowcast <- p_nowcast %>% 
  #     plotly::layout(
  #       annotations = list(yref = 'paper', xref = 'paper',
  #                          y = 1, x = .1, align = "left",
  #                          text = paste("Cumulative case notifications through:",format(max(database$Date),"%B %d, %Y"),
  #                                       display(tail(database$cum.cases,1)),"\n",
  #                                       "Cumulative deaths through:",format(max(database$Date),"%B %d, %Y"),
  #                                       display(tail(database$cum.deaths,1))
  #                          ),
  #                          showarrow = FALSE,
  #                          xanchor = "left")
  #     )
  # }
  p_nowcast
}


# plot nowcast from death reports -------------------------------------------------------------

plot_nowcast_from_death_reports <- function(database, maxy = 10^7, theme = 'default') {
  
  # VISUAL VARIABLES ***************************
  if(is.character(theme)) {
    if(theme=='default') {
      nctheme <- get0('nowcast_theme_default', envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())
    }else{
      nctheme <- get0(theme, envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())    
    }
  }else{
    if(is.list(theme)) { 
      nctheme <- nowcast_theme()
      nctheme[[names(theme)[i]]] <- theme[[i]]
    }
  }
  for(i in seq_along(nctheme)) { 
    assign(names(nctheme)[i], nctheme[[i]]) 
  }
  
  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  # END VISUAL VARIABLES ***********************
  
  
  p_nowcast <- plotly::plot_ly(data = database, x = ~Date) %>% 
    plotly::add_trace(y = ~deaths, type = 'bar', 
                      name = 'New death notifications', 
                      marker = list(color = col.cases),
                      legendgroup = 'group1'
  ) %>% 
    plotly::add_trace(y = ~I, 
                      name = 'Symptomatic cases', mode = 'lines',
                      line = list(color = col.I, width = lwd.data),
                      legendgroup = 'group2') %>%
    plotly::add_trace(y = ~I.forecast.mean, 
                      name = '(forecast average)', mode = 'lines',
                      line = list(color = col.I, width = lwd.mean, dash = 'dot'),
                      legendgroup = 'group2') %>% 
    plotly::add_ribbons(ymin = ~I.forecast.lower80, ymax = ~I.forecast.upper80,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.I, width = lwd.ci),
                        fillcolor = col.I.ci,
                        legendgroup = 'group2') %>% 
    
    plotly::add_trace(y = ~E, 
                      name = 'Latent cases', mode = 'lines',
                      line = list(color = col.E, width = lwd.data),
                      legendgroup = 'group3') %>%
    plotly::add_trace(y = ~E.forecast.mean, 
                      name = '(forecast average)', mode = 'lines',
                      line = list(color = col.E, width = lwd.mean, dash = 'dot'),
                      legendgroup = 'group3') %>% 
    plotly::add_ribbons(ymin = ~E.forecast.lower80, ymax = ~E.forecast.upper80,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.E, width = lwd.ci),
                        fillcolor = col.E.ci,
                        legendgroup = 'group3') %>% 
    
    plotly::add_trace(y = ~nowcast.mean, 
                      name = 'Total unnotified cases', mode = 'lines',
                      line = list(color = col.nowcast, width = lwd.data, dash = 'dot'),
                      legendgroup = 'group4') %>% 
    plotly::add_ribbons(ymin = ~nowcast.lower, ymax = ~nowcast.upper,
                        name = '(prediction interval)', mode='lines',
                        line = list(color = col.nowcast, width = lwd.ci),
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

# plot ascertainment --------------------------------------------------------------------------

plot_ascertainment <- function(ascertainment, logy = TRUE, theme = 'default'){

  # VISUAL VARIABLES ***************************
  if(is.character(theme)) {
    if(theme=='default') {
      nctheme <- get0('nowcast_theme_default', envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())
    }else{
      nctheme <- get0(theme, envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())    
    }
  }else{
    if(is.list(theme)) { 
      nctheme <- nowcast_theme()
      nctheme[[names(theme)[i]]] <- theme[[i]]
    }
  }
  for(i in seq_along(nctheme)) { 
    assign(names(nctheme)[i], nctheme[[i]]) 
  }
  
  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  # END VISUAL VARIABLES ***********************
  
  p_ascertainment <- plotly::plot_ly(data = ascertainment, x = ~Date , y = ~mean, type = 'scatter',
                                     name = 'Ascertainment', mode = 'lines',
                                     line = list(color = col.other, width = lwd.data)
  ) %>%
    plotly::add_trace(y = ~lower,
                      name = 'upper 80% confidence', mode = 'lines',
                      line = list(color = col.other, width = lwd.ci, dash = 'dot')) %>%
    plotly::add_trace(y = ~upper,
                      name = 'lower 80% confidence', mode = 'lines',
                      line = list(color = col.other, width = lwd.ci)) %>%
    plotly::layout(yaxis = list(type = ifelse(logy==TRUE,"log","linear"),
                                tickformat = ".2%",
                                title="Ascertainment"))
  p_ascertainment
}

# plot R effective ----------------------------------------------------------------------------

plot_R_effective <- function(database, legend=TRUE, maxy = 10, theme = 'default') {

  # VISUAL VARIABLES ***************************
  if(is.character(theme)) {
    if(theme=='default') {
      nctheme <- get0('nowcast_theme_default', envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())
    }else{
      nctheme <- get0(theme, envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())    
    }
  }else{
    if(is.list(theme)) { 
      nctheme <- nowcast_theme()
      nctheme[[names(theme)[i]]] <- theme[[i]]
    }
  }
  for(i in seq_along(nctheme)) { 
    assign(names(nctheme)[i], nctheme[[i]]) 
  }
  
  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  # END VISUAL VARIABLES ***********************
  

  p_Reff <- plotly::plot_ly(data = database, x = ~Date , y = ~R_eff.mean, type = 'scatter',
                                     name = 'R (mean)', mode = 'lines',
                                     line = list(color = col.R, width = lwd.data),
                            hoverinfo = "x+text",
                            hoverlabel = hoverlabel,
                            text = ~paste("mean R =", display(R_eff.mean),"\n",
                                          "(range:", display(R_eff.lower), "to",
                                          display(R_eff.upper),")")) %>%
    plotly::add_trace(y = ~R_eff.upper,
                      name = 'upper 80% confidence', mode = 'lines',
                      line = list(color = col.R, width = lwd.ci, dash = 'dot'),
                      hoverinfo = "none") %>%
    plotly::add_trace(y = ~R_eff.lower,
                      name = 'lower 80% confidence', mode = 'lines',
                      line = list(color = col.R, width = lwd.ci),
                      hoverinfo = "none") %>%
    plotly::layout(yaxis = list(# type = "log", 
                                # tickformat = ".2%", 
                                title="R")) %>% 
    plotly::layout(shapes = list(list(type = "line", 
                                      line = list(color = "black", width=2, dash = "dot"), 
                                      xref = "x", yref = "y",
                                      x0 = ~min(Date), x1 = ~max(Date), 
                                      y0 = 1, y1 = 1)),
                   annotations = list(yref = 'y', xref = 'paper',
                                      y = 1, x = .05, align = "left",
                                      text = "R = 1",
                                      showarrow = FALSE,
                                      xanchor = "left", yanchor = "bottom"),
                   yaxis = list(range=c(0,maxy), spikethickness = 0),
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

# plot cases ----------------------------------------------------------------------------------

plot_cases <- function(database, maxy = max(database$cases), plottrend = TRUE, 
                       logy = TRUE, legend = FALSE,
                       theme = 'default') {

  # VISUAL VARIABLES ***************************
  if(is.character(theme)) {
    if(theme=='default') {
      nctheme <- get0('nowcast_theme_default', envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())
    }else{
      nctheme <- get0(theme, envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())    
    }
  }else{
    if(is.list(theme)) { 
      nctheme <- nowcast_theme()
      nctheme[[names(theme)[i]]] <- theme[[i]]
    }
  }
  for(i in seq_along(nctheme)) { 
    assign(names(nctheme)[i], nctheme[[i]]) 
  }
  
  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  # END VISUAL VARIABLES ***********************
  
  if(plottrend==TRUE){
    database$cases.ravg <- movingAverage(database$cases, window=7,centered = FALSE)
  }
  
  p_cases <- plotly::plot_ly(data = database, x = ~Date) %>% 
    plotly::add_trace(y = ~cases, type = 'bar',
                      name = 'New case notifications',
                      marker = list(color = col.cases),
                      legendgroup = 'group_cases',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(cases),serif("new case reports")) )

  p_cases <- p_cases %>% 
    plotly::layout(yaxis = list(type = ifelse(logy==TRUE,"log","linear"),
                                range = ifelse(logy==TRUE,c(-.25,log10(maxy)),c(0,maxy)),
                                title = "Cases",
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
                   title = ""
    )
  
  if(plottrend == TRUE) {
    p_cases <- p_cases %>% 
      plotly::add_trace(y = ~cases.ravg, 
                        name = '7 day moving average', type = 'scatter', mode = 'lines',
                        line = list(color = col.cases, width = lwd.data, shape = 'hvh'),
                        legendgroup = 'group_cases',
                        hoverinfo = "x+text",
                        hoverlabel = hoverlabel,
                        text = ~paste(display(cases.ravg),serif("new cases (7 day average)")) )
  }
  
  p_cases
}


# plot deaths ---------------------------------------------------------------------------------


plot_deaths <- function(database, maxy = max(database$deaths), plottrend = TRUE, 
                        logy = TRUE, legend = FALSE,
                        theme = 'default') {
  
  # VISUAL VARIABLES ***************************
  if(is.character(theme)) {
    if(theme=='default') {
      nctheme <- get0('nowcast_theme_default', envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())
    }else{
      nctheme <- get0(theme, envir = .GlobalEnv, mode = "any", inherits = TRUE,
                      ifnotfound = nowcast_theme())    
    }
  }else{
    if(is.list(theme)) { 
      nctheme <- nowcast_theme()
      nctheme[[names(theme)[i]]] <- theme[[i]]
    }
  }
  for(i in seq_along(nctheme)) { 
    assign(names(nctheme)[i], nctheme[[i]]) 
  }
  
  plotfont <- list(family = "serif")
  hoverlabel <- list(namelength = -1,
                     font = list(family = "serif"))
  # END VISUAL VARIABLES ***********************
  
  if(plottrend==TRUE){
    database$deaths.ravg <- movingAverage(database$deaths, window=7,centered = FALSE)
  }
  
  p_deaths <- plotly::plot_ly(data = database, x = ~Date) %>% 
    plotly::add_trace(y = ~deaths, type = 'bar',
                      name = 'New death notifications',
                      marker = list(color = col.deaths),
                      legendgroup = 'group_deaths',
                      hoverinfo = "x+text",
                      hoverlabel = hoverlabel,
                      text = ~paste(display(deaths),serif("new death reports")) )
  
  p_deaths <- p_deaths %>% 
    plotly::layout(yaxis = list(type = ifelse(logy==TRUE,"log","linear"),
                                range = ifelse(logy==TRUE,c(-.25,log10(maxy)),c(0,maxy)),
                                title = "Deaths",
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
                   title = ""
    )
  
  if(plottrend == TRUE) {
    p_deaths <- p_deaths %>% 
      plotly::add_trace(y = ~deaths.ravg, 
                        name = '7 day moving average', type = 'scatter', mode = 'lines',
                        line = list(color = col.deaths, width = lwd.data, shape = 'hvh'),
                        legendgroup = 'group_deaths',
                        hoverinfo = "x+text",
                        hoverlabel = hoverlabel,
                        text = ~paste(display(deaths.ravg),serif("new deaths (7 day average)")) )
  }
  
  p_deaths
}


# plot dashboard ------------------------------------------------------------------------------


plot_dashboard <- function(database, params, 
                           daterange = c("2020-02-15", #start
                                         as.character(Sys.Date())  #end
                                         ),
                           plotcumulative = FALSE,
                           theme = 'default') {
  p_asc <- plot_ascertainment(params$q, logy = FALSE, theme = theme)
  p_reff <- plot_R_effective(database, legend = TRUE, maxy = 4, theme = theme)
  p_nc <- plot_nowcast_from_case_reports(database, legend = TRUE, logy = FALSE, 
                                         maxy = max(database$nowcast.upper),
                                         plotcases = FALSE,
                                         plotdeaths = FALSE,
                                         plotcumulative = plotcumulative,
                                         theme = theme)
  p_cases <- plot_cases(database, legend = TRUE, logy = FALSE, theme = theme)
  p_deaths <- plot_deaths(database, legend = TRUE, logy = FALSE, theme = theme)
  
  dash <- plotly::subplot(p_nc, p_cases, p_deaths, p_asc, p_reff,
                          nrows = 5,
                          heights = NULL,
                          shareX = TRUE,
                          titleY = TRUE) %>% 
    plotly::layout(
      xaxis = list(range = c(as.numeric(as.POSIXct(daterange[1], format="%Y-%m-%d"))*1000,
                             as.numeric(as.POSIXct(daterange[2], format="%Y-%m-%d"))*1000),
                   type = "date"
      )
    )
  dash
}

