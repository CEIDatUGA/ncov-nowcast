
library(glmnet)
library(signal)
library(ggplot2)
library(reshape2)

# # required for mosaic 
# install.packages(zeallot)
# install.packages(broom)

library(mosaic)
library(ModelMetrics)


#' Simulates different methods of deconvoluting desired outbreak difficulties.
#'
#' @param n The number of outbreaks to simulate. Defaults to 1000.
#' @param mean The mean incubation period.
#' @param sd The standard deviation of the incubation period.
#' @param difficulty A vector containing the desired difficulty levels. Default c("easy","medium","hard").
#' @param distribution A string containing the distribution of the incubation period ("gamma" or "lognormal"). Defaults to "gamma".
#' @param methods A vector containing each of the methods to be tested. Choices include anything in the default value of c("simple","random","ridge","rl","fourier","frequency","average")
#' @param rl0 A string containing the name of the method to be used as the initial estimate for the Richardson-Lucy method. Defaults to "random".
#' @param average_methods A vector containing boolean values for which methods to include in the average method. (simple, random, ridge, rl, fourier, frequency) Defaults to true for all values.
#' @param diagnostics A boolean value indicating whether to compute correlation, rmse, and coverage. Defaults to TRUE.
#'
#' @return A list containing all the computed pieces
#'
#' @examples
#' simulate_trials(n=500, mean=9.7, sd=5.5, distribution="gamma", rl0="random", average_methods=c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
#'
#' @export
simulate_trials = function(n=1000, 
                           mean,
                           sd,
                           difficulty=c("easy","medium","hard"),
                           distribution="gamma",
                           methods=c("simple","random","ridge","rl","fourier","frequency","average"), 
                           rl0="random", 
                           average_methods = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                           diagnostics=TRUE)
{
  if(rl0 %in% methods) print("thanks for being sensible")
  else rl0=methods[1]
  
  if(distribution == "gamma")
  {
    incubation_distribution=gamma_incubation_period
  }
  if(distribution == "lognorm")
  {
    incubation_distribution=lnorm_incubation_period
  }
  
  results=list()
  if("easy" %in% difficulty)
  {
    easy_results=list()
    easy_infection_curves = generate_easy_trials(n)
    easy_symptom_onset_curves = convolute_infection_curves(easy_infection_curves, incubation_period=incubation_distribution, parms=c(mean, sd))
    easy_results = deconvolve_trials(mean, sd, distribution, methods, rl0, average_methods, diagnostics, easy_infection_curves, easy_symptom_onset_curves)
    easy_results[["infection_curves"]]=easy_infection_curves
    easy_results[["symptom_onset_curves"]]=easy_symptom_onset_curves
    
    results[["easy"]]=easy_results
  }
  if("medium" %in% difficulty)
  {
    medium_results=list()
    medium_infection_curves = generate_medium_trials(n)
    medium_symptom_onset_curves = convolute_infection_curves(medium_infection_curves, incubation_period=incubation_distribution, parms=c(mean, sd))
    medium_results = deconvolve_trials(mean, sd, distribution, methods, rl0, average_methods, diagnostics, medium_infection_curves, medium_symptom_onset_curves)
    medium_results[["infection_curves"]]=medium_infection_curves
    medium_results[["symptom_onset_curves"]]=medium_symptom_onset_curves
    
    results[["medium"]]=medium_results
  }
  if("hard" %in% difficulty)
  {
    hard_results=list()
    hard_infection_curves = generate_hard_trials(n)
    hard_symptom_onset_curves = convolute_infection_curves(hard_infection_curves, incubation_period=incubation_distribution, parms=c(mean, sd))
    hard_results = deconvolve_trials(mean, sd, distribution, methods, rl0, average_methods, diagnostics, hard_infection_curves, hard_symptom_onset_curves)
    hard_results[["infection_curves"]]=hard_infection_curves
    hard_results[["symptom_onset_curves"]]=hard_symptom_onset_curves
    
    results[["hard"]]=hard_results
  }
  return(results)
}

#' Simulates different methods of deconvoluting outbreaks.
#'
#' @param mean The mean incubation period.
#' @param sd The standard deviation of the incubation period.
#' @param distribution A string containing the distribution of the incubation period ("gamma" or "lognormal"). Defaults to "gamma".
#' @param methods A vector containing each of the methods to be tested. Choices include anything in the default value of c("simple","random","ridge","rl","fourier","frequency","average")
#' @param rl0 A string containing the name of the method to be used as the initial estimate for the Richardson-Lucy method. Defaults to "random".
#' @param average_methods A vector containing boolean values for which methods to include in the average method. (simple, random, ridge, rl, fourier, frequency) Defaults to true for all values.
#' @param diagnostics A boolean value indicating whether to compute correlation, rmse, and coverage. Defaults to TRUE.
#' @param infection_curves A dataframe containing the original infection curves. This is used strictly for diagnostics and does not aid in the deconvolution process. Defaults to NULL.
#' @param symptom_onset_curves A dataframe containing the symptom-onset curves. Used in combination with the mean, sd, and distribution to deconvolve the curves.
#' 
#' @return A list containing all the computed pieces
#'
#' @examples
#' simulate_easy_trials(n=500, mean=9.7, sd=5.5, distribution="gamma", rl0="random", average_methods=c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
#'
#' @export
deconvolve_trials = function(mean,
                             sd,
                             distribution="gamma",
                             methods=c("simple","random","ridge","rl","fourier","frequency","average"), 
                             rl0="random", 
                             average_methods = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                             diagnostics=TRUE,
                             infection_curves=NULL,
                             symptom_onset_curves)
{
  results = list()
  
  if("simple" %in% methods)
  {
    deconvolution_simple = deconvolve_infection_curves_simple(symptom_onset_curves, distribution=distribution, parms=c(mean, sd))
    results[["deconvolution_simple"]]=deconvolution_simple
    if(diagnostics && !is.null(infection_curves))
    {
      deconvolution_simple_rmse = determine_rmse(infection_curves, deconvolution_simple)
      deconvolution_simple_correlation = determine_correlation(infection_curves, deconvolution_simple)
      deconvolution_simple_coverage = determine_coverage(infection_curves, deconvolution_simple)
      results[["deconvolution_simple_rmse"]]=deconvolution_simple_rmse
      results[["deconvolution_simple_correlation"]]=deconvolution_simple_correlation
      results[["deconvolution_simple_coverage"]]=deconvolution_simple_coverage
    }
  }
  else 
  {
    average_methods[1]=FALSE
    deconvolution_simple=0
  }
  
  if("random" %in% methods)
  {
    deconvolution_random = deconvolve_infection_curves_random(symptom_onset_curves, distribution = distribution, parms=c(mean, sd))
    results[["deconvolution_random"]]=deconvolution_random
    if(diagnostics && !is.null(infection_curves))
    {
      deconvolution_random_rmse = determine_rmse(infection_curves, deconvolution_random)
      deconvolution_random_correlation = determine_correlation(infection_curves, deconvolution_random)
      deconvolution_random_coverage = determine_coverage(infection_curves, deconvolution_random)
      results[["deconvolution_random_rmse"]]=deconvolution_random_rmse
      results[["deconvolution_random_correlation"]]=deconvolution_random_correlation
      results[["deconvolution_random_coverage"]]=deconvolution_random_coverage
    }
  }
  else 
  {
    average_methods[2]=FALSE
    deconvolution_random=0
  }
  
  if("ridge" %in% methods)
  {
    deconvolution_ridge = deconvolve_infection_curves_ridge(symptom_onset_curves, distribution = distribution, parms = c(mean, sd))
    results[["deconvolution_ridge"]]=deconvolution_ridge
    if(diagnostics && !is.null(infection_curves))
    {
      deconvolution_ridge_rmse = determine_rmse(infection_curves, deconvolution_ridge)
      deconvolution_ridge_correlation = determine_correlation(infection_curves, deconvolution_ridge)
      deconvolution_ridge_coverage = determine_coverage(infection_curves, deconvolution_ridge)
      results[["deconvolution_ridge_rmse"]]=deconvolution_ridge_rmse
      results[["deconvolution_ridge_correlation"]]=deconvolution_ridge_correlation
      results[["deconvolution_ridge_coverage"]]=deconvolution_ridge_coverage
    }
  }
  else 
  {
    average_methods[3]=FALSE
    deconvolution_ridge=0
  }
  
  if("fourier" %in% methods)
  {
    deconvolution_fourier = deconvolve_infection_curves_fourier(symptom_onset_curves, distribution = distribution, parms = c(mean, sd))
    results[["deconvolution_fourier"]]=deconvolution_fourier
    if(diagnostics && !is.null(infection_curves))
    {
      deconvolution_fourier_rmse = determine_rmse(infection_curves, deconvolution_fourier)
      deconvolution_fourier_correlation = determine_correlation(infection_curves, deconvolution_fourier)
      deconvolution_fourier_coverage = determine_coverage(infection_curves, deconvolution_fourier)
      results[["deconvolution_fourier_rmse"]]=deconvolution_fourier_rmse
      results[["deconvolution_fourier_correlation"]]=deconvolution_fourier_correlation
      results[["deconvolution_fourier_coverage"]]=deconvolution_fourier_coverage
    }
  }
  else 
  {
    average_methods[5]=FALSE
    deconvolution_fourier=0
  }
  
  if("frequency" %in% methods)
  {
    deconvolution_frequency = deconvolve_infection_curves_frequency(symptom_onset_curves, distribution = distribution, parms = c(mean, sd))
    results[["deconvolution_frequency"]]=deconvolution_frequency
    if(diagnostics && !is.null(infection_curves))
    {
      deconvolution_frequency_rmse = determine_rmse(infection_curves, deconvolution_frequency)
      deconvolution_frequency_correlation = determine_correlation(infection_curves, deconvolution_frequency)
      deconvolution_frequency_coverage = determine_coverage(infection_curves, deconvolution_frequency)
      results[["deconvolution_frequency_rmse"]]=deconvolution_frequency_rmse
      results[["deconvolution_frequency_correlation"]]=deconvolution_frequency_correlation
      results[["deconvolution_frequency_coverage"]]=deconvolution_frequency_coverage
    }
  }
  else 
  {
    average_methods[6]=FALSE
    deconvolution_frequency=0
  }
  
  if("rl" %in% methods)
  {
    if(rl0 == "simple") estimate0=deconvolution_simple
    else if(rl0 == "random") estimate0=deconvolution_random
    else if(rl0 == "ridge") estimate0=deconvolution_ridge
    else if(rl0 == "fourier") estimate0=deconvolution_fourier
    else if(rl0 == "frequency") estimate0=deconvolution_frequency
    deconvolution_rl = deconvolve_infection_curves_rl(symptom_onset_curves, estimates_matrix = estimate0)
    results[["deconvolution_rl"]]=deconvolution_rl
    if(diagnostics && !is.null(infection_curves))
    {
      deconvolution_rl_rmse = determine_rmse(infection_curves, deconvolution_rl)
      deconvolution_rl_correlation = determine_correlation(infection_curves, deconvolution_rl)
      deconvolution_rl_coverage = determine_coverage(infection_curves, deconvolution_rl)
      results[["deconvolution_rl_rmse"]]=deconvolution_rl_rmse
      results[["deconvolution_rl_correlation"]]=deconvolution_rl_correlation
      results[["deconvolution_rl_coverage"]]=deconvolution_rl_coverage
    }
  }
  else 
  {
    average_methods[4]=FALSE
    deconvolution_rl=0
  }
  
  if("average" %in% methods)
  {
    deconvolution_average_methods = list(deconvolution_simple, deconvolution_random, deconvolution_ridge, deconvolution_rl, deconvolution_fourier, deconvolution_frequency)
    deconvolution_average = deconvolve_infection_curves_average(deconvolution_average_methods[average_methods])
    results[["deconvolution_average"]]=deconvolution_average
    if(diagnostics && !is.null(infection_curves))
    {
      deconvolution_average_rmse = determine_rmse(infection_curves, deconvolution_average)
      deconvolution_average_correlation = determine_correlation(infection_curves, deconvolution_average)
      deconvolution_average_coverage = determine_coverage(infection_curves, deconvolution_average)
      results[["deconvolution_average_rmse"]]=deconvolution_average_rmse
      results[["deconvolution_average_correlation"]]=deconvolution_average_correlation
      results[["deconvolution_average_coverage"]]=deconvolution_average_coverage
    }
  }
  
  return(results)
}


#' Generates random infection curves with an easy deconvolution difficulty.
#'
#' @param n The number of trials to simulate.
#' 
#' @return A dataframe containing columns of generated infection curves
#'
#' @examples
#' generate_easy_trials(1000)
#'
#' @export
generate_easy_trials = function(n)
{
  padding = 50
  days = -100:100
  easy_infection_curves = data.frame(day = matrix(1:(length(days) + 2 * padding), nrow = (length(days) + 2 * padding)))
  for(i in 1:n)
  {
    repeat
    {
      random_variable_one = abs(rnorm(1, mean = 1, sd = 0.1))
      random_variable_two = abs(rnorm(1, mean = 250, sd = 5))
      random_variable_three = abs(rnorm(2, mean = 6, sd = 3))
      smoothed_infection_curve = -days^2/abs(random_variable_one * random_variable_two[1]) + random_variable_two/random_variable_three[2]
      infection_curve_noise = round(smoothed_infection_curve + rnorm(length(days), mean = 1, sd = 1))
      if(infection_curve_noise[1] < 0)
      {
        break
      }
    }
    infection_curve_noise[which(infection_curve_noise < 0)] = 0
    padded_infection_curve = c(rep(0, times = padding), infection_curve_noise, rep(0, times = padding))
    easy_infection_curves = cbind(easy_infection_curves, padded_infection_curve)
  }
  names(easy_infection_curves) = c("day", factor(1:n))
  return(easy_infection_curves)
}


#' Generates random infection curves with a medium deconvolution difficulty.
#'
#' @param n The number of trials to simulate.
#' 
#' @return A dataframe containing columns of generated infection curves
#'
#' @examples
#' generate_medium_trials(1000)
#'
#' @export
generate_medium_trials = function(n)
{
  padding = 50
  days = seq(-2.8, 2.8, length = 201)
  medium_infection_curves = data.frame(day = matrix(1:(length(days) + 2 * padding), nrow = (length(days) + 2 * padding)))
  for(i in 1:n)
  {
    repeat
    {
      random_variable_one = abs(rnorm(1, mean = 2.5, sd = 0.05))
      random_variable_two = abs(rnorm(1, mean = 21, sd = 0.5))
      random_variable_three = abs(rnorm(1, mean = 58, sd = 1))
      random_variable_four = abs(rnorm(1, mean = 45, sd = 4))
      smoothed_infection_curve = 0.05 * days^10 - random_variable_one * days^8 + random_variable_two * days^6 - 
        random_variable_three * days^4 + random_variable_four * days^2 + 20
      infection_curve_noise = round(smoothed_infection_curve + rnorm(length(days), mean = 1, sd = 1))
      if(infection_curve_noise[1] < 0 && max(infection_curve_noise) < 50)
      {
        break
      }
    }
    infection_curve_noise[which(infection_curve_noise < 0)] = 0
    padded_infection_curve = c(rep(0, times = padding), infection_curve_noise, rep(0, times = padding))
    medium_infection_curves = cbind(medium_infection_curves, padded_infection_curve)
  }
  names(medium_infection_curves) = c("day", factor(1:n))
  return(medium_infection_curves)
}


#' Generates random infection curves with a hard deconvolution difficulty.
#'
#' @param n The number of trials to simulate.
#' 
#' @return A dataframe containing columns of generated infection curves
#'
#' @examples
#' generate_hard_trials(1000)
#'
#' @export
generate_hard_trials = function(n)
{
  padding = 50
  days = seq(-5, 5, length = 201)
  hard_infection_curves = data.frame(day = matrix(1:(length(days) + 2 * padding), nrow = (length(days) + 2 * padding)))
  for(i in 1:n)
  {
    repeat
    {
      random_variable_one = abs(rnorm(1, mean = 1, sd = 0.03))
      random_variable_two = abs(rnorm(1, mean = 5, sd = 1))
      random_variable_three = abs(rnorm(1, mean = 2, sd = 0.4))
      smoothed_infection_curve = -random_variable_one * days^2 + random_variable_two * sin(random_variable_three * pi * days) + 20
      infection_curve_noise = round(smoothed_infection_curve + rnorm(length(days), mean = 1, sd = 1))
      if(infection_curve_noise[1] < 0)
      {
        break
      }
    }
    infection_curve_noise[which(infection_curve_noise < 0)] = 0
    padded_infection_curve = c(rep(0, times = padding), infection_curve_noise, rep(0, times = padding))
    hard_infection_curves = cbind(hard_infection_curves, padded_infection_curve)
  }
  names(hard_infection_curves) = c("day", factor(1:n))
  return(hard_infection_curves)
}


#' Generates random incubation according to the gamma distribution and the given parameters
#'
#' @param parms A vector containing the mean and standard deviation for the incubation period.
#' 
#' @return A random incubation period.
#'
#' @export
gamma_incubation_period = function(parms)
{
  shape=parms[1]^2/parms[2]^2
  scale=parms[2]^2/parms[1]
  return(round(rgamma(n = 1, shape = shape, scale = scale)))
}


#' Generates random incubation according to the lognormal distribution and the given parameters
#'
#' @param parms A vector containing the mean and standard deviation for the incubation period.
#' 
#' @return A random incubation period.
#'
#' @export
lnorm_incubation_period = function(parms)
{
  mean=log(parms[1])-0.5*log((parms[2]/parms[1])^2 + 1)
  sd=sqrt(log((parms[2]/parms[1])^2 + 1))
  return(round(rlnorm(n=1, meanlog=mean, sdlog=sd)))
}

#' Convolves a dataframe filled with infection curves.
#'
#' @param infection_curves A data frame containing columns of infection curves to be convolved, preferably from one of the generate trial functions. The first column understood to be a day index.
#' @param incubation_period The name of a function for randomly generating incubation periods using given parameters. Defaults to a gamma distribution.
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing columns of convolved symptom-onset curves
#'
#' @examples
#' convolute_infection_curves(infection curves = easy_infection_curves, incubation_period = lnorm_incubation_period)
#'
#' @export
convolute_infection_curves = function(infection_curves, incubation_period = gamma_incubation_period, parms=c(9.7, 5.5))
{
  convoluted_trials = as.data.frame(infection_curves[,1])
  for(i in 2:length(infection_curves[1,]))
  {
    symptom_onset_curve = generate_symptom_onset_curve(infection_curves[,i], incubation_period, parms)
    convoluted_trials = cbind(convoluted_trials, symptom_onset_curve)
  }
  return(convoluted_trials)
}

#' Generates a single symptom-onset curve from a single infection curve.
#'
#' @param infection_curve A vector containing an infection curve to be convolved.
#' @param incubation_period The name of a function for randomly generating incubation periods using given parameters. Defaults to a gamma distribution.
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A vector containing a single symptom-onset curve.
#'
#' @examples
#' generate_symptom_onset_curve(infection_curve = c(0,0,0,...,0,0,0), incubation_period = lnorm_incubation_period, parms = c(15.3, 9.1))
#'
#' @export
generate_symptom_onset_curve = function(infection_curve, incubation_period = gamma_incubation_period, parms = c(9.7, 5.5))
{
  outbreak_duration = length(infection_curve)
  symptom_onset_curve = rep(0, outbreak_duration)
  for(day in 1:outbreak_duration)
  {
    if(infection_curve[day] > 0)
    {
      for(patient in 1:infection_curve[day])
      {
        symptom_onset_day = day + incubation_period(parms)
        if(symptom_onset_day > outbreak_duration)
        {
          symptom_onset_day = outbreak_duration
        }
        symptom_onset_curve[symptom_onset_day] = symptom_onset_curve[symptom_onset_day] + 1
      }
    }
  }
  return(symptom_onset_curve)
}

#' Determine the median incubation period.
#'
#' @param distribution A string containing the desired distribution, whether "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period.
#' 
#' @return The median incubation period.
#'
#' @examples
#' median_incubation_period(parms = c(15.3, 9.1))
median_incubation_period = function(distribution="gamma", parms)
{
  if(distribution == "gamma")
  {
    shape=parms[1]^2/parms[2]^2
    scale=parms[2]^2/parms[1]
    return(round(qgamma(0.5, shape = shape, scale = scale)))
  }
  else if(distribution == "lognorm") 
  {
    mean=log(parms[1])-0.5*log((parms[2]/parms[1])^2 + 1)
    sd=sqrt(log((parms[2]/parms[1])^2 + 1))
    return(round(qlnorm(0.5, meanlog=mean, sdlog=sd)))
  }
}


#' Determines the number of days before 99.9% of all infected people are estimated to have begun showing symptoms
#'
#' @param distribution A string containing the desired distribution, whether "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period.
#' 
#' @return The length of the incubation period.
#'
#' @examples
#' determine_incubation_length(parms = c(15.3, 9.1))
determine_incubation_length = function(distribution, parms)
{
  if(distribution == "gamma")
  {
    shape=parms[1]^2/parms[2]^2
    scale=parms[2]^2/parms[1]
    return(round(qgamma(0.999, shape = shape, scale = scale)))
  }
  else if(distribution == "lognorm")
  {
    mean=log(parms[1])-0.5*log((parms[2]/parms[1])^2 + 1)
    sd=sqrt(log((parms[2]/parms[1])^2 + 1))
    return(round(qlnorm(0.999, meanlog=mean, sdlog=sd)))
  }
}


#' Generates a vector containing the probability that someone infected on the day of the first index begins showing symptoms on subsequent days.
#'
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A vector containing the probability that someone begins showing symptoms i days after infection. This vector is zero based.
#'
#' @examples
#' incubation_period_daily_probability(distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
incubation_period_daily_probability = function(distribution = "gamma", parms = c(9.7, 5.5))
{
  if(distribution == "gamma")
  {
    shape=parms[1]^2/parms[2]^2
    scale=parms[2]^2/parms[1]
    
    incubation_period_length = round(qgamma(0.9999, shape = shape, scale = scale))
    probability = c(pgamma(0, shape = shape, scale = scale), 
                    pgamma(seq(0.5, incubation_period_length - 0.5, by = 1), shape = shape, scale = scale),
                    pgamma(incubation_period_length, shape = shape, scale = scale))
  }
  else if(distribution == "lognorm")
  {
    mean=log(parms[1])-0.5*log((parms[2]/parms[1])^2 + 1)
    sd=sqrt(log((parms[2]/parms[1])^2 + 1))
    incubation_period_length = round(qlnorm(0.9999, meanlog=mean, sdlog=sd))
    probability = c(plnorm(0, meanlog=mean, sdlog=sd), 
                    plnorm(seq(0.5, incubation_period_length - 0.5, by = 1), meanlog=mean, sdlog=sd),
                    plnorm(incubation_period_length, meanlog=mean, sdlog=sd))
    
  }
  daily_probability = c()
  for(i in 2:length(probability))
  {
    daily_probability = c(daily_probability, (probability[i] - probability[i - 1]))
  }
  return(daily_probability)
}


#' This function returns a matrix with each entry containing the probability that someone infected on day i begins showing symptoms on day j, where i is the column number and j is the row number. The format is slightly different when using this matrix for the Richardson-Lucy method, where the first rows are 0 filled according to the length of the incubation period.
#'
#' @param outbreak_duration A value indicating the 99.9th percentile for the incubation period
#' @param incubation_probability The name of a function that returns the probability of showing symptoms on i-1 days after becoming infected, where i is the index of the vector.
#' @param rl A boolean value indicating whether this matrix is intended to be used for the Richardson-Lucy method. Defaults to FALSE.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe with indices indicating the probability that someone infected on a given day begins showing symptoms on subsequent days.
#'
#' @examples
#' incubation_period_distribution_matrix(outbreak_duration = 60, rl = TRUE, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
incubation_period_distribution_matrix = function(outbreak_duration, 
                                                 incubation_probability = incubation_period_daily_probability, 
                                                 rl = FALSE,
                                                 distribution = "gamma",
                                                 parms = c(9.7,5.5))
{
  daily_probability = incubation_probability(distribution, parms)
  daily_probability_transition = c(daily_probability, rep(0, times = outbreak_duration))
  daily_probability_transition = daily_probability_transition[1:outbreak_duration]
  incubation_distribution_matrix = matrix(nrow = outbreak_duration)
  for(i in 1:outbreak_duration)
  {
    incubation_distribution_matrix = cbind(incubation_distribution_matrix, daily_probability_transition)
    daily_probability_transition = c(0, daily_probability_transition[1:outbreak_duration - 1])
  }
  if(rl)
  {
    for(i in 1:(length(daily_probability) - 1))
    {
      incubation_distribution_matrix[i,] = rep(0, times = length(incubation_distribution_matrix[i,]))
    }
  }
  return(incubation_distribution_matrix[,2:(outbreak_duration + 1)])
}


#' Rounds the estimated infection curve while maintaining the correct number of infections.
#'
#' @param symptom_curve The original symptom-onset curve.
#' @param infection_curve The estimated infection curve.
#' 
#' @return A properly rounded infection curve.
#'
#' @export
round_infection_curve = function(symptom_curve, infection_curve)
{
  infection_curve_rounded = round(infection_curve)
  symptom_curve_sum = sum(symptom_curve)
  infection_curve_sum = sum(infection_curve)
  rounding_error = infection_curve_rounded - infection_curve
  repeat
  {
    if(symptom_curve_sum > infection_curve_sum)
    {
      index_to_increase = which.min(rounding_error)
      infection_curve_rounded[index_to_increase] = infection_curve_rounded[index_to_increase] + 1
      infection_curve_sum = sum(infection_curve_rounded)
      rounding_error[index_to_increase] = 0
    }
    else if(symptom_curve_sum < infection_curve_sum)
    {
      index_to_decrease = which.max(rounding_error)
      infection_curve_rounded[index_to_decrease] = infection_curve_rounded[index_to_decrease] - 1
      infection_curve_sum = sum(infection_curve_rounded)
      rounding_error[index_to_decrease] = 0
    }
    if(symptom_curve_sum == infection_curve_sum)
    {
      break
    }
  }
  return(infection_curve_rounded)
}


#' Takes a weekly symptom-onset curve and converts it to an average daily count.
#'
#' @param symptom_onset_curve The weekly symptom-onset curve.
into_days <- function(symptom_onset_curve)
{
  final = c()
  for(i in symptom_onset_curve)
  {
    final = c(final, rep(i, times = 7)/7)
  }
  return(final)
}


#' Deconvolves a dataframe of symptom-onset curves using simple deconvolution.
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param median_incubation_period_function The name of a function that returns the median incubation period according to the distribution and paramaters given.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing estimated infection curves using simple deconvolution.
#'
#' @examples
#' deconvolve_infection_curves_simple(easy_symptom_onset_curves, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
deconvolve_infection_curves_simple = function(symptom_onset_curves, 
                                              median_incubation_period_function = median_incubation_period,
                                              distribution="gamma",
                                              parms=c(9.7, 5.5))
{
  median_incubation_period = median_incubation_period_function(distribution, parms)
  deconvolved_infection_curves_simple = as.data.frame(symptom_onset_curves[,1])
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_simple(symptom_onset_curves[,i], median_incubation_period)
    deconvolved_infection_curves_simple = cbind(deconvolved_infection_curves_simple, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_simple)
}


#' Deconvolves a single symptom-onset curve using simple deconvolution.
#'
#' @param symptom_onset_curve A single symptom-onset curve.
#' @param median_incubation_period The median incubation period according to the distribution and paramaters for the incubation period.
#' 
#' @return A single estimated infection curve using simple deconvolution.
#'
#' @examples
#' deconvolve_infection_curve_simple(easy_symptom_onset_curve, 10)
#'
#' @export
deconvolve_infection_curve_simple = function(symptom_onset_curve, median_incubation_period)
{
  median_incubation_period=round(median_incubation_period)
  outbreak_duration = length(symptom_onset_curve)
  deconvolved_infection_curve_simple = c(symptom_onset_curve[(median_incubation_period + 1):outbreak_duration], 
                                         rep(0, times = median_incubation_period))
  #deconvolved_infection_curve_simple_rounded = round_infection_curve(symptom_onset_curve, deconvolved_infection_curve_simple)
  return(deconvolved_infection_curve_simple)
}


#' Generates a random incubation period under the defined distribution, mean, and standard deviation.
#'
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A single randomly generated incubation period.
#'
#' @examples
#' generate_incubation_period(distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
generate_incubation_period = function(distribution = "gamma", parms = c(9.7,5.5))
{
  if(distribution == "gamma") 
  {
    shape=parms[1]^2/parms[2]^2
    scale=parms[2]^2/parms[1]
    return(round(rgamma(n = 1, shape = shape, scale = scale)))
  }
  else if(distribution == "lognorm")
  {
    mean=log(parms[1])-0.5*log((parms[2]/parms[1])^2 + 1)
    sd=sqrt(log((parms[2]/parms[1])^2 + 1))
    return(round(rlnorm(n = 1, meanlog = mean, sdlog = sd)))
  }
}


#' Deconvolves a dataframe of symptom-onset curves using Random deconvolution.
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param incubation_period The name of a function that returns a randomly generated incubation period according to the distribution and paramaters given. Defaults to generate_incubation_period.
#' @param trials The number of simulated infection curves to generate. Defaults to 35.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing estimated infection curves using Random deconvolution.
#'
#' @examples
#' deconvolve_infection_curves_random(easy_symptom_onset_curves, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
deconvolve_infection_curves_random = function(symptom_onset_curves, 
                                              incubation_period = generate_incubation_period, 
                                              trials = 35,
                                              distribution="gamma",
                                              parms=c(9.7,5.5))
{
  deconvolved_infection_curves_random = as.data.frame(symptom_onset_curves[,1])
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_random(symptom_onset_curves[,i], incubation_period, trials, distribution, parms)
    deconvolved_infection_curves_random = cbind(deconvolved_infection_curves_random, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_random)
}


#' Deconvolves a single symptom-onset curve using Random deconvolution.
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param incubation_period The name of a function that returns a randomly generated incubation period according to the distribution and paramaters given.
#' @param trials The number of simulated infection curves to generate.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm".
#' @param parms A vector containing the mean and standard deviation for the incubation period.
#' 
#' @return A single estimated infection curve using Random deconvolution.
#'
#' @examples
#' deconvolve_infection_curve_random(easy_symptom_onset_curve, trials = 50, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
deconvolve_infection_curve_random = function(symptom_onset_curve, incubation_period, trials, distribution, parms, curves=FALSE)
{
  duration = length(symptom_onset_curve)
  possible_infection_curves = data.frame(matrix(nrow = duration, ncol = 0))
  for(trial in 1:trials){
    possible_infection_curve = rep(0, times = duration)
    for(day in 1:duration){
      if(symptom_onset_curve[day] > 0)
      {
        for(patient in 1:symptom_onset_curve[day])
        {
          estimated_infection_day = day - incubation_period(distribution, parms)
          if(estimated_infection_day <= 0)
          {
            estimated_infection_day = 1
          }
          possible_infection_curve[estimated_infection_day] = possible_infection_curve[estimated_infection_day] + 1
        }
      }
    }
    possible_infection_curves = cbind(possible_infection_curves, possible_infection_curve)
  }
  estimated_infection_curve = average_possible_infection_curves(possible_infection_curves)
  #estimated_infection_curve_rounded = round_infection_curve(symptom_onset_curve, estimated_infection_curve)
  if(curves==TRUE){
    return(list(result=estimated_infection_curve, estimates=possible_infection_curves))
  }
  return(estimated_infection_curve)
}

#' Finds the average estimated infection curve given a dataframe of possible infection curves.
#'
#' @param possible_infection_curves A dataframe containing possible infection curves.
#' 
#' @return A single estimated infection curve using Random deconvolution.
average_possible_infection_curves = function(possible_infection_curves)
{
  estimated_infection_curve = rowMeans(possible_infection_curves)
  return(estimated_infection_curve)
}


#' Deconvolves a dataframe of symptom-onset curves using a Ridge regression. This ridge regression  assumes that there is a linear equation that determines the number of patients who begin showing symptoms on each day of the outbreak. This equation can be solved easly, but because of the noise from the incubation period, the solution will be overfitted. The ridge regression is used to minimize overfitting by setting a penalty on large error values. 
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param incubation_probability The name of a function that returns a data frame containing daily symptom-onset probabilities according to the distribution and paramaters given. Defaults to incubation_period_distribution_matrix.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing estimated infection curves using a Ridge regression.
#'
#' @examples
#' deconvolve_infection_curves_ridge(easy_symptom_onset_curves, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
deconvolve_infection_curves_ridge = function(symptom_onset_curves, 
                                             incubation_probability = incubation_period_distribution_matrix,
                                             distribution = "gamma",
                                             parms = c(9.7,5.5))
{
  incubation_matrix = incubation_probability(length(symptom_onset_curves[,1]), distribution = distribution, parms = parms)
  deconvolved_infection_curves_ridge = as.data.frame(symptom_onset_curves[,1])
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_ridge(symptom_onset_curves[,i], incubation_matrix)
    deconvolved_infection_curves_ridge = cbind(deconvolved_infection_curves_ridge, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_ridge)
}


#' Deconvolves a single symptom-onset curve using a Ridge regression. This ridge regression  assumes that there is a linear equation that determines the number of patients who begin showing symptoms on each day of the outbreak. This equation can be solved easly, but because of the noise from the incubation period, the solution will be overfitted. The ridge regression is used to minimize overfitting by setting a penalty on large error values. 
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param probability_matrix A data frame containing daily symptom-onset probabilities according to the incubation period.
#' 
#' @return A single estimated infection curve using a Ridge regression.
#'
#' @examples
#' deconvolve_infection_curve_ridge(easy_symptom_onset_curve, probability_matrix)
#'
#' @export
deconvolve_infection_curve_ridge = function(symptom_onset_curve, probability_matrix)
{
  glmnet_fit = cv.glmnet(x = probability_matrix, y = symptom_onset_curve, lower.limits = 0, lambda.min = 0.00001, alpha = 0, intercept = TRUE, nlambda = 100)
  deconvolved_coefficients = coef(glmnet_fit)
  deconvolved_infection_curve_ridge = deconvolved_coefficients[-1]
  deconvolved_infection_curve_ridge = deconvolved_infection_curve_ridge/sum(deconvolved_infection_curve_ridge)
  deconvolved_infection_curve_ridge = deconvolved_infection_curve_ridge * sum(symptom_onset_curve)
  #deconvolved_infection_curve_ridge_rounded = round_infection_curve(symptom_onset_curve, deconvolved_infection_curve_ridge)
  return(deconvolved_infection_curve_ridge)
}


#' Creates a filter kernal for deconvolving a symptom-onset curve using a fourier filter. The filter takes the frequencies present in the incubation period and determines what frequency changes are needed to make the frequencies resemble a hanning window of size 20. These frequency changes are then turned back into the time domain as the filter kernel.
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param incubation_probability The name of a function that returns a data frame containing daily symptom-onset probabilities according to the distribution and paramaters given. Defaults to incubation_period_distribution_matrix.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing estimated infection curves using a Ridge regression.
#'
#' @examples
#' deconvolve_infection_curves_ridge(easy_symptom_onset_curves, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
construct_fourier_filter_kernel = function(daily_probability = incubation_period_daily_probability, 
                                           distribution = "gamma", 
                                           parms = c(9.7,5.5))
{
  daily_probability = c(daily_probability(distribution, parms), rep(0, times = 100))
  daily_probability = daily_probability[1:100]
  incubation_frequencies = fft(daily_probability)
  incubation_frequencies_max = max(Re(incubation_frequencies)^2 + Im(incubation_frequencies)^2)
  
  average = median_incubation_period(distribution, parms)
  hanning_window = c(rep(0, times = average), hanning(20), rep(0, times = 100))
  hanning_window = hanning_window[1:100]
  hanning_window_frequencies = fft(hanning_window)
  hanning_window_max = max(Re(hanning_window_frequencies)^2 + Im(hanning_window_frequencies)^2)
  
  filter_kernel = Re(fft((hanning_window_frequencies * incubation_frequencies_max)/
                           (hanning_window_max * incubation_frequencies), inverse = TRUE))
  return(filter_kernel)
}


#' Deconvolves a dataframe of symptom-onset curves by applying a filter kernel to the symptom onset curves.
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param construct_filter_kernel The name of a function that returns a fourier filter kernel according to the distribution and paramaters given. Defaults to construct_fourier_filter_kernel.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing estimated infection curves using a Fourier filter.
#'
#' @examples
#' deconvolve_infection_curves_fourier(easy_symptom_onset_curves, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
deconvolve_infection_curves_fourier = function(symptom_onset_curves, 
                                               construct_filter_kernel = construct_fourier_filter_kernel,
                                               distribution = "gamma",
                                               parms = c(9.7,5.5))
{
  fourier_filter_kernel = construct_filter_kernel(distribution = distribution, parms = parms)
  deconvolved_infection_curves_fourier = as.data.frame(symptom_onset_curves[,1])
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_fourier(symptom_onset_curve = symptom_onset_curves[,i], filter_kernel = fourier_filter_kernel)
    deconvolved_infection_curves_fourier = cbind(deconvolved_infection_curves_fourier, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_fourier)
}


#' Deconvolves a single symptom-onset curve by applying a filter kernel.
#'
#' @param symptom_onset_curve A single symptom-onset curves
#' @param filter_kernel A vector containing the filter kernel according to the distribution, mean, and standard deviation for the incubation period.
#' 
#' @return A single estimated infection curves using a Fourier filter.
#'
#' @examples
#' deconvolve_infection_curve_fourier(easy_symptom_onset_curve, filter_kernel)
#'
#' @export
deconvolve_infection_curve_fourier = function(symptom_onset_curve, filter_kernel)
{
  deconvolved_infection_curve_fourier = rep(0, times = length(symptom_onset_curve))
  for(day in 1:length(deconvolved_infection_curve_fourier))
  {
    for(filter_index in 1:length(filter_kernel))
    {
      if(day + filter_index < length(symptom_onset_curve))
      {
        deconvolved_infection_curve_fourier[day] = deconvolved_infection_curve_fourier[day] + 
          symptom_onset_curve[day + filter_index - 1] * filter_kernel[filter_index]
      }
    }
  }
  deconvolved_infection_curve_fourier[which(deconvolved_infection_curve_fourier < 0)] = 0
  deconvolved_infection_curve_fourier = deconvolved_infection_curve_fourier/sum(deconvolved_infection_curve_fourier)
  deconvolved_infection_curve_fourier = deconvolved_infection_curve_fourier * sum(symptom_onset_curve)
  #deconvolved_infection_curve_fourier_rounded = round_infection_curve(symptom_onset_curve, deconvolved_infection_curve_fourier)
  return(deconvolved_infection_curve_fourier)
}


#' Returns the discrete fourier transform of a vector.
dft <- function(z, inverse=FALSE)
{
  n <- length(z)
  if(n == 0) return(z)
  k <- 0:(n-1)
  ff <- (if(inverse) 1 else -1) * 2*pi * 1i * k/n
  vapply(1:n, function(h) sum(z * exp(ff*(h-1))), complex(1))
}


#' Returns the hermitian of a matrix.
hermitian = function(matrix)
{
  hermitian_matrix = t(matrix)
  for(i in 1:length(matrix[,1]))
  {
    for(j in 1:length(matrix[1,]))
    {
      hermitian_matrix[i,j] = complex(real = Re(hermitian_matrix[i,j]), imaginary = -1 * Im(hermitian_matrix[i,j]))
    }
  }
  return(hermitian_matrix)
}


#' Returns a matrix used for a regularized frequency filter.
#'
#' @param outbreak_length The length of the outbreak.
#' @param daily_probability A function that takes in the given distribution and paramaters and returns the daily probability of showing symptoms i days after infection.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A regularized frequency filter matrix.
#'
#' @examples
#' find_frequency_matrix(301, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
find_frequency_matrix = function(outbreak_length,
                                 daily_probability = incubation_period_daily_probability,
                                 distribution = "gamma",
                                 parms = c(9.7, 5.5))
{
  incubation_probabilities = daily_probability(distribution, parms)
  incubation_probability_length = length(incubation_probabilities)
  l_c = incubation_probability_length
  l_e = outbreak_length
  l_e.tilde = 2
  repeat
  {
    l_e.tilde = l_e.tilde * 2
    if(l_e.tilde > l_e) break
  }
  sigma_z = (l_e + l_c)/l_e.tilde
  sigma_n = sqrt(((l_c - 1)/l_e.tilde * sum(incubation_probabilities^2)))
  extension = 2
  reps=l_e.tilde - incubation_probability_length
  if(reps < 0) reps = 0
  incubation_probabilities_extended = c(incubation_probabilities, rep(0, times = reps))
  lambda = diag(dft(incubation_probabilities_extended))
  lambda_inverse = (solve(hermitian(lambda) %*% lambda + sigma_n/sigma_z * diag(ncol(lambda))) %*% hermitian(lambda))
  z_matrix = diag(l_e.tilde)[,1:l_e]
  omega = exp(-2 * pi * 1i / l_e.tilde)
  f_matrix = outer(0:(l_e.tilde - 1), 0:(l_e.tilde - 1), function(i,j) omega^(i * j))/sqrt(l_e.tilde)
  w_tilde = f_matrix %*% lambda_inverse %*% hermitian(f_matrix) %*% z_matrix
  return(w_tilde)
}


#' Deconvolves a dataframe of symptom-onset curves using a regularized frequency filter.
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves
#' @param construct_frequency_matrix The name of a function that returns a frequency matrix. Defaults to find_frequency_matrix.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing estimated infection curves using a Fourier filter.
#'
#' @examples
#' deconvolve_infection_curves_frequency(easy_symptom_onset_curves, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
deconvolve_infection_curves_frequency = function(symptom_onset_curves, 
                                                 construct_frequency_matrix = find_frequency_matrix,
                                                 distribution = "gamma",
                                                 parms = c(9.7, 5.5))
{
  frequency_matrix = find_frequency_matrix(length(symptom_onset_curves[,2]), distribution = distribution, parms = parms)
  deconvolved_infection_curves_frequency = as.data.frame(symptom_onset_curves[,1])
  incubation_length = determine_incubation_length(distribution, parms)
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_frequency(symptom_onset_curve = symptom_onset_curves[,i], matrix = frequency_matrix, incubation_length)
    deconvolved_infection_curves_frequency = cbind(deconvolved_infection_curves_frequency, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_frequency)
}


#' Deconvolves a single symptom-onset curve using a regularized frequency filter.
#'
#' @param symptom_onset_curve A single symptom-onset curve
#' @param matrix A regularized frequency filter matrix.
#' @param incubation_length The length of the incubation period
#' 
#' @return An estimated infection using a regularized frequency filter.
#'
#' @examples
#' deconvolve_infection_curve_frequency(easy_symptom_onset_curve, frequency_matrix, 60)
#'
#' @export
deconvolve_infection_curve_frequency = function(symptom_onset_curve, matrix, incubation_length)
{
  identity_matrix = diag(nrow(matrix))
  deconvolved_infection_curve_frequency = c()
  for(i in round((incubation_length/2)):min(round(length(symptom_onset_curve) + (incubation_length/2) - 1),length(matrix[,1])))
  {
    patients_infected = identity_matrix[i,] %*% matrix %*% symptom_onset_curve
    patients_infected = sqrt(Re(patients_infected)^2 + Im(patients_infected)^2)
    deconvolved_infection_curve_frequency = c(deconvolved_infection_curve_frequency, patients_infected)
  }
  deconvolved_infection_curve_frequency[which(deconvolved_infection_curve_frequency < 0)] = 0
  deconvolved_infection_curve_frequency = deconvolved_infection_curve_frequency/sum(deconvolved_infection_curve_frequency)
  deconvolved_infection_curve_frequency = deconvolved_infection_curve_frequency * sum(symptom_onset_curve)
  deconvolved_infection_curve_frequency = c(deconvolved_infection_curve_frequency, rep(0, times=length(symptom_onset_curve)-length(deconvolved_infection_curve_frequency)))
  #deconvolved_infection_curve_frequency_rounded = round_infection_curve(symptom_onset_curve, deconvolved_infection_curve_frequency)
  return(deconvolved_infection_curve_frequency)
}


#' Deconvolves a dataframe of symptom-onset curves using the Richardson-Lucy algorithm.
#'
#' @param symptom_onset_curves A dataframe containing symptom-onset curves.
#' @param estimates_matrix A dataframe containing preliminary estimates for infection curves. 
#' @param incubation_probability The name of a function that returns a data frame containing daily symptom-onset probabilities according to the distribution and paramaters given. Defaults to incubation_period_distribution_matrix.
#' @param distribution A string containing what distribution the incubation period follows, "gamma" or "lognorm". Defaults to "gamma".
#' @param parms A vector containing the mean and standard deviation for the incubation period. Defaults to c(9.7, 5.5), the mean and standard deviation for Ebola.
#' 
#' @return A dataframe containing estimated infection curves using the Richardson-Lucy algorithm.
#'
#' @examples
#' deconvolve_infection_curves_rl(easy_symptom_onset_curves, estimates_matrix = deconvolution_random, distribution = "lognorm", parms = c(15.3, 9.1))
#'
#' @export
deconvolve_infection_curves_rl = function(symptom_onset_curves,
                                          estimates_matrix,
                                          incubation_probability = incubation_period_distribution_matrix,
                                          distribution = "gamma",
                                          parms = c(9.7, 5.5))
{
  incubation_matrix = incubation_period_distribution_matrix(length(symptom_onset_curves[,1]), rl = TRUE, distribution = distribution, parms = parms)
  deconvolved_infection_curves_rl = as.data.frame(symptom_onset_curves[,1])
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_rl(symptom_onset_curve = symptom_onset_curves[,i], matrix = incubation_matrix, estimate = estimates_matrix[,i])
    deconvolved_infection_curves_rl = cbind(deconvolved_infection_curves_rl, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_rl)
}

#' A diagnostic tool used to assist in executing the Richardson-Lucy algorithm.
#'
#' @param symptom_onset_curve The original symptom-onset curve.
#' @param estimated_symptom_onset_curve An estimate of the original symptom-onset curve using a preliminary estimated infection curve. 
determine_chi2 = function(symptom_onset_curve, estimated_symptom_onset_curve)
{
  error = symptom_onset_curve - estimated_symptom_onset_curve
  error = error^2
  error_over_total = error/symptom_onset_curve
  error_over_total[is.nan(error_over_total)] = 0
  error_over_total[is.infinite(error_over_total)] = 0
  sum_difference = sum(error_over_total)
  chi2_score = sum_difference/sum(symptom_onset_curve)
  return(chi2_score)
}


#' Deconvolves a single symptom-onset curve using the Richardson-Lucy algorithm.
#'
#' @param symptom_onset_curve A single symptom-onset curve.
#' @param estimate A preliminary estimate for the infection curve. 
#' @param matrix A data frame containing daily symptom-onset probabilities.
#' 
#' @return A single estimated infection curve using the Richardson-Lucy algorithm.
#'
#' @examples
#' deconvolve_infection_curve_rl(easy_symptom_onset_curve, matrix = incubation_matrix, estimates = deconvolution_random)
#'
#' @export
deconvolve_infection_curve_rl = function(symptom_onset_curve, matrix, estimate)
{
  temporary_deconvolved_infection_curve_rl = estimate
  deconvolved_infection_curve_rl = c()
  old_chi2 = 20
  i = 1
  repeat
  {
    i = i + 1
    estimated_symptom_onset_curve = as.matrix(matrix) %*% as.matrix(temporary_deconvolved_infection_curve_rl, ncol = 1)
    new_chi2 = determine_chi2(symptom_onset_curve, estimated_symptom_onset_curve)
    new_estimate_proportions = as.vector(symptom_onset_curve)/as.vector(estimated_symptom_onset_curve)
    new_estimate_proportions[is.nan(new_estimate_proportions)] = 0
    new_estimate_proportions[is.infinite(new_estimate_proportions)] = 0
    transition_vector = t(as.matrix(new_estimate_proportions, ncol = 1)) %*% as.matrix(matrix)
    temporary_deconvolved_infection_curve_rl = as.vector(transition_vector * temporary_deconvolved_infection_curve_rl)
    ratio = new_chi2/old_chi2
    if(is.na(ratio)) break
    else if(is.infinite(ratio)) break
    else if(i > 10) break
    else if(ratio > 0.98) break
    old_chi2 = new_chi2
    deconvolved_infection_curve_rl = temporary_deconvolved_infection_curve_rl
  }
  deconvolved_infection_curve_rl[is.infinite(deconvolved_infection_curve_rl)] = 0
  deconvolved_infection_curve_rl[which(deconvolved_infection_curve_rl < 0)] = 0
  deconvolved_infection_curve_rl = deconvolved_infection_curve_rl/sum(deconvolved_infection_curve_rl)
  deconvolved_infection_curve_rl = deconvolved_infection_curve_rl * sum(symptom_onset_curve)
  #deconvolved_infection_curve_rl_rounded = round_infection_curve(symptom_onset_curve, deconvolved_infection_curve_rl)
  return(deconvolved_infection_curve_rl)
}


#' Deconvolves a dataframe of symptom-onset curves using the Average method.
#'
#' @param deconvoluted_curves A list containing dataframes of symptom-onset curves.
#' 
#' @return A dataframe containing estimated infection curves using the Average method.
#'
#' @export
deconvolve_infection_curves_average = function(deconvoluted_curves)
{
  sum_of_curves = data.frame(matrix(0, nrow = nrow(deconvoluted_curves[[1]]), ncol = ncol(deconvoluted_curves[[1]])))
  for(method in deconvoluted_curves)
  {
    sum_of_curves = sum_of_curves + method
  }
  deconvoluted_infection_curves_average = sum_of_curves/length(deconvoluted_curves)
  return(deconvoluted_infection_curves_average)
}


#' Deconvolves a single symptom-onset curve using the Average method.
#'
#' @param deconvoluted_curves A list containing deconvolved symptom-onset curves.
#' 
#' @return A single estimated infection curve using the Average method.
#'
#' @export
deconvolve_infection_curve_average = function(deconvoluted_curves)
{
  sum_of_curves = rep(0, times = length(deconvoluted_curves[1]))
  for(method in deconvoluted_curves)
  {
    sum_of_curves = sum_of_curves + method
  }
  deconvoluted_infection_curve_average = sum_of_curves/length(deconvoluted_curves)
  return(deconvoluted_infection_curve_average)
}

#' A diagnostics tool to determine the RMSE of estimated infection curves.
#'
#' @param original_curve A dataframe containing the original infection curves. 
#' @param estimated_curve A dataframe containing the estimated infection curves. 
#' 
#' @return A vector containing the RMSE for each comparison.
#'
#' @export
determine_rmse = function(original_curve, estimated_curve)
{
  difference = original_curve - estimated_curve
  squared_difference = difference^2
  mean_squared_difference = colMeans(squared_difference)
  rmse = sqrt(mean_squared_difference)
  rmse = rmse[2:length(rmse)]
  return(rmse)
}


#' A diagnostics tool to determine the correlation of estimated infection curves.
#'
#' @param original_curve A dataframe containing the original infection curves. 
#' @param estimated_curve A dataframe containing the estimated infection curves. 
#' 
#' @return A vector containing the correlation for each comparison.
#'
#' @export
determine_correlation = function(original_curve, estimated_curve){
  correlations = c()
  for(i in 2:length(original_curve[1,])){
    correlation = cor.test(original_curve[,i], estimated_curve[,i])$estimate
    correlations = c(correlations, correlation)
  }
  return(correlations)
}


#' A diagnostics tool to determine the coverage of estimated infection curves.
#'
#' @param original_curve A dataframe containing the original infection curves. 
#' @param estimated_curve A dataframe containing the estimated infection curves. 
#' 
#' @return A vector containing the coverage for each comparison.
#'
#' @export
determine_coverage = function(original_curve, estimated_curve){
  coverage = c()
  for(i in 2:length(original_curve[1,]))
  {
    number_of_samples = 10000
    residuals = original_curve[,i] - estimated_curve[,i]
    trials = resample(residuals, size = number_of_samples, invisibly.return=TRUE)
    upper = suppressMessages(confint(trials, level = 0.95, method = "quantile")[[2]])
    lower = suppressMessages(confint(trials, level = 0.95, method = "quantile")[[1]])
    trials1 = trials[trials<upper]
    trials2 = trials1[trials1>lower]
    trial_coverage = length(trials2)/number_of_samples
    coverage = c(coverage, trial_coverage)
  }
  return(coverage)
}



















#' Simulates different methods of deconvoluting easy outbreaks.
#'
#' @param n The number of outbreaks to simulate. Defaults to 1000.
#' @param mean The mean incubation period.
#' @param sd The standard deviation of the incubation period.
#' @param distribution A string containing the distribution of the incubation period ("gamma" or "lognormal"). Defaults to "gamma".
#' @param methods A vector containing each of the methods to be tested. Choices include anything in the default value of c("simple","random","ridge","rl","fourier","frequency","average")
#' @param rl0 A string containing the name of the method to be used as the initial estimate for the Richardson-Lucy method. Defaults to "random".
#' @param average_methods A vector containing boolean values for which methods to include in the average method. (simple, random, ridge, rl, fourier, frequency) Defaults to true for all values.
#' @param diagnostics A boolean value indicating whether to compute correlation, rmse, and coverage. Defaults to TRUE.
#'
#' @return A list containing all the computed pieces
#'
#' @examples
#' simulate_easy_trials(n=500, mean=9.7, sd=5.5, distribution="gamma", rl0="random", average_methods=c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
#'
#' @export
simulate_easy_trials = function(n=1000, 
                                mean,
                                sd,
                                distribution="gamma",
                                methods=c("simple","random","ridge","rl","fourier","frequency","average"), 
                                rl0="random", 
                                average_methods = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                                diagnostics=TRUE)
{
  if(rl0 %in% methods) print("thanks for being sensible")
  else rl0=methods[1]
  
  if(distribution == "gamma")
  {
    incubation_distribution=gamma_incubation_period
  }
  if(distribution == "lognorm")
  {
    incubation_distribution=lnorm_incubation_period
  }
  
  results=list()
  easy_infection_curves = generate_easy_trials(n)
  easy_symptom_onset_curves = convolute_infection_curves(easy_infection_curves, incubation_period=incubation_distribution, parms=c(mean, sd))
  results[["easy_infection_curves"]]=easy_infection_curves
  results[["easy_symptom_onset_curves"]]=easy_symptom_onset_curves
  
  if("simple" %in% methods)
  {
    easy_deconvolution_simple = deconvolve_infection_curves_simple(easy_symptom_onset_curves, distribution=distribution, parms=c(mean, sd))
    results[["easy_deconvolution_simple"]]=easy_deconvolution_simple
    if(diagnostics && is.null(easy_infection_curves))
    {
      easy_deconvolution_simple_rmse = determine_rmse(easy_infection_curves, easy_deconvolution_simple)
      easy_deconvolution_simple_correlation = determine_correlation(easy_infection_curves, easy_deconvolution_simple)
      easy_deconvolution_simple_coverage = determine_coverage(easy_infection_curves, easy_deconvolution_simple)
      results[["easy_deconvolution_simple_rmse"]]=easy_deconvolution_simple_rmse
      results[["easy_deconvolution_simple_correlation"]]=easy_deconvolution_simple_correlation
      results[["easy_deconvolution_simple_coverage"]]=easy_deconvolution_simple_coverage
    }
  }
  else 
  {
    average_methods[1]=FALSE
    easy_deconvolution_simple=0
  }
  
  if("random" %in% methods)
  {
    easy_deconvolution_random = deconvolve_infection_curves_random(easy_symptom_onset_curves, distribution = distribution, parms=c(mean, sd))
    results[["easy_deconvolution_random"]]=easy_deconvolution_random
    if(diagnostics && is.null(easy_infection_curves))
    {
      easy_deconvolution_random_rmse = determine_rmse(easy_infection_curves, easy_deconvolution_random)
      easy_deconvolution_random_correlation = determine_correlation(easy_infection_curves, easy_deconvolution_random)
      easy_deconvolution_random_coverage = determine_coverage(easy_infection_curves, easy_deconvolution_random)
      results[["easy_deconvolution_random_rmse"]]=easy_deconvolution_random_rmse
      results[["easy_deconvolution_random_correlation"]]=easy_deconvolution_random_correlation
      results[["easy_deconvolution_random_coverage"]]=easy_deconvolution_random_coverage
    }
  }
  else 
  {
    average_methods[2]=FALSE
    easy_deconvolution_random=0
  }
  
  if("ridge" %in% methods)
  {
    easy_deconvolution_ridge = deconvolve_infection_curves_ridge(easy_symptom_onset_curves, distribution = distribution, parms = c(mean, sd))
    results[["easy_deconvolution_ridge"]]=easy_deconvolution_ridge
    if(diagnostics && is.null(easy_infection_curves))
    {
      easy_deconvolution_ridge_rmse = determine_rmse(easy_infection_curves, easy_deconvolution_ridge)
      easy_deconvolution_ridge_correlation = determine_correlation(easy_infection_curves, easy_deconvolution_ridge)
      easy_deconvolution_ridge_coverage = determine_coverage(easy_infection_curves, easy_deconvolution_ridge)
      results[["easy_deconvolution_ridge_rmse"]]=easy_deconvolution_ridge_rmse
      results[["easy_deconvolution_ridge_correlation"]]=easy_deconvolution_ridge_correlation
      results[["easy_deconvolution_ridge_coverage"]]=easy_deconvolution_ridge_coverage
    }
  }
  else 
  {
    average_methods[3]=FALSE
    easy_deconvolution_ridge=0
  }
  
  if("fourier" %in% methods)
  {
    easy_deconvolution_fourier = deconvolve_infection_curves_fourier(easy_symptom_onset_curves, distribution = distribution, parms = c(mean, sd))
    results[["easy_deconvolution_fourier"]]=easy_deconvolution_fourier
    if(diagnostics && is.null(easy_infection_curves))
    {
      easy_deconvolution_fourier_rmse = determine_rmse(easy_infection_curves, easy_deconvolution_fourier)
      easy_deconvolution_fourier_correlation = determine_correlation(easy_infection_curves, easy_deconvolution_fourier)
      easy_deconvolution_fourier_coverage = determine_coverage(easy_infection_curves, easy_deconvolution_fourier)
      results[["easy_deconvolution_fourier_rmse"]]=easy_deconvolution_fourier_rmse
      results[["easy_deconvolution_fourier_correlation"]]=easy_deconvolution_fourier_correlation
      results[["easy_deconvolution_fourier_coverage"]]=easy_deconvolution_fourier_coverage
    }
  }
  else 
  {
    average_methods[5]=FALSE
    easy_deconvolution_fourier=0
  }
  
  if("frequency" %in% methods)
  {
    easy_deconvolution_frequency = deconvolve_infection_curves_frequency(easy_symptom_onset_curves, distribution = distribution, parms = c(mean, sd))
    results[["easy_deconvolution_frequency"]]=easy_deconvolution_frequency
    if(diagnostics && is.null(easy_infection_curves))
    {
      easy_deconvolution_frequency_rmse = determine_rmse(easy_infection_curves, easy_deconvolution_frequency)
      easy_deconvolution_frequency_correlation = determine_correlation(easy_infection_curves, easy_deconvolution_frequency)
      easy_deconvolution_frequency_coverage = determine_coverage(easy_infection_curves, easy_deconvolution_frequency)
      results[["easy_deconvolution_frequency_rmse"]]=easy_deconvolution_frequency_rmse
      results[["easy_deconvolution_frequency_correlation"]]=easy_deconvolution_frequency_correlation
      results[["easy_deconvolution_frequency_coverage"]]=easy_deconvolution_frequency_coverage
    }
  }
  else 
  {
    average_methods[6]=FALSE
    easy_deconvolution_frequency=0
  }
  
  if("rl" %in% methods)
  {
    if(rl0 == "simple") estimate0=easy_deconvolution_simple
    else if(rl0 == "random") estimate0=easy_deconvolution_random
    else if(rl0 == "ridge") estimate0=easy_deconvolution_ridge
    else if(rl0 == "fourier") estimate0=easy_deconvolution_fourier
    else if(rl0 == "frequency") estimate0=easy_deconvolution_frequency
    easy_deconvolution_rl = deconvolve_infection_curves_rl(easy_symptom_onset_curves, estimates_matrix = estimate0)
    results[["easy_deconvolution_rl"]]=easy_deconvolution_rl
    if(diagnostics && is.null(easy_infection_curves))
    {
      easy_deconvolution_rl_rmse = determine_rmse(easy_infection_curves, easy_deconvolution_rl)
      easy_deconvolution_rl_correlation = determine_correlation(easy_infection_curves, easy_deconvolution_rl)
      easy_deconvolution_rl_coverage = determine_coverage(easy_infection_curves, easy_deconvolution_rl)
      results[["easy_deconvolution_rl_rmse"]]=easy_deconvolution_rl_rmse
      results[["easy_deconvolution_rl_correlation"]]=easy_deconvolution_rl_correlation
      results[["easy_deconvolution_rl_coverage"]]=easy_deconvolution_rl_coverage
    }
  }
  else 
  {
    average_methods[4]=FALSE
    easy_deconvolution_rl=0
  }
  
  if("average" %in% methods)
  {
    easy_deconvolution_average_methods = list(easy_deconvolution_simple, easy_deconvolution_random, easy_deconvolution_ridge, easy_deconvolution_rl, easy_deconvolution_fourier, easy_deconvolution_frequency)
    easy_deconvolution_average = deconvolve_infection_curves_average(easy_deconvolution_average_methods[average_methods])
    results[["easy_deconvolution_average"]]=easy_deconvolution_average
    if(diagnostics && is.null(easy_infection_curves))
    {
      easy_deconvolution_average_rmse = determine_rmse(easy_infection_curves, easy_deconvolution_average)
      easy_deconvolution_average_correlation = determine_correlation(easy_infection_curves, easy_deconvolution_average)
      easy_deconvolution_average_coverage = determine_coverage(easy_infection_curves, easy_deconvolution_average)
      results[["easy_deconvolution_average_rmse"]]=easy_deconvolution_average_rmse
      results[["easy_deconvolution_average_correlation"]]=easy_deconvolution_average_correlation
      results[["easy_deconvolution_average_coverage"]]=easy_deconvolution_average_coverage
    }
  }
  
  return(results)
}