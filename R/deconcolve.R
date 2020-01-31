deconvolve_single_curve <- function(curve, parms) {

  # simple
  simple <- deconvolve_infection_curve_simple(curve, parms[1])
  
  # random
  random <- deconvolve_infection_curve_random(
    curve,
    generate_incubation_period,
    trials = 100,
    distribution = "gamma",
    parms = parms
  )

  # ridge
  incubation_matrix <- incubation_period_distribution_matrix(length(curve), 
                                                             distribution = distribution, 
                                                             parms = parms
                                                             )
  ridge <- deconvolve_infection_curve_ridge(curve, incubation_matrix)
  
  # RL
  incubation_matrix <- incubation_period_distribution_matrix(length(curve),
                                                             distribution = distribution,
                                                             parms = parms,
                                                             rl = TRUE
                                                             )
  rl <- deconvolve_infection_curve_rl(curve, incubation_matrix, random)
  
  # fourier filter - doesn't do well with shorter outbreaks. 
  # testing on simulated outbreaks still has strong correlation and low RMSE
  fourier_filter_kernel <- construct_fourier_filter_kernel(distribution = distribution, 
                                                           parms = parms
                                                           )
  fourier <- deconvolve_infection_curve_fourier(curve, fourier_filter_kernel)
  
  # frequency filter
  frequency_matrix <- find_frequency_matrix(length(curve), 
                                            distribution = distribution, 
                                            parms = parms
                                            )
  matrix <- frequency_matrix
  incubation_length <- determine_incubation_length(distribution, parms)
  frequency <- deconvolve_infection_curve_frequency(curve, frequency_matrix, incubation_length)
  
  #list of estimates to average
  deconvolutions = list(simple
                        , random
                        #, ridge
                        , rl
                        , frequency
                        , c(frequency, 0, 0, 0))
  
  average <- deconvolve_infection_curve_average(deconvolutions)
  
  curves <- list(curve=curve, simple=simple, random=random, ridge=ridge, rl=rl, fourier=fourier, frequency=frequency, average=average)
  return(curves)
}