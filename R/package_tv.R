
incubation_period_daily_probability_tv = function(distribution = "gamma", parms = c(9.7, 5.5))
{
  if(min(parms) == 0){
    daily_probability = c(1)
  }
  else {
    if(distribution == "gamma")
    {
      shape = parms[1]^2/parms[2]^2
      scale = parms[2]^2/parms[1]
      
      incubation_period_length = round(qgamma(0.9999, shape = shape, scale = scale))
      probability = c(pgamma(0, shape = shape, scale = scale)
                      , pgamma(seq(0.5, incubation_period_length - 0.5, by = 1), shape = shape, scale = scale)
                      , pgamma(incubation_period_length, shape = shape, scale = scale))
    }
    else if (distribution == "lognorm")
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
  }
  return(daily_probability)
}



incubation_period_distribution_matrix_tv = function(outbreak_parms, 
                                                    incubation_probability = incubation_period_daily_probability_tv, 
                                                    rl = FALSE,
                                                    distribution = "gamma")
{
  outbreak_duration = length(outbreak_parms[,1])
  #daily_probability = incubation_probability(distribution, parms)
  #daily_probability_transition = c(daily_probability, rep(0, times = outbreak_duration))
  #daily_probability_transition = daily_probability_transition[1:outbreak_duration]
  incubation_distribution_matrix = matrix(nrow = outbreak_duration)
  for(i in 1:outbreak_duration)
  {
    daily_probability_transition = incubation_probability( outbreak_parms[i,"distribution"]
                                                         , c(outbreak_parms[i,"mean"]
                                                         , outbreak_parms[i,"sd"]))
    daily_probability_transition = c(rep(0, times = i-1)
                                     , daily_probability_transition
                                     , rep(0, times = outbreak_duration))
    daily_probability_transition = daily_probability_transition[1:outbreak_duration]
    incubation_distribution_matrix = cbind(incubation_distribution_matrix
                                           , daily_probability_transition)
    #daily_probability_transition = c(0, daily_probability_transition[1:outbreak_duration - 1])
    
  }
  if(rl)
  {
    for(i in 1:(length(incubation_probability(distribution
                                              , c(outbreak_parms[1,"mean"]
                                                  ,outbreak_parms[1,"sd"]))) - 1))
    {
      incubation_distribution_matrix[i,] = rep(0, times = length(incubation_distribution_matrix[i,]))
    }
  }
  return(incubation_distribution_matrix[,2:(outbreak_duration + 1)])
}


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


deconvolve_infection_curve_random_tv = function(casemix){
  matrix = incubation_period_distribution_matrix_tv(outbreak_parms = casemix)
  estimates = rep(0, times=length(matrix[1,]))
  
  incubation_period = generate_incubation_period
  distribution = "lognorm"
  for (distro in 1:length(matrix[1,])) {
    for (day in 1:length(matrix[,1])) {
      if(casemix$cases[day]!=0 & matrix[day, distro]!=0){
        for (patient in 1:casemix$cases[day]) {
          parms = c(casemix$mean[distro], casemix$sd[distro])
          if(sum(parms)!=0){
            estimated_infection_day = day - incubation_period(distribution = casemix$distribution[day], parms = c(casemix$mean[distro], casemix$sd[distro]))
          }
          else{
            estimated_infection_day = day
          }
          if(estimated_infection_day <= 0)
          {
            estimated_infection_day = 1
          }
          estimates[estimated_infection_day] = estimates[estimated_infection_day] + matrix[day, distro]/sum(matrix[day,])
        }
      }
    }
  }
  return(estimates)
}


deconvolve_infection_curves_ridge_tv = function(symptom_onset_curves, 
                                                incubation_probability = incubation_period_distribution_matrix_tv,
                                                distros)
{
  incubation_matrix = incubation_probability(outbreak_parms = distros)
  deconvolved_infection_curves_ridge = as.data.frame(symptom_onset_curves[,1])
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_ridge(symptom_onset_curves[,i], incubation_matrix)
    deconvolved_infection_curves_ridge = cbind(deconvolved_infection_curves_ridge, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_ridge)
}



deconvolve_infection_curves_rl_tv = function(symptom_onset_curves,
                                             estimates_matrix,
                                             incubation_probability = incubation_period_distribution_matrix_tv,
                                             distros)
{
  incubation_matrix = incubation_probability(outbreak_parms = distros, rl = TRUE)
  deconvolved_infection_curves_rl = as.data.frame(symptom_onset_curves[,1])
  for(i in 2:length(symptom_onset_curves[1,]))
  {
    deconvolved_infection_curve = deconvolve_infection_curve_rl(symptom_onset_curve = symptom_onset_curves[,i], matrix = incubation_matrix, estimate = estimates_matrix[,i])
    deconvolved_infection_curves_rl = cbind(deconvolved_infection_curves_rl, deconvolved_infection_curve)
  }
  return(deconvolved_infection_curves_rl)
}



