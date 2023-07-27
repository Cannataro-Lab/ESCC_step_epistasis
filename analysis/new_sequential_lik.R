#' Development of new sequential selection method
#' Besides, sequential_mut_prop, arguments are the same as the old sequential model likelihood function,
#' and they'll be supplied automatically by ces_variant(). Use custom_lik_args to supply
#' sequential_mut_prop.
#' 
#' @param rates_tumors_with named vector of site-specific mutation rates for all tumors
#'   with variant
#' @param rates_tumors_without named vector of site-specific mutation rates for all
#'   eligible tumors without variant
#' @param sample_index data.table with columns Unique_Patient_Identifier, group_name, group_index
#' @param sequential_mut_prop proportion of mutations estimated to occur during each sequential state
sequential_lik_dev <- function(rates_tumors_with, rates_tumors_without, sample_index, sequential_mut_prop) {
  stages_tumors_with = sample_index[names(rates_tumors_with), group_index, on = "Unique_Patient_Identifier"]
  stages_tumors_without = sample_index[names(rates_tumors_without), group_index, on = "Unique_Patient_Identifier"]
  num_pars = sample_index[, uniqueN(group_index)]
  
  fn = function(gamma) {
    gamma = unname(gamma)
    sums = cumsum(gamma)
    gamma_sums = sums[stages_tumors_without]
    
    # samples without mutation
    sum_log_lik = -1 * sum(mapply(
        function(rate, stage) {
          flux = gamma[1:stage] * rate * sequential_mut_prop[1:stage]
          return(sum(flux))
        }, rates_tumors_without, stages_tumors_without))

    # samples with mutation
    if(length(rates_tumors_with) > 0) {
      sum_log_lik = sum_log_lik + sum(mapply(
        function(rate, stage) {
          # stage-specific likelihoods of mutation
          lik_no_mutation = exp(-1 * gamma * rate * sequential_mut_prop)
          lik_mutation = 1 - lik_no_mutation
          
          cum_lik_no_mut = c(1, cumprod(lik_no_mutation))
          return(log(sum(cum_lik_no_mut[1:stage] * lik_mutation[1:stage])))
        }, rates_tumors_with, stages_tumors_with))
    }
    
    # in case it tried all the max at once
    if(!is.finite(sum_log_lik)){
      return(-1e200)
    }
    return(-1 * sum_log_lik) # actually returning negative loglikelihood
  }
  
  # Set default values for all parameters, which ces_variant will use to set starting values of optimization
  formals(fn)[["gamma"]] = rep.int(1000, num_pars)
  
  # Optimization tool, bbmle::mle, requires that vector of parameters to optimize have named elements
  group_names = unique(sample_index[, .(group_index, group_name)], by = 'group_index')[order(group_index), group_name]
  bbmle::parnames(fn) = paste0("si_", group_names)
  return(fn)
}
