wrap_raschtree <- function (n, I, nvars, seed, 
                            betavalues = rep(0, I), 
                            deltavalues = matrix(c(-0.5, 0, 0, 0, 0.5), 
                                                 nvars, I),
                            thetavalues = rep(0,n),
                            covariatevalues = matrix(rnorm(n *nvars), n, nvars), 
                            no_of_cores_used = "unknown",
                            number_of_conditions = 6) 
{ 
  set.seed(seed)
  scen_1 <- sim_scenario_DIF_data(n, I, nvars, seed, betavalues = betavalues, 
                                  deltavalues = deltavalues, thetavalues = thetavalues, 
                                  covariatevalues = covariatevalues)
  rt_data <- as.data.frame(scen_1$covariatevalues)
  rt_data$resp <- scen_1$data
  
  time_lbfgs <- system.time(result <- raschtree(resp ~ . , data = rt_data))
  
  raschtree_extract_variables <- function(resultobj) sort(unique(do.call(c,lapply(1:length(resultobj) , function(x) resultobj[[x]]$node$split$varid))))
  selected_variables <- rep(FALSE,I)

    if ( !is.null(raschtree_extract_variables(result)))
          {selected_variables[raschtree_extract_variables(result) - 1] = TRUE}   # (first variable is resp)
  
  scen_2 <- list(n=nrow(scen_1$data), I=ncol(scen_1$data), nvars=nrow(deltavalues),
                 effect_strength = unique(sort(-abs(deltavalues)))[1],
                 effect_covariates = which(colSums(abs(deltavalues)) > 0),
                 condition = (seed - 1998)%%number_of_conditions)
  
  list(result = result, selected_variables = selected_variables, time_lbfgs = time_lbfgs, scenario = c(scen_1,scen_2), 
       session = session_all(no_of_cores_used = no_of_cores_used))
}
