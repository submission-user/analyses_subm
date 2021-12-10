wrap_cml_lasso_lbfgs_optimized_temp <- function (n, I, nvars, seed, nlambdas, 
                                                 betavalues = rep(0, I), 
                                                 deltavalues = matrix(c(-0.5, 0, 0, 0, 0.5), 
                                                                      nvars, I),
                                                 thetavalues = rep(0,n),
                                                 covariatevalues = matrix(rnorm(n *nvars), n, nvars), 
                                                 no_of_cores_used = "unknown",
                                                 tolerance_df = 1e-08) 
{ 
  set.seed(seed)
  scen_1 <- sim_scenario_DIF_data(n, I, nvars, seed, betavalues = betavalues, 
                                  deltavalues = deltavalues, thetavalues = thetavalues, 
                                  covariatevalues = covariatevalues)
  time_lbfgs <- system.time(result <- 
                              cml_lasso_win_lbfgs_optimized_temp(data = scen_1$data, 
                                                            covariates = scen_1$covariatevalues, nlambdas = nlambdas,
                                                            tolerance_df = tolerance_df))
  list(result = result, time_lbfgs = time_lbfgs, scenario = scen_1, 
       lambdas = result$lambdas, session = session_all(no_of_cores_used = no_of_cores_used))
}
