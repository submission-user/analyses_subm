#####
#
#
# check accuracy and computing times of cml_lasso_win_lbfgs
#
#
####



library(doSNOW) # load first for progress bar
library(cmlDIFlasso2)



########
## wrap function for simulating binary covariates
##
######

wrap_cml_lasso_lbfgs_optimized_binarycov <- function (n, I, nvars, seed, nlambdas, betavalues = rep(0, I), 
                                                      deltavalues = matrix(c(-0.5, 0, 0, 0, 0.5), nvars, I), 
                                                      thetavalues = rep(0, n), covariatevalues = matrix(rnorm(n * nvars), n, nvars), 
                                                      no_of_cores_used = "unknown") {
  
  covariatevalues <- (covariatevalues > 0) *1
  scen_1 <- sim_scenario_DIF_data(n, I, nvars, seed, betavalues = betavalues, 
                                  deltavalues = deltavalues, thetavalues = thetavalues, 
                                  covariatevalues = covariatevalues)
  time_lbfgs <- system.time(result <- cml_lasso_win_lbfgs_optimized(data = scen_1$data, 
                                                                    covariates = scen_1$covariatevalues, nlambdas = nlambdas))
  list(result = result, time_lbfgs = time_lbfgs, scenario = scen_1, 
       lambdas = result$lambdas, session = session_all(no_of_cores_used = no_of_cores_used))
}








repetitions <- 100

pars_help <- expand.grid(n = c(50,200,500),
                    I = c(5,8,15,25), 
                    nvars = c(1,2,6,11))
pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))
seeds <- 1999:(1999+nrow(pars)-1)
no_of_cores = 53
nlambda_wrap = 21


# library(doSNOW)

cl<-makeCluster(no_of_cores,type="SOCK")  # Cluster 
registerDoSNOW(cl)   # Cluster 


# make progress bar:
pb <- txtProgressBar(max = nrow(pars), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

time_full <- system.time(res_foreach <- foreach(n = pars$n, 
                    I = pars$I, 
                    nvars =pars$nvars, 
                    seed = seeds,
                    .packages = c("cmlDIFlasso2","psychotools","lbfgs"),
                    .options.snow = opts) %dopar% {
                      wrap_cml_lasso_lbfgs_optimized_binarycov(n, I, nvars, seed, nlambda_wrap,
                                           no_of_cores_used = no_of_cores)
                    })
close(pb)
stopCluster(cl)

save(res_foreach,file = "results_cmlDIFlasso_opt_100x100scenarios.Rdata")

t_and_scen <- lapply(1:nrow(pars),FUN = function(x) unlist(list(scenario = pars[x,],
                                                                time = res_foreach[[x]]$time_lbfgs, 
                                                                extract_obs_alpha_beta_power(x,res_foreach))))

# dataframe with comp times of all scenarios
df_t_and_scen <- as.data.frame(do.call(rbind,t_and_scen))

df_t_and_scen_cmlDIFlasso_opt <- df_t_and_scen
save(df_t_and_scen_cmlDIFlasso_opt,
     file="df_t_100x100cases_cmldif_opt_tolerancedf_0.Rdata")
######
#
# inspecting the data
#
#######
library(dplyr)

# means for all cases:
tbl1 <- df_t_and_scen_cmlDIFlasso_opt %>% 
  group_by(scenario.n,scenario.I, scenario.nvars) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))

df_all_reduced_cmlDIFlasso_DIF_100samples <- tbl1
save(df_all_reduced_cmlDIFlasso_DIF_100samples,
     file = "df_all_reduced_cmlDIFlasso_DIF_100samples.Rdata")

