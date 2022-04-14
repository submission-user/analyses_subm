####
#

####

library(DIFlasso)
library(cmlDIFlasso2) # for simulation and dosnow
require(benchmarkme) # for RAM
source("session_all.R") ## for infos about R and computer etc.
library(pushoverr)

########
## function for simulating binary covariates
##
######

sim_scenario_DIF_data_binarycov <- function (npersons, nitems, ncovs, 
                                             seed, betavalues = rep(0,nitems), 
                                             deltavalues = matrix(c(-0.5, 0, 0, 0, 0.5), ncovs, nitems), 
                                             thetavalues = rep(0, npersons), 
                                             covariatevalues = matrix(rnorm(npersons * ncovs), npersons, ncovs)) {
  
  covariatevalues <- (covariatevalues > 0) *1
  set.seed(seed)
  return(list(data = simulate_DIF_data(betavalues, deltavalues, 
                                       covariatevalues, thetavalues), betavalues = betavalues, 
              deltavalues = deltavalues, thetavalues = thetavalues, 
              covariatevalues = covariatevalues, seed = seed))
}


########
## function for wrapping DIFlasso
##
## and calculate alphaobs, powerobs from wrap_DIFlasso resultlist
######

wrap_DIFlasso <-  function(n, I, nvars, seed, no_of_cores_used, l.lambda = 21){
  
  scen_1 <- sim_scenario_DIF_data_binarycov(n, I, nvars, seed,
                                            deltavalues = matrix(c(0, 0, 0, 0, 0), nvars, I ))
  
  time_DIFlasso <- system.time(result <- DIFlasso(Y = as.data.frame(scen_1$data),
                                                  X = as.data.frame(scale(scen_1$covariatevalues)),
                                                  l.lambda = l.lambda,
                                                  grouped = FALSE))
  
  out <- list(result = result[-11],
       time_DIFlasso = time_DIFlasso,
       scenario = scen_1,
       session = session_all(no_of_cores_used = no_of_cores_used))
  
  save(out, file=paste("result_files/",seed,"_DIFlasso_null_n_",n,"_I_",I,"_nvars_",nvars,".Rdata", sep=""))
  
  return(out)
  
}


extract_obs_alpha_beta_power_DIFlasso <- function (x, resultlist) 
{
  tb1 <- table(abs(resultlist[[x]]$result$dif.mat) > 
                 1e-07, abs(resultlist[[x]]$scenario$deltavalues) != 0)
  tb <- matrix(0, 2, 2)
  if ("FALSE" %in% rownames(tb1) & "FALSE" %in% 
      colnames(tb1)) {
    tb[1, 1] <- tb1["FALSE", "FALSE"]
  }
  if ("FALSE" %in% rownames(tb1) & "TRUE" %in% 
      colnames(tb1)) {
    tb[1, 2] <- tb1["FALSE", "TRUE"]
  }
  if ("TRUE" %in% rownames(tb1) & "FALSE" %in% 
      colnames(tb1)) {
    tb[2, 1] <- tb1["TRUE", "FALSE"]
  }
  if ("TRUE" %in% rownames(tb1) & "TRUE" %in% colnames(tb1)) {
    tb[2, 2] <- tb1["TRUE", "TRUE"]
  }
  list(obs.alpha = round(prop.table(tb, 2)[2, 1], 3), obs.beta = round(prop.table(tb, 
                                                                                  2)[1, 2], 3), obs.pwr = round(prop.table(tb, 2)[2, 2], 
                                                                                                                3))
}


################
#
# DIFlasso results (40 lambdas to better compare with 2step approach )
#
################



repetitions <- 100

pars_help <- expand.grid(n = c(50,200,500),
                    I = c(5,8,15,25), 
                    nvars = c(1,2,6,11))
pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))
seeds <- 1999:(1999+nrow(pars)-1)
no_of_cores = 53
nlambda_wrap = 21


library(doSNOW)
cl<-makeCluster(no_of_cores,type="SOCK", outfile="")  # Cluster 
registerDoSNOW(cl)   # Cluster start

# make progress bar:
pb <- txtProgressBar(max = nrow(pars), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

### calculate

time_full <- system.time(res_foreach <- foreach(n = pars$n, 
                                                I = pars$I, 
                                                nvars =pars$nvars, 
                                                seed = seeds,
                                                .packages = c("DIFlasso","cmlDIFlasso2"),
                                                .options.snow = opts) %dopar% {
                                                  wrap_DIFlasso(n, I, nvars, seed,
                                                                no_of_cores_used = no_of_cores,
                                                                l.lambda = nlambda_wrap)
                                                })
close(pb)
stopCluster(cl)

pushover("Simulation done!")

save(res_foreach, file = "results_DIFlasso_100nullcases_scenarios2_100t.Rdata")

t_and_scen <- lapply(1:nrow(pars),FUN = function(x) unlist(list(scenario = pars[x,], 
                                                                time = res_foreach[[x]]$time_DIFlasso,
                                                                extract_obs_alpha_beta_power_DIFlasso(x, res_foreach))))

df_t_and_scen <- as.data.frame(do.call(rbind,t_and_scen))

df_t_and_scen_DIFlasso_null <- df_t_and_scen
save(df_t_and_scen_DIFlasso_null,
     file="df_t_null_diflasso_100each.Rdata")

###
#
# add alpha and power
#
#####


extract_obs_alpha_beta_power_DIFlasso(21,res_foreach)


library(dplyr)
tbl1 <- df_t_and_scen %>% 
  group_by(scenario.n,scenario.I, scenario.nvars) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))

df_all_reduced_DIFlasso_null_100samples <- tbl1
save(df_all_reduced_DIFlasso_null_100samples,
     file = "df_all_reduced_DIFlasso_null_100samples.Rdata")
