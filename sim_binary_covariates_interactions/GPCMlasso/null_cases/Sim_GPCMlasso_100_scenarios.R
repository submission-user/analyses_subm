####
#
####

library(GPCMlasso)
library(cmlDIFlasso2) # for simulation and dosnow
require(benchmarkme) # for RAM
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

wrap_GPCMlasso <-  function(n, I, nvars, seed, no_of_cores_used, l.lambda = 21, main.effects = TRUE){
  
  scen_1 <- sim_scenario_DIF_data_binarycov(n, I, nvars, seed,
                                  deltavalues = matrix(c(0, 0, 0, 0, 0), nvars, I ))
  alldata <- as.data.frame(cbind(scen_1$data, scen_1$covariatevalues))
  names(alldata) <- c(paste("item", 1:I, sep=""),paste("cov", 1:nvars, sep=""))
  
  form <- as.formula(paste("cbind(",paste(colnames(alldata)[1:I],collapse=","),")~."))
  
  time_GPCMlasso <- system.time(result <- GPCMlasso(form, alldata, model = "RM",
                                                   control = ctrl_GPCMlasso(l.lambda = l.lambda,
                                                                            cores = 1),
                                                   main.effects = main.effects))
  
  out <- list(result = result,
       time_GPCMlasso = time_GPCMlasso,
       scenario = scen_1,
       session = session_all(no_of_cores_used = no_of_cores_used),
       dataforGPCM = alldata)
  
  save(out, file=paste("result_files/",seed,"_GPCMlasso_null_n_",n,"_I_",I,"_nvars_",nvars,".Rdata", sep=""))
  
  return(out)
  
}


extract_obs_alpha_beta_power_GPCMlasso <- function (x, resultlist) 
{ 
  nitems <- dim(resultlist[[x]]$scenario$deltavalues)[2]
  nvars <- dim(resultlist[[x]]$scenario$deltavalues)[1]
  nmaineffects <- ifelse(resultlist[[x]]$result$main.effects,nvars,0)
  dif.mat <- matrix(resultlist[[x]]$result$coefficients[which.min(resultlist[[x]]$result$BIC),(nitems+nmaineffects+1):(nitems+nmaineffects+nitems*nvars)], nvars, nitems)
  
  tb1 <- table(abs(dif.mat) > 
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
# GPCMlasso results (21 lambdas)
#
################



repetitions <- 100

pars_help <- expand.grid(n = c(50,200,500),
                    I = c(5,8,15,25), 
                    nvars = c(1,2,6,11))
pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))
seeds <- 1999:(1999+nrow(pars)-1)
no_of_cores = 31
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
                                                .packages = c("GPCMlasso","cmlDIFlasso2"),
                                                .options.snow = opts) %dopar% {
                                                  wrap_GPCMlasso(n, I, nvars, seed,
                                                                no_of_cores_used = no_of_cores,
                                                                l.lambda = nlambda_wrap,
                                                                main.effects = TRUE)
                                                })
close(pb)
stopCluster(cl)


# function to replace all designX bei NULL to save memory space and fir results into memory
for (i in 1:4800){
  res_foreach[[i]]$result$design_list$designX <- NULL
}

# res_foreach <- lapply(1:10000, function(x) res_foreach[[x]][[1]])

save(res_foreach, file = "results_GPCMlasso_100nullcases_100x100.Rdata")
pushover("Loading null cases done!")


t_and_scen <- lapply(1:nrow(pars),FUN = function(x) unlist(list(scenario = pars[x,], 
                                                                time = res_foreach[[x]]$time_GPCMlasso,
                                                                extract_obs_alpha_beta_power_GPCMlasso(x, res_foreach))))
# unlist(t_and_scen[[1]])
df_t_and_scen <- as.data.frame(do.call(rbind,t_and_scen))


df_t_and_scen_GPCMlasso_null <- df_t_and_scen
save(df_t_and_scen_GPCMlasso_null,
     file="df_t_GPCM_diflasso_null_100each.Rdata")

###
#
# add alpha and power
#
#####

library(dplyr)
tbl1 <- df_t_and_scen %>% 
  group_by(scenario.n,scenario.I, scenario.nvars) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))

df_all_reduced_GPCMlasso_null_100x100_samples <- tbl1
save(df_all_reduced_GPCMlasso_null_100x100_samples,
     file = "df_all_reduced_GPCMlasso_null_100x100_samples.Rdata")


