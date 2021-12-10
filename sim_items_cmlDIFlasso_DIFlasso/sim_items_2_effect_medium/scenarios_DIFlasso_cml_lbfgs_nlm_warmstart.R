#####
####
#
# DIFlasso-scenarios
# sim 100 per condition
# compare DIFlasso, cmlDIFlasso nlm, 
# weak effects (half of the large effects)
# medium .75 times effect.
# (no nlm warm start because of unstable results)
#
#####
# load all additional functions first
# lbfgs
source("cml_lasso_win_lbfgs_optimized_temp.R")  # tolerance_df = 1e-8
source("wrap_cml_lasso_lbfgs_optimized_temp.R")  # tolerance_df = 1e-8
# nlm
source("function_cml_lasso_win_nlm_optimized.R")
source("wrap_cml_lasso_nlm_optimized.R")# grouplasso
# write wrapper
# nlm_warm_start
source("function_cml_lasso_win_nlm_optimized_warm_start.R")
source("wrap_cml_lasso_nlm_warm_start.R")

library(doSNOW) # load first for progress bar
library(cmlDIFlasso2)
library(pushoverr)


##### 
# make simulation scenarios
#####

repetitions <- 100
n <- c(250,rep(500,4))
set.seed(1999)
I <- c(20,20,20,40,20)
beta_help <- rnorm(40)
beta_vals <- lapply(I, function(x) beta_help[1:x])

# make delta values (with minus for easiness-parameters)
delta_help <- -.75*(matrix(c( -0.8,0.6,0,0,0.8, 
                         0,0.8,-0.7,0,0.7,
                         0.6,0,0.8,-0.8,0,
                         0,0,0.8,0.7,-0.5),5,4))
delta_vals <- lapply(I, function(x) matrix(0,5,x))
## all setting item 1 to 4 with DIF as above. Setting 2 and 3 also ites 5 to 8 as above
for (i in 1:5) {delta_vals[[i]][1:5,1:4] <-  delta_help  }
for (i in 3:4) {delta_vals[[i]][1:5,5:8] <-  delta_help  }

covariate_help <- matrix(rnorm(500*5),500,5)
covariate_help[,1:2] <- (covariate_help[,1:2] > 0)*1 
covariate_vals <- lapply(n, function(x) covariate_help[1:x,])

theta_help <- rnorm(500)
theta_vals <- lapply(n, function(x) theta_help[1:x])
theta_vals[[5]] <- theta_vals[[5]] + covariate_vals[[5]][,1] ## correlated first variable

####
# set simulation and estimation procedure properties

# repetitions <- 2 for testing
seeds <- 2499:(2499+5*repetitions-1)
no_of_cores = 59
nlambda_wrap = 21


####
# start parallel execution separately for DIFlasso, cmlDIFlasso nlm and nlm warm start 
############



#######
# nlm
######

cl<-makeCluster(no_of_cores,type="SOCK") 
registerDoSNOW(cl) 

# make progress bar:
pb <- txtProgressBar(max = 5*repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

time_full <- system.time(res_foreach <-  
                           foreach (i=1:5,n=n, I=I,
                                    betas = beta_vals,
                                    deltas = delta_vals,
                                    covariates = covariate_vals,
                                    thetas = theta_vals,
                                    .combine = 'c') %:%
                           foreach(seed = seeds[(i-1)*repetitions + (1:repetitions)],
                                   .packages = c("cmlDIFlasso2","psychotools","lbfgs"),
                                   .options.snow = opts,
                                   .inorder=TRUE) %dopar% {
                                     wrap_cml_lasso_nlm_optimized(n, I, 5, seed, nlambda_wrap,
                                                                          betas, deltas, thetas, covariates,
                                                                          no_of_cores_used = no_of_cores,
                                                                          DIF_type = "items")})

close(pb)
stopCluster(cl)

pushover("nlm done!")

do.call(c,lapply(1:500,function(x)dim(res_foreach[[x]]$scenario$data)))
do.call(c,lapply(1:500,function(x)res_foreach[[x]]$time_lbfgs[3]))
do.call(c,lapply(1:500,function(x)res_foreach[[x]]$scenario$seed))

res_nlm_medium_effects <- res_foreach
save(res_foreach,file = "results_nlm_medium_effects.Rdata")

time_full_nlm <- time_full

# #######
# # nlm_warm_start
# ######
# 
# cl<-makeCluster(no_of_cores,type="SOCK")  
# registerDoSNOW(cl)
# 
# # make progress bar:
# pb <- txtProgressBar(max = 5*repetitions, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# 
# time_full <- system.time(res_foreach <-  
#                            foreach (i=1:5,n=n, I=I,
#                                     betas = beta_vals,
#                                     deltas = delta_vals,
#                                     covariates = covariate_vals,
#                                     thetas = theta_vals,
#                                     .combine = 'c') %:%
#                            foreach(seed = seeds[(i-1)*repetitions + (1:repetitions)],
#                                    .packages = c("cmlDIFlasso2","psychotools","lbfgs"),
#                                    .options.snow = opts,
#                                    .inorder=FALSE) %dopar% {
#                                      wrap_cml_lasso_nlm_warm_start(n, I, 5, seed, nlambda_wrap,
#                                                                   betas, deltas, thetas, covariates,
#                                                                   no_of_cores_used = no_of_cores,
#                                                                   DIF_type = "items")})
# 
# close(pb)
# stopCluster(cl)
# 
# pushover("nlm warm start done!")
# 
# do.call(c,lapply(1:500,function(x)dim(res_foreach[[x]]$scenario$data)))
# do.call(c,lapply(1:500,function(x)res_foreach[[x]]$time_lbfgs[3]))
# do.call(c,lapply(1:500,function(x)res_foreach[[x]]$scenario$seed))
# 
# res_nlm_warm_start_medium_effects <- res_foreach
# save(res_foreach,file = "results_nlm_warm_start_medium_effects.Rdata")
# 
# time_full_nlm_warm_start <- time_full


#######
######
## DIFlasso
######

library(DIFlasso)
library(cmlDIFlasso2) # for simulation and dosnow
require(benchmarkme) # for RAM
# source("session_all.R") ## for infos about R and computer etc.


########
## function for wrapping DIFlasso
##
## and calculate alphaobs, powerobs from wrap_DIFlasso resultlist
######

wrap_DIFlasso <-  function(n, I, nvars, seed, nlambdas = 21, 
                           betavalues = rep(0, I), 
                           deltavalues = matrix(c(-0.5, 0, 0, 0, 0.5), 
                                                nvars, I),
                           thetavalues = rep(0,n),
                           covariatevalues = matrix(rnorm(n *nvars), n, nvars), 
                           no_of_cores_used = "unknown"){
  
  set.seed(seed)
  scen_1 <- sim_scenario_DIF_data(n, I, nvars, seed, betavalues = betavalues, 
                                  deltavalues = deltavalues, thetavalues = thetavalues, 
                                  covariatevalues = covariatevalues)
  time_DIFlasso <- system.time(result <- DIFlasso(Y = as.data.frame(scen_1$data),
                                                  X = as.data.frame(scale(scen_1$covariatevalues)),
                                                  l.lambda = nlambdas,
                                                  grouped = TRUE))
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



######
## DIFlasso ordered:
##########
##
# simulation
##

cl<-makeCluster(no_of_cores,type="SOCK")  
registerDoSNOW(cl)

# make progress bar:
pb <- txtProgressBar(max = 5*repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

time_full <- system.time(res_foreach <-  
                           foreach (i=1:5,n=n, I=I,
                                    betas = beta_vals,
                                    deltas = delta_vals,
                                    covariates = covariate_vals,
                                    thetas = theta_vals,
                                    .combine = 'c') %:%
                           foreach(seed = seeds[(i-1)*repetitions + (1:repetitions)],
                                   .packages = c("DIFlasso","cmlDIFlasso2"),
                                   .options.snow = opts,
                                   .inorder=TRUE) %dopar% {
                                     wrap_DIFlasso(n, I, 5, seed, nlambda_wrap,
                                                   betas, deltas, thetas, covariates,
                                                   no_of_cores_used = no_of_cores)})

close(pb)
stopCluster(cl)

pushover("DIFlasso done!")

do.call(c,lapply(1:500,function(x)dim(res_foreach[[x]]$scenario$data)))
do.call(c,lapply(1:500,function(x)res_foreach[[x]]$time_DIFlasso[3]))
do.call(c,lapply(1:500,function(x)res_foreach[[x]]$scenario$seed))

res_DIFlasso_grouped_ordered_medium_effects <- res_foreach
save(res_foreach,file = "results_DIFlasso_grouped_ordered_medium_effects.Rdata")

time_full_DIFlasso_grouped_ordered <- time_full

