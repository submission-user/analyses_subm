####
#
# first run with this file file. Other runs for missing cases with "...-missing cases.R"
#
####

library(DIFlasso)
library(cmlDIFlasso2) # for simulation and dosnow
require(benchmarkme) # for RAM
source("session_all.R") ## for infos about R and computer etc.
library(pushoverr)

########
## function for wrapping DIFlasso
##
## and calculate alphaobs, powerobs from wrap_DIFlasso resultlist
######

wrap_DIFlasso <-  function(n, I, nvars, seed, no_of_cores_used, l.lambda = 21){
  
  scen_1 <- sim_scenario_DIF_data(n, I, nvars, seed,
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
# DIFlasso results (21 lambdas )
#
################

# n = c(50,200,500,1000)
# I = c(5,8,15,25,40) 
# nvars = c(1,2,5,10, 20) 


repetitions <- 100

pars_help <- expand.grid(n = c(50,200,500,1000),
                    I = c(5,8,15,25,40), 
                    nvars = c(1,2,6,11, 21))
pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))
seeds <- 1999:(1999+nrow(pars)-1)
no_of_cores = 53
nlambda_wrap = 21


library(doSNOW)
cl<-makeCluster(no_of_cores,type="SOCK", outfile="")
registerDoSNOW(cl)


########
# Split1
############################
# make progress bar:
pb <- txtProgressBar(max = nrow(pars), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

pars1 <- pars[1:2500,]
seeds1 <- seeds[1:2500]

### calculate

time_full1 <- system.time(res_foreach1 <- foreach(n = pars1$n, 
                                                I = pars1$I, 
                                                nvars = pars1$nvars, 
                                                seed = seeds1,
                                                .packages = c("DIFlasso","cmlDIFlasso2"),
                                                .options.snow = opts) %dopar% {
                                                  wrap_DIFlasso(n, I, nvars, seed,
                                                                no_of_cores_used = no_of_cores,
                                                                l.lambda = nlambda_wrap)
                                                })
close(pb)
###################################
pushover("Simulation 1 done!")
save(res_foreach1, file = "results_DIFlasso_100nullcases_100t_1.Rdata")

########
# Split2
############################
# make progress bar:
pb <- txtProgressBar(max = nrow(pars), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

pars2 <- pars[2500 + (1:2500), ]
seeds2 <- seeds[2500 + (1:2500)]

### calculate

time_full2 <- system.time(res_foreach2 <- foreach(n = pars2$n, 
                                                I = pars2$I, 
                                                nvars = pars2$nvars, 
                                                seed = seeds2,
                                                .packages = c("DIFlasso","cmlDIFlasso2"),
                                                .options.snow = opts) %dopar% {
                                                  wrap_DIFlasso(n, I, nvars, seed,
                                                                no_of_cores_used = no_of_cores,
                                                                l.lambda = nlambda_wrap)
                                                })
close(pb)
###################################
pushover("Simulation 2 done!")
save(res_foreach2, file = "results_DIFlasso_100nullcases_100t_2.Rdata")


########
# Split3
############################
# make progress bar:
pb <- txtProgressBar(max = nrow(pars), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

pars3 <- pars[5000 + (1:2500), ]
seeds3 <- seeds[5000 + (1:2500)]

### calculate

time_full3 <- system.time(res_foreach3 <- foreach(n = pars3$n, 
                                                  I = pars3$I, 
                                                  nvars = pars3$nvars, 
                                                  seed = seeds3,
                                                  .packages = c("DIFlasso","cmlDIFlasso2"),
                                                  .options.snow = opts) %dopar% {
                                                    wrap_DIFlasso(n, I, nvars, seed,
                                                                  no_of_cores_used = no_of_cores,
                                                                  l.lambda = nlambda_wrap)
                                                  })
close(pb)
###################################
pushover("Simulation 3 done!")
save(res_foreach3, file = "results_DIFlasso_100nullcases_100t_3.Rdata")


########
# Split4
############################
# make progress bar:
pb <- txtProgressBar(max = nrow(pars), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

pars4 <- pars[7500 + (1:2500), ]
seeds4 <- seeds[7500 + (1:2500)]

### calculate

time_full4 <- system.time(res_foreach4 <- foreach(n = pars4$n, 
                                                  I = pars4$I, 
                                                  nvars = pars4$nvars, 
                                                  seed = seeds4,
                                                  .packages = c("DIFlasso","cmlDIFlasso2"),
                                                  .options.snow = opts) %dopar% {
                                                    wrap_DIFlasso(n, I, nvars, seed,
                                                                  no_of_cores_used = no_of_cores,
                                                                  l.lambda = nlambda_wrap)
                                                  })
close(pb)
###################################
save(res_foreach4, file = "results_DIFlasso_100nullcases_100t_4.Rdata")

stopCluster(cl)

pushover("Simulation 4 done!")

## combine all result lists
res_foreach <- c(res_foreach1, res_foreach2, res_foreach3, res_foreach4)

save(res_foreach, file = "results_DIFlasso_100nullcases_100t_all.Rdata")

t_and_scen <- lapply(1:nrow(pars),FUN = function(x) unlist(list(scenario = pars[x,], 
                                                                time = res_foreach[[x]]$time_DIFlasso,
                                                                extract_obs_alpha_beta_power_DIFlasso(x, res_foreach))))

df_t_and_scen <- as.data.frame(do.call(rbind,t_and_scen))

#  plot times vs. n, I and nvars. colours changing by n, I or nvars
library(GGally)
ggpairs(df_t_and_scen[,c(1:3,6)], aes(colour = as.factor(scenario.n), alpha = 0.4))
ggpairs(df_t_and_scen[,c(1:3,6)], aes(colour = as.factor(scenario.I), alpha = 0.4))
ggpairs(df_t_and_scen[,c(1:3,6)], aes(colour = as.factor(scenario.nvars), alpha = 0.4))

# influencing computation times
#
summary(lm(time.elapsed~scenario.n+scenario.I+scenario.nvars, data = df_t_and_scen))


## interactions
summary(lm(time.elapsed~scenario.n*scenario.I*scenario.nvars, data = df_t_and_scen))


# polynomials:
summary(lm(time.elapsed ~ poly(scenario.n, 2, raw=TRUE)+
             poly(scenario.I, 2, raw=TRUE)+
             poly(scenario.nvars, 2, raw=TRUE), data = df_t_and_scen))

summary(lm(time.elapsed ~ poly(scenario.n, 3, raw=TRUE)+
             poly(scenario.I, 3, raw=TRUE)+
             poly(scenario.nvars, 3, raw=TRUE), data = df_t_and_scen))




df_t_and_scen_DIFlasso_null <- df_t_and_scen
save(df_t_and_scen_DIFlasso_null,
     file="df_t_null_diflasso_100each.Rdata")

###
#
# add alpha and power
#
#####


extract_obs_alpha_beta_power_DIFlasso(21,res_foreach)

View(df_t_and_scen[,c(1:3,6,9,11)])
plot(density(df_t_and_scen$obs.alpha[1:100]),xlim=c(-.1,1))

# power:
plot(density(df_t_and_scen$obs.pwr),xlim=c(-.1,1))

###
#
###

library(dplyr)
tbl1 <- df_t_and_scen %>% 
  group_by(scenario.n,scenario.I, scenario.nvars) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))

df_all_reduced_DIFlasso_null_50samples <- tbl1
save(df_all_reduced_DIFlasso_null_50samples,
     file = "df_all_reduced_DIFlasso_null_50samples.Rdata")

####
#
# inspect certain results, regarding detection of DIF
#
#####


# plot all 100 BIC trajectories
for (i in 100:1) plot(x=res_foreach[[i]]$result$lambda,y=res_foreach[[i]]$result$BIC, 
                      main = paste("n =", df_t_and_scen$scenario.n[i],
                                   ", I =", df_t_and_scen$scenario.I[i],
                                   ", nvars =", df_t_and_scen$scenario.nvars[i],
                                   ",       i =",i))

## check how many of the 100 scenarios deliver no selected parameters
# via delta estimates

cases_noselection_deltahat <- unlist(lapply(1:length(res_foreach),
                                            function(x) sum(res_foreach[[x]]$result$dif.mat) == 0))
sum(cases_noselection_deltahat)

# cases with identified DIF
which(!cases_noselection_deltahat)  


# check if there is pattern for selection according to n, I and nvars
plot(1:100 %in% which(!cases_noselection_deltahat)%in% 1:100)

dput(unlist(lapply(1:100,function(x)which.min(res_foreach[[x]]$result$BIC))))

