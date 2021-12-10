####
#
####

library(GPCMlasso)
library(cmlDIFlasso2) # for simulation and dosnow
require(benchmarkme) # for RAM
# source("session_all.R") ## for infos about R and computer etc.
library(pushoverr)

########
## function for wrapping DIFlasso
##
## and calculate alphaobs, powerobs from wrap_DIFlasso resultlist
######

wrap_GPCMlasso <-  function(n, I, nvars, seed, no_of_cores_used, l.lambda = 21, main.effects = TRUE){
  
  scen_1 <- sim_scenario_DIF_data(n, I, nvars, seed,
                                  deltavalues = matrix(c(-0.5, 0, 0, 0, 0.5), nvars, I ))
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
  
  save(out, file=paste("result_files/",seed,"_GPCMlasso_DIF_n_",n,"_I_",I,"_nvars_",nvars,".Rdata", sep=""))
  
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

# n = c(50,200,500,1000)
# I = c(5,8,15,25,40) 
# nvars = c(1,2,5,10, 20) 


repetitions <- 100

pars_help <- expand.grid(n = c(50,200,500,1000),
                    I = c(5,8,15,25,40), 
                    nvars = c(1,2,6,11, 21))
pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))
seeds <- 1999:(1999+nrow(pars)-1)
no_of_cores = 31
nlambda_wrap = 21


library(doSNOW)
cl<-makeCluster(no_of_cores,type="SOCK", outfile="") 
registerDoSNOW(cl)   

# make progress bar:
pb <- txtProgressBar(max = nrow(pars), style = 3)
opts <- list(progress = progress)

progress <- function(n) setTxtProgressBar(pb, n)
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

pushover("Simulation done!")

## combine all result lists into res_foreach
## prior to this: copy all files into folder "result_files_all"

all_result_files_names <- list.files("result_files_all/")
# sort by number:
all_result_files_names_sorted <- all_result_files_names[order(as.numeric(sapply(all_result_files_names, function(x) strsplit(x,"_")[[1]][1])))]
# make full path, attaching folder name in front
all_result_files_names_sorted_fullpath <- sapply(all_result_files_names_sorted, function(x) paste("result_files_all/",x,sep=""))

# function to not only load object into env, but return it (and not write to env)
load_and_return <- function(x) {
  load(x)
  return(out)}

load_and_return2 <- function(x) {
  load(x)
  out$result$design_list$designX <- NULL
  out$result$design_list$acoefs <- NULL
  out$result$formula <- c()
  return(out)}



res_foreach <- lapply(1:10000,function(x) {print(paste(x,"done"))
  lapply(all_result_files_names_sorted_fullpath[x], 
         load_and_return2)})

# # function to replace all designX bei NULL to save memory space and fir results into memory
# for (i in 1:25){
#   for (j in 1:100) {
#     res_foreach_help1[[i]][[j]]$result$design_list$designX <- NULL
#   }
# }

res_foreach <- lapply(1:10000, function(x) res_foreach[[x]][[1]])

save(res_foreach, file = "results_GPCMlasso_100DIFcases_100x100.Rdata")
pushover("Loading done!")


# before running this return to initial pars-matrix (run first commands to produce pars)
t_and_scen <- lapply(1:nrow(pars),FUN = function(x) unlist(list(scenario = pars[x,], 
                                                                time = res_foreach[[x]]$time_GPCMlasso,
                                                                extract_obs_alpha_beta_power_GPCMlasso(x, res_foreach))))

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




df_t_and_scen_GPCMlasso_DIF <- df_t_and_scen
save(df_t_and_scen_GPCMlasso_DIF,
     file="df_t_GPCMlasso_DIF_100each.Rdata")

###
#
# add alpha and power
#
#####


extract_obs_alpha_beta_power_GPCMlasso(21,res_foreach)

View(df_t_and_scen[,c(1:3,6,9,11)])
plot(density(df_t_and_scen$obs.alpha[1:100]),xlim=c(-.1,1))

# power
plot(density(na.omit(df_t_and_scen$obs.pwr)),xlim=c(-.1,1))

# # no estimates (NaN) for scenarios:
# which(is.na(df_t_and_scen$obs.pwr))
# # [1] 1481 2281 6681 6985 7581 8581
# # inspect further

excluded_cases <- which(is.na(df_t_and_scen$obs.pwr))

###
#
###

library(dplyr)
tbl1 <- df_t_and_scen[-excluded_cases, ] %>% 
  group_by(scenario.n,scenario.I, scenario.nvars) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))

df_all_reduced_GPCMlasso_DIF_100x100_samples <- tbl1
save(df_all_reduced_GPCMlasso_DIF_100x100_samples,
     file = "df_all_reduced_GPCMlasso_DIF_100x100_samples.Rdata")


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

