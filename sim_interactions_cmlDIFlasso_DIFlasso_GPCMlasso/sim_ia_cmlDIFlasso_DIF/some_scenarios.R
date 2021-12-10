#####
#
#
# check accuracy and computing times of cml_lasso_win_lbfgs
#
#
####



#####
####
#
# cml lasso wrap lbfgs
#
#####

library(doSNOW) # load first for progress bar
library(cmlDIFlasso2)


repetitions <- 100

pars_help <- expand.grid(n = c(50,200,500,1000),
                    I = c(5,8,15,25,40), 
                    nvars = c(1,2,6,11, 21))
pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))
seeds <- 1999:(1999+nrow(pars)-1)
no_of_cores = 53
nlambda_wrap = 21




cl<-makeCluster(no_of_cores,type="SOCK")  
registerDoSNOW(cl)  

## for parallel package:
# cl<-makeCluster(no_of_cores) 
# registerDoParallel(cl)

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
                      wrap_cml_lasso_lbfgs_optimized(n, I, nvars, seed, nlambda_wrap,
                                           no_of_cores_used = no_of_cores)
                    })
close(pb)
stopCluster(cl)

save(res_foreach,file = "results_cmlDIFlasso_opt_100x100scenarios.Rdata")


table(abs(res_foreach[[1]]$result$result_picked$delta_estimates)>.0000001,
      abs(res_foreach[[1]]$scenario$deltavalues) != 0)

t_and_scen <- lapply(1:nrow(pars),FUN = function(x) unlist(list(scenario = pars[x,],
                                                                time = res_foreach[[x]]$time_lbfgs, 
                                                                extract_obs_alpha_beta_power(x,res_foreach))))

# dataframe with comp times
df_t_and_scen <- as.data.frame(do.call(rbind,t_and_scen))


for(i in 1:4) plot(density(df_t_and_scen$time.elapsed[df_t_and_scen$scenario.n==c(50,200,500,1000)[i]]),
                   main = paste("times for n =",c(50,200,500,1000)[i]), xlim = c(0,350),ylim = c(0, .003))

for(i in 1:5) plot(density(df_t_and_scen$time.elapsed[df_t_and_scen$scenario.I==c(5,8,15,25,40)[i]]),
                   main = paste("times for nitems =",c(5,8,15,25,40)[i]), xlim = c(0,350),ylim = c(0, .0015))

for(i in 1:5) plot(density(df_t_and_scen$time.elapsed[df_t_and_scen$scenario.nvars==c(1,2,6,11, 21)[i]]),
                   main = paste("times for nvars =",c(1,2,6,11, 21)[i]), xlim = c(0,350),ylim = c(0, .0008))

df_t_and_scen[which(df_t_and_scen$time.elapsed>5000),]
### all times under 300s


checkdeltaresults <- function(zeileinpars, round = 7) list(round(res_foreach[[zeileinpars]]$result$result_picked$delta_estimates,round),
                                                           res_foreach[[zeileinpars]]$scenario$deltavalues)

pars[3,]
checkdeltaresults(48,3)


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

df_t_and_scen$time.elapsed[df_t_and_scen$scenario.n==1000]
max(df_t_and_scen$time.elapsed)
min(df_t_and_scen$time.elapsed[df_t_and_scen$scenario.n==1000])


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

tbl1
View(tbl1)
plot(tbl1$avgtime)
plot(tbl1$se_time)
spgroesse = 200
## Average power by number of deltas for n 200
plot((tbl1$scenario.nvars*tbl1$scenario.I)[tbl1$scenario.n==spgroesse], 
     tbl1$avgalpha[tbl1$scenario.n==spgroesse])

# combined for n
tbl_n <- df_t_and_scen_cmlDIFlasso_opt %>% 
  group_by(scenario.n) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))
# combined for I
tbl_I <- df_t_and_scen_cmlDIFlasso_opt %>% 
  group_by(scenario.I) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))
# combined for nvars
tbl_nvars <- df_t_and_scen_cmlDIFlasso_opt %>% 
  group_by(scenario.nvars) %>% 
  summarise(avgtime = mean(time.elapsed),
            se_time = sd(time.elapsed)/99, 
            avgalpha = mean(obs.alpha),
            avgbeta = mean(obs.beta),
            avgpwr = mean(obs.pwr))

# means for each combination of I and nvars for  n = 1000

tbl_n1000 <- summarise(group_by(df_t_and_scen_cmlDIFlasso_opt[df_t_and_scen_cmlDIFlasso_opt$scenario.n==1000,], 
                                scenario.I, scenario.nvars), avgtime = mean(time.elapsed))
View(tbl_n1000)


# inspect first 100 results
View(df_t_and_scen[1:100,c(1:3,6,9,11)])
plot(density(df_t_and_scen$obs.alpha[1:100]),xlim=c(-.1,1))

### cmlDIFlasso power:
plot(density(df_t_and_scen$obs.pwr),xlim=c(-.1,1))


