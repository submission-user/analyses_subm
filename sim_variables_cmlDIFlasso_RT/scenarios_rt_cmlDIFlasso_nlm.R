#####
####
#
# variable selection
# sim 100 per condition
# compare cmlDIFlasso and rt
#
#####
# load all additional functions first

source("function_cml_lasso_win_nlm_optimized.R")
source("wrap_cml_lasso_nlm_optimized.R")# grouplasso
# write wrapper
# raschtree
source("wrap_raschtree.R")
raschtree_extract_variables <- function(resultobj) sort(unique(do.call(c,lapply(1:length(resultobj) , function(x) resultobj[[x]]$node$split$varid))))
# nlm_warm_start
# source("function_cml_lasso_win_nlm_optimized_warm_start.R")
# source("wrap_cml_lasso_nlm_warm_start.R")

library(doSNOW) # load first for progress bar
library(cmlDIFlasso2)
library(pushoverr)
library(psychotree)


##### 
# make simulation scenarios
#####

##
# 18times n=250 then 18*n=500
# 2 times (n) 6 (scenarios) times  each effect = 1.5, 6x .75, 6x 0
# 2 times 6 times scenarios
# (1) 1 binary covariate
# (2) binary and continuous covariate (effect on binary)
# (3) binary and continuous covariate (effect on continuous)
# (4) binary and continuous covariate (effect on both)
# (5) binary and 4 continuous covariate (effect on first two covariates)
# (6) binary and 4 continuous covariate (effect on first two covariates) with two groups differing in ability along the binary covariate (+- half the effect)
# effects if present are are on the second item for covariate 1 and on the fifth item for covariate 2

# added two more scenarios later, as described in submission (7 and 8)

# scenarios 3 4 5 (1,2,3 in sub) and 7 and 8 (4 in sub)  used in submission.


repetitions <- 100
n <- rep(c(250,500),each=24)
set.seed(1999)
I <- rep(20,48)
ncovs <- rep(c(1,2,2,2,5,5,5,5),6)
beta_help <- rnorm(20)
beta_vals <- lapply(1:48,function(x)beta_help)

# make delta values (with minus for easiness-parameters)
delta_vals <- lapply(ncovs, function(x) matrix(0,x,20))
## all setting item 1 to 4 with DIF as above. Setting 2 and 3 also items 5 to 8 as above
##
for (i in (1:length(delta_vals))[-(rep(0:5,each=3)*8+c(3,7,8))]) {delta_vals[[i]][1,2] = -1.5}
for (i in (1:length(delta_vals))[-(rep(0:5,each=4)*8+c(1:2,7:8))]) {delta_vals[[i]][2,5] = -1.5}
for (i in rep(0:5,each=2)*8+7 ) {delta_vals[[i]][1:2,] = matrix(c(-.2,0,.2),2,20)}
for (i in rep(0:5,each=2)*8+8 ) {delta_vals[[i]][1:2,] = matrix(c(-.05,0,.05),2,20)}
for (i in c(1:8,25:32)+8) {delta_vals[[i]] <- delta_vals[[i]]*.5}
for (i in c(1:8,25:32)+16) {delta_vals[[i]] <- delta_vals[[i]]*0}


covariate_help <- matrix(rnorm(500*5),500,5)
covariate_help[,1] <- (covariate_help[,1] > 0)*1 
covariate_vals <- lapply(1:48, function(x) covariate_help[1:n[x],1:ncovs[x],drop=FALSE])

theta_help <- rnorm(500)
theta_vals <- lapply(n, function(x) theta_help[1:x])

# for (i in 1:48) theta_vals[[i]] <- theta_vals[[i]] + (covariate_vals[[i]][,1] - 0.5) * rep(c(c(rep(1.5,6),.2,.05), 
#                                                                                            .5*c(rep(1.5,6),.2,.05),
#                                                                                           0*c(rep(1.5,6),.2,.05)),2)[i] ## correlated first variable

####
# set simulation and estimation procedure properties

# repetitions <- 2 for testing
seeds <- 1999:(1999+48*repetitions-1)
no_of_cores = 59
nlambda_wrap = 21




#######
# nlm
######

cl<-makeCluster(no_of_cores,type="SOCK")  # Corecluster erzeugen, für snow
registerDoSNOW(cl)   # Cluster starten

# make progress bar:
pb <- txtProgressBar(max = 48*repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

time_full <- system.time(res_foreach <-  
                           foreach (i=1:48,n=n, I=I,
                                    betas = beta_vals,
                                    deltas = delta_vals,
                                    nvars = ncovs,
                                    covariates = covariate_vals,
                                    thetas = theta_vals,
                                    .combine = 'c') %:%
                           foreach(seed = seeds[(i-1)*repetitions + (1:repetitions)],
                                   .packages = c("cmlDIFlasso2","psychotools","lbfgs"),
                                   .options.snow = opts,
                                   .inorder=FALSE) %dopar% {
                                     wrap_cml_lasso_nlm_optimized(n, I, nvars, seed, nlambda_wrap,
                                                                          betas, deltas, thetas, covariates,
                                                                          no_of_cores_used = no_of_cores,
                                                                          DIF_type = "variables")})

close(pb)
stopCluster(cl)

pushover("nlm done!")

hist(do.call(c,lapply(1:4800,function(x)dim(res_foreach[[x]]$scenario$data))))
hist(do.call(c,lapply(1:4800,function(x)res_foreach[[x]]$time_lbfgs[3])))
table(do.call(c,lapply(1:4800,function(x)res_foreach[[x]]$scenario$seed)) - 1999:(1999+4800-1))

seed_list <- do.call(c,lapply(1:length(res_foreach), function(x) res_foreach[[x]]$scenario$seed))
res_foreach<- lapply(order(seed_list), function(x) res_foreach[[x]])

res_nlm <- res_foreach
save(res_foreach,file = "results_nlm.Rdata")
pushover("nlm save done!")

time_full_nlm <- time_full



#######################
#######################
# rt
#######################
#######################

cl<-makeCluster(no_of_cores,type="SOCK")  # Corecluster for snow
registerDoSNOW(cl)   # Cluster start

# make progress bar:
pb <- txtProgressBar(max = 48*repetitions, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

time_full <- system.time(res_foreach <-  
                           foreach (i=1:48,n=n, I=I, ncovs = ncovs,
                                    betas = beta_vals,
                                    deltas = delta_vals,
                                    covariates = covariate_vals,
                                    thetas = theta_vals,
                                    .combine = 'c') %:%
                           foreach(seed = seeds[(i-1)*repetitions + (1:repetitions)],
                                   .packages = c("cmlDIFlasso2","psychotree"),
                                   .export = c("wrap_raschtree"),
                                   .options.snow = opts,
                                   .inorder=TRUE) %dopar% {
                                     wrap_raschtree(n, I, ncovs, seed,
                                                    betas, deltas, thetas, covariates,
                                                    no_of_cores_used = no_of_cores,
                                                    number_of_conditions = 8)})

close(pb)
stopCluster(cl)

pushover("raschtree done!")

hist(do.call(c,lapply(1:4800,function(x)dim(res_foreach[[x]]$scenario$data))))
hist(do.call(c,lapply(1:4800,function(x)res_foreach[[x]]$time_lbfgs[3])))
table(do.call(c,lapply(1:4800,function(x)res_foreach[[x]]$scenario$seed)) - 1999:(1999+4800-1))

# # if not in order
# seed_list <- do.call(c,lapply(1:length(res_foreach), function(x) res_foreach[[x]]$scenario$seed))
# res_foreach<- lapply(order(seed_list), function(x) res_foreach[[x]])
res_raschtree <- res_foreach
save(res_foreach,file = "results_raschtree.Rdata")

time_full_raschtree <- time_full

