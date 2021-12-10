###
#
# Math as example for all cases
# variable selection: cml and RTs
# item selection: cml and DIFlasso
# all interactions: cml, DIFlasso, GPCM
# 
###


library("cmlDIFlasso2")
library("DIFlasso")
library("psychotree")
library("DIFtree")
library("GPCMlasso")
# library("eRm")
library("pushoverr")
source("function_cml_lasso_win_nlm_optimized.R")

# load("Example_data/datasets.RData") # so packages with data don´t have to be installed


prepare_data_DIFlasso <- function(data, covariates){
  ifelse(is.data.frame(data)== FALSE, 
         dfd <-  as.data.frame(data), 
         dfd <-  data )
  ifelse(is.data.frame(covariates)== FALSE, 
         dfc <- as.data.frame(covariates), 
         dfc <- covariates )
  
  dfc_num <- sapply(dfc,as.numeric)
  dfc_stand <- scale(dfc_num)
  return(list(df_data = dfd,
              df_cov = dfc_stand))
}

cores_to_use = 1



###
# examples sirt-Package
###

####
# get data from packages (not necessary if data loaded via .Rdata-file)

data(data.pisaMath, package = "sirt")


# prepare data
math <- prepare_data_DIFlasso(data.pisaMath$data[,6:16], data.pisaMath$data[,3:5])

pushover("Examples started!")

###
# selecting variables
###

# cmlDIFlasso grouplasso nlm

time_math_cmlDIFlasso_variables <- system.time(
  result_math_cmlDIFlasso_variables <- 
    cml_lasso_win_nlm_optimized_foreach(as.matrix(math[[1]]),math[[2]],
                                        DIF_type = "variables",20, tolerance_df = .0001,
                                        no_of_cores_to_use = cores_to_use))

# rt

math_rt <- data.pisaMath$data[,3:5]
math_rt$resp <- as.matrix(data.pisaMath$data[,6:16])

time_math_rt_variables <- system.time(
  result_math_raschtree_variables <- 
    raschtree(resp~female + hisei + migra, data = math_rt))



###
# selecting items
###

# cmlDIFlasso grouplasso nlm

time_math_cmlDIFlasso_items <- system.time(
  result_math_cmlDIFlasso_items <- 
    cml_lasso_win_nlm_optimized_foreach(as.matrix(math[[1]]),math[[2]],
                                        DIF_type = "items",20, tolerance_df = .0001,
                                        no_of_cores_to_use = cores_to_use))

# DIFlasso grouplasso

time_math_DIFlasso_items <- system.time(
  result_math_DIFlasso_items <- 
    DIFlasso(Y = math$df_data,
             X = as.data.frame(math$df_cov),
             trace = TRUE,
             grouped = TRUE))

###
# selecting interactions
###

# cmlDIFlasso L1 lbfgs

time_math_cmlDIFlasso_all_int_lbfgs <- system.time(
  result_math_cmlDIFlasso <- 
    cml_lasso_win_lbfgs_optimized_foreach(data = as.matrix(math[[1]]),
                                          covariates = math[[2]],
                                          DIF_type = "all.interactions",
                                          nlambdas = 20,
                                          no_of_cores = cores_to_use))

time_math_cmlDIFlasso_all_int_nlm <- system.time(
  result_math_cmlDIFlasso_all_int_nlm <- 
    cml_lasso_win_nlm_optimized_foreach(as.matrix(math[[1]]),math[[2]],
                                        DIF_type = "all.interactions",20, tolerance_df = .0001,
                                        no_of_cores_to_use = cores_to_use))



# DIFlasso L1

time_math_DIFlasso_all_int <- system.time(
  result_math_DIFlasso_all_int <- 
    DIFlasso(Y = math$df_data,
             X = as.data.frame(math$df_cov),
             trace = TRUE,
             grouped = FALSE))
# GPCMlasso

df_math_all <- as.data.frame(cbind(math[[1]],math[[2]])) 
names(df_math_all)[1:11] <- paste("item", 1:11, sep="")

form <- as.formula(paste("cbind(",paste(colnames(df_math_all)[1:11],collapse=","),")~."))

time_math_GPCMlasso_all_int <- system.time(
  result_math_GPCMlasso_all_int <- 
    GPCMlasso(form, df_math_all, model = "RM",
              control = ctrl_GPCMlasso(l.lambda = 20,
                                       log.lambda = FALSE,
                                       cores = 1),
              main.effects = TRUE))


# IFTs

time_math_DIFtree_all_int <- system.time(
  result_math_DIFtree_all_int<- 
    DIFtree(Y=data.pisaMath$data[,6:16],X=data.pisaMath$data[,3:5], 
            model="Logistic",type="udif",
            alpha=0.05,nperm=1000,trace=TRUE))

print(result_math_DIFtree_all_int)

# save.image(file='all_example_results_math.RData')

pushover("Examples done!")

(ls())[grep("time",ls())]
lapply( (ls())[grep("time_math_cml",ls())], get)
lapply( (ls())[grep("time_math_DIFl",ls())], get)
lapply( (ls())[grep("time_math_DIFt",ls())], get)
lapply( (ls())[grep("time_math_rt",ls())], get)


