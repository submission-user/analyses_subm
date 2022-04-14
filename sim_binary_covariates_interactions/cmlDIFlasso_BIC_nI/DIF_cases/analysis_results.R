####
#
# comp.times (21 lambdas each)
#
#####

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


#####
#
# results
#
#####

# extract par_estimates
deltapars <- lapply(1:length(res_foreach),function(x) unlist(res_foreach[[x]]$result$result_picked$delta_estimates))

## boxplots for certain scenario
make_box <- function(n1,I1,ncovs1) {
  deltas_picked <- deltapars[which(apply(pars,1,function(x) all(x == c(n1,I1,ncovs1))))]
  deltas_picked2 <- lapply(deltas_picked, function(x) as.vector(x))
  deltas_matrix <- do.call(rbind,deltas_picked2)
  boxplot(deltas_matrix, xlab = "deltas", ylab = "estimates", las = 1,
          main = paste("n=",n1,"I=",I1,"ncovs=",ncovs1))
}
make_box_parscombined <- 
  function(n1,I1,ncovs1) {
    deltas_picked <- deltapars[which(apply(pars,1,function(x) all(x == c(n1,I1,ncovs1))))]
    deltas_picked2 <- lapply(deltas_picked, function(x) as.vector(x))
    deltas_matrix <- do.call(rbind,deltas_picked2)
    boxplot(as.vector(deltas_matrix[,c(TRUE,FALSE,FALSE,FALSE,FALSE)]),
            as.vector(deltas_matrix[,c(FALSE,TRUE,TRUE,TRUE,FALSE)]),
            as.vector(deltas_matrix[,c(FALSE,FALSE,FALSE,FALSE,TRUE)]),
            xlab = "delta negative, positive, null", ylab = "estimates", las = 1,
            main = paste("n=",n1,"I=",I1,"ncovs=",ncovs1))
  }
make_box_parscombined(500,8,6)

# calc_MSE_var_bias_exp for paameters
calc_MSE_var_bias_exp <- 
  function(n1,I1,ncovs1) {
    deltas_picked <- deltapars[which(apply(pars,1,function(x) all(x == c(n1,I1,ncovs1))))]
    deltas_picked2 <- lapply(deltas_picked, function(x) as.vector(x))
    deltas_matrix <- do.call(rbind,deltas_picked2)
    par_est <- list(deltas_neg = as.vector(deltas_matrix[,c(TRUE,FALSE,FALSE,FALSE,FALSE)]),
            deltas_null = as.vector(deltas_matrix[,c(FALSE,TRUE,TRUE,TRUE,FALSE)]),
            deltas_pos =as.vector(deltas_matrix[,c(FALSE,FALSE,FALSE,FALSE,TRUE)]))
    unlist(list(MSEneg = mean((par_est[[1]] + .5)^2),
         MSEnull = mean((par_est[[2]])^2),
         MSEpos = mean((par_est[[3]] - .5)^2),
         varneg = var(par_est[[1]]),
         varnull = var(par_est[[2]]),
         varpos = var(par_est[[3]]),
         biasneg = mean(par_est[[1]]) + .5,
         biasnull = mean(par_est[[2]]),
         biaspos = mean(par_est[[3]]) - .5,
         expneg = mean(par_est[[1]]),
         expnull = mean(par_est[[2]]),
         exppos = mean(par_est[[3]])
         ))
  }

test <- calc_MSE_var_bias_exp(500,8,6)
# scenarios:
repetitions <- 100
pars_help <- expand.grid(n = c(50,200,500),
                         I = c(5,8,15,25), 
                         nvars = c(1,2,6,11))
pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))

mse_etc <- lapply(1:48,function(x)calc_MSE_var_bias_exp(pars_help[x,1],
                                                         pars_help[x,2],
                                                         pars_help[x,3]))
library(dplyr)
df_mse_etc <- do.call(rbind, mse_etc)
df_mse_etc2 <- cbind(pars_help, df_mse_etc)
df_all <- cbind(tbl1,dplyr::arrange(df_mse_etc2,n,I,nvars)[,-(1:3)])
# combine MSE, var and exp for po


summary(lm(MSEpos~scenario.n+scenario.I+scenario.nvars, data = df_all))
boxplot(as.list(df_all[,9:20]))

## combine positive and negative deltas
# keep only n,I,K, time, T.I err., Mean, Variance, Power, Mean, Variance
library(dplyr)
names(df_all)
df_all_helper <- df_all[,c(1,2,3,4,6,8,12,13,14,18,19,20)]
df_all_helper$Mean_DIF_items <- (df_all_helper$exppos - df_all_helper$expneg)/2
df_all_helper$Variance_DIF_items <-   (df_all_helper$varpos + df_all_helper$varneg)/2
df_all_helper$Mean_noDIF_items <- df_all_helper$expnull
df_all_helper$Variance_noDIF_items <- df_all_helper$varnull

df_all_reduced <- df_all_helper[,c(1:6,13:16)]
names(df_all_reduced) <- c("n","I","K", "Comp.time","Type I error", "Power", "Mean_DIF", "Variance_DIF", "Mean_noDIF","Variance_noDIF") 

df_all_reduced_DIF <- df_all_reduced
save(df_all_reduced_DIF,
     file = "df_all_reducedDIF_100x100DIFcases_cmldif.Rdata")

# plot alpha by K, separate lines per I, sep. plots per n
# separately for time, power DIF/noDIF mean and variance
# two plots each (except means), second with adjusted scales per plot
library(ggplot2)
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Comp.time, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n,nrow = 2) + labs(x= "Number of covariates", y="Computation time",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Comp.time, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free',  nrow = 2) + labs(x= "number of covariates", y="Computation time",color ="Number \n of items")


# DIF-deltas
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Power, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Power",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Power, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Power",color ="Number \n of items")

ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Mean_DIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Mean of estimates (Items with DIF)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Mean_DIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Mean of estimates (Items with DIF)",color ="Number \n of items")

ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Variance_DIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Variance of estimates (Items with DIF)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Variance_DIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n,scales = 'free', nrow =2) + labs(x= "number of covariates", y="Variance of estimates (Items with DIF)",color ="Number \n of items")


# noDIF-delta
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=get('Type I error'), group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Type I error",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=get('Type I error'), group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Type I error",color ="Number \n of items")

ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Mean_noDIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Mean of estimates (Items without DIF)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Mean_noDIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Mean of estimates (Items without DIF)",color ="Number \n of items")

ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Variance_noDIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Variance of estimates (Items without DIF)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Variance_noDIF, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Variance of estimates (Items without DIF)",color ="Number \n of items")




