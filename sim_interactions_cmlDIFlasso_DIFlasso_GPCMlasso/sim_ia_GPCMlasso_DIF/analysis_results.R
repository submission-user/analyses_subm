####
#
# comp.times (21 lambdas each)
#
#####

#  plot times vs. n, I and nvars. colours changing by n, I or nvars
library(GGally)
ggpairs(df_t_and_scen[,c(3,6:8)], aes(colour = as.factor(scenario.n), alpha = 0.4))
ggpairs(df_t_and_scen[,c(3,6:8)], aes(colour = as.factor(scenario.I), alpha = 0.4))
ggpairs(df_t_and_scen[,c(3,6:8)], aes(colour = as.factor(scenario.nvars), alpha = 0.4))

# influencing computation times
#
summary(lm(time.elapsed~scenario.n+scenario.I+scenario.nvars, data = df_t_and_scen))

## interactions
summary(lm(time.elapsed~scenario.n*scenario.I*scenario.nvars, data = df_t_and_scen))

# View(tbl_n1000)
# View(tbl_n)
# View(tbl_I)
# View(tbl_nvars)
View(tbl1)

#####
#
# results
#
#####



deltapars <- lapply((1:length(res_foreach))#[-excluded_cases]
                    ,
                    function(x) matrix(res_foreach[[x]]$result$coefficients[which.min(res_foreach[[x]]$result$BIC),
                                                                            (pars$I[x]+pars$nvars[x]+1):(pars$I[x]+pars$nvars[x]+pars$I[x]*pars$nvars[x])],
                                       pars$nvars[x],
                                       pars$I[x]))

## boxplots for certain scenario
make_box <- function(n1,I1,ncovs1) {
  deltas_picked <- deltapars[which(apply(pars,1,function(x) all(x == c(n1,I1,ncovs1))))]
  deltas_picked2 <- lapply(deltas_picked, function(x) as.vector(x))
  deltas_matrix <- do.call(rbind,deltas_picked2)
  boxplot(deltas_matrix, xlab = "deltas", ylab = "estimates", las = 1,
          main = paste("n=",n1,"I=",I1,"ncovs=",ncovs1))
}
####
# not necessary for noDIF case:
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
    par_est <- as.vector(deltas_matrix)
    unlist(list(MSE = mean((par_est)^2, na.rm = TRUE),
                var = var(par_est, na.rm = TRUE),
                 bias = mean(par_est, na.rm = TRUE)
                ))
  }

(test <- calc_MSE_var_bias_exp(50,8,6))
# scenarios:
# already calculated
# repetitions <- 100
# pars_help <- expand.grid(n = c(50,200,500,1000),
#                          I = c(5,8,15,25,40),
#                          nvars = c(1,2,6,11, 21))
# pars <- do.call(rbind,replicate(repetitions,pars_help, simplify = FALSE))

mse_etc <- lapply(1:100,function(x)calc_MSE_var_bias_exp(pars_help[x,1],
                                                         pars_help[x,2],
                                                         pars_help[x,3]))
library(dplyr)
df_mse_etc <- do.call(rbind, mse_etc)
df_mse_etc2 <- cbind(pars_help, df_mse_etc)
df_all <- cbind(tbl1,dplyr::arrange(df_mse_etc2,n,I,nvars)[,-(1:3)])
# combine MSE, var and exp


summary(lm(MSE~scenario.n+scenario.I+scenario.nvars, data = df_all))

boxplot(as.list(df_all[,9:11]))

## combine positive and negative deltas
# keep only n,I,K, time, T.I err., Mean, Variance, Power, Mean, Variance

names(df_all)
df_all_reduced <- df_all[, -c(5,7)]

names(df_all_reduced) <- c("n","I","K", "Comp.time","Type I error", "Power", "MSE", "Variance", "Bias") 

df_all_reduced_DIF <- df_all_reduced
save(df_all_reduced_DIF,
     file = "df_all_reducedDIF_100x100DIFcases_GPCMlasso.Rdata")


# plot alpha by K, separate lines per I, sep. plots per n
# separately for time, power DIF/noDIF mean and variance
# two plots each (except means), second with adjusted scales per plot
library(ggplot2)
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Comp.time, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n,nrow = 2) + labs(x= "Number of covariates", y="Computation time (DIF-samples)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Comp.time, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free',  nrow = 2) + labs(x= "number of covariates", y="Computation time (DIF-samples)",color ="Number \n of items")


# DIF-delta
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=get('Type I error'), group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Type I error (DIF-samples)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=get('Type I error'), group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Type I error (DIF-samples)",color ="Number \n of items")

ggplot(data = df_all_reduced, aes(x=as.factor(K), y=get('Power'), group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Power (DIF-samples)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=get('Power'), group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Power (DIF-samples)",color ="Number \n of items")


ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Bias, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Mean of estimates (DIF-samples)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Bias, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Mean of estimates (DIF-samples)",color ="Number \n of items")

ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Variance, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, nrow =2) + labs(x= "number of covariates", y="Variance of estimates (DIF-samples)",color ="Number \n of items")
ggplot(data = df_all_reduced, aes(x=as.factor(K), y=Variance, group = I, color = as.factor(I))) + geom_line(size =.8) +
  geom_point() + facet_wrap(~ n, scales = 'free', nrow =2) + labs(x= "number of covariates", y="Variance of estimates (DIF-samples)",color ="Number \n of items")


