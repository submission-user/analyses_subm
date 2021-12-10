## Wrapper for RM_DIF_nlm
## selection of items or variables possible (DIF_type = "items" or "variables")

cml_lasso_win_nlm_optimized <- function (data, covariates, DIF_type = "all.interactions", nlambdas = 21, 
                                         tolerance_df = 0.0001, ...) 
{
  npersons <- nrow(data)
  nitems <- ncol(data)
  ncov <- ncol(covariates)
  time_unres <- system.time(res_unrestricted <- RM_DIF_lbfgs(data = data, 
                                                             covariates = covariates, lambda = 0, invisible = 1))
  time_res <- system.time(beta_estimates_restricted <- -itempar(raschmodel(data)))
  lambda_max <- calc_lambda_max_cml(beta_intercepts = beta_estimates_restricted, 
                                    covariates = covariates, Y = data)[[DIF_type]]
  lambdas <- seq(10^-5, lambda_max, length.out = nlambdas - 
                   1)
  lambdas <- c(lambdas, lambdas[nlambdas - 1] + (lambdas[nlambdas - 
                                                           1] - lambdas[nlambdas - 2]))
  time_nlm <- system.time(results <- lapply(1:length(lambdas), 
                                              function(x) RM_DIF_nlm(data = data, covariates = covariates, 
                                                                       lambda = lambdas[x], starting_betas = res_unrestricted$beta_estimates, 
                                                                       starting_deltas = res_unrestricted$delta_estimates, 
                                                                       DIF_type = DIF_type,
                                                                       lambda_restr = 0, restr_pen = FALSE,
                                                                       ...)))
  delta_est <- array(unlist(lapply(1:length(lambdas), function(x) results[[x]]$estimatevalues$deltas)), 
                     c(ncov, nitems, length(lambdas)))
  beta_est <- t(sapply(1:length(lambdas), function(x) results[[x]]$estimatevalues$betas))
  beta_est <- beta_est - rowMeans(beta_est)
  BIC_list <- NULL
  for (j in 1:length(lambdas)) {
    BIC_list[j] <- results[[j]]$estimatevalues$BIC <- BIC_lasso(Y = data, 
                                                                betas = results[[j]]$estimatevalues$betas, deltas = results[[j]]$estimatevalues$deltas, 
                                                                covariates = covariates, DIF_type = DIF_type, delta_unres = res_unrestricted$delta_estimates, 
                                                                tolerance = tolerance_df)
  }
  model_picked <- ifelse(1 %in% which(BIC_list == min(BIC_list)), 
                         max(which(BIC_list == min(BIC_list))), which.min(BIC_list))
  res_picked <- results[[model_picked]]
  result <- list(results_per_lambda = results, BICs = BIC_list, 
                 beta_estimates = beta_est, delta_estimates = delta_est, 
                 lambdas = lambdas, times = list(time_unrestricted_model = time_unres, 
                                                 time_penalized_part = time_nlm, time_restricted_model = time_res), 
                 model_picked = model_picked, lambda_picked = lambdas[model_picked], 
                 result_picked = res_picked, result_unrestricted = res_unrestricted, 
                 beta_estimates_RM = beta_estimates_restricted)
}


cml_lasso_win_nlm_optimized_foreach <- function (data, covariates, DIF_type = "all.interactions", nlambdas = 21, 
                                         tolerance_df = 0.0001,
                                         no_of_cores_to_use = ifelse(detectCores() == 
                                                                       1, 1, detectCores() - 1),...) 
{
  npersons <- nrow(data)
  nitems <- ncol(data)
  ncov <- ncol(covariates)
  time_unres <- system.time(res_unrestricted <- RM_DIF_lbfgs(data = data, 
                                                             covariates = covariates, lambda = 0, invisible = 1))
  time_res <- system.time(beta_estimates_restricted <- -itempar(raschmodel(data)))
  lambda_max <- calc_lambda_max_cml(beta_intercepts = beta_estimates_restricted, 
                                    covariates = covariates, Y = data)[[DIF_type]]
  lambdas <- seq(10^-5, lambda_max, length.out = nlambdas - 
                   1)
  lambdas <- c(lambdas, lambdas[nlambdas - 1] + (lambdas[nlambdas - 
                                                           1] - lambdas[nlambdas - 2]))
  time_nlm <- system.time({
    cl = makeCluster(no_of_cores_to_use)
    registerDoParallel(cl)
    results <- foreach(i = lambdas, .packages = c("cmlDIFlasso2", 
                                                  "psychotools", "lbfgs")) %dopar% {
      RM_DIF_nlm(data = data, covariates = covariates, 
                 lambda = i, starting_betas = res_unrestricted$beta_estimates, 
                 starting_deltas = res_unrestricted$delta_estimates, 
                 DIF_type = DIF_type,
                 lambda_restr = 0, restr_pen = FALSE,
                 ...)
    }
    stopCluster(cl)
  })
  delta_est <- array(unlist(lapply(1:length(lambdas), function(x) results[[x]]$estimatevalues$deltas)), 
                     c(ncov, nitems, length(lambdas)))
  beta_est <- t(sapply(1:length(lambdas), function(x) results[[x]]$estimatevalues$betas))
  beta_est <- beta_est - rowMeans(beta_est)
  BIC_list <- NULL
  for (j in 1:length(lambdas)) {
    BIC_list[j] <- results[[j]]$estimatevalues$BIC <- BIC_lasso(Y = data, 
                                                                betas = results[[j]]$estimatevalues$betas, deltas = results[[j]]$estimatevalues$deltas, 
                                                                covariates = covariates, DIF_type = DIF_type, delta_unres = res_unrestricted$delta_estimates, 
                                                                tolerance = tolerance_df)
  }
  model_picked <- ifelse(1 %in% which(BIC_list == min(BIC_list)), 
                         max(which(BIC_list == min(BIC_list))), which.min(BIC_list))
  res_picked <- results[[model_picked]]
  result <- list(results_per_lambda = results, BICs = BIC_list, 
                 beta_estimates = beta_est, delta_estimates = delta_est, 
                 lambdas = lambdas, times = list(time_unrestricted_model = time_unres, 
                                                 time_penalized_part = time_nlm, time_restricted_model = time_res), 
                 model_picked = model_picked, lambda_picked = lambdas[model_picked], 
                 result_picked = res_picked, result_unrestricted = res_unrestricted, 
                 beta_estimates_RM = beta_estimates_restricted)
}

