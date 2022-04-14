###
#
# function for all information criteria
#
###


calc_all_IC_RM_lasso <- function (Y, betas, deltas, covariates, DIF_type, delta_unres, tolerance) {
  
  cloglik <- cl_sample_2(Y = Y, betas = calc_itemvecs(betas, 
                                                      deltas, covariates))
  # cloglik <- cl_sample_PCM(Y = Y, 
  #                          sigma_mat = calc_sigma_p_sample_mat(sigmas, deltas, covariates, catmax), 
  #                          catmax = catmax)
  df_lambda_result <- df_lambda(deltas, DIF_type, delta_unres, 
                                tolerance)
  n <- nrow(Y)
  I <- ncol(Y)
  BIC <- -2 * cloglik + df_lambda_result * log(n * I)
  BIC_n <- -2 * cloglik + df_lambda_result * log(n)
  AIC <- -2 * cloglik + 2 * df_lambda_result
  CAIC <- -2 * cloglik + df_lambda_result * (log(n * I) + 1)
  CAIC_n <- -2 * cloglik + df_lambda_result * (log(n) + 1)
  AICc <- AIC + ((2 * df_lambda_result * (df_lambda_result + 
                                            1))/(n * I - df_lambda_result - 1))
  AICc_n <- AIC + ((2 * df_lambda_result * (df_lambda_result + 
                                              1))/(n - df_lambda_result - 1))
  return(list(BIC = BIC, BIC_n = BIC_n, AIC = AIC, CAIC = CAIC, 
              CAIC_n = CAIC_n, AICc = AICc, AICc_n = AICc_n))
}

# calculate selected models
# add to result
# replace $BIC
# replace $model_picked (number) $lambda_picked and 
# $result_picked depending on new_ic

new_ic = "BIC_n" # one of "BIC"    "BIC_n"  "AIC"    "CAIC"   "CAIC_n" "AICc"   "AICc_n"

# i=10000
for (i in 1:length(res_foreach)) {
  # calculate all inf crit
  if ((i %% 100) == 0) cat(i,"  ")
  criterion_list <- NULL
  for (j in 1:length(res_foreach[[i]]$lambdas)) {
    
    criterion_list[[j]] <- 
      res_foreach[[i]]$result$results_per_lambda[[j]]$estimatevalues <- 
      calc_all_IC_RM_lasso(res_foreach[[i]]$scenario$data, res_foreach[[i]]$result$results_per_lambda[[j]]$beta_estimates,
                           res_foreach[[i]]$result$results_per_lambda[[j]]$delta_estimates, res_foreach[[i]]$scenario$covariatevalues,
                           "all.interactions", res_foreach[[i]]$result$result_unrestricted$delta_estimates, 0)
    
  }
  # make inf criterion matrix
  IC_matrix <- do.call(rbind, lapply(criterion_list, unlist))
  # calc model picked numbers per criterion
  res_foreach[[i]]$result$models_picked_all_IC <- models_picked_all_IC <- apply(do.call(rbind, lapply(criterion_list, 
                                                                                                      unlist)), 2, which.min)
  # model picked for IC selected
  res_foreach[[i]]$result$model_picked <- model_picked <- models_picked_all_IC[[new_ic]]
  
  res_foreach[[i]]$result$BIC_and_other_IC <- IC_matrix # store IC_matrix
  
  res_foreach[[i]]$result$BICs <- IC_matrix[,new_ic] # new IC list (still named BICs, but new IC)
  # change IC values per lambda result
  for (j in 1:length(res_foreach[[i]]$lambdas))  res_foreach[[i]]$result$results_per_lambda[[j]]$estimatevalues$BIC <- IC_matrix[j,new_ic]
  # change model_picked and result_picked
  res_foreach[[i]]$result$lambda_picked <- res_foreach[[i]]$result$lambdas[model_picked]
  res_foreach[[i]]$result$result_picked <- res_foreach[[i]]$result$results_per_lambda[[model_picked]]
  res_foreach[[i]]$result$IC_chosen <- new_ic
  
}
