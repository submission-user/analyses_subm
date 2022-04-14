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

