#' Computing posterior distribution for r.squared in compositional data analysis
#'
#' `bayes_R2.CoDa` Function to compute posterior distribution of r.squared in compositional data analysis
#'
#' @param model compositional model fitted using bamlss
#' @param nsamples number of samples to be used when the r.squared is computed
#' @param type if like residuals based residual variance, if not, based on residuals
#' @param resp_names names of the response variable transformed in the ilr coordinates. If NULL, the variables which contains the word ilr are used.
#' @return simulations of the posterior distribution of r.squared
#'
#' @export
#' @import brms
#' @import berryFunctions
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
bayes_R2.CoDa <- function(model1, type = "like", resp_names = NULL) 
{
  y_pred <- posterior_epred(model1)
  
  #For each linear predictor
  var_fit <- matrix(nrow = dim(y_pred)[1], ncol = dim(y_pred)[3])
  for (i in 1:dim(y_pred)[3])
  {
    var_fit[,i] <- apply(y_pred[,,i], 1, var)
  }
  var_fit_tot <- apply(var_fit, 1, sum)
  
  if(type == "like")
  {
    # Residual fit
    sigma2 <- as.matrix(model1, pars = c("sigma"))^2
    # sigma2 <- posterior::as_draws_matrix(model1)[,"sigma"]^2
    # 
    # sigma2 <- as.matrix(model1, variable = "sigma")^2
    
    var_res <- apply(sigma2, 1, sum)
  } else{
    dim(y_pred)
    error <- array(NA, dim = dim(y_pred))
    if(is.null(resp_names))
    {
      resp_names <- model1$data %>% colnames(.) 
      resp_names <- resp_names[grep("ilr", resp_names)]
    }
    y <- model1$data[,resp_names]
    error <- 1:dim(y_pred)[1] %>%
      lapply(., function(x) y_pred[x,,] - y) %>%
      berryFunctions::l2array(.) %>%
      aperm(., c(3,1,2))
    
    var_error <- matrix(nrow = dim(y_pred)[1], ncol = dim(y_pred)[3])
    for (i in 1:dim(y_pred)[3])
    {
      var_error[,i] <- apply(error[,,i], 1, var)
    }
    var_res <- apply(var_error, 1, sum)
  }
  r.squared <- var_fit_tot/(var_fit_tot + var_res)
  r.squared
}
