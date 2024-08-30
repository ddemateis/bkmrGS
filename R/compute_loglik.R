#' @importFrom mvtnorm dmvnorm
log_norm_density <- function(x, y, X1, Z, data.comps){
  
  #extract parameters from input vector x
  sigma2 <- as.numeric(x[grep("sigma2", names(x))])
  r <- as.numeric(x[grep("r", names(x))])
  lambda <- as.numeric(x[grep("lambda", names(x))])
  beta_mat <- as.matrix(x[grep("beta", names(x))])
  
  #computations
  K <- exp(-makeKpart(r, Z))
  V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-K) #I + lambda*K
  if (data.comps$nlambda == 2) {
    V <- V + lambda[2]*data.comps$crossTT
  }
  cov_mat <- as.matrix(sigma2*V)
  
  #loglikelihood matrix
  ll_mat <- dmvnorm(x = y, 
                    mean = X1%*%t(beta_mat), #X * beta
                    sigma = cov_mat, #sigma2 * V
                    log = T)
  
  return(ll_mat)
}

#' Compute log likelihood matrix
#'
#' Compute log likelihood matrix for WAIC computation
#'
#' @export
#' 
#' @param fit an object of class "bkmrfit"
#' @return a matrix with a row for each posterior sample and column for each data point

compute_loglik <- function(fit){
  
  X1 <- as.matrix(fit$X)
  Z <- as.matrix(fit$Z)
  y <- fit$y
  modifier <- fit$modifier
  Z <- cbind(Z, modifier)
  X1 <- as.matrix(cbind(X1, modifier))
  data.comps <- fit$data.comps
  
  posterior_mat <- data.frame(sigma2 = fit$sigsq.eps,
                              r = fit$r,
                              beta = fit$beta,
                              lambda = fit$lambda)
  
  loglik_samps <- matrix(NA, nrow = nrow(posterior_mat), ncol = length(y))
  for(i in 1:nrow(posterior_mat)){
    loglik_samps[i,] <- log_norm_density(x = posterior_mat[i,], 
                                         y = y, 
                                         X1 = X1, 
                                         Z = Z, 
                                         data.comps = data.comps)
  }
  
  return(loglik_samps)
  
  #error: Error in X1 %*% t(beta_mat) : non-conformable arguments
  # loglik_samps <- apply(X = posterior_mat, 
  #                       MARGIN = 1, 
  #                       FUN = log_norm_density, 
  #                       y = y, 
  #                       X1 = X1, 
  #                       Z = Z, 
  #                       data.comps = data.comps)
}
