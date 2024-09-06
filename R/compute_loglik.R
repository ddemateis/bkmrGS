#' Compute log likelihood matrix
#'
#' Compute log likelihood matrix for WAIC computation
#'
#' @export
#' 
#' @importFrom stats dnorm
#' 
#' @param fit an object of class "bkmrfit"
#' @param sel A vector selecting which iterations of the BKMR fit should be retained for inference. If not specified, will default to keeping every 10 iterations after dropping the first half of samples, or if this results in fewer than 100 iterations, than 100 iterations are kept
#' @return a matrix with a row for each posterior sample and column for each data point
compute_loglik <- function(fit, sel = NULL){
  
  #add condition for default if sel not provided
  if(is.null(sel)){
    
  }
  
  loglik_samps <- matrix(NA, nrow = length(sel), ncol = length(y))
  for(s in sel){
    sigma2 <- fit$sigsq.eps[s]
    HXB <- SamplePred(fit, sel = s)
    
    #point-wise loglikelihood vector for n observations
    loglik_samps[which(s == sel),] <- dnorm(x = fit$y, #length n
                                            mean = HXB, #h + X * beta, length n
                                            sd = sigma2, #sigma2, length 1
                                            log = T) #ll_vec is length n
  }
  
  return(loglik_samps)
  
}
