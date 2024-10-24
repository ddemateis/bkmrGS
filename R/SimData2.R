HFun <- function (z, opt = 1){ 
  
  if(opt==1){
   h <- 4*plogis(1/4 * (z[1] + z[2] + 1/2 * z[1] * z[2]), 0, 0.3)
  }else if(opt == 2){
   h <- 4 - 4 * plogis(1/4 * (z[1] + z[2] + 1/2 * z[1] * z[2]), 0, 0.3)
  }else if(opt == 3){
    h <- 0
  }else if(opt == 4){
    h <- 2*plogis(1/4 * (z[1] + z[2] + 1/2 * z[1] * z[2]), 0, 0.3)
  }
  
  return(h)
  
}

#' Simulate dataset with modifier
#'
#' Simulate predictor, covariate, and continuous outcome data
#'
#' @export
#'
#' @inheritParams kmbayes
#' @param opt Simulation scenario option: 1 for no modification, 2 for opposite effects, 3 for an effect in one group, and 4 for a scaled effect between groups
#' @param SNR signal-to-noise ratio
#' @param bin_mod 0 for group 1, 1 for group 2
#' @param mod_DGM indicator for including modifier in data generating mechanism
#' @param sim_exp indicator for using simulated exposure values and covariates. Not recommended for simulation
SimData2 <- function (opt = 1,
                      SNR = 10,
                      bin_mod = 0,
                      mod_DGM = T,
                      sim_exp = F){
  
  #data_standardized is lazy loaded
  med_vita <- median(data_standardized$vita)
  if(bin_mod == 0){
    dta <- data_standardized[data_standardized$vita < med_vita,]
  }else if(bin_mod == 1){
    dta <- data_standardized[data_standardized$vita >= med_vita,]
  }else{
    stop("wrong bin_mod value. either 0 or 1")
  }
  n <- nrow(dta)
  
  #generate exposures
  M <- 2
  if(sim_exp){
    Z <- matrix(rnorm(n * M), n, M)
  }else{
    Z <- cbind(scale(dta$pb_ln), scale(dta$mn_ln))
  }
  colnames(Z) <- paste0("z", 1:M)
  
  #generate covariates
  if(sim_exp){
    X <- cbind(3 * cos(Z[,1]) + 2 * rnorm(n),
               3 * cos(Z[,2]) + 2 * rnorm(n))
  }else{
    X <- cbind(scale(dta$X1),
               scale(dta$X2),
               scale(dta$X3))
  }
  
  #randomly generate covariate coefficients
  beta.true <- rnorm(ncol(X))
  
  #construct exposure-response curve
  h <- apply(Z, 1, HFun, opt = opt)
  
  
  #mean response
  if(mod_DGM){
    mu <- X %*% beta.true + h + bin_mod #bin_mod is modifier, defaulting main effect to 1 for binary value 1
    
  }else{
    mu <- X %*% beta.true + h
  }
  
  #set noise for specified SNR
  signal <- sd(mu)
  noise <- signal/SNR
  
  #add noise to obtain response
  eps <- rnorm(n, sd = noise)
  y <- mu + eps
  
  #data structure
  dat <- list(n = n, 
              M = M, 
              beta.true = beta.true, 
              Z = Z, 
              h = h, 
              X = X, 
              y = y, 
              HFun = HFun, 
              opt = opt,
              noise2 = noise^2)
  return(dat)
  
}
