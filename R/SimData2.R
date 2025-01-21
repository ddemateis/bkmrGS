HFun <- function (z, opt = 1){ 
  
  if(length(z) == 2){
    z_sum <- (z[1] + z[2] + 1/2 * z[1] * z[2])
  }else if(length(z == 3)){
    z_sum <- (z[1] + z[2] + z[3] + 1/2 * z[1] * z[2] * z[3])
  }
  
  if(opt==1){
   h <- 4*plogis(1/4 * z_sum, 0, 0.3)
  }else if(opt == 2){
   h <- 4 - 4 * plogis(1/4 * z_sum, 0, 0.3)
  }else if(opt == 3){
    h <- 0
  }else if(opt == 4){
    h <- 2*plogis(1/4 * z_sum, 0, 0.3)
  }else if(opt == 5){
    h <- plogis(1/4 * z_sum, 0, 0.3)
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
#' @param scenario Simulation scenario option: "none" for no modification between two groups, "oneGroup" for one group effect out of two groups, "scaled2" for one groups having a scaled effect of the other group with 2 exposures, "scaled3" is the same but for 3 exposures, and "multi" for one group with no effect and the other two groups having a scaled effect
#' @param SNR signal-to-noise ratio
#' @param mod_DGM indicator for including modifier in data generating mechanism
#' @param sim_exp indicator for using simulated exposure values and covariates. Not recommended for simulation
SimData2 <- function (scenario = "none",
                      SNR = 10,
                      mod_DGM = T,
                      sim_exp = F){
  
  
  #data_standardized is lazy loaded
  dta <- data_standardized
  n <- nrow(dta)
  
  if(scenario == "scaled3"){
    M <- 3
  }else{
    M <- 2
  }
  
  #sort by modifier
  if(scenario != "multi"){
    dta <- dta[order(dta$vita),]
  }else{
    dta <- dta[order(dta$X16),]
  }
  
  #generate exposures
  if(sim_exp){
    Z <- matrix(rnorm(n * M), n, M)
  }else{
    if(M==1){
      Z <- scale(dta$pb_ln)
    }else if(M==2){
      Z <- cbind(scale(dta$pb_ln), scale(dta$mn_ln))
    }else if(M==3){
      Z <- cbind(scale(dta$pb_ln), scale(dta$mn_ln), scale(dta$as_ln))
    }else{
      stop("M not supported.")
    }
    
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
  h <- c()
  if(scenario != "multi"){
    med_vita <- median(data_standardized$vita)
    modifier <- ifelse(data_standardized$vita < med_vita, "low", "high")
    for(i in 1:n){
      if(scenario == "none"){
        opt = 1
      }else if(scenario == "oneGroup"){
        opt = 3
      }else if(scenario == "scaled2" | scenario == "scaled3"){
        opt = 4
      }
      h[i] <- HFun(Z[i,], opt = ifelse(modifier[i] == "low", 1, opt))
    }
    mod_effect <- as.matrix(model.matrix(~modifier)[,-1]) * 1
    modifier <- factor(modifier, levels = c("low", "high"))
  }else{
    modifier <- vector("character", n)
    ref <- factor(round(dta$X16, 2))
    modifier[ref == "0.02"] <- "low" #13% of observations
    modifier[ref == "0.05"] <- "medium" #66% of observations
    modifier[ref == "0.07"] <- "high" #21% of observations
    for(i in 1:n){
      if(modifier[i] == "low"){
        opt = 1
      }else if(modifier[i] == "medium"){
        opt = 4
      }else if(modifier[i] == "high"){
        opt = 3
      }
      h[i] <- HFun(Z[i,], opt = opt)
    }
    mod_effect <- model.matrix(~modifier)[,-1] %*% c(1,1)
    modifier <- factor(modifier, levels = c("low", "medium", "high"))
  }
  
  #mean response
  if(mod_DGM){
    mu <- X %*% beta.true + h + mod_effect
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
              modifier = modifier,
              HFun = HFun, 
              scenario = scenario,
              noise2 = noise^2)
  return(dat)
  
}
