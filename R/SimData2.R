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
#' @param scenario Simulation scenario option: "none" for no modification between two groups, "oneGroup" for one group effect out of two groups, "scaled2" for one groups having a scaled effect of the other group with 2 exposures, "scaled3" is the same but for 3 exposures, and "multi" for one group with no effect and the other two groups having a scaled effect
#' @param SNR signal-to-noise ratio
#' @param mod_DGM indicator for including modifier in data generating mechanism
#' @param covariates "old" for the 10 used in the first submission, "new" for the 14 used in the analysis and new submission
SimData2 <- function (scenario = "none",
                      SNR = 10,
                      mod_DGM = TRUE,
                      covariates = "new"){
  
  #ex_data is lazy loaded
  dta <- cbind(ex_data$Z, ex_data$X)
  n <- nrow(dta) #350

  #sort by modifier
  if(scenario != "multi"){
    dta <- dta[order(dta$Sex),] #child sex
  }else{ 
    dta <- dta[order(dta$HOME_play),] #HOME score (play)
  }
  
  #set exposures
  Z <- cbind(scale(dta$`Log Lead`), scale(dta$`Log Manganese`))
  
  #set covariates
  if(covariates == "old"){ #covariates from original submission
    #data_standardized needs to be loaded externally from the package
    #this is only here to check an issue, remove after found
    load("data_standardized.rda")
    data_standardized <- data_standardized[-c(50),]
    X <- cbind(scale(data_standardized$X1),
               scale(data_standardized$X3),
               scale(data_standardized$X4),
               scale(data_standardized$X5),
               scale(data_standardized$X6),
               scale(data_standardized$X7),
               scale(data_standardized$X8),
               scale(data_standardized$X9),
               scale(data_standardized$X10),
               scale(data_standardized$X11))
  }else if(covariates == "new"){ #covariates from analysis and new submission
    if(scenario != "multi"){
      X_full <- data.frame(dta$Age,
                           dta$Gestation,
                           dta$Delivery,
                           dta$Birth_order,
                           dta$Education_parent1,
                           dta$Education_parent2,
                           dta$Smoking,
                           dta$HOME_emotional,
                           dta$HOME_avoid,
                           dta$HOME_careg,
                           dta$HOME_env,
                           dta$HOME_play,
                           dta$HOME_stim,
                           dta$Energy) 
    }else{
      X_full <- data.frame(dta$Age,
                           dta$Sex,
                           dta$Gestation,
                           dta$Delivery,
                           dta$Birth_order,
                           dta$Education_parent1,
                           dta$Education_parent2,
                           dta$Smoking,
                           #dta$HOME_emotional, #dropping home covariates due to inversion issue in MCMC
                           #dta$HOME_avoid,
                           #dta$HOME_careg,
                           #dta$HOME_env,
                           #dta$HOME_stim,
                           dta$Energy)
    }
    X <- model.matrix(~., data=X_full)[,-1]
  }
  
  #randomly generate covariate coefficients
  beta.true <- rnorm(ncol(X))
  
  #construct exposure-response curve
  h <- c()
  if(scenario != "multi"){
    modifier <- ifelse(dta$Sex == "male", 
                       "low", 
                       "high") 
    for(i in 1:n){
      if(scenario == "none"){
        opt = 1
      }else if(scenario == "oneGroup"){
        opt = 3
      }else if(scenario == "scaled2"){
        opt = 4
      }
      h[i] <- HFun(Z[i,], opt = ifelse(modifier[i] == "low", 1, opt))
    }
    mod_effect <- as.matrix(model.matrix(~modifier)[,-1]) * 1
    modifier <- factor(modifier, levels = c("low", "high"))
  }else{
    modifier <- vector("character", n)
    ref <- factor(round(dta$HOME_play, 2))
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
  
  #including other null exposures
  Z <- cbind(Z, scale(dta$`Log Arsenic`))
  colnames(Z) <- paste0("z", 1:ncol(Z))
  
  #data structure
  dat <- list(n = n, 
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
