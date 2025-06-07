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
  dta <- data_standardized[-c(50),] #50 is delivery type, 253 is pica
  n <- nrow(dta)
  
  if(scenario == "scaled3"){
    M <- 3
  }else{
    M <- 2
  }
  
  #sort by modifier
  if(scenario != "multi"){
    dta <- dta[order(dta$X2),]
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
  
  #generate covariates
  if(sim_exp){
    X <- cbind(3 * cos(Z[,1]) + 2 * rnorm(n),
               3 * cos(Z[,2]) + 2 * rnorm(n))
  }else{
    X_full <- data.frame(X1 = dta$X1, #continuous: visit age
                         X2 = as.factor(dta$X2), #categorical: child sex 
                         X3 = dta$X3, #continuous: gestation weeks, 
                         X4 = as.factor(dta$X4), #categorical: delivery type, 2 categories
                         X5 = as.factor(case_when(round(dta$X5, 2) == 0.03 ~ "birth_order_1",
                                                  round(dta$X5, 2) == 0.05 ~ "birth_order_2",
                                                  round(dta$X5, 2) >= 0.08 ~ "birth_order_3")), #categorical: birth order, originally 6 categories, we made 3 by grouping the last 4 together
                         #X6 = data_standardized$X6, #continuous: drink water cups, 
                         #X7 = as.factor(data_standardized$X7), #categorical: hospital child (y/n),
                         #X8 = as.factor(data_standardized$X8), #categorical: pica (?) (y/n),
                         X9 = as.factor(case_when(round(dta$X9, 2) <= 0.02 ~ "ed_1", #ed_1 contains the 0 group as well as the next lowest level of ed
                                                  round(dta$X9, 2) == 0.04 ~ "ed_2",
                                                  round(dta$X9, 2) >= 0.06 ~ "ed_3")), #categorical: education for birthing parent,
                         X10 = as.factor(case_when(round(dta$X10, 2) <= 0.02 ~ "ed_1", #ed_1 contains the 0 group as well as the next lowest level of ed
                                                   round(dta$X10, 2) == 0.04 ~ "ed_2",
                                                   round(dta$X10, 2) >= 0.07 ~ "ed_3")), #categorical: education for non-birthing parent,
                         X11 = as.factor(dta$X11), #categorical: smoke environment,
                         X12 = dta$X12, #continuous: HOME score (emotional), 
                         X13 = dta$X13, #continuous: HOME score (avoid), 
                         X14 = dta$X14, #continuous: HOME score (careg), 
                         X15 = dta$X15, #continuous: HOME score (env), 
                         X16 = dta$X16, #continuous: HOME score (play), 
                         X17 = dta$X17, #continuous: HOME score (stim), 
                         X18 = dta$X18 #continuous: child estimate daily energy intake, 
                         #X19 = data_standardized$X19, #continuous: baseline lead, 
                         #X20 = data_standardized$X20, #continuous: baseline manganese, 
                         #X21 = data_standardized$X21 #continuous: baseline arsenic, 
    )
    if(scenario != "multi"){
      X_full <- X_full[,-c(2)] #remove X12 because it is the modifier and is included as a covariate in the model by default
    }else{
      X_full <- X_full[,-c(13)] #remove X16 because it is the modifier and is included as a covariate in the model by default
    }
    
    X <- model.matrix(~., data=X_full)[,-1]
  }
  
  #randomly generate covariate coefficients
  beta.true <- rnorm(ncol(X))
  
  #construct exposure-response curve
  h <- c()
  if(scenario != "multi"){
    #med_vita <- median(data_standardized$vita)
    modifier <- ifelse(round(dta$X2,1) == 0, 
                       "low", 
                       "high") #ifelse(data_standardized$vita < med_vita, "low", "high")
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
  
  #including other null exposures
  Z <- cbind(Z, scale(dta$as_ln))
  colnames(Z) <- paste0("z", 1:ncol(Z))
  
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
