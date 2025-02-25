#used in non-interactive effect functions (overall, single-exposure)
riskSummary.approx <- function(point1, point2, modnew = NULL, preds.fun, ...) {
  cc <- c(-1, 1)
  newz <- rbind(point1, point2)
  # if(!is.null(modnew)){
  #   modnew <- matrix(rep(modnew, nrow(newz)), ncol=1)
  # }
  preds <- preds.fun(newz, modnew, ...)
  diff <- drop(cc %*% preds$postmean) #E[AX] = A*E[X]
  diff.sd <- drop(sqrt(cc %*% preds$postvar %*% cc)) #V[AX] = A * V[X] *A'
  c(est = diff, sd = diff.sd)
}

#not used
riskSummary.samp <- function(point1, point2, preds.fun, ...) {
  cc <- c(-1, 1)
  newz <- rbind(point1, point2)
  preds <- preds.fun(newz, ...)
  diff.preds <- drop(preds %*% cc)
  c(est = mean(diff.preds), sd = sd(diff.preds))
}

#used in interaction effect functions
interactionSummary.approx <- function(newz.q1, newz.q2, modnew.1, modnew.2, preds.fun, ...) {
  newz <- rbind(newz.q1, newz.q2)
  modnew <- c(modnew.1, modnew.2)
  preds <- preds.fun(newz, modnew, ...)
  cc <- c(-1*c(-1, 1), c(-1, 1))
  if("matrix" %in% class(preds)){
    post_samp <- preds %*% matrix(cc, ncol=1)
    int <- drop(mean(post_samp))
    int.lb <- drop(quantile(post_samp, probs = 0.025))
    int.ub <- drop(quantile(post_samp, probs = 0.975))
    c(est = int, lb = int.lb, ub = int.ub)
  }else{
    int <- drop(cc %*% preds$postmean) #E[AX] = A*E[X]
    int.se <- drop(sqrt(cc %*% preds$postvar %*% cc)) #V[AX] = A * V[X] *A'
    c(est = int, sd = int.se)
  }
  
}

#not used
interactionSummary.samp <- function(newz.q1, newz.q2, preds.fun, ...) {
  cc <- c(-1*c(-1, 1), c(-1, 1))
  newz <- rbind(newz.q1, newz.q2)
  preds <- preds.fun(newz, ...)
  int.preds <- drop(preds %*% cc)
  c(est = mean(int.preds), sd = sd(int.preds))
}



#' Calculate overall risk summaries
#' 
#' Compare estimated \code{h} function when all predictors are at a particular quantile to when all are at a second fixed quantile
#' @inheritParams kmbayes
#' @inheritParams ComputePostmeanHnew
#' @inherit ComputePostmeanHnew details
#' @param qs vector of quantiles at which to calculate the overall risk summary 
#' @param q.fixed a second quantile at which to compare the estimated \code{h} function
#' @param m.fixed the modifier value at which to compare the estimated \code{h} function
#' @export
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the overall risk measures
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' risks.overall <- OverallRiskSummaries(fit = fitkm, qs = seq(0.25, 0.75, by = 0.05), 
#' q.fixed = 0.5, method = "exact")
OverallRiskSummaries <- function(fit, y = NULL, Z = NULL, X = NULL, 
                                 modifier = NULL, 
                                 qs = seq(0.25, 0.75, by = 0.05), 
                                 q.fixed = 0.5, m.fixed = NULL,
                                 method = "approx", sel = NULL) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) 
      y <- fit$y
    if (is.null(Z)) 
      Z <- fit$Z
    if (is.null(X)) 
      X <- fit$X
    if (is.null(modifier)) 
      modifier <- fit$modifier   
  }
  
  # happens in computepostmeanhnew
  # #convert modifier to factor and construct contrast matrix
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  # }
  
  if(!is.null(m.fixed) & kernel.method == "one"){#BKMR-mod
    modnew <- matrix(rep(m.fixed, 2), ncol=1)
    Z_for_quants <- Z
  }else if(!is.null(m.fixed) & kernel.method == "two"){#GS-BKMR
    modnew <- matrix(rep(m.fixed, 2), ncol=1)
    Z_for_quants <- Z[modifier == m.fixed,]
  }else{
    modnew <- NULL
    Z_for_quants <- Z
  }
  
  #get the quantiles of the exposures for the first point of comparison specified by q.fixed
  point1 <- apply(Z_for_quants, 2, quantile, q.fixed)

  if (method %in% c("approx", "exact")) {
    
    preds.fun <- function(znew, modnew=NULL) {
    
      ComputePostmeanHnew(fit = fit,
                          y = y,
                          Z = Z,
                          X = X,
                          modifier = modifier,
                          Znew = znew,
                          mod_new = modnew,
                          sel = sel,
                          method = method)
      
      }
    riskSummary <- riskSummary.approx
  }else {
    stop("method must be one of c('approx', 'exact')")
  }
  
  tmp_fn <- function(quant) { 
    riskSummary(point1 = point1, 
                point2 = apply(Z_for_quants, 2, quantile, quant),
                modnew = modnew,
                preds.fun = preds.fun) 
  }
  
  #for each quantile in qs, call call riskSummary function
  risks.overall <- t(sapply(qs, tmp_fn))
  risks.overall <- data.frame(quantile = qs, risks.overall)
}

#' Calculate modifier interaction summaries
#' 
#' Compare estimated \code{h} function when all predictors are at a particular quantile to when all are at a second fixed quantile
#' @inheritParams kmbayes
#' @inheritParams ComputePostmeanHnew
#' @inherit ComputePostmeanHnew details
#' @param q.fixed vector of quantiles at which to calculate the summary 
#' @param m.diff vector of two modifier values at which to compare the estimated \code{h} function
#' @export
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the overall risk measures

ModifierIntSummaries <- function(fit, y = NULL, Z = NULL, X = NULL, 
                                 modifier = NULL, m.diff,
                                 q.fixed = seq(0.25, 0.75, by = 0.05), 
                                 method = "approx", sel = NULL) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) 
      y <- fit$y
    if (is.null(Z)) 
      Z <- fit$Z
    if (is.null(X)) 
      X <- fit$X
    if (is.null(modifier)) 
      modifier <- fit$modifier   
  }
  
  # happens in computepostmeanhnew
  # #convert modifier to factor and construct contrast matrix
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  # }
  
  if (method %in% c("approx", "exact")) {
    
    preds.fun <- function(znew, modnew=NULL) {
      
      ComputePostmeanHnew(fit = fit,
                          y = y,
                          Z = Z,
                          X = X,
                          modifier = modifier,
                          Znew = znew,
                          mod_new = modnew,
                          sel = sel,
                          method = method)
      
    }
    riskSummary <- riskSummary.approx
  }else {
    stop("method must be one of c('approx', 'exact')")
  }
  
  tmp_fn <- function(quant) { 
    riskSummary(point1 = apply(Z, 2, quantile, quant), 
                point2 = apply(Z, 2, quantile, quant),
                modnew = m.diff,
                preds.fun = preds.fun) 
  }
  
  #for each quantile in qs, call call riskSummary function
  risks.overall <- t(sapply(q.fixed, tmp_fn))
  risks.overall <- data.frame(quantile = q.fixed, risks.overall)
}

#used in SingVarRiskSummaries() below
#Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile
VarRiskSummary <- function (whichz = 1, fit, y = NULL, Z = NULL, X = NULL,
                            modifier = NULL, qs.diff = c(0.25, 0.75), 
                            q.fixed = 0.5, method = "approx", m.fixed = NULL, 
                            sel = NULL, ...){ 
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) 
      y <- fit$y
    if (is.null(Z)) 
      Z <- fit$Z
    if (is.null(X)) 
      X <- fit$X
    if(is.null(modifier))
      modifier <- fit$modifier
  }
  
  # #convert modifier to factor and construct contrast matrix
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  #   m.fixed <- factor(m.fixed, levels = levels(factor(orig_modifier)))
  #   m.fixed <- as.matrix(model.matrix(~m.fixed)[,-1])
  # }
  
  if(!is.null(m.fixed) & kernel.method == "one"){#BKMR-mod
    modnew <- matrix(rep(m.fixed, 2), ncol=1)
    Z_for_quants <- Z
  }else if(!is.null(m.fixed) & kernel.method == "two"){#GS-BKMR
    modnew <- matrix(rep(m.fixed, 2), ncol=1)
    Z_for_quants <- Z[modifier == m.fixed,]
  }else{
    modnew <- NULL
    Z_for_quants <- Z
  }
  
  point2 <- point1 <- apply(Z_for_quants, 2, quantile, q.fixed) #get exposures at fixed quantile q.fixed, along with the fixed modifier value
  point2[whichz] <- apply(Z_for_quants[, whichz, drop = FALSE], 2, quantile, 
                          qs.diff[2]) #for given exposure whichz, change that exposure to the larger quantile in the difference, qs.diff[2]
  point1[whichz] <- apply(Z_for_quants[, whichz, drop = FALSE], 2, quantile, 
                          qs.diff[1]) #for given exposure whichz, change that exposure to the larger quantile in the difference, qs.diff[1]
  
  if (method %in% c("approx", "exact")) {
    
    preds.fun <- function(znew, modnew=NULL){
      
      ComputePostmeanHnew(fit = fit, 
                          y = y, 
                          Z = Z, 
                          X = X,
                          modifier = modifier, 
                          Znew = znew, 
                          mod_new = modnew,
                          sel = sel, 
                          method = method)
      }
    riskSummary <- riskSummary.approx
  }
  else {
    stop("method must be one of c('approx', 'exact')")
  }
  
  riskSummary(point1 = point1,
              point2 = point2,
              modnew = modnew,
              preds.fun = preds.fun, 
              ...)
}

#' Single Variable Risk Summaries
#' 
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
#' 
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @inheritParams OverallRiskSummaries
#' @inherit ComputePostmeanHnew details
#' @param qs.diff vector indicating the two quantiles \code{q_1} and \code{q_2} at which to compute \code{h(z_{q2}) - h(z_{q1})}
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in \code{Z}
#' @param z.names optional vector of names for the columns of \code{z}
#' @param ... other arguments to pass on to the prediction function
#' @param which.z vector indicating which variables (columns of \code{Z}) for which the summary should be computed
#' @export
#' 
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the single-predictor risk measures
#' 
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' risks.singvar <- SingVarRiskSummaries(fit = fitkm, method = "exact")
SingVarRiskSummaries <- function(fit, y = NULL, Z = NULL, X = NULL, 
                                 modifier = NULL, which.z = 1:ncol(Z),
                                 qs.diff = c(0.25, 0.75), 
                                 q.fixed = c(0.25, 0.50, 0.75), 
                                 m.fixed = NULL, method = "approx", 
                                 sel = NULL, z.names = colnames(Z), ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) 
      y <- fit$y
    if (is.null(Z)) 
      Z <- fit$Z
    if (is.null(X)) 
      X <- fit$X
    if(is.null(modifier)) 
      modifier <- fit$modifier 
  }
  
  # #convert modifier to factor and construct contrast matrix - not needed because happens in VarRiskSummary()
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  #   m.fixed <- factor(m.fixed, levels = levels(factor(orig_modifier)))
  #   m.fixed <- as.matrix(model.matrix(~m.fixed)[,-1])
  # }
  
  if(!is.null(m.fixed) & is.null(modifier)){
    stop("Cannot specify modifier values for model without modification.")
  }
  
  if (is.null(z.names)) 
    z.names <- paste0("z", 1:ncol(Z))
  df <- dplyr::tibble()
  
  for (i in seq_along(q.fixed)) { #for each quantile (for fixing other exposures) in q.fixed
    for (j in seq_along(which.z)) { #for each exposure in which.z
      risk <- VarRiskSummary(whichz = which.z[j], 
                             fit = fit, 
                             y = y, 
                             Z = Z, 
                             X = X, 
                             qs.diff = qs.diff, 
                             q.fixed = q.fixed[i], 
                             method = method, 
                             m.fixed = m.fixed, 
                             sel = sel) 
      df0 <- dplyr::tibble(q.fixed = q.fixed[i], variable = z.names[j], 
                           est = risk["est"], sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)
    }
  }
  
  df <- dplyr::mutate_at(df, "variable", function(x) factor(x, 
                                                            levels = z.names[which.z]))
  df <- dplyr::mutate_at(df, "q.fixed", function(x) as.factor(x))
  attr(df, "qs.diff") <- qs.diff
  df
}

#used in SingVarIntSummaries() below
SingVarIntSummary <- function(whichz = 1, fit, y = NULL, Z = NULL, 
                              X = NULL, modifier = NULL, 
                              qs.diff = c(0.25, 0.75), 
                              qs.fixed = c(0.25, 0.75), mod.diff = NULL,
                              method = "approx", sel = NULL, ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
    if (is.null(modifier)) modifier <- fit$modifier
  }
  
  # #convert modifier to factor and construct contrast matrix
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  #   mod.diff <- factor(mod.diff, levels = levels(factor(orig_modifier)))
  #   mod.diff <- as.matrix(model.matrix(~mod.diff)[,-1])
  # }
  
  #modifier for difference, if modification 
  if(!is.null(modifier)){#need modnew even for two-kernel model, SamplePred handles it
    modnew.1 <- rep(mod.diff[1], 2)
    modnew.2 <-rep(mod.diff[2], 2)
    
    #subset exposures so each exposure has the same min and max per group
    Z_for_quants <- Z_support(Z = Z, 
                              mod.diff = mod.diff, 
                              modifier = modifier)
  }else{
    modnew.1 <- rep(NULL, 2)
    modnew.2 <- rep(NULL, 2)
    Z_for_quants <- Z
  }

  #second difference: h(z_1^{qs2}, z_2^{q2}, \dots, z_M^{q2}) - h(z_1^{qs1}, z_2^{q2}, \dots, z_M^{q2})
  q.fixed <- qs.fixed[1] #first quantile to fix remaining exposures
  point2 <- point1 <- apply(Z_for_quants, 2, quantile, q.fixed) #get desired quantiles for all exposures, for group 1 if there is a group 1
  point2[whichz] <- quantile(Z_for_quants[, whichz], qs.diff[2]) #h(z_1^{qs2}, z_2^{q2}, \dots, z_M^{q2})
  point1[whichz] <- quantile(Z_for_quants[, whichz], qs.diff[1]) #h(z_1^{qs1}, z_2^{q2}, \dots, z_M^{q2})
  newz.q1 <- rbind(point1, point2) 
  
  #first difference: h(z_1^{qs2}, z_2^{q1}, \dots, z_M^{q1}) - h(z_1^{qs1}, z_2^{q1}, \dots, z_M^{q1})
  q.fixed <- qs.fixed[2] #first quantile to fix remaining exposures
  point2 <- point1 <- apply(Z_for_quants, 2, quantile, q.fixed) #get desired quantiles for all exposures, for group 2 if there is a group 2
  point2[whichz] <- quantile(Z_for_quants[, whichz], qs.diff[2]) #h(z_1^{qs2}, z_2^{q1}, \dots, z_M^{q1})
  point1[whichz] <- quantile(Z_for_quants[, whichz], qs.diff[1]) #h(z_1^{qs1}, z_2^{q1}, \dots, z_M^{q1})
  newz.q2 <- rbind(point1, point2)
  
  if (method %in% c("approx", "exact", "fullpost")) {
    preds.fun <- function(znew, modnew){
      if(method %in% c("approx", "exact")){
        ComputePostmeanHnew(fit = fit, 
                            y = y, 
                            Z = Z, 
                            X = X, 
                            Znew = znew, 
                            mod_new = modnew,
                            sel = sel, 
                            method = method)
      }else{
        SamplePred(fit = fit,
                   y = y,
                   Z = Z,
                   X = X,
                   modifier = modifier,
                   Znew = znew,
                   Xnew = matrix(0, nrow = nrow(znew), ncol = ncol(X)),
                   mod_new = modnew,
                   sel = sel)
      }
    }
    interactionSummary <- interactionSummary.approx
  } else {
    stop("method must be one of c('approx', 'exact', 'fullpost')")
  }
  interactionSummary(newz.q1, newz.q2, modnew.1, modnew.2, preds.fun, ...)
}

#' Single Variable Interaction Summaries
#' 
#' Compare the single-predictor health risks when all of the other predictors in Z are fixed to their a specific quantile to when all of the other predictors in Z are fixed to their a second specific quantile. For interaction overall effect, set \code{mod.diff} to two separate levels of the modifier, and make \code{qs.diff} and \code{qs.fixed} the same vector. 
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @inheritParams SingVarRiskSummaries
#' @inherit ComputePostmeanHnew details
#' @param qs.diff vector indicating the two quantiles at which to compute the single-predictor risk summary
#' @param qs.fixed vector indicating the two quantiles at which to fix all of the remaining exposures in \code{Z}
#' @param mod.diff vector of two modifier values at which to compute the single-predictor risk summary (only for models with modification.) 
#' @export
#' 
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the single-predictor risk measures
#' 
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' risks.int <- SingVarIntSummaries(fit = fitkm, method = "exact")
SingVarIntSummaries <- function(fit, y = NULL, Z = NULL, X = NULL,
                                modifier = NULL, which.z = 1:ncol(Z), 
                                qs.diff = c(0.25, 0.75), 
                                qs.fixed = c(0.25, 0.75), 
                                mod.diff = NULL, method = "approx", 
                                sel = NULL, z.names = colnames(Z), ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
    if (is.null(modifier)) modifier <- fit$modifier
  }
  
  #convert modifier to factor and construct contrast matrix
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  #   mod.diff <- factor(mod.diff, levels = levels(factor(orig_modifier)))
  #   mod.diff <- as.matrix(model.matrix(~mod.diff)[,-1])
  # }
  
  if(!is.null(mod.diff) & is.null(modifier)){
    stop("Cannot specify modifier values for model without modification.")
  }
  
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))
  
  #for each exposure, apply SingVarIntSummary()
  ints <- sapply(which.z, function(whichz)
    SingVarIntSummary(whichz = whichz, 
                      fit = fit, 
                      Z = Z, 
                      X = X, 
                      y = y, 
                      modifier = modifier, 
                      qs.diff = qs.diff, 
                      qs.fixed = qs.fixed, 
                      mod.diff = mod.diff, 
                      method = method, 
                      sel = sel,
                      ...)
  )
  
  if("sd" %in% rownames(ints)){#is sd returned
    df <- dplyr::tibble(variable = factor(z.names[which.z], levels = z.names), est = ints["est", ], sd = ints["sd", ])
  }else{#if lb ub returned
    df <- dplyr::tibble(variable = factor(z.names[which.z], levels = z.names), est = ints["est", ], lb = ints["lb.2.5%", ], ub = ints["ub.97.5%", ])
  }
}

#used in OverallIntSummaries() below
OverallIntSummary <- function(whichz = 1, fit, y = NULL, Z = NULL, 
                              X = NULL, modifier = NULL, 
                              qs = 0.25, q.fixed = 0.5, mod.diff = NULL,
                              method = "approx", sel = NULL, ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
    if (is.null(modifier)) modifier <- fit$modifier
  }
  
  # #convert modifier to factor and construct contrast matrix
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  # }
  
  #modifier
  modnew.1 <-rep(mod.diff[1], 2)
  modnew.2 <-rep(mod.diff[2], 2)
  
  #subset exposures so each exposure has the same min and max per group
  Z_for_quants <- Z_support(Z = Z,
                            mod.diff = mod.diff,
                            modifier = modifier)

  #differences, regardless of modifier: h(z_1^{qs}, z_2^{qs}, \dots, z_M^{qs}) - h(z_1^{q}, z_2^{q}, \dots, z_M^{q})
  point2 <- apply(Z_for_quants, 2, quantile, q.fixed)
  point1 <- apply(Z_for_quants, 2, quantile, qs)
  newz.q1 <- newz.q2 <- rbind(point1, point2) 
  
  if (method %in% c("approx", "exact", "fullpost")) {
    preds.fun <- function(znew, modnew){
      if(method %in% c("approx", "exact")){
        ComputePostmeanHnew(fit = fit, 
                            y = y, 
                            Z = Z, 
                            X = X, 
                            Znew = znew, 
                            mod_new = modnew,
                            sel = sel, 
                            method = method)
      }else{
        SamplePred(fit = fit,
                   y = y,
                   Z = Z,
                   X = X,
                   modifier = modifier,
                   Znew = znew,
                   Xnew = matrix(0, nrow = nrow(znew), ncol = ncol(X)),
                   mod_new = modnew,
                   sel = sel)
      }
    }
    interactionSummary <- interactionSummary.approx
  } else {
    stop("method must be one of c('approx', 'exact', 'fullpost')")
  }

  interactionSummary(newz.q1, newz.q2, modnew.1, modnew.2, preds.fun, ...)
}

#' Overall Effect Interaction Summaries
#' 
#' Compare the overall effect difference between groups for models with modification 
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @inheritParams SingVarRiskSummaries
#' @inherit ComputePostmeanHnew details
#' @param qs vector of quantiles at which to calculate the overall risk summary 
#' @param q.fixed a second quantile at which to compare the estimated \code{h} function
#' @param mod.diff a vector of two unique modifier values to compare the estimated \code{h} function
#' @export
#' 
#' @return a data frame containing the (posterior mean) estimate and posterior standard deviation of the single-predictor risk measures
#' 
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' 
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
#' 
#' risks.int <- SingVarIntSummaries(fit = fitkm, method = "exact")
OverallIntSummaries <- function(fit, y = NULL, Z = NULL, X = NULL,
                                modifier = NULL, which.z = 1:ncol(Z), 
                                qs = seq(0.25, 0.75,0.05), mod.diff,
                                q.fixed = 0.5, method = "approx", 
                                sel = NULL, z.names = colnames(Z), ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
    if (is.null(modifier)) modifier <- fit$modifier
  }
  
  # #convert modifier to factor and construct contrast matrix
  # if(!is.null(modifier)){
  #   orig_modifier <- modifier
  #   modifier <- as.factor(modifier)
  #   modifier <- as.matrix(model.matrix(~modifier)[,-1])
  # }
  
  if(is.null(modifier)){
    stop("Cannot compute overall interactive effect between groups for model without modification.")
  }
  
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(Z))
  
  #for each quantile, apply OverallIntSummary()
  risks.overall <- t(sapply(qs, function(q_s)
    OverallIntSummary(fit = fit, 
                      Z = Z, 
                      X = X, 
                      y = y, 
                      modifier = modifier, 
                      qs = q_s, 
                      q.fixed = q.fixed, 
                      mod.diff = mod.diff,
                      method = method, 
                      sel = sel,
                      ...)))
  risks.overall <- data.frame(quantile = qs, risks.overall)
}

