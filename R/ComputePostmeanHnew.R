#' Compute the posterior mean and variance of \code{h} at a new predictor values
#'
#' @inheritParams kmbayes
#' @param fit An object containing the results returned by a the \code{kmbayes} function 
#' @param Znew matrix of new predictor values at which to predict new \code{h}, where each row represents a new observation. If set to NULL then will default to using the observed exposures Z.
#' @param mod_new vector of new modifier values at which to predict new \code{h}. If set to NULL then will default to using the observed modifiers.
#' @param method method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets; see details
#' @param sel selects which iterations of the MCMC sampler to use for inference; see details
#' @details
#' \itemize{
#'   \item If \code{method == "approx"}, the argument \code{sel} defaults to the second half of the MCMC iterations.
#'   \item If \code{method == "exact"}, the argument \code{sel} defaults to keeping every 10 iterations after dropping the first 50\% of samples, or if this results in fewer than 100 iterations, than 100 iterations are kept
#' }
#' For guided examples and additional information, go to \url{https://jenfb.github.io/bkmr/overview.html}
#' @export
#' 
#' @return a list of length two containing the posterior mean vector and posterior variance matrix 
#' 
#' @examples
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
#' med_vals <- apply(Z, 2, median)
#' Znew <- matrix(med_vals, nrow = 1)
#' h_true <- dat$HFun(Znew)
#' h_est1 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "approx")
#' h_est2 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "exact")
ComputePostmeanHnew <- function(fit, y = NULL, Z = NULL, X = NULL, modifier = NULL, Znew = NULL, mod_new = NULL, sel = NULL, method = "approx") {

  if (method == "approx") {
    res <- ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z, X = X, modifier = modifier, Znew = Znew, mod_new = mod_new, sel = sel)
  } else if (method == "exact") {
    res <- ComputePostmeanHnew.exact(fit = fit, y = y, Z = Z, X = X, modifier = modifier, Znew = Znew, mod_new = mod_new, sel = sel)
  }
  res
}

#' Compute the posterior mean and variance of \code{h} at a new predictor values
#'
#' Function to approximate the posterior mean and variance as a function of the estimated model parameters (e.g., tau, lambda, beta, and sigsq.eps)
#' @param Znew matrix of new predictor values at which to predict new \code{h}, where each row represents a new observation. If set to NULL then will default to using the observed exposures Z.
#' @param mod_new vector of new modifier values at which to predict new \code{h}. If set to NULL then will default to using the observed modifiers.
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @noRd
ComputePostmeanHnew.approx <- function(fit, y = NULL, Z = NULL, 
                                       X = NULL, modifier = NULL,
                                       Znew = NULL,  mod_new = NULL,
                                       sel = NULL) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
    if (is.null(modifier)) modifier <- fit$modifier
  }
  
  #convert modifier to factor and construct contrast matrix
  if(!is.null(modifier)){
    orig_modifier <- modifier
    modifier <- as.factor(modifier)
    modifier <- as.matrix(model.matrix(~modifier)[,-1])
    if(!is.null(mod_new)){
      orig_mod_new <- mod_new
      mod_new <- factor(mod_new, levels = levels(factor(orig_modifier)))
      mod_new <- as.matrix(model.matrix(~mod_new)[,-1])
    }
  }
  
  kernel.method <- fit$kernel.method
  if(kernel.method == "one"){
    kern_modifier <- NULL
  }else if (kernel.method == "two"){
    kern_modifier <- modifier
  }
  
  if (!is.null(Znew)) {
    if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
    if(inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
  }
  if(is.null(dim(X))) X <- matrix(X, ncol=1)
  
  Znew <- cbind(Znew, mod_new)
  Z <- cbind(Z, modifier)
  X <- cbind(X, modifier)
  
  ests <- ExtractEsts(fit, sel = sel)
  sigsq.eps <- ests$sigsq.eps[, "mean"]
  r <- ests$r[, "mean"]
  beta <- ests$beta[, "mean"]
  lambda <- ests$lambda[, "mean"]
  if (fit$family == "gaussian") {
    ycont <- y
  } else if (fit$family == "binomial") {
    ycont <- ests$ystar[, "mean"]
  }
  
  Kpart <- makeKpart(r, Z)
  K <- exp(-Kpart)
  if(kernel.method == "two"){
    K <- block_kernel(mod_vec1 = modifier,
                      mod_vec2 = modifier,
                      K = K)
  }
  scalar_lambda <- lambda_scalar(mod1 = modifier,
                                 mod2 = modifier,
                                 lambda = lambda,
                                 data.comps = data.comps)
  V <- diag(1, nrow(Z), nrow(Z)) + scalar_lambda*K
  cholV <- chol(V)
  Vinv <- chol2inv(cholV)
  
  if (!is.null(Znew)) {
    # if(is.null(data.comps$knots)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    Kpartall <- makeKpart(r, rbind(Z, Znew))
    Kmat <- exp(-Kpartall)
    if(kernel.method == "two"){
      Kmat <- block_kernel(mod_vec1 = modifier,
                           mod_vec2 = mod_new,#this was modifier, but should be mod_new to match in Kpartall
                           K = Kmat)
    }
    Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
    Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
    
    scalar_lambda <- lambda_scalar(mod1 = modifier,
                                   mod2 = mod_new,
                                   lambda = lambda,
                                   data.comps = data.comps)
    if(data.comps$gs.tau){ #scalar_lambda returns a matrix that needs to be adjusted to match the Ks
      scalar_lambda_10 <- scalar_lambda[(n0+1):nall,1:n0 ,drop=FALSE]
      scalar_lambda_1 <- scalar_lambda[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    }else{ #scalar_lambda returns a scalar, no adjustment needed
      scalar_lambda_10 <- scalar_lambda 
      scalar_lambda_1 <- scalar_lambda
    }
    lamK10Vinv <- scalar_lambda_10*Kmat10 %*% Vinv
    postvar <- scalar_lambda_1*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
    postmean <- lamK10Vinv %*% (ycont - X%*%beta)
    # } else {
    # stop("GPP not yet implemented")
    # }
  } else {
    scalar_lambda <- lambda_scalar(mod1 = modifier,
                                   mod2 = modifier,
                                   lambda = lambda,
                                   data.comps = data.comps)
    lamKVinv <- scalar_lambda*K%*%Vinv
    postvar <- scalar_lambda*sigsq.eps*(K - lamKVinv%*%K)
    postmean <- lamKVinv %*% (ycont - X%*%beta)
  }
  ret <- list(postmean = drop(postmean), postvar = postvar)
  ret
}

#' Compute the posterior mean and variance of \code{h} at a new predictor values
#'
#' Function to estimate the posterior mean and variance by obtaining the posterior mean and variance at particular iterations and then using the iterated mean and variance formulas
#' 
#' @inheritParams kmbayes
#' @inheritParams SamplePred
#' @inheritParams ExtractEsts
#' 
#' @noRd
ComputePostmeanHnew.exact <- function(fit, y = NULL, Z = NULL, X = NULL, modifier = NULL, Znew = NULL, mod_new = NULL, sel = NULL) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
    if (is.null(modifier)) modifier <- fit$modifier
  }
  
  #convert modifier to factor and construct contrast matrix
  if(!is.null(modifier)){
    orig_modifier <- modifier
    modifier <- as.factor(modifier)
    modifier <- as.matrix(model.matrix(~modifier)[,-1])
    if(!is.null(mod_new)){
      orig_mod_new <- mod_new
      mod_new <- factor(mod_new, levels = levels(factor(orig_modifier)))
      mod_new <- as.matrix(model.matrix(~mod_new)[,-1])
    }
  }
  
  kernel.method <- fit$kernel.method
  if(kernel.method == "one"){
    kern_modifier <- NULL
  }else if (kernel.method == "two"){
    kern_modifier <- modifier
  }
  
  if (!is.null(Znew)) {
    if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow = 1)
    if (inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
    if (ncol(Z) != ncol(Znew)) {
      stop("Znew must have the same number of columns as Z")
    }
  }
  
  if (is.null(dim(X))) X <- matrix(X, ncol=1)
  
  if(kernel.method == "one"){
    Znew <- cbind(Znew, mod_new)
    Z <- cbind(Z, modifier)
  }else if(kernel.method == "two"){
    Znew <- Znew
    Z <- Z
  }
  X <- cbind(X, modifier)
  
  # if (!is.null(fit$Vinv)) {
  #   sel <- attr(fit$Vinv, "sel")
  # }
  
  if (is.null(sel)) {
    sel <- with(fit, seq(floor(iter/2) + 1, iter, 10))
    if (length(sel) < 100) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter, length.out = 100))
    }
    sel <- unique(floor(sel))
  }
  
  family <- fit$family
  data.comps <- fit$data.comps
  post.comps.store <- list(postmean = vector("list", length(sel)),
                           postvar = vector("list", length(sel))
  )
  
  for (i in seq_along(sel)) {
    s <- sel[i]
    beta <- fit$beta[s, ]
    lambda <- fit$lambda[s, ]
    sigsq.eps <- fit$sigsq.eps[s]
    r <- fit$r[s, ]
    
    if (family == "gaussian") {
      ycont <- y
    } else if (family == "binomial") {
      ycont <- fit$ystar[s, ]
    }
    
    Kpart <- makeKpart(r, Z)
    K <- exp(-Kpart)
    if(kernel.method == "two"){
      K <- block_kernel(mod_vec1 = modifier,
                        mod_vec2 = modifier,
                        K = K)
    }
    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps, modifier = kern_modifier)
    Vinv <- Vcomps$Vinv
    # if (is.null(fit$Vinv)) {
    # V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*K
    # cholV <- chol(V)
    # Vinv <- chol2inv(cholV)
    # } else {
    #   Vinv <- fit$Vinv[[i]]
    # }
    
    if (!is.null(Znew)) {
      # if(is.null(data.comps$knots)) {
      n0 <- nrow(Z)
      n1 <- nrow(Znew)
      nall <- n0 + n1
      Kpartall <- makeKpart(r, rbind(Z, Znew))
      Kmat <- exp(-Kpartall)
      if(kernel.method == "two"){
        Kmat <- block_kernel(mod_vec1 = rbind(modifier, mod_new),
                             mod_vec2 = rbind(modifier, mod_new),
                             K = Kmat)
      }
      Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
      Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
      Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
      
      scalar_lambda <- lambda_scalar(mod1 = rbind(modifier, mod_new),
                                     mod2 = rbind(modifier, mod_new),
                                     lambda = lambda,
                                     data.comps = data.comps)
      if(data.comps$gs.tau){ #scalar_lambda returns a matrix that needs to be adjusted to match the Ks
        scalar_lambda_10 <- scalar_lambda[(n0+1):nall,1:n0 ,drop=FALSE]
        scalar_lambda_1 <- scalar_lambda[(n0+1):nall,(n0+1):nall ,drop=FALSE]
      }else{ #scalar_lambda returns a scalar, no adjustment needed
        scalar_lambda_10 <- scalar_lambda 
        scalar_lambda_1 <- scalar_lambda
      }
      lamK10Vinv <- scalar_lambda_10*Kmat10 %*% Vinv
      postvar <- scalar_lambda_1*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
      postmean <- lamK10Vinv %*% (ycont - X%*%beta)
      # } else {
      # stop("GPP not yet implemented")
      # }
    } else {
      scalar_lambda <- lambda_scalar(mod1 = modifier,
                                     mod2 = modifier,
                                     lambda = lambda,
                                     data.comps = data.comps)
      lamKVinv <- scalar_lambda*K%*%Vinv
      postvar <- scalar_lambda*sigsq.eps*(K - lamKVinv%*%K)
      postmean <- lamKVinv %*% (ycont - X%*%beta)
    }
    
    post.comps.store$postmean[[i]] <- postmean
    post.comps.store$postvar[[i]] <- postvar
    
  }
  
  postmean_mat <- t(do.call("cbind", post.comps.store$postmean))
  m <- colMeans(postmean_mat)
  postvar_arr <- with(post.comps.store, 
                      array(unlist(postvar), 
                            dim = c(nrow(postvar[[1]]), ncol(postvar[[1]]), length(postvar)))
  )
  ve <- var(postmean_mat)
  ev <- apply(postvar_arr, c(1, 2), mean)
  v <- ve + ev
  ret <- list(postmean = m, postvar = v)
  
  ret
}

