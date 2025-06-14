#' Obtain posterior samples of predictions at new points
#'
#' Obtains posterior samples of \code{E(Y) = h(Znew) + beta*Xnew} or of \code{g^{-1}[E(y)]}
#' 
#' @param sel A vector selecting which iterations of the BKMR fit should be retained for inference. If not specified, will default to keeping every 10 iterations after dropping the first 50\% of samples, or if this results in fewer than 100 iterations, than 100 iterations are kept
#' @param Znew optional matrix of new predictor values at which to predict new \code{h}, where each row represents a new observation. If not specified, defaults to using observed Z values
#' @param mod_new optional vector of new modifier values at which to predict new \code{h}, If not specified, defaults to using observed modifier values
#' @param Xnew optional matrix of new covariate values at which to obtain predictions. If not specified, defaults to using observed X values
#' @param type whether to make predictions on the scale of the link or of the response; only relevant for the binomial outcome family
#' @param ... other arguments; not currently used
#' @inheritParams kmbayes
#' @inheritParams ExtractEsts
#' @details For guided examples, go to \url{https://jenfb.github.io/bkmr/overview.html}
#' @export
#' 
#' @return a matrix with the posterior samples at the new points
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
#' set.seed(111)
#' samps3 <- SamplePred(fitkm, Znew = Znew, Xnew = cbind(0))
#' head(samps3)
SamplePred <- function(fit, Znew = NULL, Xnew = NULL, mod_new = NULL, 
                       Z = NULL, X = NULL, modifier = NULL, y = NULL, 
                       sel = NULL, type = c("link", "response"), ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
    if (is.null(modifier)) modifier <- fit$modifier #added by DD
  }
  if (length(type) > 1) type <- type[1]
  
  #convert modifier to factor and construct contrast matrix
  if(!is.null(modifier)){
    orig_modifier <- modifier
    modifier <- as.factor(modifier)
    modifier <- as.matrix(model.matrix(~modifier)[,-1])
    if(!is.null(mod_new)){
      orig_mod_new <- mod_new
      mod_new <- factor(mod_new, levels = levels(factor(orig_modifier)))
      mod_new <- matrix(model.matrix(~mod_new)[,-1], ncol = ncol(modifier))
    }
  }
  
  kernel.method <- fit$kernel.method
  if(kernel.method == "one"){
    Z <- cbind(Z, modifier)
  }else if (kernel.method == "two"){
    Z <- Z
  }

  if (!is.null(Znew)) { #what is Znew is NULL?
    if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow = 1)
    if (inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
    
    if(kernel.method == "one"){
      Znew <- cbind(Znew, mod_new)
    }else if(kernel.method == "two"){
      Znew <- Znew
    }
    
    if (ncol(Z) != ncol(Znew)) {
      stop("Znew must have the same number of columns as Z")
    }
  }

  if (is.null(Xnew)){
    Xnew <- X
  }
  if(is.null(mod_new)){
    mod_new <- modifier #added by DD
  }
  if (!inherits(Xnew, "matrix")) Xnew <- matrix(Xnew, nrow = 1)
  X <- cbind(X, modifier) #added by DD
  Xnew <- cbind(Xnew, mod_new) #added by DD
  if (ncol(X) != ncol(Xnew)) {
    stop("Xnew must have the same number of columns as X")
  }
  
  if (is.null(sel)) {
    sel <- with(fit, seq(floor(iter/2) + 1, iter, 10))
    if (length(sel) < 100) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter, length.out = 100))
    }
    sel <- unique(floor(sel))
  }
  
  family <- fit$family
  data.comps <- fit$data.comps
  lambda <- fit$lambda
  sigsq.eps <- fit$sigsq.eps
  beta <- fit$beta
  r <- fit$r
  
  if (!is.null(Znew)) {
    preds <- matrix(NA, length(sel), nrow(Znew))
    colnames(preds) <- paste0("znew", 1:nrow(Znew))
  } else {
    preds <- matrix(NA, length(sel), nrow(Z))
    colnames(preds) <- paste0("z", 1:nrow(Z))
  }
  rownames(preds) <- paste0("iter", sel)
  for (s in sel) {
    beta.samp <- beta[s, ]
    
    if (family == "gaussian") {
      ycont <- y
    } else if (family == "binomial") {
      ycont <- fit$ystar[s, ]
    }
    if (!is.null(Znew)) { #Z and X d
      hsamp <- newh.update(Z = Z, Znew = Znew, mod_new = mod_new, Vcomps = NULL, lambda = lambda[s, ], sigsq.eps = sigsq.eps[s], r = r[s, ], y = ycont, X = X, beta = beta.samp, data.comps = data.comps, modifier = modifier, kernel.method = kernel.method)  
    } else { 
      hsamp <- h.update(lambda = lambda[s, ], Vcomps = NULL, sigsq.eps = sigsq.eps[s], y = ycont, X = X, beta = beta.samp, r = r[s, ], Z = Z, data.comps = data.comps, modifier = modifier, kernel.method = kernel.method)$hsamp
    }
    
    Xbeta <- drop(Xnew %*% beta.samp)
    linpred <- hsamp + Xbeta
    
    if (type == "link") {
      pred <- linpred
    } else if (type == "response") {
      if (family == "gaussian") {
        pred <- linpred
      } else if (family == "binomial") {
        pred <- pnorm(linpred)
      }
    }
    preds[paste0("iter", s), ] <- pred
  }
  attr(preds, "type") <- type
  attr(preds, "family") <- family
  preds
  
}