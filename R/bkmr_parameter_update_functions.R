beta.update <- function(X, Vinv, y, sigsq.eps) {
	XVinv <- crossprod(X, Vinv)
	Vbeta <- chol2inv(chol(XVinv %*% diag((sigsq.eps)^(-1)) %*% X)) #when sigsq.eps could not vary by group: 	Vbeta <- chol2inv(chol(XVinv %*% X))
	cholVbeta <- chol(Vbeta)
	betahat <- Vbeta %*% XVinv %*% diag((sigsq.eps)^(-1)) %*% y #when sigsq.eps could not vary by group: 	betahat <- Vbeta %*% XVinv %*% y
	n01 <- rnorm(ncol(X))
	betahat + crossprod(cholVbeta, n01) #when sigsq.eps could not vary by group: 	betahat + crossprod(sqrt(sigsq.eps)*cholVbeta, n01)
}

sigsq.eps.update <- function(y, X, beta, Vinv, a.eps=1e-3, b.eps=1e-3) {
	mu <- y - X%*%beta
	prec.y <- rgamma(1, shape=a.eps + nrow(X)/2, rate=b.eps + 1/2*crossprod(mu, Vinv)%*%mu)
	1/prec.y
}

ystar.update <- function(y, X, beta, h) {
  mu <-  drop(h + X %*% beta)
  lower <- ifelse(y == 1, 0, -Inf)
  upper <- ifelse(y == 0, 0,  Inf)
  samp <- truncnorm::rtruncnorm(1, a = lower, b = upper, mean = mu, sd = 1)
  drop(samp)
}
#' @importFrom tmvtnorm rtmvnorm
ystar.update.noh <- function(y, X, beta, Vinv, ystar) {
  mu <-  drop(X %*% beta)
  lower <- ifelse(y == 1, 0, -Inf)
  upper <- ifelse(y == 0, 0,  Inf)
  samp <- tmvtnorm::rtmvnorm(1, mean = mu, H = Vinv, lower = lower, upper = upper, algorithm = "gibbs", start.value = ystar)
  #samp <- truncnorm::rtruncnorm(1, a = lower, b = upper, mean = mu, sd = 1)
  drop(samp)
}

r.update <- function(r, whichcomp, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, data.comps, control.params, rprop.gen, rprop.logdens, rprior.logdens, modifier = NULL, ...) {
	# r.params <- set.r.params(r.prior = control.params$r.prior, comp = whichcomp, r.params = control.params$r.params)
	r.params <- make_r_params_comp(control.params$r.params, whichcomp)
	rcomp <- unique(r[whichcomp])
	if(length(rcomp) > 1) stop("rcomp should only be 1-dimensional")
	
	## generate a proposal
	rcomp.star <- rprop.gen(current = rcomp, r.params = r.params)
	lambda.star <- lambda
	delta.star <- delta
	move.type <- NA

	## part of M-H ratio that depends on the proposal distribution
	negdifflogproposal <- -rprop.logdens(rcomp.star, rcomp, r.params = r.params) + rprop.logdens(rcomp, rcomp.star, r.params = r.params)

	## prior distribution
	diffpriors <- rprior.logdens(rcomp.star, r.params = r.params) - rprior.logdens(rcomp, r.params = r.params)

	r.star <- r
	r.star[whichcomp] <- rcomp.star

	## M-H step
	r.return <- MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps, modifier = modifier)
	return(r.return)
}

rdelta.comp.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, ztest, data.comps, control.params, rprop.gen2, rprop.logdens1, rprior.logdens, rprior.logdens2, rprop.logdens2, rprop.gen1, modifier = NULL, ...) { ## individual variable selection
	r.params <- control.params$r.params
	a.p0 <- control.params$a.p0
	b.p0 <- control.params$b.p0
	delta.star <- delta
	r.star <- r

	move.type <- ifelse(all(delta[ztest] == 0), 1, sample(c(1,2),1))
	move.prob <- ifelse(all(delta[ztest] == 0), 1, 1/2)
	if(move.type == 1) {
		comp <- ifelse(length(ztest) == 1, ztest, sample(ztest, 1))
		r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

		delta.star[comp] <- 1 - delta[comp]
		move.prob.star <- ifelse(all(delta.star[ztest] == 0), 1, 1/2)
		r.star[comp] <- ifelse(delta.star[comp] == 0, 0, rprop.gen1(r.params = r.params))

		diffpriors <- (lgamma(sum(delta.star[ztest]) + a.p0) + lgamma(length(ztest) - sum(delta.star[ztest]) + b.p0) - lgamma(sum(delta[ztest]) + a.p0) - lgamma(length(ztest) - sum(delta[ztest]) + b.p0)) + ifelse(delta[comp] == 1, -1, 1)*with(list(r.sel = ifelse(delta[comp] == 1, r[comp], r.star[comp])), rprior.logdens(x = r.sel, r.params = r.params))

 		negdifflogproposal <- -log(move.prob.star) + log(move.prob) - ifelse(delta[comp] == 1, -1, 1)*with(list(r.sel = ifelse(delta[comp] == 1, r[comp], r.star[comp])), rprop.logdens1(x = r.sel, r.params = r.params))

	} else if(move.type == 2) {
		comp <- ifelse(length(which(delta == 1)) == 1, which(delta == 1), sample(which(delta == 1), 1))
		r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

		r.star[comp] <- rprop.gen2(current = r[comp], r.params = r.params)

		diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[comp], r.params = r.params)

		negdifflogproposal <- -rprop.logdens2(r.star[comp], r[comp], r.params = r.params) + rprop.logdens2(r[comp], r.star[comp], r.params = r.params)
	}

	lambda.star <- lambda

	## M-H step
	return(MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps, modifier = modifier))
}

rdelta.group.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, ztest, data.comps, control.params, rprop.gen1, rprior.logdens, rprop.logdens1, rprop.gen2, rprop.logdens2, modifier = NULL, ...) { ## grouped variable selection
	r.params <- control.params$r.params
	a.p0 <- control.params$a.p0
	b.p0 <- control.params$b.p0
	groups <- control.params$group.params$groups
	sel.groups <- control.params$group.params$sel.groups
	neach.group <- control.params$group.params$neach.group
	delta.star <- delta
	r.star <- r

	# if(length(mu.r) == 1) mu.r <- rep(mu.r, nz)
	# if(length(sigma.r) == 1) sigma.r <- rep(sigma.r, nz)

	delta.source <- sapply(sel.groups, function(x) ifelse(any(delta[which(groups == groups[x])] == 1), 1, 0))
	delta.source.star <- delta.source

	## randomly select move type
	if(all(delta.source == 0)) {
		move.type <- 1
		move.prob <- 1
	} else if(length(which(neach.group > 1 & delta.source == 1)) == 0) {
		move.type <- sample(c(1, 3), 1)
		move.prob <- 1/2
	} else {
		move.type <- sample(1:3, 1)
		move.prob <- 1/3
	}
	# move.type <- ifelse(all(delta.source == 0), 1, ifelse(length(which(neach.group > 1 & delta.source == 1)) == 0, sample(c(1, 3), 1), sample(1:3, 1)))

	# print(move.type)

	if(move.type == 1) { ## randomly select a source and change its state (e.g., from being in the model to not being in the model)

		source <- sample(seq_along(delta.source), 1)
		source.comps <- which(groups == source)

		# r.params <- set.r.params(r.prior = control.params$r.prior, comp = source.comps, r.params = r.params)

		delta.source.star[source] <- 1 - delta.source[source]
		delta.star[source.comps] <- rmultinom(1, delta.source.star[source], rep(1/length(source.comps), length(source.comps)))
		move.prob.star <- ifelse(all(delta.source.star == 0), 1, ifelse(length(which(neach.group > 1 & delta.source.star == 1)) == 0, 1/2, 1/3))

		## which component got switched
		comp <- ifelse(delta.source[source] == 1, source.comps[which(delta[source.comps] == 1)], source.comps[which(delta.star[source.comps] == 1)])
		r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

		r.star[comp] <- ifelse(delta.star[comp] == 0, 0, rprop.gen1(r.params = r.params))

		# diffpriors <- ifelse(delta.source[source] == 1, log(length(sel.groups) - sum(delta.source) + b.p0) - log(sum(delta.source.star) + a.p0), log(sum(delta.source) + a.p0) - log(length(sel.groups) - sum(delta.source.star) + b.p0)) + ifelse(delta.source[source] == 1, 1, -1)*log(length(source.comps)) + ifelse(delta.source[source] == 1, -1, 1)*with(list(r.sel = ifelse(delta.source[source] == 1, r[source.comps][which(delta[source.comps] == 1)], r.star[source.comps][which(delta.star[source.comps] == 1)])), rprior.logdens(x = r.sel, r.params = r.params))
		diffpriors <- ifelse(delta.source[source] == 1, log(length(sel.groups) - sum(delta.source) + b.p0) - log(sum(delta.source.star) + a.p0), log(sum(delta.source) + a.p0) - log(length(sel.groups) - sum(delta.source.star) + b.p0)) + ifelse(delta.source[source] == 1, 1, -1)*log(length(source.comps)) + ifelse(delta.source[source] == 1, -1, 1)*with(list(r.sel = ifelse(delta.source[source] == 1, r[comp], r.star[comp])), rprior.logdens(x = r.sel, r.params = r.params))

		# negdifflogproposal <- -log(move.prob.star) + log(move.prob) -ifelse(delta.source[source] == 1, 1, -1)*(log(length(source.comps)) - with(list(r.sel = ifelse(delta.source[source] == 1, r[source.comps][which(delta[source.comps] == 1)], r.star[source.comps][which(delta.star[source.comps] == 1)])), rprop.logdens1(x = r.sel, r.params = r.params)))
		negdifflogproposal <- -log(move.prob.star) + log(move.prob) -ifelse(delta.source[source] == 1, 1, -1)*(log(length(source.comps)) - with(list(r.sel = ifelse(delta.source[source] == 1, r[comp], r.star[comp])), rprop.logdens1(x = r.sel, r.params = r.params)))

	} else if(move.type == 2) { ## randomly select a multi-component source that is in the model and change which component is included

		tmp <- which(neach.group > 1 & delta.source == 1)
		source <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
		source.comps <- which(groups == source)

		oldcomp <- source.comps[delta[source.comps] == 1]
		tmp <- source.comps[delta[source.comps] == 0]
		comp <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))

		r.params.oldcomp <- set.r.params(r.prior = control.params$r.prior, comp = oldcomp, r.params = r.params)
		r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

		delta.star[oldcomp] <- 0
		delta.star[comp] <- 1

		r.star[oldcomp] <- 0
		r.star[comp] <- rprop.gen1(r.params = r.params)

		diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[oldcomp], r.params = r.params.oldcomp)

		negdifflogproposal <- -rprop.logdens1(r.star[comp], r.params = r.params) + rprop.logdens1(r[oldcomp], r.params = r.params.oldcomp)

	} else if(move.type == 3) { ## randomly select a component that is in the model and update it
		tmp <- which(delta == 1)
		comp <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))

		r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

		r.star[comp] <- rprop.gen2(current = r[comp], r.params = r.params)

		diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[comp], r.params = r.params)

		negdifflogproposal <- -rprop.logdens2(r.star[comp], r[comp], r.params = r.params) + rprop.logdens2(r[comp], r.star[comp], r.params = r.params)
	}

	lambda.star <- lambda

	## M-H step
	return(MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps, modifier = modifier))
}

lambda.update <- function(r, delta, lambda, whichcomp=1, y, X, Z = Z, beta, sigsq.eps, Vcomps, data.comps, control.params, modifier = NULL) {
	lambda.jump <- control.params$lambda.jump[whichcomp]
	mu.lambda <- control.params$mu.lambda[whichcomp]
	sigma.lambda <- control.params$sigma.lambda[whichcomp]
	lambdacomp <- lambda[whichcomp]

	## generate a proposal
	lambdacomp.star <- rgamma(1, shape=lambdacomp^2/lambda.jump^2, rate=lambdacomp/lambda.jump^2)
	r.star <- r
	delta.star <- delta
	move.type <- NA

	## part of M-H ratio that depends on the proposal distribution
	negdifflogproposal <- -dgamma(lambdacomp.star, shape=lambdacomp^2/lambda.jump^2, rate=lambdacomp/lambda.jump^2, log=TRUE) + dgamma(lambdacomp, shape=lambdacomp.star^2/lambda.jump^2, rate=lambdacomp.star/lambda.jump^2, log=TRUE)

	## prior distribution
	diffpriors <- dgamma(lambdacomp.star, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE) - dgamma(lambdacomp, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE)

	lambda.star <- lambda
	lambda.star[whichcomp] <- lambdacomp.star

	## M-H step
	return(MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps, modifier = modifier))
}

MHstep <- function(r, lambda, lambda.star, r.star, delta, delta.star, y, X, Z, beta, sigsq.eps, diffpriors, negdifflogproposal, Vcomps, move.type, data.comps, modifier = NULL) {
	## compute log M-H ratio
	Vcomps.star <- makeVcomps(r.star, lambda.star, Z, data.comps, modifier = modifier)
	mu <- y - X%*%beta
	if(data.comps$gs.sig){
	  sigs <- sigsq.eps
	  names(sigs) <- data.comps$levels
	  sigsq.eps.mh <- sigs[apply(modifier, 1, paste, collapse = "")]
	}else{
	  sigsq.eps.mh <- rep(sigsq.eps, nrow(X)) #if single sigsq.eps, make vector of the same value
	}
	diffliks <- 1/2*Vcomps.star$logdetVinv - 1/2*Vcomps$logdetVinv - 1/2*crossprod(mu, Vcomps.star$Vinv - Vcomps$Vinv)%*% diag((sigsq.eps.mh)^(-1))%*%mu #added %*% diag((sigsq.eps)^(-1)) to generalize for both group-specific and one error variance
	logMHratio <- diffliks + diffpriors + negdifflogproposal
	logalpha <- min(0,logMHratio)

	## return value
	acc <- FALSE
	if( log(runif(1)) <= logalpha ) {
		r <- r.star
		delta <- delta.star
		lambda <- lambda.star
		Vcomps <- Vcomps.star
		acc <- TRUE
	}
	return(list(r=r, lambda=lambda, delta=delta, acc=acc, Vcomps=Vcomps, move.type=move.type))
}

h.update <- function(lambda, Vcomps, sigsq.eps, y, X, beta, r, Z, data.comps, modifier=NULL, kernel.method) {
  
  #set the modifier that is passed to Vcomps to be null if using one-kernel method, because it is not needed to block the kernel
  #and set modifier passed to Vcomps to be the modifier if using two-kernel method, because it is needed to block the kernel
  if(kernel.method == "one"){
    kern_modifier <- NULL
  }else if(kernel.method == "two"){
    kern_modifier <- modifier
  }
  
  if (is.null(Vcomps)) {
    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps, modifier = kern_modifier)
  }
	if(is.null(Vcomps$Q)) {
		Kpart <- makeKpart(r, Z)
		K <- exp(-Kpart)
		if(kernel.method == "two"){
		  K <- block_kernel(mod_vec1 = modifier,
		                    mod_vec2 = modifier,
		                    K = K)
		}
		Vinv <- Vcomps$Vinv
		scalar_lambda <- lambda_scalar(mod1 = modifier,
		                               mod2 = modifier,
		                               lambda = lambda, 
		                               data.comps = data.comps) ## in case with random intercept (randint==TRUE), don't include lambda for random intercept matrix
		lamKVinv <- scalar_lambda*K%*%Vinv
		h.postmean <- lamKVinv%*%(y-X%*%beta)
		##h.postvar <- sigsq.eps*lamKVinv
		if(gs.sig){
		  sigs <- sigsq.eps
		  names(sigs) <- lvls
		  sigsq.eps.h <- sigs[apply(modifier, 1, paste, collapse = "")]
		}else{
		  sigsq.eps.h <- rep(sigsq.eps, nrow(X)) #if single sigsq.eps, make vector of the same value
		}
		h.postvar <- scalar_lambda* (sigsq.eps.h*(K - lamKVinv%*%K))
		h.postvar.sqrt <- try(chol(h.postvar), silent=TRUE)
		if(inherits(h.postvar.sqrt, "try-error")) {
			sigsvd <- svd(h.postvar)
			h.postvar.sqrt <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
		}
		hsamp <- h.postmean + crossprod(h.postvar.sqrt, rnorm(length(h.postmean)))
		hcomps <- list(hsamp = hsamp)
	} else {## predictive process approach
	  if(data.comps$gs.tau){
	    stop("Predictive process approach not available for group-separable tau model.")
	  }
		h.star.postvar.sqrt <- sqrt(sigsq.eps*lambda)*forwardsolve(t(Vcomps$cholR), Vcomps$Q)
		h.star.postmean <- lambda[1]*Vcomps$Q %*% Vcomps$Rinv %*% Vcomps$K10 %*% (y - X %*% beta)
		hsamp.star <- h.star.postmean + crossprod(h.star.postvar.sqrt, rnorm(length(h.star.postmean)))
		hsamp <- t(Vcomps$K10) %*% Vcomps$Qinv %*% hsamp.star
		hcomps <- list(hsamp = hsamp, hsamp.star = hsamp.star)
	}
	hcomps
}

newh.update <- function(Z, Znew, mod_new, Vcomps, lambda, sigsq.eps, r, y, X, beta, data.comps, modifier = NULL, kernel.method) {

  #set the modifier that is passed to Vcomps to be null if using two-kernel method (because modifier is not in kernel)
  if(kernel.method == "one"){
    kern_modifier <- NULL
  }else if(kernel.method == "two"){
    kern_modifier <- modifier
  }
  
	if(is.null(data.comps$knots)) {
		n0 <- nrow(Z)
		n1 <- nrow(Znew)
		nall <- n0 + n1
		# Kpartall <- makeKpart(r, rbind(Z, Znew))
		# Kmat <- exp(-Kpartall)
		# Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
		# Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
		# Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
		Kmat1 <- exp(-makeKpart(r, Znew))
		if(kernel.method == "two"){
		  Kmat1 <- block_kernel(mod_vec1 = mod_new,
		                        mod_vec2 = mod_new,
		                        K = Kmat1)
		}
		Kmat10 <- exp(-makeKpart(r, Znew, Z))
		if(kernel.method == "two"){
		  Kmat10 <- block_kernel(mod_vec1 = mod_new,
		                         mod_vec2 = modifier,
		                         K = Kmat10)
		}

		if(is.null(Vcomps)) {
			Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps, modifier = kern_modifier)
		}
		Vinv <- Vcomps$Vinv
    
		scalar_lambda <- lambda_scalar(mod1 = mod_new, 
		                               mod2 = modifier,
		                               lambda = lambda,
		                               data.comps = data.comps)
		lamK10Vinv <- scalar_lambda*Kmat10 %*% Vinv
		scalar_lambda <- lambda_scalar(mod1 = mod_new,
		                               mod2 = mod_new,
		                               lambda = lambda,
		                               data.comps = data.comps)
		if(data.comps$gs.sig){
		  sigs <- sigsq.eps
		  names(sigs) <- data.comps$levels
		  sigsq.eps.h <- sigs[apply(mod_new, 1, paste, collapse = "")]
		}else{
		  sigsq.eps.h <- rep(sigsq.eps, nrow(Znew)) #if single sigsq.eps, make vector of the same value
		}
		Sigma.hnew <- scalar_lambda*(sigsq.eps.h*(Kmat1 - lamK10Vinv %*% t(Kmat10)))
		mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
		root.Sigma.hnew <- try(chol(Sigma.hnew), silent=TRUE)
		if(inherits(root.Sigma.hnew, "try-error")) {
			sigsvd <- svd(Sigma.hnew)
			root.Sigma.hnew <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
		}
		hsamp <- mu.hnew + crossprod(root.Sigma.hnew, rnorm(n1))
	} else {
	  if(data.comps$gs.tau){
	    stop("Predictive process approach not available for group-separable tau model.")
	  }
		n0 <- nrow(data.comps$knots)
		n1 <- nrow(Znew)
		nall <- n0 + n1
		# Kpartall <- makeKpart(r, rbind(data.comps$knots, Znew))
		# Kmat <- exp(-Kpartall)
		# Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
		# Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
		# Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
		Kmat10 <- exp(-makeKpart(r, Znew, data.comps$knots))
		if(!is.null(modifier)){
		  Kmat10 <- block_kernel(mod_vec1 = modifier,
		                         mod_vec2 = modifier,
		                         K = Kmat10)
		}

		if(is.null(Vcomps)) {
			Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps, modifier = kern_modifier)
			h.star.postvar.sqrt <- sqrt(sigsq.eps*lambda[1])*forwardsolve(t(Vcomps$cholR), Vcomps$Q)
			h.star.postmean <- lambda[1]*Vcomps$Q %*% Vcomps$Rinv %*% Vcomps$K10 %*% (y - X %*% beta)
			Vcomps$hsamp.star <- h.star.postmean + crossprod(h.star.postvar.sqrt, rnorm(length(h.star.postmean)))
		}
		hsamp <- Kmat10 %*% Vcomps$Qinv %*% Vcomps$hsamp.star
	}

	hsamp
}