#' Trace plot
#'
#' @inheritParams ExtractEsts
#' @param par which parameter to plot
#' @param comp which component of the parameter vector to plot
#' @param main title
#' @param xlab x axis label
#' @param ylab y axis label
#' @param ... other arguments to pass onto the plotting function
#' @export
#' @import graphics
#' @details For guided examples, see vignette(bkmrGSOverview)
#' 
#' @return No return value, generates plot
#' 
#' @examples
#' ## First generate data set
#' y <- ex_data$y
#' Z <- ex_data$Z
#' modifier <- ex_data$X$Sex
#' X_full <- ex_data$X[,-2] #remove Sex from the covariate matrix because it is the modifier
#' #create design matrix to account for factor variables, remove the intercept column
#' X <- model.matrix(~., data=X_full)[,-1] 
#' 
#' ## Fit model 
#' ## Using only 10 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, modifier = modifier, X = X, iter = 10, verbose = FALSE) 
#' 
#' TracePlot(fit = fitkm, par = "beta")
#' TracePlot(fit = fitkm, par = "lambda")
#' TracePlot(fit = fitkm, par = "sigsq.eps")
#' TracePlot(fit = fitkm, par = "r", comp = 1)
TracePlot <- function(fit, par, comp = 1, sel = NULL, main = "", xlab = "iteration", ylab = "parameter value", ...) {
    samps <- ExtractSamps(fit, sel = sel)[[par]]
    if (!is.null(ncol(samps))) {
        nm <- colnames(samps)[comp]
        samps <- samps[, comp]
    } else {
        nm <- par
    }
    main <- paste0(main, "\n(", nm, " = ", format(mean(samps), digits = 2), ")")
    plot(samps, type = "l", main = main, xlab = xlab,  ylab = ylab, ...)
    abline(h = mean(samps), col = "blue", lwd = 2)
}
