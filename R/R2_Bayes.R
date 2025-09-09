#' @title Bayesian R-squared for \code{flexreg} Objects
#'
#' @description Bayesian version of R-squared for flexible regression models for bounded continuous and discrete responses.
#'
#' @param model an object (or a list of objects) of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#'
#' @return A list with the same length as the number of objects of class \code{`flexreg`} passed in the \code{model} argument.
#' Each element of the list contains a vector of  Bayesian R-squared values with the same length as the Markov Chain(s) after warmup.
#'
#' @details  The function provides a Bayesian version of the R-squared measure, defined as the variance of the predicted values divided by itself plus the expected variance of the errors.
#' @references {
#' Gelman, A., Goodrich, B., Gabry, J., Vehtari, A. (2019). R-squared for Bayesian Regression Models, The American Statistician, 73:3, 307--309. doi: 10.1080/00031305.2018.1549100
#' }
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, type = "FB", n.iter=1000)
#' hist(R2_bayes(FB)[[1]])
#' }
#'
#' @export
#'
#' @import rstantools
#'

R2_bayes <- function(model){
  x <- model
  if(!inherits(x, "flexreg")){
    if(any(unlist(lapply(x, function(x) !inherits(x, "flexreg")))))
      stop("The argument must be an object (or a list of objects) of class `flexreg`")
  }

  if(inherits(x, "flexreg")){
    out <- list(R2Bayes_internal(x)) #aggiungo list perchè dichiaro che restituisco una lista
  } else {
    out <- lapply(x, R2Bayes_internal) #restituisce una lista di vettori
  }
  return(out)
}

#' internal function
#' @keywords internal

R2Bayes_internal <- function(x){
  y <- x$response
  X <- x$design.X
  posterior <- x$model
  aug <- x$aug

  if(length(posterior)>1){
      mu.chain <-  predict_mu.chain(x, X)
      if(aug == "01"){
      q.chain <- rstan::extract(posterior[[2]], pars="q", permuted=T)[[1]]
      y_tilde <- q.chain[,,2]+ q.chain[,,3] * mu.chain
      } else {
        if(aug == "1"){
        q1.chain <- rstan::extract(posterior[[2]], pars="mu", permuted=T)[[1]]
        q0.chain <- 0
        }
        if(aug == "0"){
          q0.chain <- rstan::extract(posterior[[2]], pars="mu", permuted=T)[[1]]
          q1.chain <- 0
        }
        y_tilde <- q1.chain+ (1-q0.chain-q1.chain) * mu.chain
      }
  } else {
    y_tilde <- rstan::extract(posterior[[1]], pars="mu", permuted=T)[[1]]
  }

  #bayes_R2(y_tilde, y) #si potrebbe semplicemente usare la funzione di rstantools
  var_fit <- apply(y_tilde, 1, var)

  # Calcolo la distribuzione dei residui per ogni unita':
  res <- apply(rbind(y,y_tilde), 2, function(x) x[-1]-x[1])
  # Varianza dei residui condizionatamente al parametro simulato:
  var_res <- apply(res, 1, var)
  R2 <- var_fit/(var_fit + var_res)
  return(R2)
}





