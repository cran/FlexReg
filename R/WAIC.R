#' @title WAIC and LOO
#'
#' @description The function computes widely applicable information criterion (WAIC) and efficient approximate leave-one-out cross-validation (LOO) from fitted  regression model objects of class \code{`flexreg`}.
#'
#' @param model an object (or a list of objects) of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param ... additional arguments.
#'
#' @returns A named list with components from \code{\link[loo]{loo}} and \code{\link[loo]{waic}}.
#'
#' @details This function takes advantage of the \pkg{loo} package to compute the widely applicable information criterion (WAIC) and leave-one-out cross-validation (LOO) for objects of class \code{`flexreg`}.
#' If a list of two or more objects of class \code{`flexreg`} is provided, the function returns the difference in their expected predictive accuracy (see \code{\link[loo]{loo_compare}} for further details).
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, type="FB", n.iter=1000)
#' WAIC(FB)
#'}
#'
#' @references {
#' Vehtari, A., Gelman, A., Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. \bold{27}(5), 1413--1432. doi:10.1007/s11222-016-9696-4 \cr
#' \cr
#'
#' }
#'
#' @import loo
#'
#' @export
#'


WAIC <-  function(model,...){
  x <- model
  if(!inherits(x, "flexreg")){
    if(any(unlist(lapply(x, function(x) !inherits(x, "flexreg")))))
      stop("The argument must be an object (or a list of objects) of class `flexreg`")
  }
  if(inherits(x, "flexreg")){
    loglik <- WAIC_internal(model = x$model,  response = x$response)
    waic_out <- suppressWarnings(loo::waic(loglik))
    loo_out <- suppressWarnings(loo::loo(loglik))
  } else {
    #per quando si ha lista di modelli
    loglik <- lapply(x, function(x) WAIC_internal(x$model, x$response))
    loos <- lapply(loglik, function(x) suppressWarnings(loo::loo(x)))
    waics <- lapply(loglik, function(x) suppressWarnings(loo::waic(x)))
    loo_out <- loo::loo_compare(loos)
    waic_out <- loo::loo_compare(waics)
  }
  output <- list(loo_out = loo_out, waic_out = waic_out)
  class(output) <- "WAIC.flexreg"

  return(suppressWarnings(output))
}


#' internal function
#' @keywords internal
#'
WAIC_internal <- function(model, response,...){
  if(length(model)>1){

    ll_cont <- loo::extract_log_lik(model[[1]])
    ll_aug <- loo::extract_log_lik(model[[2]])

    index01 <- which(response==0 | response ==1)
    N <- length(response)
    index_cont <- (1:N)[-index01]
    loglik <- ll_aug
    loglik[,index_cont] <- ll_aug[,index_cont] + ll_cont
  } else {
    loglik <- loo::extract_log_lik(model[[1]])
  }
  return(loglik=loglik)
}


#' Print method for WAIC.flexreg objects
#'
#' @param x an object of class \code{`WAIC.flexreg`}, usually the result of \code{\link{WAIC}}.
#'
#' @rdname WAIC
#'
#' @export
#'

print.WAIC.flexreg <- function(x, ...){
  cat("Waic method:\n")
  print(suppressWarnings(x$waic_out))

  cat("\nLoo method:\n")
  print(suppressWarnings(x$loo_out))
}
