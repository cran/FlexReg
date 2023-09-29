#' Coefficient Methods for flexreg Objects
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param ... additional arguments. Currently not used.
#'
#' @rdname summary.flexreg
#'
#' @export
#'
coef.flexreg <- function(object, ...){
  summ <- summary.flexreg(object)
  mu.model <- summ$Summary.mu[,1]
  theta.model <- summ$Summary.theta[,1]
  phi.model <- summ$Summary.phi[,1]
  q0.model <- summ$Summary.q0[,1]
  q1.model <- summ$Summary.q1[,1]

  names(mu.model) <- rownames(summ$Summary.mu)
  names(theta.model) <- rownames(summ$Summary.theta)
  names(phi.model) <- rownames(summ$Summary.phi)
  names(q0.model) <- rownames(summ$Summary.q0)
  names(q1.model) <- rownames(summ$Summary.q1)

  l <- list("mean_model" = mu.model,
            "overdispersion_model" = theta.model,
            "precision_model" = phi.model,
            "zero_augmentation" = q0.model,
            "one_augmentation" = q1.model)

  if(!is.null(summ$Summary.add)){
    additional.par <- summ$Summary.add[,1]
    l[[length(l)+1]] <- additional.par
    names(l)[length(l)] <- "additional_par"
  }
  l <- l[which(!sapply(l, is.null))]#elimina eventuali null dalla lista (se non c'è zero/one augmentation)
  return(l)
}

