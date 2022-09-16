#' @title var.fun
#' @export
#' @keywords internal
#'

var.fun <- function(model, mu.chain, phi.chain,q0.chain, q1.chain){
  posterior <- model$model
  model.name <- posterior@model_name

  if("Beta" %in% substr(model.name,1,4)){
    var1 <- NULL; var2 <- NULL
    if (is.na(dim(phi.chain)[2])){
      cond.variance <- apply(mu.chain*(1-mu.chain),2, function(x) x/(1+phi.chain))
    } else  {cond.variance <- (mu.chain*(1-mu.chain))/(1+phi.chain)}
  }

  if("VIB" %in% substr(model.name,1,3)){
    p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
    k.chain <- rstan::extract(posterior, pars="k", permuted=T)[[1]]

    if (is.na(dim(phi.chain)[2])){
      var1 <- apply(mu.chain*(1-mu.chain),2, function(x) p.chain*x/(1+phi.chain*k.chain))
      var2 <- apply(mu.chain*(1-mu.chain),2, function(x) (1-p.chain)*x/(1+phi.chain))
      cond.variance  <- var1+var2
    } else {
      var1 <- (mu.chain*(1-mu.chain))/apply(phi.chain, 2, function(x) (1+x*k.chain)/p.chain)
      var2 <- apply((mu.chain*(1-mu.chain))/(1+phi.chain), 2, function(x) (1-p.chain)*x)
      cond.variance  <- var1+var2
    }
  }

  if("FB" %in% substr(model.name,1,2)){
    p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
    w.chain <- rstan::extract(posterior, pars="w", permuted=T)[[1]]
    wtilde.chain <- apply(mu.chain, 2, function(x) w.chain*pmin(x/p.chain, (1-x)/(1-p.chain)))
    lambda1.chain <- mu.chain+ apply(wtilde.chain, 2, function(x) (1-p.chain)*x)
    lambda2.chain <- mu.chain+ apply(wtilde.chain, 2, function(x) p.chain*x)
    if (is.na(dim(phi.chain)[2])){
      var1 <- apply(lambda1.chain*(1-lambda1.chain),2, function(x) p.chain*x/(1+phi.chain))
      var2 <- apply(lambda2.chain*(1-lambda2.chain),2, function(x) (1-p.chain)*x/(1+phi.chain))
      cond.variance  <- apply(mu.chain*(1-mu.chain)+ apply(wtilde.chain, 2, function(x) x^2*phi.chain*p.chain*(1-p.chain)), 2, function(x) x/(1+phi.chain))
    } else {
      var1 <- (lambda1.chain*(1-lambda1.chain))/(1+phi.chain)
      var2 <- (lambda2.chain*(1-lambda2.chain))/(1+phi.chain)
      cond.variance <- (mu.chain*(1-mu.chain)+apply(wtilde.chain^2*phi.chain,2, function(x) x*p.chain*(1-p.chain)))/(1+phi.chain)
    }
  }
  q2.chain <- (1-q0.chain-q1.chain)
  variance <- q2.chain*cond.variance+q1.chain+q2.chain*mu.chain^2-(q1.chain+q2.chain*mu.chain)^2
  return(list(variance=variance, var1=var1, var2=var2))
}

#' @title mu.chain.nd
#' @export
#' @keywords internal
#'
mu.chain.nd <- function(posterior, newdata, link.mu){
  beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
  eta.chain <- beta.chain %*% t(newdata)
  if(link.mu == "logit") mu.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x))) else
    if(link.mu == "probit") mu.chain <- apply(eta.chain,c(1,2), function(x) pnorm(x)) else
      if(link.mu == "cloglog") mu.chain <- apply(eta.chain,c(1,2), function(x) 1-exp(-exp(x))) else
        if(link.mu == "loglog") mu.chain <- apply(eta.chain,c(1,2), function(x) exp(-exp(x)))
  return(mu.chain)
}

#' @title q.chain.nd
#' @export
#' @keywords internal
#'
q.chain.nd <- function(model, newdata.q0, newdata.q1){
  posterior <- model$model
  if(is.null(model$call$zero.formula) & is.null(model$call$one.formula)){
    q0.chain <- NULL
    q1.chain <- NULL
  } else if(!is.null(model$call$zero.formula) & is.null(model$call$one.formula)) {
    q0.chain <- q0.chain.nd(posterior, newdata.q0)
    q1.chain <- NULL
  } else if(is.null(model$call$zero.formula) & !is.null(model$call$one.formula)){
    q0.chain <- NULL
    q1.chain <- q1.chain.nd(posterior, newdata.q1)
  } else if(!is.null(model$call$zero.formula) & !is.null(model$call$one.formula)){
    q0.chain <- q01.chain.nd(posterior, newdata.q0, newdata.q1)[[1]]
    q1.chain <- q01.chain.nd(posterior, newdata.q0, newdata.q1)[[2]]
  }
  return(list(q0.chain, q1.chain))
}


#' @title q0.chain.nd
#' @export
#' @keywords internal
#'
q0.chain.nd <- function(posterior, newdata){
  omega0.chain <- rstan::extract(posterior, pars="omega0", permuted=T)[[1]]
  eta.chain <- omega0.chain %*% t(newdata)
  q0.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x)))
    return(q0.chain)
}

#' @title q1.chain.nd
#' @export
#' @keywords internal
#'
q1.chain.nd <- function(posterior, newdata){
  omega1.chain <- rstan::extract(posterior, pars="omega1", permuted=T)[[1]]
  eta.chain <- omega1.chain %*% t(newdata)
  q1.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x)))
    return(q1.chain)
}

#' @title q1.chain.nd
#' @export
#' @keywords internal
#'
q01.chain.nd <- function(posterior, newdata.q0, newdata.q1){
  omega0.chain <- rstan::extract(posterior, pars="omega0", permuted=T)[[1]]
  omega1.chain <- rstan::extract(posterior, pars="omega1", permuted=T)[[1]]
  eta0.chain <- exp(omega0.chain %*% t(newdata.q0))
  eta1.chain <- exp(omega1.chain %*% t(newdata.q1))
  eta.chain <- cbind(eta0.chain,eta1.chain)
  q0.chain <- eta0.chain/(1+rowSums(eta.chain))
  q1.chain <- eta1.chain/(1+rowSums(eta.chain))
  return(list(q0.chain, q1.chain))
}

#' @title phi.chain.nd
#' @export
#' @keywords internal
#'
phi.chain.nd <- function(posterior, newdata, link.phi){
  if(link.phi == "identity") {
    phi.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]} else {
      psi.chain <- rstan::extract(posterior, pars="psi", permuted=T)[[1]]
      eta.chain <- psi.chain %*% t(newdata)
      if(link.phi == "log") phi.chain <- apply(eta.chain,c(1,2), function(x) exp(x)) else
        if(link.phi == "sqrt") phi.chain <- apply(eta.chain,c(1,2), function(x) x^2)
    }
  return(phi.chain)
}


#' @title theta.chain.nd
#' @export
#' @keywords internal
#'
theta.chain.nd <- function(posterior, newdata, link.theta){
  if(is.null(link.theta)) {
    theta.chain <- rstan::extract(posterior, pars="theta", permuted=T)[[1]]
    } else {
      psi.chain <- rstan::extract(posterior, pars="psi", permuted=T)[[1]]
      eta.chain <- psi.chain %*% t(newdata)
      if(link.theta == "logit") theta.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x))) else
        if(link.theta == "probit") theta.chain <- apply(eta.chain,c(1,2), function(x) pnorm(x)) else
          if(link.theta == "cloglog") theta.chain <- apply(eta.chain,c(1,2), function(x) 1-exp(-exp(x))) else
            if(link.theta == "loglog") theta.chain <- apply(eta.chain,c(1,2), function(x) exp(-exp(x)))
    }
  return(theta.chain)
}
