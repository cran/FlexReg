#' @title newdata.adjust
#' @keywords internal
#'

newdata.adjust <- function(newdata, formula){
  newdata <- as.data.frame(newdata)
  Terms <- delete.response(terms(formula))
  n <- newdata$n
  newdata <- (model.frame(Terms, newdata))#if a variable is missing from newdata it returns an error
  #that is clear enough.

  #add the intercept if missing
  if(attr(Terms, "intercept")==1 & ("(Intercept)" %in% colnames(newdata)) ==F){
    newdata$`(Intercept)` <- rep(1, nrow(newdata)) # .. add the intercept
  }
  newdata$n <- n
  return(newdata)
}


#' @title predict_mu.chain
#' @keywords internal
#'
#'
predict_mu.chain <- function(model, newdata){
  posterior <- model$model[[1]]
  if(is.null(newdata)){
    mu.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
  } else {
    link.mu <- model$link.mu
    newdata.X <- newdata[,match(colnames(model$design.X),colnames(newdata))]
    mu.chain <- mu.chain.nd(posterior, newdata.X, link.mu)
  }
  return(mu.chain)
}


#' @title predict_response
#' @keywords internal
#'

predict_response <- function(model, newdata, cluster, n){


  if(is.null(newdata)){
    newdata.temp <- model$design.X
  } else newdata.temp <- newdata

  mu.chain <- predict_mu.chain(model, newdata.temp)

  if("flexreg_bound" %in% class(model)){
    response.binom <- 0

    if(model$aug == "No") {
      response <- 0
      q.chain <- list(q0.chain = 0, q1.chain = 0)
    } else{
      q.chain <- predict_q.chain(model, newdata)
      q2.chain <-  1-q.chain$q0.chain-q.chain$q1.chain
      q.chain <- append(q.chain, list(q2.chain=q2.chain))

      response <- q.chain$q1.chain + q.chain$q2.chain*mu.chain
    }

    if(cluster == T){
      posterior <- model$model[[1]]
      lambda.chain <- predict_lambda.chain(posterior, mu.chain, newdata.temp)
    } else {
      lambda.chain <- list(l1.chain = 0, l2.chain = 0)
    }

  } else { # The model is for binomial data
    response <- 0
    q.chain <- list(q0.chain = 0, q1.chain = 0)

    response.binom <- t(apply(mu.chain, 1, function(x) x*n))

    if(cluster == T){
      posterior <- model$model[[1]]
      lambda.chain <- predict_lambda.chain(posterior, mu.chain, newdata)
    } else {
      lambda.chain <- list(l1.chain = 0, l2.chain = 0)
    }
  }

  return(pred.chain = list(response.binom = response.binom,
                           response.aug = response,
                           response = mu.chain,
                           q0 = q.chain$q0.chain,
                           q1 = q.chain$q1.chain,
                           l1 = lambda.chain$l1.chain,
                           l2 = lambda.chain$l2.chain))
}


#' @title predict_link
#' @keywords internal
#'

predict_link <- function(model,  newdata){
  posterior <- model$model[[1]]
  beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
  X <- model$design.X

  if(!is.null(newdata))  X <- newdata[,match(colnames(X),colnames(newdata))]

  pred.chain  <- list(link=beta.chain %*% t(X))

  return(pred.chain)
}

#' @title predict_precision
#' @keywords internal
#'

predict_precision <- function(model, newdata){
  # n <- length(model$response)
  posterior <- model$model[[1]]
  link.phi <- model$link.phi
  if(is.null(link.phi)) link.phi <- "identity"
  if(link.phi=="identity"){
    pred.chain <- list(precision = as.matrix(rstan::extract(posterior, pars="phi", permuted=T)[[1]]))
   }else{
    if(is.null(newdata)) {
      newdata.Z <- model$design.Z
    }else{
    newdata.Z <- newdata[,match(colnames(model$design.Z),colnames(newdata))]
    }
    pred.chain <- list(precision=as.matrix(phi.chain.nd(posterior, newdata.Z, link.phi)))
  }

  #if(is.na(dim(pred.chain[[1]])[2])) pred.chain <- list(precision=matrix(rep(pred.chain[[1]], n),ncol=n))
  return(pred.chain)
}

#' @title phi.chain.nd
#' @keywords internal
#'
phi.chain.nd <- function(posterior, newdata.Z, link.phi){
  #if(link.phi == "identity") {
   # phi.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]
  #} else {
    psi.chain <- rstan::extract(posterior, pars="psi", permuted=T)[[1]]
    eta.chain <- psi.chain %*% t(newdata.Z)
    if(link.phi == "log") phi.chain <- apply(eta.chain,c(1,2), function(x) exp(x)) else
      if(link.phi == "sqrt") phi.chain <- apply(eta.chain,c(1,2), function(x) x^2)
 # }
  return(phi.chain)
}


#' @title predict_variance
#' @keywords internal
#'

predict_variance <- function(model, newdata, cluster, cluster.var = T, model.type, model.class, n){

  posterior <- model$model[[1]]
  response.chain <- predict_response(model, newdata, cluster, n)
  if("flexreg_binom" %in% class(model)){
    theta.chain <- predict_over(model, newdata)[[1]]
    phi.chain <- NULL
  } else{ # if the model is for continuous data

    phi.chain <- predict_precision(model, newdata)[[1]]
    theta.chain <- NULL
  }

  variance.chain <- var.fun(model.class = model.class,
                            model.type = model.type,
                            posterior = posterior,
                            mu.chain = response.chain$response,
                            q0.chain = response.chain$q0,
                            q1.chain = response.chain$q1,
                            l1.chain = response.chain$l1,
                            l2.chain = response.chain$l2,
                            phi.chain = phi.chain,
                            theta.chain = theta.chain, n = n)
  pred.chain <- list(variance = variance.chain$variance,
                     cluster1 = variance.chain$var1,
                     cluster2 = variance.chain$var2)
  if(cluster.var == F){
    pred.chain <- list(variance = pred.chain$variance)
  }
  return(pred.chain)
}

#' @title predict_over
#' @keywords internal
#'

predict_over <- function(model, newdata){
  posterior <- model$model[[1]]
  #n <- length(model$response)
  if(model$type == "Bin"){
    pred.chain <- NULL
  } else {
    if(is.null(newdata)){
      pred.chain <- list(overdispersion = as.matrix(rstan::extract(posterior, pars="theta", permuted=T)[[1]]))
    } else {
      link.theta <- model$link.theta
      newdata.theta <- newdata[,match(colnames(model$design.Z),colnames(newdata))]
      pred.chain <- list(overdispersion=as.matrix(theta.chain.nd(posterior, newdata.theta, link.theta)))
    }
  }
  return(pred.chain)
}

#' @title mu.chain.nd
#' @keywords internal
#'
mu.chain.nd <- function(posterior, newdata.X, link.mu){
  beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
  eta.chain <- beta.chain %*% t(as.matrix(newdata.X))
  if(link.mu == "logit") mu.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x))) else
    if(link.mu == "probit") mu.chain <- apply(eta.chain,c(1,2), function(x) pnorm(x)) else
      if(link.mu == "cloglog") mu.chain <- apply(eta.chain,c(1,2), function(x) 1-exp(-exp(x))) else
        if(link.mu == "loglog") mu.chain <- apply(eta.chain,c(1,2), function(x) exp(-exp(x)))
  return(mu.chain)
}


#' @title predict_q.chain
#' @keywords internal
#'
predict_q.chain <- function(model, newdata = NULL){
  posterior.aug <- model$model[[2]]

  if(is.null(model$call$zero.formula) & is.null(model$call$one.formula)){
    q0.chain <- 0
    q1.chain <- 0
  } else if(!is.null(model$call$zero.formula) & is.null(model$call$one.formula)) {
    q1.chain <- 0

    if(is.null(newdata)){
      q0.chain <- rstan::extract(posterior.aug, pars="mu", permuted=T)[[1]]
    } else{
      newdata.X0 <- newdata[,match(colnames(model$design.X0),colnames(newdata))]
      q0.chain <- q0.chain.nd(posterior.aug, newdata.X0)
    }

  } else if(is.null(model$call$zero.formula) & !is.null(model$call$one.formula)){
    q0.chain <- 0

    if(is.null(newdata)){
      q1.chain <- rstan::extract(posterior.aug, pars="mu", permuted=T)[[1]]
    } else{
      newdata.X1 <- newdata[,match(colnames(model$design.X1),colnames(newdata))]
      q1.chain <- q1.chain.nd(posterior.aug, newdata.X1)
    }

  } else if(!is.null(model$call$zero.formula) & !is.null(model$call$one.formula)){
    if(is.null(newdata)){
      q.chain <- rstan::extract(posterior.aug, pars="q", permuted=T)[[1]]
      q0.chain <- q.chain[,,1]
      q1.chain <- q.chain[,,2]
    }else{
      newdata.X0 <- newdata[,match(colnames(model$design.X0),colnames(newdata))]
      newdata.X1 <- newdata[,match(colnames(model$design.X1),colnames(newdata))]
      q.chain <- q01.chain.nd(posterior.aug, newdata.X0, newdata.X1)
      q0.chain <- q.chain$q0.chain
      q1.chain <- q.chain$q1.chain
    }
  }

  return(list(q0.chain = q0.chain, q1.chain = q1.chain))
}



#' @title q0.chain.nd
#' @keywords internal
#'
q0.chain.nd <- function(posterior.aug, newdata.X0){
  omega0.chain <- rstan::extract(posterior.aug, pars="beta", permuted=T)[[1]]
  eta.chain <- omega0.chain %*% t(newdata.X0)
  q0.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x)))
  return(q0.chain)
}

#' @title q1.chain.nd
#' @keywords internal
#'
q1.chain.nd <- function(posterior.aug, newdata.X1){
  omega1.chain <- rstan::extract(posterior.aug, pars="beta", permuted=T)[[1]]
  eta.chain <- omega1.chain %*% t(newdata.X1)
  q1.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x)))
  return(q1.chain)
}

#' @title q01.chain.nd
#' @keywords internal
#'
q01.chain.nd <- function(posterior.aug, newdata.X0, newdata.X1){
  omega0.chain <- rstan::extract(posterior.aug, pars="omega0", permuted=T)[[1]]
  omega1.chain <- rstan::extract(posterior.aug, pars="omega1", permuted=T)[[1]]
  eta0.chain <- exp(omega0.chain %*% t(newdata.X0))
  eta1.chain <- exp(omega1.chain %*% t(newdata.X1))
  eta.chain <- cbind(eta0.chain,eta1.chain)
  q0.chain <- eta0.chain/(1+rowSums(eta.chain))
  q1.chain <- eta1.chain/(1+rowSums(eta.chain))
  return(list(q0.chain=q0.chain, q1.chain=q1.chain))
}


#' @title var.fun
#' @keywords internal
#'

var.fun <- function(model.class, model.type, posterior, mu.chain, phi.chain, theta.chain, q0.chain, q1.chain, l1.chain, l2.chain, n){

  if("flexreg_binom" %in% model.class){
    if(model.type == "Bin"){
      var1 <- var2 <- NULL
      cond.variance <- t(apply(mu.chain, 1, function(x) n*(x*(1-x))))
    } else{
      if(ncol(theta.chain) == 1)  theta.chain <- matrix(rep(theta.chain, length(n)), ncol=length(n))#theta.chain <- as.vector(theta.chain)
      if(model.type == "BetaBin"){
        var1 <- var2 <- NULL
        cond.variance <- t(apply(mu.chain,1, function(x) n*(x*(1-x)))) *(1+
                                                                           t(apply(theta.chain, 1, function(x) x*(n-1))))
      } else   if(model.type == "FBB"){
        p.chain <- as.vector(rstan::extract(posterior, pars="p", permuted=T)[[1]])
        m.mu.p <- apply(mu.chain, 2, function(x) pmin((x*(1-p.chain))/(p.chain*(1-x)),
                                                      (p.chain*(1-x))/(x*(1-p.chain))))

        w.chain <- as.vector(rstan::extract(posterior, pars="w", permuted=T)[[1]])

        if(is.null(dim(l1.chain ))) {
          var1 <- var2 <- NULL
        } else {
          var1 <- t(apply(l1.chain,1, function(x) n*(x*(1-x)))) *(1+t(apply(theta.chain, 1, function(x) x*(n-1))))
          var2 <- t(apply(l2.chain,1, function(x) n*(x*(1-x)))) *(1+t(apply(theta.chain, 1, function(x) x*(n-1))))
        }

        cond.variance <- t(apply(mu.chain,1, function(x) n*(x*(1-x))))*
          (1+t(apply(theta.chain, 1, function(x) x*(n-1)))+t(apply(theta.chain, 1, function(x) x*(n-1)))*w.chain^2*m.mu.p)
      }
    }
  } else {

    #check for collapsing phi.chain, if necessary
    if(ncol(phi.chain) == 1) phi.chain <- as.vector(phi.chain)

    if(model.type == "Beta"){
      var1 <- var2 <- NULL
      cond.variance <- (mu.chain*(1-mu.chain))/(1+phi.chain)
    }

    if(model.type == "VIB"){
      p.chain <- as.vector(rstan::extract(posterior, pars="p", permuted=T)[[1]])
      k.chain <- as.vector(rstan::extract(posterior, pars="k", permuted=T)[[1]])

      var1 <- (mu.chain * (1-mu.chain)) / (1+phi.chain * k.chain)
      var2 <- (mu.chain * (1-mu.chain)) / (1+phi.chain)
      cond.variance  <- p.chain * var1 + (1 - p.chain) * var2
    }

    if(model.type == "FB"){

      var1 <- (l1.chain*(1-l1.chain))/(1+phi.chain)
      var2 <- (l2.chain*(1-l2.chain))/(1+phi.chain)

      p.chain <- as.vector(rstan::extract(posterior, pars="p", permuted=T)[[1]])
      w.chain <- as.vector(rstan::extract(posterior, pars="w", permuted=T)[[1]])

      wtilde.chain <- apply(mu.chain, 2, function(x) w.chain*pmin(x/p.chain, (1-x)/(1-p.chain)))
      cond.variance <- (mu.chain*(1-mu.chain)+apply(wtilde.chain^2*phi.chain,2, function(x) x*p.chain*(1-p.chain)))/(1+phi.chain)
    }
  }
  q2.chain <- (1-q0.chain-q1.chain)
  variance <- q2.chain*cond.variance+q1.chain+q2.chain*mu.chain^2-(q1.chain+q2.chain*mu.chain)^2
  return(list(variance = variance, cond.variance = cond.variance, var1 = var1, var2 = var2))
}


#' @title predict_lambda.chain
#' @keywords internal
#'
predict_lambda.chain <- function(posterior, mu.chain, newdata){
  if(is.null(newdata)){
    l1.chain <- rstan::extract(posterior, pars="lambda1", permuted=T)[[1]]
    l2.chain <- rstan::extract(posterior, pars="lambda2", permuted=T)[[1]]
  } else{

    p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
    w.chain <- rstan::extract(posterior, pars="w", permuted=T)[[1]]
    parz.min <- pmin(apply(mu.chain, 2, function(x) x/p.chain) , apply(1-mu.chain, 2, function(x) x/(1-p.chain)))

    l1.chain <- mu.chain + apply(parz.min,2, function(x) x*(1-p.chain)*w.chain)
    l2.chain <- mu.chain - apply(parz.min,2, function(x) x*p.chain*w.chain)
  }
  return(list(l1.chain = l1.chain, l2.chain = l2.chain))
}


#' @title theta.chain.nd
#' @keywords internal
#'
theta.chain.nd <- function(posterior, newdata.theta, link.theta){
  if(link.theta == "identity") {
    theta.chain <- rstan::extract(posterior, pars="theta", permuted=T)[[1]]
  } else {
    psi.chain <- rstan::extract(posterior, pars="psi", permuted=T)[[1]]
    eta.chain <- psi.chain %*% t(newdata.theta)
    if(link.theta == "logit") theta.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x))) else
      if(link.theta == "probit") theta.chain <- apply(eta.chain,c(1,2), function(x) pnorm(x)) else
        if(link.theta == "cloglog") theta.chain <- apply(eta.chain,c(1,2), function(x) 1-exp(-exp(x))) else
          if(link.theta == "loglog") theta.chain <- apply(eta.chain,c(1,2), function(x) exp(-exp(x)))
  }
  return(theta.chain)
}

#' @title extract.pars
#' @keywords internal
#'
extract.pars <- function(posterior){
  pars.full <- names(posterior)
  pars <- c()
  pars <- c(pars, pars.full[grep("beta",pars.full)])
  pars <- c(pars, pars.full[grep("psi",pars.full)])#if model.phi is TRUE or if model.theta is TRUE
  pars <- c(pars, pars.full[which(pars.full=="phi")])#if model.phi or model.theta is FALSE
  pars <- c(pars, pars.full[grep("omega0",pars.full)])#for 0 augmentation
  pars <- c(pars, pars.full[grep("omega1",pars.full)])#for 1 augmentation
  pars <- c(pars, pars.full[which(pars.full=="theta")])#if model.theta is FALSE
  pars <- c(pars, pars.full[which(pars.full=="p")])#if type is FB, VIB, or FBB
  pars <- c(pars, pars.full[which(pars.full=="w")])#if type is FB or FBB
  pars <- c(pars, pars.full[which(pars.full=="k")])#if type is VIB
  return(pars)
}


#' @title rate_plot
#' @keywords internal
#'
#plot for rate of convergence
rate_plot <- function(chains, pars, n.warmup = n.warmup){
  S <- dim(chains)[1]#n.iter
  n.chain <- dim(chains)[2]#n.chain
  sum.parz <- apply(chains, c(2,3), cumsum)
  #dim(sum.parz) <- dim(chains)
  mean.parz <- apply(sum.parz, 3, function(x) x/(1:S))
  data.plot <- as.data.frame(mean.parz)
  names(data.plot) <- pars

  if (n.chain > 1) {
    data.plot$Chain <- as.factor(rep(1:n.chain, each=S))
    data.plot$iter  <- rep( 1:S, n.chain)
  } else  {
    data.plot$iter  <- 1:S
    data.plot$Chain <- as.factor(rep(1, S))
  }

  plot.out <- lapply(pars, function(g)
    ggplot(data.plot,
           aes(x=data.plot$iter, y=!!as.name(g), color= data.plot$Chain))+
      annotate("rect",xmin = 0,xmax = n.warmup,
               ymin = -Inf, ymax = Inf, fill =  "grey20", alpha = 0.1)+
      geom_line()+ theme(legend.position = 'none',panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y=element_text(angle=0,hjust=1),
                         axis.line = element_line(colour = "black"))+
      ggtitle(paste("Rate plot of ", g, sep=" ")))
  return(plot.out)
}
