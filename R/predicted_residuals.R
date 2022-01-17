#' @title Prediction Method for flexreg Objects
#'
#' @description Method that computes various types of prediction from objects of class \code{`flexreg`}. If the model type is \code{FB} or \code{FBB} and \code{cluster = T}, the function returns also cluster means.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}}.
#' @param newdata an optional data frame containing variables with which to predict. If omitted, the fitted values are used.
#' @param cluster logical. If the model is \code{"FB"} or \code{"FBB"}, \code{cluster = T} returns the cluster means. By default, \code{cluster = F}.
#' @param type a character indicating the type of predictions. Available options are the  fitted means of response/relative response (\code{response}), the linear predictor (\code{link}),
#' the fitted precision parameter phi (\code{precision}), the fitted overdispersion parameter theta (\code{overdispersion}), and the fitted variances of the response (\code{variance}).
#' @param estimate the type of estimate: \code{mean} (default), \code{median} or \code{quantile}.
#' @param q if estimate is \code{quantile}, numeric value of probability in (0, 1).
#' @param ... additional arguments. Currently not used.
#'
#' @examples{
#' data("Reading")
#' FB <- flexreg(accuracy ~ iq, Reading, type="FB", n.iter=1000)
#' predict(FB, type="response", cluster=TRUE)
#' }
#'
#' @references {
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018) A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' \cr
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' }
#'
#' @import  stats
#'
#' @export
#'
#'

predict.flexreg <- function(object, newdata=NULL, cluster=F,
                            type = "response",
                            estimate="mean", q=NULL, ...){
  model <- object
  model.name <- model$model@model_name
  n <- length(model$response)

  if(!(model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta")) & cluster==T)
    stop("Cluster means are available only for FB and FBB models")
  if(!(estimate %in% c("mean", "median", "quantile")))
    stop("Argument `estimate` must be set equal to `mean`, `median` or `quantile`")
  if(estimate=="quantile" & is.null(q))
    stop("You have to specify the order of the quantile")
  if(!(type %in% c("response", "link", "precision", "variance", "overdispersion")))
    stop("Argument `type` must be set equal to `response`, `link`, `precision`, `variance`, or `overdispersion`")

  if(type=="variance" & !(model.name %in% c("Beta", "FB", "VIB", "Beta_phi", "FB_phi", "VIB_phi")))
    stop("predicted variances are available only for Beta, FB, and VIB models")
  if(type=="overdispersion" & !(model.name %in% c("BetaBin", "BetaBin_theta", "FBB", "FBB_theta")))
    stop("predicted ovedispersion is available only for BetaBin and FBB models")
  if(type=="precision" & !(model.name %in% c("Beta", "FB", "VIB", "Beta_phi", "FB_phi", "VIB_phi")))
    stop("predicted precision is available only for Beta, FB, and VIB models")

  posterior <- model$model
  # If we do not have new data:
  if(is.null(newdata)){
    if( type == "response") pred.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
    if( type == "link") {
      beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
      X <- model$design.X
      pred.chain <- beta.chain %*% t(X)
    }

    if( type == "precision" & model.name %in% c("Beta", "FB", "VIB",
                                                "Beta_phi", "FB_phi", "VIB_phi")) {
      pred.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }
    if( type == "variance" & model.name %in% c("Beta", "FB", "VIB",
                                               "Beta_phi", "FB_phi", "VIB_phi")) {
      mu.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
      phi.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]
      pred.chain <- var.fun(model, mu.chain, phi.chain)
    }
    if( type == "overdispersion" & model.name %in% c("BetaBin", "BetaBin_theta", "FBB", "FBB_theta")) {
      pred.chain <- rstan::extract(posterior, pars="theta", permuted=T)[[1]]
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }

    # If cluster == T, report the predicted lambda without considering the type argument
    if ((model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta")) & cluster==T){
      l1.chain <- rstan::extract(posterior, pars="lambda1", permuted=T)[[1]]
      l2.chain <- rstan::extract(posterior, pars="lambda2", permuted=T)[[1]]
    }
  } else { # cioè se newdata non è NULL --> Se ho dati da prevedere
    newdata <- as.data.frame(newdata)
    if((c("(Intercept)") %in% colnames(newdata)) ==F & # If newdata does not have the intercept..
       ((c("(Intercept)") %in% colnames(model$design.X)) ==T) |
       (c("(Intercept)") %in% colnames(model$design.Z) ==T)) # .. and the intercept is required..
      newdata$`(Intercept)` <- rep(1, nrow(newdata)) # .. add the intercept

    if(!all(c(colnames(model$design.X),colnames(model$design.Z))  %in% colnames(newdata))){
      stop("newdata must be a data frame containing the same predictors as the ones in formula")
    }
    n <- nrow(newdata)
    if(type == "response") {
      link.mu <- model$link.mu
      newdata.mu <- newdata[,which(colnames(newdata) %in% colnames(model$design.X) )]
      pred.chain <- mu.chain.nd(posterior, newdata.mu, link.mu)
    }
    if(type == "link") {
      beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
      newdata.mu <- newdata[,which(colnames(newdata) %in% colnames(model$design.X) )]
      pred.chain  <- beta.chain %*% t(newdata.mu)
    }
    if(type == "precision" & model.name %in% c("Beta", "FB", "VIB",
                                               "Beta_phi", "FB_phi", "VIB_phi")){
      link.phi <- model$link.phi
      newdata.phi <- newdata[,which(colnames(newdata) %in% colnames(model$design.Z) )]
      pred.chain <- phi.chain.nd(posterior, newdata.phi, link.phi)
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }
    if(type == "overdispersion" & model.name %in% c("BetaBin", "BetaBin_theta", "FBB", "FBB_theta")) {
      link.theta <- model$link.theta
      newdata.theta <- newdata[,which(colnames(newdata) %in% colnames(model$design.Z) )]
      pred.chain <- theta.chain.nd(posterior, newdata.theta, link.theta)
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }
    if(type == "variance"& model.name %in% c("Beta", "FB", "VIB",
                                             "Beta_phi", "FB_phi", "VIB_phi")){
      link.mu <- model$link.mu
      newdata.mu <- newdata[,which(colnames(newdata) %in% colnames(model$design.X) )]
      mu.chain <- mu.chain.nd(posterior, newdata.mu, link.mu)
      link.phi <- model$link.phi
      newdata.phi <- newdata[,which(colnames(newdata) %in% colnames(model$design.Z) )]
      phi.chain <- phi.chain.nd(posterior, newdata.phi, link.phi)
      if(is.na(dim(phi.chain)[2])) phi.chain <- matrix(rep(phi.chain, n),ncol=n)
      pred.chain <- var.fun(model, mu.chain, phi.chain)$variance
    }

    if((model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta")) & cluster==T) {
      link.mu <- model$link.mu
      newdata.mu <- newdata[,which(colnames(newdata) %in% colnames(model$design.X) )]
      mu.chain <- mu.chain.nd(posterior, newdata.mu, link.mu)

      p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
      w.chain <- rstan::extract(posterior, pars="w", permuted=T)[[1]]
      parz.min <- pmin(apply(mu.chain, 2, function(x) x/p.chain) , apply(1-mu.chain, 2, function(x) x/(1-p.chain)))
      l1.chain <-  mu.chain + apply(parz.min,2, function(x) x*(1-p.chain)*w.chain)
      l2.chain <-  mu.chain - apply(parz.min,2, function(x) x*p.chain*w.chain)
    }
  }

  if(estimate == "quantile"){
    predicted <- apply(pred.chain, 2, quantile, probs=q)
    if ((model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta")) & cluster==T){
      predicted.l1 <- apply(l1.chain, 2, quantile, probs=q)
      predicted.l2 <- apply(l2.chain, 2, quantile, probs=q)}
  } else {
    predicted <- apply(pred.chain, 2, eval(str2expression(estimate)))
    if ((model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta")) & cluster==T){
      predicted.l1 <- apply(l1.chain, 2, eval(str2expression(estimate)))
      predicted.l2 <- apply(l2.chain, 2, eval(str2expression(estimate)))
    }
  }
  if ((model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta")) & cluster==T){
    predicted <- cbind(predicted, predicted.l1, predicted.l2)
    colnames(predicted) <- c(type, "lambda1", "lambda2")
  } else {
    predicted <- cbind(predicted)
    colnames(predicted) <- type
  }
  return(predicted)
}

#' @title Residuals Method for flexreg Objects
#'
#' @description Method that computes various types of residuals from objects of class \code{`flexreg`}. If the model type is \code{FB} or \code{FBB} and \code{cluster = T}, the method returns also residuals with respect to cluster means.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}}.
#' @param type a character indicating type of residuals (\code{raw} or \code{standardized}).
#' @param cluster logical. If the model is \code{"FB"} or \code{"FBB"}, \code{cluster=T} returns the cluster means. By default \code{cluster = F}.
#' @param estimate a character indicating the type of estimate: \code{mean} (default), \code{median}, or \code{quantile}.
#' @param q if estimate is \code{quantile}, a numeric value of probability in (0, 1).
#' @param ... additional arguments. Currently not used.
#'
#' @details Raw residuals are defined as \eqn{r_i=y_i-\hat{\mu}_i} (or \eqn{r_i= y_i/n_i-\hat{\mu}_i} for binomial data)
#' for \eqn{i=1, \dots, n}. The values \eqn{y_i} for \eqn{i,\dots,n} are referred to the observed
#' response variable and they are specified on the left-hand side of \code{formula} in the
#' \code{\link{flexreg}} function.
#' \eqn{\hat{\mu}_i} for \eqn{i=1, \dots, n} is the predicted value. It can be computed separately
#'  through the \code{\link{predict}} function by setting \code{type=response}.
#' Standardized residuals are defined as \eqn{\frac{r_i}{\widehat{Var}(y_i)}} where
#'  \eqn{\widehat{Var}(y_i)}
#' is the variance of the dependent variable evaluated at the posterior means
#' (default, otherwise quantile of order q) of the parameters.
#' If the model is \code{"FB"} or \code{"FBB"} and \code{cluster=T}, the cluster residuals are computed as
#' the difference between the observed response/relative response  and the cluster means
#' \eqn{\hat{\lambda}_{1i}} and \eqn{\hat{\lambda}_{2i}} for \eqn{i=1, \dots, n}.
#'
#' @examples{
#' data("Reading")
#' FB <- flexreg(accuracy ~ iq, Reading, type="FB", n.iter=1000)
#' residuals(FB, type="raw", cluster=TRUE)
#' }
#'
#' @references {
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018) A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' \cr
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' }
#'
#' @import stats
#'
#' @export
#'
#'

residuals.flexreg <- function(object, type = "raw",
                              cluster=FALSE, estimate="mean", q=NULL, ...){
  model <- object
  model.name <- model$model@model_name
  n <- model$n

  if(model.name %in% c("Beta", "Beta_phi", "VIB", "VIB_phi", "Bin", "BetaBin", "BetaBin_theta") & cluster==T)
    stop("Cluster residuals are available only for FB models")
  if(!(estimate %in% c("mean", "median", "quantile")))
    stop("Argument `estimate` must be set equal to `mean`, `median` or `quantile`")
  if(estimate=="quantile" & is.null(q))
    stop("You have to specify the order of the quantile")
  if(!(type %in% c("raw", "standardized")))
    stop("Argument `type` must be set equal to `raw` or `standardized`")

  if(model.name %in% c("Bin", "BetaBin", "BetaBin_theta", "FBB", "FBB_theta")){
    response <- model$response/n
  } else {
    response <- model$response
  }

  predicted <- predict.flexreg(model, type = "response", newdata=NULL, cluster=cluster, estimate=estimate, q=q)

  if(type == "raw") {
    if(model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta") & cluster==T){
      residuals <- t(apply(cbind(y=response, predicted), 1, function(x) c(x[1]-x[2],x[1]-x[3], x[1]-x[4])))
      colnames(residuals) <- c("response", "cluster1", "cluster2")
      label <- apply(residuals[,-1], 1, function(x) ifelse(abs(x[1])<=abs(x[2]), 1, 2))
      residuals <- cbind(residuals, label)} else {
        residuals <- response-predicted
      }
  } else {

    posterior <- model$model
    mu.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
    phi.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]

    if(model.name %in% c("FB", "FB_phi", "FBB", "FBB_theta") & cluster==T){
      residuals <- t(apply(cbind(y=response, predicted), 1, function(x) c(x[1]-x[2],x[1]-x[3], x[1]-x[4])))
      colnames(residuals) <- c("response", "cluster1", "cluster2")
      variance <- matrix(unlist(lapply((var.fun(model, mu.chain, phi.chain)), colMeans)), ncol=3)
      residuals  <- residuals/sqrt(variance)
      label <- apply(residuals[,-1], 1, function(x) ifelse(abs(x[1])<=abs(x[2]), 1, 2))
      residuals <- cbind(residuals, label)} else {
        variance <- colMeans(var.fun(model, mu.chain, phi.chain)$variance)
        residuals <- (response-predicted)/sqrt(variance)
      }
  }
  return(residuals)
}



