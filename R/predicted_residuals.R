#' @title Prediction Method for flexreg Objects
#'
#' @description Method that computes various types of prediction from objects of class \code{`flexreg`}. If the model type is \code{"FB"} without augmentation or \code{"FBB"} and \code{cluster = T}, the function returns also cluster means.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}}.
#' @param newdata an optional data frame containing variables with which to predict. If omitted, the fitted values are used.
#' @param cluster logical. If the model is \code{"FB"} without augmentation or \code{"FBB"}, \code{cluster = T} returns the cluster means. By default, \code{cluster = F}.
#' @param type a character indicating the type of predictions. Available options are: \code{"response"}, that returns the marginal fitted means of response/relative response;
#' \code{"link"}, the linear predictor of the mean model;
#' \code{"precision"}, the fitted precision parameter \eqn{phi};
#' \code{"overdispersion"}, the fitted overdispersion parameter \eqn{theta};
#'  \code{"variance"}, the fitted variance of the response.
#' @param estimate the type of estimate: \code{"mean"} (default), \code{"median"} or \code{"quantile"}.
#' @param q if \code{estimate = "quantile"}, numeric value of probability in (0, 1).
#' @param ... additional arguments. Currently not used.
#'
#' @details If \code{type="response"} the function returns the marginal mean that is \eqn{\mu} in case of no augmentation and
#' \eqn{q_1+(1-q_0-q_1)\mu} in case of augmentation. If \code{type="variance"} the function returns \eqn{Var(Y|0<Y<1)} in case of no augmentation and
#' \eqn{(1-q_0-q_1)Var(Y|0<Y<1)+q_1^2+(1-q_0-q_1)\mu^2-(q_1+(1-q_0-q_1)\mu)^2} in case of augmentation. See Di Brisco and Migliorati (2020) for details.
#' The option \code{type = "overdispersion"}  is available only for beta-binomial and flexible beta-binomial models.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data=Reading, type="FB")
#' predict(FB, type="response", cluster=TRUE)
#' }
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
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
  posterior <- model$model
  model.name <- posterior@model_name
  model.type <- model$type

  n <- length(model$response)

  if(!(model.name %in% c("FBNo", "FBNo_phi", "FBB", "FBB_theta")) & cluster==T)
    stop("Cluster means are available only for FBB and FB models without augmentation")
  if(!(estimate %in% c("mean", "median", "quantile")))
    stop("Argument `estimate` must be set equal to `mean`, `median` or `quantile`")
  if(estimate=="quantile" & is.null(q))
    stop("You have to specify the order of the quantile")
  if(!(type %in% c("response", "link", "precision", "variance", "overdispersion")))
    stop("Argument `type` must be set equal to `response`, `link`, `precision`, `variance`, or `overdispersion`")

  if(type=="variance" & (model.name %in% c("BetaBin", "BetaBin_theta","Bin", "FBB", "FBB_theta" )))
    stop("predicted variances are available only for Beta, FB, and VIB models")
  if(type=="overdispersion" & !(model.name %in% c("BetaBin", "BetaBin_theta", "FBB", "FBB_theta")))
    stop("predicted ovedispersion is available only for BetaBin and FBB models")
  if(type=="precision" & (model.name %in% c("BetaBin", "BetaBin_theta","Bin", "FBB", "FBB_theta" )))
    stop("predicted precision is available only for Beta, FB, and VIB models")

  # If we do not have new data:
  if(is.null(newdata)){
    if(is.null(model$call$zero.formula)){q0.chain <- 0
    }else {
      q0.chain <- rstan::extract(posterior, pars="q0", permuted=T)[[1]]
    }
    if(is.null(model$call$one.formula)){q1.chain <- 0
    }else {
      q1.chain <- rstan::extract(posterior, pars="q1", permuted=T)[[1]]
    }
    if( type == "response"){
      mu.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
      pred.chain <- q1.chain+(1-q0.chain-q1.chain)*mu.chain
    }
    if( type == "link") {
      beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
      X <- model$design.X
      pred.chain <- beta.chain %*% t(X)
    }

    if( type == "precision") {#non serve ulteriore controllo su type, perchè ci sono gli stop precedenti
      pred.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }
    if( type == "variance") {
      mu.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
      phi.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]
      pred.chain <- var.fun(model, mu.chain, phi.chain, q0.chain, q1.chain)$variance
    }
    if( type == "overdispersion" & model.name %in% c("BetaBin", "BetaBin_theta", "FBB", "FBB_theta")) {
      pred.chain <- rstan::extract(posterior, pars="theta", permuted=T)[[1]]
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }

    # non serve controllare che il modello sia FB perchè ho già imposto il blocco prima
    if (cluster==T){#(substr(model.name,1,2)=="FB") &
      l1.chain <- rstan::extract(posterior, pars="lambda1", permuted=T)[[1]]
      l2.chain <- rstan::extract(posterior, pars="lambda2", permuted=T)[[1]]
    }
  } else { # cioè se newdata non è NULL --> Se ho dati da prevedere
    newdata <- as.data.frame(newdata)
    if((c("(Intercept)") %in% colnames(newdata)) ==F & # If newdata does not have the intercept..
       ((c("(Intercept)") %in% colnames(model$design.X)) ==T) |
       (c("(Intercept)") %in% colnames(model$design.Z) ==T)|
       (c("(Intercept)") %in% colnames(model$design.X0) ==T)|
       (c("(Intercept)") %in% colnames(model$design.X1) ==T)) # .. and the intercept is required..
      newdata$`(Intercept)` <- rep(1, nrow(newdata)) # .. add the intercept

    if(!all(c(colnames(model$design.X),colnames(model$design.Z),
              colnames(model$design.X0),colnames(model$design.X1))  %in% colnames(newdata))){
      stop("newdata must be a data frame containing the same predictors as the ones in formula")
    }
    n <- nrow(newdata)

    newdata.q0 <- newdata[,match(colnames(model$design.X0),colnames(newdata))]
    newdata.q1 <- newdata[,match(colnames(model$design.X1),colnames(newdata))]
    q0.chain <- q.chain.nd(model, newdata.q0, newdata.q1)[[1]]
    if(is.null(q0.chain)) q0.chain <- 0
    q1.chain <- q.chain.nd(model, newdata.q0, newdata.q1)[[2]]
    if(is.null(q1.chain)) q1.chain <- 0

    if(type == "response") {
      link.mu <- model$link.mu
      newdata.mu <- newdata[,match(colnames(model$design.X),colnames(newdata))]
      mu.chain <- mu.chain.nd(posterior, newdata.mu, link.mu)
      pred.chain <- q1.chain+(1-q0.chain-q1.chain)*mu.chain
    }
    if(type == "link") {
      beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
      newdata.mu <- newdata[,match(colnames(model$design.X),colnames(newdata))]
      pred.chain  <- beta.chain %*% t(newdata.mu)
    }
    if(type == "precision" & model.type %in% c("Beta", "FB", "VIB")){
      link.phi <- model$link.phi
      newdata.phi <- newdata[,match(colnames(model$design.Z),colnames(newdata))]
      pred.chain <- phi.chain.nd(posterior, newdata.phi, link.phi)
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }
    if(type == "overdispersion" & model.name %in% c("BetaBin", "BetaBin_theta", "FBB", "FBB_theta")) {
      link.theta <- model$link.theta
      newdata.theta <- newdata[,match(colnames(model$design.Z),colnames(newdata))]
      pred.chain <- theta.chain.nd(posterior, newdata.theta, link.theta)
      if(is.na(dim(pred.chain)[2])) pred.chain <- matrix(rep(pred.chain, n),ncol=n)
    }
    if(type == "variance"& model.type %in% c("Beta", "FB", "VIB")){
      link.mu <- model$link.mu
      newdata.mu <- newdata[,match(colnames(model$design.X),colnames(newdata))]
      mu.chain <- mu.chain.nd(posterior, newdata.mu, link.mu)
      link.phi <- model$link.phi
      newdata.phi <- newdata[,match(colnames(model$design.Z),colnames(newdata))]
      phi.chain <- phi.chain.nd(posterior, newdata.phi, link.phi)
      if(is.na(dim(phi.chain)[2])) phi.chain <- matrix(rep(phi.chain, n),ncol=n)
      pred.chain <- var.fun(model, mu.chain, phi.chain, q0.chain,q1.chain)$variance
    }

    if(cluster==T) {
      link.mu <- model$link.mu
      newdata.mu <- newdata[,match(colnames(model$design.X),colnames(newdata))]
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
    if (cluster==T){
      predicted.l1 <- apply(l1.chain, 2, quantile, probs=q)
      predicted.l2 <- apply(l2.chain, 2, quantile, probs=q)
      }
  } else {
    predicted <- apply(pred.chain, 2, eval(str2expression(estimate)))
    if (cluster==T){
      predicted.l1 <- apply(l1.chain, 2, eval(str2expression(estimate)))
      predicted.l2 <- apply(l2.chain, 2, eval(str2expression(estimate)))
    }
  }

  if (cluster==T){
    predicted <- cbind(predicted, predicted.l1, predicted.l2)
    colnames(predicted) <- c(type,"lambda1", "lambda2")
  } else {
    predicted <- cbind(predicted)
    colnames(predicted) <- type
  }
  return(predicted)
}

#' @title Residuals Method for flexreg Objects
#'
#' @description Method that computes various types of residuals from objects of class \code{`flexreg`}. If the model type is  \code{FB} without augmentation or \code{FBB} and \code{cluster = T}, the method returns also residuals with respect to cluster means.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}}.
#' @param type a character indicating type of residuals (\code{"raw"} or \code{"standardized"}).
#' @param cluster logical. If the model is \code{"FB"} without augmentation or \code{"FBB"}, \code{cluster = T} returns the cluster means. By default \code{cluster = F}.
#' @param estimate a character indicating the type of estimate: \code{"mean"} (default), \code{"median"}, or \code{"quantile"}.
#' @param q if \code{estimate = "quantile"}, a numeric value of probability in (0, 1).
#' @param ... additional arguments. Currently not used.
#'
#' @details Raw residuals are defined as \eqn{r_i=y_i-\hat{\mu}_i} (or \eqn{r_i= y_i/n_i-\hat{\mu}_i} for binomial data).
#' The values \eqn{y_i} or \eqn{y_i/n_i} are  the observed
#' responses which are specified on the left-hand side of \code{formula} in the
#' \code{\link{flexreg}} or \code{\link{flexreg_binom}} function, respectively.
#'  \eqn{\hat{\mu}_i}  is the predicted value, the result of
#'   the \code{\link{predict}} function with \code{type = "response"}.
#' Standardized residuals are defined as \eqn{\frac{r_i}{\widehat{Var}(y_i)}} where
#'  \eqn{\widehat{Var}(y_i)}
#' is the variance of the dependent variable evaluated at the posterior means
#' (default, otherwise quantile of order q) of the parameters.
#' If the model is \code{"FB"} without augmentation or \code{"FBB"} and \code{cluster = T}, the cluster residuals are computed as
#' the difference between the observed response/relative response  and the cluster means
#' \eqn{\hat{\lambda}_{1i}} and \eqn{\hat{\lambda}_{2i}}.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data=Reading, type="FB")
#' residuals(FB, type="raw", cluster=TRUE)
#' }
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
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
  posterior <- model$model
  model.name <- posterior@model_name
  n <- model$n

  if(model.name %in% c("BetaNo", "BetaNo_phi", "Beta0", "Beta0_phi",
                       "Beta1", "Beta1_phi", "Beta01", "Beta01_phi",
                       "VIBNo", "VIBNo_phi", "VIB0", "VIB0_phi",
                       "VIB1", "VIB1_phi", "VIB01", "VIB01_phi",
                       "Bin", "BetaBin", "BetaBin_theta") & cluster==T)
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



