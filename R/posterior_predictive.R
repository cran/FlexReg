#' Posterior Predictive Method for \code{`flexreg`} objects
#'
#' @description The function takes an object of class \code{`flexreg`} and generates values from the posterior predictive distribution.
#'
#' @param model an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param newdata an optional  \code{data.frame} containing variables with which to predict. If omitted, the fitted values are used.
#' @param n.new an optional vector containing the total number of trials with which to predict. It must be specified if \code{newdata} is not \code{NULL} and the
#' \code{\link{flexreg}} object is the result of the \code{\link{flexreg_binom}} function (i.e., the fitted model is binomial, beta-binomial, or flexible beta-binomial). The vector must have the same length as \code{nrow(newdata)}.
#'
#' @details The function generates values from the posterior predictive distribution, which is the distribution of a  future outcome given the observed data.
#' The posterior predictive distribution is computed for \eqn{y} in case of bounded continuous  responses and
#' for \eqn{y/n} in case of bounded discrete responses.
#' @return An object of class \code{`flexreg_postpred`} containing a matrix with the simulated posterior predictions. Each column refers to a statistical unit to predict.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, n.iter=1000)
#' pp <- posterior_predict(FB)
#' plot(pp)
#' }
#'
#' @import stats Formula
#' @method posterior_predict flexreg
#' @importFrom rstantools posterior_predict
#' @export
#'
#' @references{
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020). Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309. doi:10.1177/1471082X18821213 \cr
#' \cr
#' Gelman, A., Carlin, J. B., Stern, H. S., Rubin, D. B. (2014). Bayesian Data Analysis, 3th edition. Chapman and Hall/CRC. doi:10.1201/b16018 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' }
#'
#'

posterior_predict.flexreg <- function(model, newdata = NULL, n.new = NULL)
{
  posterior <- model$model[[1]]
  model.type <- model$type
  nsim <- dim(posterior)[1]*dim(posterior)[2]

  if("flexreg_binom" %in% class(model)){
    if(is.null(newdata))  {
      n <- model$n
    } else {
      n <- n.new
      if(is.null(n)) stop("Please specify an additional variable in newdata containing the values of n.")
    }
  } else {
    n <- NULL
  }

  formula <- model$formula

  #adjust the newdata
  if(!is.null(newdata)) newdata <- newdata.adjust(newdata, formula)
  N <- ifelse(is.null(newdata), length(model$response), nrow(newdata))

  mu.chain <- predict_mu.chain(model, newdata)

  if(length(model$model) >1){
    q.chain <- predict_q.chain(model, newdata)
    q0.chain <- q.chain$q0.chain
    q1.chain <- q.chain$q1.chain
  } else{
    q0.chain <- NULL
    q1.chain <- NULL
  }

  #if bound model is phi, otherwise is a theta, but saved in same object
  if(("flexreg_bound" %in% class(model))){
    phi.theta.chain <- predict_precision(model, newdata)[[1]]
  } else  if( model.type != "Bin" & ("flexreg_binom" %in% class(model))){
    phi.theta.chain <- predict_over(model, newdata)[[1]]
  } else{#if model is binomial
    phi.theta.chain <- matrix(0)
  }

  #additional parameters
  if(model.type %in% c("FB", "FBB")){
    p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
    w.chain <- rstan::extract(posterior, pars="w", permuted=T)[[1]]
  }

  if(model.type == "VIB"){
    p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
    k.chain <- rstan::extract(posterior, pars="k", permuted=T)[[1]]
  }

  post.pred <- matrix(NA, ncol=N, nrow=nsim)
  for(l in 1:N){
    mu.post <- mu.chain[,l]
    q0.post <- q0.chain[,l]
    q1.post <- q1.chain[,l]

    if(ncol(phi.theta.chain) == 1){
      phi.theta.post <- as.vector(phi.theta.chain)
    } else {
      phi.theta.post <- phi.theta.chain[,l]
    }


    if(model.type == "Beta") {
      param <- as.data.frame(cbind(mu.post, phi.theta.post, q0.post, q1.post))
      post.pred[,l] <- unlist(lapply(1:nsim, function(z) rBeta(n=1, mu=param$mu.post[z], phi=param$phi.theta.post[z],
                                                               q0=param$q0.post[z], q1=param$q1.post[z])))
    } else if(model.type == "FB") {
      param <- as.data.frame(cbind(mu.post, phi.theta.post, p.chain, w.chain, q0.post, q1.post))
      post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rFB(n=1, mu=param$mu.post[z], phi=param$phi.theta.post[z],
                                                              p=param$p.chain[z], w=param$w.chain[z],
                                                              q0=param$q0.post[z], q1=param$q1.post[z])))
    } else if(model.type == "VIB") {
      param <- as.data.frame(cbind(mu.post, phi.theta.post, p.chain, k.chain, q0.post, q1.post))
      post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rVIB(n=1, mu=param$mu.post[z], phi=param$phi.theta.post[z],
                                                               p=param$p.chain[z], k=param$k.chain[z],
                                                               q0=param$q0.post[z], q1=param$q1.post[z])))
    } else if(model.type == "Bin"){
      param <- cbind(mu.post)
      post.pred[,l] <- unlist(lapply(1:nsim, function(z) rbinom(n=1, size=n[l], prob=param[z,1])))
      post.pred[,l] <- post.pred[,l]/n[l]
    } else if(model.type == "BetaBin"){
      param <- cbind(mu.post, phi.theta.post)
      post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rBetaBin(n=1, size=n[l], mu=param[z,1], theta=param[z,2])))
      post.pred[,l] <- post.pred[,l]/n[l]
    } else if(model.type == "FBB"){
      param <- cbind(mu.post, phi.theta.post, p.chain, w.chain)
      post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rFBB(n=1, size=n[l], mu=param[z,1], theta=param[z,2], p=param[z,3], w=param[z,4])))
      post.pred[,l] <- post.pred[,l]/n[l]
    }
  }
  class(post.pred) <- "flexreg_postpred"
  return(post.pred)
}


#' @title posterior_predict
#' @export
#' @keywords internal
#'
posterior_predict <- function(model, newdata = NULL, n.new = NULL)
{
  UseMethod("posterior_predict")
}

#' Plot Method for \code{`flexreg_postpred`} objects
#'
#' Method for an object of class \code{`flexreg_postpred`} containing the simulated posterior predictive distribution, usually the result of \code{\link{posterior_predict}} function. The plot shows the posterior predictive interval for each statistical unit.
#' Additionally, the mean of the posterior predictives and the values of the observed response (either \eqn{y} or \eqn{y/n} for bounded continuous  or discrete responses, respectively) can be added.
#'
#'
#' @param x an object of class \code{`flexreg_postpred`} containing the simulated posterior predictives, usually the result of \code{\link{posterior_predict}}.
#' @param prob the interval probability for the posterior predictives (default is 0.9).
#' @param p_mean a logical value indicating whether the posterior predictives' mean should be plotted.
#' @param response a numerical vector containing the response (either \eqn{y} or \eqn{y/n} for bounded continuous or discrete responses, respectively) to be added to the plot. If \code{NULL}, observed values are not plotted.
#' @param ... additional arguments. Currently not used.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, n.iter=1000)
#' pp <- posterior_predict(FB)
#' plot(pp)
#' }
#'
#' @import stats graphics
#'
#' @method plot flexreg_postpred
#'
#' @export
#'
#'
#'
plot.flexreg_postpred <- function(x, prob=0.9, p_mean=F, response=NULL, ...){
  id <- pp_inf <- pp_sup <- y <- NULL
  post.pred <- x
  N <- ncol(post.pred)
  type = c()
  # limits computation (quantile):
  pred.int <- t(apply(post.pred, 2, function(x) quantile(x, p=c((1-prob)/2, prob+(1-prob)/2))))

  #Definition of the data.frame involved by ggplot:
  dd <- data.frame(pp_inf=pred.int[,1],
                   pp_sup=pred.int[,2],
                   id=factor(1:N))

  ll <-
    ggplot(dd) + theme_minimal() +
    scale_y_continuous(name="Posterior Predictive", breaks=seq(0,1,by=.25), limits=c(0,1))+
    geom_point(aes(x=id, y=pp_inf), shape="-", size=5)+
    geom_point(aes(x=id, y=pp_sup), shape="-", size=5)+
    labs(x="Unit ID")
  # Add segments:
  for(e in 1:N) {
    #e <- as.numeric(e)
    ll <- ll +
      geom_segment(x=e, y=dd$pp_inf[e], xend=e, yend=dd$pp_sup[e], linetype="dashed")
  }

  # If p_mean = T, add the posterior predictive mean
  dd.points <- data.frame(x=NULL, y=NULL, type=NULL)

  if(p_mean){
    pred.mean <- as.numeric(t(apply(post.pred, 2, function(x) mean(x))))
    dd.points <- rbind(dd.points, data.frame(x=c(1:N),y=pred.mean, type = "Mean"))
    #ll <- ll + geom_point(data=data.frame(x=c(1:N),y=pred.mean),aes(x=x,y=y), color="#0072B2", size=1.5)
  }
  # Add the response:
  if(!is.null(response)) {
    dd.points <- rbind(dd.points,data.frame(x=c(1:N),y=response, type = "Response"))
    #ll <- ll + geom_point(data=data.frame(x=c(1:N),y=response),aes(x=x,y=y), color="#D55E00", size=1.5)
  }

  if(nrow(dd.points) > 0){
    ll <- ll +
      geom_point(data=dd.points,aes(x=x, y=y, color=type), size=1.5) +
      theme(legend.position="bottom")+
      theme(axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.title = element_blank(),
            legend.text = element_text(size = 15))
  }

  return(ll)
}


#' @title Summary Method for \code{`flexreg_postpred`} objects
#'
#' @description Summary method for an object of class \code{`flexreg_postpred`}, containing the simulated posterior predictive distribution.
#'
#'
#' @param object an object of class \code{`flexreg_postpred`} containing the simulated posterior predictives, usually the result of \code{\link{posterior_predict}}.
#' @param ... additional arguments.
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, n.iter=1000)
#' pp <- posterior_predict(FB)
#' summary(pp)
#' }
#'
#' @import stats graphics
#'
#' @method summary flexreg_postpred
#'
#' @export
#'
#' @return The function \code{\link{summary.flexreg_postpred}} returns an array with the statistical units by row. The number of rows of the array is equal to the number of columns of the object of class \code{`flexreg_postpred`} that is given to the function.
#' By column there are some synthesis values that are the minimum, the first quartile, the media, the mean, the third quartile, and the maximum.
#'
#'
summary.flexreg_postpred <- function(object, ...){
  n <- ncol(object)
  x <- t(object)
  summa <- cbind(apply(x, 1, min),
                 apply(x, 1, quantile, 0.25),
                 apply(x,1,median),
                 rowMeans(x),
                 apply(x,1, quantile, 0.75),
                 apply(x,1, max))
  colnames(summa) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")
  rownames(summa) <- 1:n

  return(summa)
}
