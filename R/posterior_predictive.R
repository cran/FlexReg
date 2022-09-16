#' Posterior Predictive
#'
#' The function takes an object of class \code{`flexreg`} and generates values from the posterior predictive distribution.
#'
#' @param model an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}}.
#' @param newdata an optional data frame containing variables with which to predict. If omitted, the fitted values are used.
#'
#' @details The function generates values from the posterior predictive distribution, which is the distribution of a  future outcome given the observed data.
#' The posterior predictive distribution is computed for \eqn{y} in case of bounded responses and
#' for \code{y/n} in case of binomial responses.
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
#' @export
#'
#'

posterior_predict.flexreg <- function(model, newdata=NULL)
{
  posterior <- model$model
  model.name <- posterior@model_name
  model.type <- model$type
  N <- length(model$response)
  nsim <- dim(posterior)[1]*dim(posterior)[2]
  n <- model$n

  link.mu <- model$link.mu
  link.phi <- model$link.phi
  link.theta <- model$link.theta

  covariate.names.mu <- colnames(model$design.X)
  covariate.names.phi <- colnames(model$design.Z)
  covariate.names.q0 <- colnames(model$design.X0)
  covariate.names.q1 <- colnames(model$design.X1)
  covariate.names.theta <- colnames(model$design.Z)

  ##############################
  # Estraggo e salvo le catene:
  ##############################


  if(model.type %in% c("Beta", "FB", "VIB")){
    # Se non ho dei nuovi dati da prevedere, predico i dati osservati. Le catene sono giÃ  salvate nel modello
    # quindi estraggo la catena ed applico la funzione di stima scelta (di default la media a posteriori)
    if(is.null(newdata)){
      #Ncol <- N
      mu.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
      phi.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]
      if(!is.null(model$call$zero.formula)){
        q0.chain <- rstan::extract(posterior, pars="q0", permuted=T)[[1]]
      } else{
        q0.chain <- NULL
      }
      if(!is.null(model$call$one.formula)){
        q1.chain <- rstan::extract(posterior, pars="q1", permuted=T)[[1]]
      } else{
        q1.chain <- NULL#matrix(rep(NULL,N), ncol=N)
      }
    } else { # --> Se ho dati da prevedere
      if(!all(colnames(newdata) %in% unique(c(covariate.names.mu, covariate.names.phi,
                                              covariate.names.q0, covariate.names.q1))))#check if correct
        stop("`newdata` must contain the same predictors as the one ones in formula")
      N <- nrow(newdata)
      if((c("(Intercept)") %in% colnames(newdata)) ==F & # If newdata does not have the intercept..
         ((c("(Intercept)") %in% colnames(model$design.X)) ==T) |
         (c("(Intercept)") %in% colnames(model$design.Z) ==T)|
         (c("(Intercept)") %in% colnames(model$design.X0) ==T)|
         (c("(Intercept)") %in% colnames(model$design.X1) ==T)) # .. and the intercept is required..
        newdata$`(Intercept)` <- rep(1, nrow(newdata)) # .. add the intercept

      newdata.mu  <- newdata[,match(colnames(model$design.X),colnames(newdata))]
      newdata.phi <- newdata[,match(colnames(model$design.Z),colnames(newdata))]
      newdata.q0 <- newdata[,match(colnames(model$design.X0),colnames(newdata))]
      newdata.q1 <- newdata[,match(colnames(model$design.X1),colnames(newdata))]
      mu.chain <-   mu.chain.nd(posterior, newdata.mu , link.mu)
      phi.chain <-   phi.chain.nd(posterior, newdata.phi , link.phi)
      q0.chain <- q.chain.nd(model, newdata.q0, newdata.q1)[[1]]
      q1.chain <- q.chain.nd(model, newdata.q0, newdata.q1)[[2]]
    }

    #additional parameters
    if(model.type == "FB"){
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

      if(grepl("phi", model.name)){
        phi.post <- phi.chain[,l]
      } else {
        phi.post <- phi.chain
      }
      if(model.type=="Beta") {
        param <- as.data.frame(cbind(mu.post, phi.post, q0.post, q1.post))
        post.pred[,l] <- unlist(lapply(1:nsim, function(z) rBeta_mu(n=1, mu=param$mu.post[z], phi=param$phi.post[z],
                                                                    q0=param$q0.post[z], q1=param$q1.post[z])))
      } else if(model.type=="FB") {
        param <- as.data.frame(cbind(mu.post, phi.post, p.chain, w.chain, q0.post, q1.post))
        post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rFB(n=1, mu=param$mu.post[z], phi=param$phi.post[z],
                                                                p=param$p.chain[z], w=param$w.chain[z],
                                                                q0=param$q0.post[z], q1=param$q1.post[z])))
      } else if(model.type=="VIB") {
        param <- as.data.frame(cbind(mu.post, phi.post, p.chain, k.chain, q0.post, q1.post))
        post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rVIB(n=1, mu=param$mu.post[z], phi=param$phi.post[z],
                                                                 p=param$p.chain[z], k=param$k.chain[z],
                                                                 q0=param$q0.post[z], q1=param$q1.post[z])))
      } #else stop("Model not defined")#penso si possa togliere questo stop perchè inutile

    }
  } else {
    # ELSE, if the model is for binomial data:

    # Se non ho dei nuovi dati da prevedere, predico i dati osservati. Le catene sono gia' salvate nel modello
    # quindi estraggo la catena ed applico la funzione di stima scelta (di default la media a posteriori)
    if(is.null(newdata)){
      #Ncol <- N
      mu.chain <- rstan::extract(posterior, pars="mu", permuted=T)[[1]]
      if(model.name != "Bin") theta.chain <- rstan::extract(posterior, pars="theta", permuted=T)[[1]]
    } else {
      # --> Se ho dati da prevedere
      if(!all(unique(c(covariate.names.mu, covariate.names.phi)) %in% colnames(newdata)))
        stop("`newdata` must contain the same predictors as the one ones in formula")
      N <- nrow(newdata)
      if((c("(Intercept)") %in% colnames(newdata)) ==F & # If newdata does not have the intercept..
         ((c("(Intercept)") %in% colnames(model$design.X)) ==T) |
         (c("(Intercept)") %in% colnames(model$design.Z) ==T)) # .. and the intercept is required..
        newdata$`(Intercept)` <- rep(1, nrow(newdata)) # .. add the intercept

      newdata.mu  <- newdata[,match(colnames(model$design.X),colnames(newdata))]
      newdata.theta <- newdata[,match(colnames(model$design.Z),colnames(newdata))]
      mu.chain <-   mu.chain.nd(posterior, newdata.mu , link.mu)
      theta.chain <-   theta.chain.nd(posterior, newdata.phi , link.phi)
    }

    #additional parameters
    if(model.name %in% c("FBB", "FBB_theta")){
      p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
      w.chain <- rstan::extract(posterior, pars="w", permuted=T)[[1]]
    }

    post.pred <- matrix(NA, ncol=N, nrow=nsim)
    for(l in 1:N){
      mu.post <- mu.chain[,l]
      if(model.name %in% c("BetaBin_theta", "FBB_theta")){
        theta.post <- theta.chain[,l]
      } else if(model.name != "Bin") {
        theta.post <- theta.chain
      }
      if(model.name %in% c("Bin")) {
        param <- cbind(mu.post)
        type <- "Bin"
      } else if(model.name %in% c("FBB", "FBB_theta")) {
        param <- cbind(mu.post, theta.post,p.chain, w.chain)
        type <- "FBB"
      } else if(model.name %in% c("BetaBin", "BetaBin_theta")) {
        param <- cbind(mu.post, theta.post)
        type <- "BetaBin"
      } else stop("Model not defined")

      if(type == "Bin"){
        post.pred[,l] <- unlist(lapply(1:nsim, function(z) rbinom(n=1, size=n[l], prob=param[z,1])))
      } else if(type == "FBB"){
        post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rFBB(n=1, size=n[l], mu=param[z,1], theta=param[z,2], p=param[z,3], w=param[z,4])))
      } else if(type=="BetaBin") {
        post.pred[,l] <-  unlist(lapply(1:nsim, function(z) rBetaBin(n=1, size=n[l], mu=param[z,1], theta=param[z,2])))
      }

      # Propriotion:
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
# sistemare qui...gli argomenti della funzione devono matchare con quelli della precedente
posterior_predict <- function(model, newdata=NULL)
{
  UseMethod("posterior_predict")
}


#' Posterior Predictives Plot
#'
#' Method for plotting the simulated posterior predictive distribution from an object of class \code{`flexreg_postpred`}. The plot shows the posterior predictive interval for each statistical unit.
#' Additionally, the mean of the posterior predictives and the values of the observed response (either \eqn{y} or \eqn{y/n} for bounded or binomial responses, respectively) can be added.
#'
#'
#' @param x an object of class \code{`flexreg_postpred`} containing the simulated posterior predictives, usually the result of \code{\link{posterior_predict}}.
#' @param prob the interval probability for the posterior predictives (default is 0.9).
#' @param p_mean a logical value indicating whether the posterior predictives' mean should be plotted.
#' @param response a numerical vector containing the response (either \eqn{y} or \eqn{y/n} for bounded or binomial responses, respectively) to be added to the plot. If \code{NULL}, observed values are not plotted.
#' @param ... additional arguments. Currently not used.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy ~ iq, data = Reading)
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
  # limits computation (quantile):
  pred.int <- t(apply(post.pred, 2, function(x) quantile(x, p=c((1-prob)/2, prob+(1-prob)/2))))

  #Definition of the data.frame involved by ggplot:
  dd <- data.frame(pp_inf=pred.int[,1],
                   pp_sup=pred.int[,2],
                   id=c(1:N))

  ll <-
    ggplot(dd) + theme_minimal() +
    scale_y_continuous(name="Posterior Predictive", breaks=seq(0,1,by=.25), limits=c(0,1.1))+
    geom_point(aes(x=id, y=pp_inf), shape="-",size=5)+
    geom_point(aes(x=id, y=pp_sup), shape="-",size=5)+
    labs(x="Unit ID")
  # Add segments:
  for(e in dd$id) {
    ll <- ll + geom_segment(x=e, y=dd$pp_inf[e], xend=e, yend=dd$pp_sup[e], linetype="dashed",size=.1)
  }

  # If p_mean = T, add the posterior predictive mean
  if(p_mean){
    pred.mean <- as.numeric(t(apply(post.pred, 2, function(x) mean(x))))
    ll <- ll + geom_point(data=data.frame(x=c(1:N),y=pred.mean),aes(x=x,y=y), color="#0072B2", size=1.5)
  }
  # Add the response:
  if(!is.null(response)) {
    ll <- ll + geom_point(data=data.frame(x=c(1:N),y=response),aes(x=x,y=y), color="#D55E00", size=1.5)
  }

  return(ll)
}
# plot.flexreg_postpred <- function(object, prob=0.9, p_mean=F, response=NULL, ...){
#   post.pred <- object
#   N <- ncol(post.pred)
#   pred.int <- t(apply(post.pred, 2, function(x) quantile(x, p=c((1-prob)/2, prob+(1-prob)/2))))
#
#
#   plot(pred.int[,1] ~ c(1:N), pch=-95, ylim=c(0,1.1), ylab="Posterior Predictive", xlab="Unit ID", main="Posterior Predictive")
#   points(c(1:N),pred.int[,2], pch=-95)
#   if(p_mean){
#     pred.mean <- t(apply(post.pred, 2, function(x) mean(x)))
#     points(c(1:N),pred.mean, pch=20, col=4)
#   }
#   if(!is.null(response)) points(c(1:N), response, pch=20, col=2)
#   for(i in 1:N) lines(c(i,i), c(pred.int[i,1],pred.int[i,2]), lty=3)
# }
