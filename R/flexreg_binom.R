#' Flexible Regression Models for Bounded Discrete Responses
#'
#' @description The function fits some flexible regression models for bounded discrete responses via a Bayesian approach to inference based on Hamiltonian Monte Carlo algorithm.
#' Available regression models are the flexible beta-binomial (\code{type = "FBB"}, default), the beta-binomial (\code{type = "BetaBin"}), and the binomial one (\code{type = "Bin"}).
#'
#'
#' @param formula an object of class "\code{\link{formula}}": a symbolic description of the model to be fitted (\code{y ~ x} or \code{y ~ x | z}, see Details).
#' @param data an optional  \code{data.frame}, list, or object that is coercible to a  \code{data.frame} through \code{\link{as.data.frame}} containing the variables in the model. If not found in \code{data}, the variables in \code{formula} are taken from the environment from which the function \code{\link{flexreg_binom}} is called.
#' @param type a character specifying the type of regression model. Current options are  \code{"FBB"} (flexible beta-binomial, default), \code{"BetaBin"} (beta-binomial), and  \code{"Bin"} (binomial).
#' @param n a character specifying the name of the variable containing the total number of trials.
#' @param link.mu a character specifying the link function for the mean model. Currently, \code{"logit"} (default), \code{"probit"}, \code{"cloglog"}, and \code{"loglog"} are supported.
#' @param prior.beta a character specifying the prior distribution for the  regression coefficients of the mean model, \code{beta}. Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.beta a positive numeric (vector of length 1) specifying the hyperprior scale parameter for the prior distribution of \code{beta} regression coefficients. The default is 100 if the prior is \code{"normal"}, 2.5 if it is \code{"cauchy"}.
#' @param link.theta a character specifying the link function for the overdispersion model. Currently, \code{"identity"} (default), \code{"logit"}, \code{"probit"}, \code{"cloglog"}, and \code{"loglog"} are supported. If \code{link.theta = "identity"}, the prior distribution for \code{theta} is a beta.
#' @param hyper.theta.a a numeric (vector of length 1) specifying the first shape parameter for the beta prior distribution of \code{theta}.
#' @param hyper.theta.b a numeric (vector of length 1) specifying the second shape parameter for the beta prior distribution of \code{theta}.
#' @param prior.psi a character specifying the prior distribution for the regression coefficients of the overdispersion model,\code{psi}. Not supported if \code{link.theta="identity"}. Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.psi a positive numeric (vector of length 1) specifying the hyperprior scale parameter for the prior distribution of \code{psi} regression coefficients. The default is 100 if the prior is \code{"normal"}, 2.5 if it is \code{"cauchy"}.
#' @param n.chain a positive integer specifying the number of Markov chains. The default is 1.
#' @param n.iter 	a positive integer specifying the number of iterations for each chain (including warm-up). The default is 5000.
#' @param warmup.perc the percentage of iterations per chain to discard.
#' @param thin a positive integer specifying the period for saving samples. The default is 1.
#' @param verbose a logical (with default \code{TRUE}) indicating whether to print intermediate output.
#' @param ... additional arguments from \code{\link[rstan]{sampling}}.
#'
#' @return The \code{\link{flexreg_binom}} function returns an object of class \code{`flexreg`}, i.e. a list with the following elements:
#' \item{\code{call}}{the function call.}
#' \item{\code{type}}{the type of regression model.}
#' \item{\code{formula}}{the original formula.}
#' \item{\code{link.mu}}{a character specifing the link function in the mean model.}
#' \item{\code{link.theta}}{a character specifing the link function in the overdispersion model.}
#' \item{\code{model}}{an object of class \code{`stanfit`} containing the fitted model.}
#' \item{\code{response}}{the response variable, assuming values in (0, 1).}
#' \item{\code{design.X}}{the design matrix for the mean model.}
#' \item{\code{design.Z}}{the design matrix for the overdispersion model (if defined).}
#'
#' @details Let Y be a random variable whose distribution can be specified in the \code{type} argument and \eqn{\mu} be the mean of Y/n.
#' The \code{\link{flexreg_binom}} function links the parameter \eqn{\mu} to a linear predictor through a function  \eqn{g_1(\cdot)} specified in \code{link.mu}:
#' \deqn{g_1(\mu) = x^t \bold{\beta},} where \eqn{\bold{\beta}} is the vector of regression coefficients for the mean model.
#' The prior distribution and the related hyperparameter of \eqn{\bold{\beta}} can be specified in \code{prior.beta} and \code{hyperparam.beta}.
#' By default, \code{link.theta="identity"}, meaning that the overdispersion parameter \eqn{\theta} is assumed to be constant.
#' In that case, the prior distribution for \eqn{\theta} is a beta with shape hyperparameters \eqn{a} and \eqn{b} that can be specified in \code{hyper.theta.a} and \code{hyper.theta.b}.
#' If not specified, \eqn{a=b=1}, otherwise if only one hyperparameter is specified, the other  is set equal.
#' It is possible to extend the model by linking \eqn{\theta} to an additional (possibly overlapping) set of covariates through a proper link
#' function \eqn{g_2(\cdot)}  specified in the \code{link.theta} argument: \deqn{g_2(\theta) = z^t \bold{\psi},} where \eqn{\bold{\psi}} is the vector of regression coefficients for the overdispersion model.
#' The prior distribution and the related hyperparameter of \eqn{\bold{\psi}} can be specified in \code{prior.psi} and \code{hyperparam.psi}.

#' In \code{\link{flexreg_binom}}, the regression model for the mean and, where appropriate, for the overdispersion parameter can be specified in the
#' \code{formula} argument with a formula of type \code{y ~ x1 + x2 | z1 + z2} where covariates on the left of "|" are included in the regression model
#' for the mean, whereas covariates on the right of "|" are included in the regression model for the overdispersion.
#'
#' If the second part is omitted, i.e., \code{y ~ x1 + x2}, the overdispersion is assumed constant for each observation.
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @examples
#' \dontrun{
#' data(Bacteria)
#' fbb <- flexreg_binom(y ~ females,  n = "n", data = Bacteria, type = "FBB")
#' }
#' @import Rcpp methods rstan
#'
#' @export

flexreg_binom <- function(formula, data, type="FBB", n,
                     link.mu="logit", prior.beta = "normal", hyperparam.beta = 100,
                     hyper.theta.a=NULL, hyper.theta.b=NULL,
                     link.theta=NULL, prior.psi = NULL, hyperparam.psi = NULL,
                     n.chain=1, n.iter=5000, warmup.perc=.5, thin=1, verbose=TRUE, ...)
{
  cl <- match.call()
  n_arg <- n

  if(is.na(match(type, c("FBB", "Bin", "BetaBin")))) stop("Please specify the type of model correctly.")

  if (missing(data)) data <- environment(formula) else  data <- as.data.frame(data)

  ###############################################
  ###############################################

  formula <- Formula::as.Formula(formula)

  if(is.character(link.theta)) link_code_theta <-
                  pmatch(link.theta, c("identity", "logit", "probit", "cloglog", "loglog"))

  if(missing(n)) stop("Please specify the argument n.")

  if(length(formula)[2] >= 2){
    model.theta <- TRUE

    if(type=="Bin") stop("parameter theta is not defined in the binomial regression.")

    # Check on link functions for theta:
    if(is.null(link.theta)) link_code_theta <- 2
    if(link_code_theta <2) stop("invalid link function for regression model for theta")
    if(!is.null(hyper.theta.a) | !is.null(hyper.theta.b)) stop("hyperparameters chosen for theta but regression model for theta in formula. Please define priors for coefficients psi instead.")
    if(length(formula)[2] > 2)  warning("formula must not have more than two RHS parts")

    } else if(length(formula)[2] < 2){
      model.theta <- FALSE
      link_prior_psi <- NULL
      hyperparam.psi <- NULL
      prior.psi <- NULL

      if(is.null(link.theta)) link_code_theta <- 1

      if(link_code_theta >=2) stop("a link function has been selected for theta but no covariates appear on the RHS part. Please specify covariates for a regression model on theta in the RHS part.")

      # Check on hyperparameters:
      if(is.null(hyper.theta.a) & is.null(hyper.theta.b)) {
        hyper.theta.a <- 1
        hyper.theta.b <- 1
      } else if(is.null(hyper.theta.a) & !is.null(hyper.theta.b)){
        hyper.theta.a <- hyper.theta.b
      } else if(!is.null(hyper.theta.a) & is.null(hyper.theta.b)){
        hyper.theta.b <- hyper.theta.a
      }

      if(!is.null(prior.psi) | !is.null(hyperparam.psi)) stop("Error: prior for psi coefficients defined but no regression model for phi. Please define covariates for regression model for theta")
    }

  # Check on link functions for mu:
  if(is.character(link.mu)) {
    link_code_mu <- pmatch(link.mu, c("logit", "probit", "cloglog", "loglog"))
  } else {
    stop("Invalid link function for mu.")
  }

  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  N <- length(y)
  n <- get(n, as.environment(data))
  miss <- as.numeric(attr(attr(mf, "na.action"), "names"))
  if(length(miss)>0) n <-  n[-miss]


  if(N < 1)
    stop("Empty model")
  if(min(y) < 0)
    stop("Invalid dependent variable: y cannot assume negative values.")
  if(max(y/n) > 1)
    stop("Invalid dependent variable: y must be less or equal to n for every observation.")

  # Design matrix:
  if(n_arg %in% colnames(mf)){
    X <- model.matrix(formula, data = mf[,-which(colnames(mf)==n_arg)], rhs = 1)
  } else {
    X <- model.matrix(formula, data = mf, rhs = 1)
  }

  if((model.theta == F)){
    Z <- NULL
  } else {
    if(n_arg %in% colnames(mf)){
      Z  <- model.matrix(formula, data = mf[,-which(colnames(mf)==n_arg)], rhs = 2)
    } else {
      Z <- model.matrix(formula, data = mf, rhs = 2)
    }
  }

  #check for deleting elements in y, X, corresponding to NA in n
  if(anyNA(n)){
  miss.n <- !is.na(n)
    y <- y[miss.n]
    X <- X[miss.n,]
    Z <- Z[miss.n,]
  }

  if(is.null(prior.beta)) {
    link_prior_beta <- 1
    } else if(is.character(prior.beta)) {
        link_prior_beta  <- pmatch(prior.beta, c("normal", "cauchy"))
        if(is.na(link_prior_beta)) stop("Invalid prior for beta parameters.")
      } else stop("Invalid prior for beta parameters.")



  # Checks on hyperparameters:
  if(hyperparam.beta < 0) stop("Hyperprior for beta coefficients must be positive")
  if(link_prior_beta == 1 & hyperparam.beta < 10) warning("An hyperparameter of at least 10 for normal prior sd for beta coefficients is recommended")
  if(link_prior_beta == 2 & hyperparam.beta != 2.5) warning("An hyperparameter of 2.5 for cauchy prior for beta coefficients is recommended")

  # If no regression is specified for theta....
  if(model.theta == F){
    if(hyper.theta.a <= 0 | hyper.theta.b <= 0) stop("Hyperprior for theta must be positive")
  } else { # else, if a regression model is specified for theta...
    hyper.theta.a <- NULL
    hyper.theta.b <- NULL

    if(is.null(prior.psi)) {
      link_prior_psi <- 1
      } else if (is.character(prior.psi)) {
        link_prior_psi  <- pmatch(prior.psi, c("normal", "cauchy"))
        if(is.na(link_prior_psi)) stop("Invalid prior for psi parameters.")
      } else stop("Invalid prior for psi parameters.")


    if(is.null(hyperparam.psi)) {
      hyperprior_psi <- ifelse(link_prior_psi==1, 100, 2.5)
      } else if(hyperparam.psi < 0) stop("Hyperparameter for psi coefficients must be positive") else {
        if(link_prior_psi == 1 & hyperparam.psi < 10) warning("An hyperparameter of at least 10 for normal for psi coefficients is recommended")
        if(link_prior_psi == 2 & hyperparam.psi != 2.5) warning("An hyperparameter of 2.5 for cauchy prior for psi coefficients is recommended")
      }
  }

  model <- fit.model_binom(model.theta = model.theta, type = type, N = N, n=n,
                     y = y,  X = X, Z = Z,
                     link_code_mu = link_code_mu,
                     link_prior_beta = link_prior_beta,
                     hyperparam.beta = hyperparam.beta,

                     link_code_theta = link_code_theta,
                     hyper.theta.a = hyper.theta.a,
                     hyper.theta.b = hyper.theta.b,

                     link_prior_psi = link_prior_psi,
                     hyperparam.psi = hyperparam.psi,

                     n.iter = n.iter, warmup.perc = warmup.perc,
                     n.chain = n.chain, thin = thin, verbose=verbose, ...)

#questa riga serve per summary
link.theta <- c("identity", "logit", "probit", "cloglog", "loglog")[link_code_theta]
  output <- list(call=cl, type=type, formula=formula, link.mu=link.mu, link.theta=link.theta,
                 model=model, response=y, design.X=X, design.Z=Z, n=n)
  class(output) <- c("flexreg", "flexreg_binom")
  #invisible(output)
  return(output)
}


fit.model_binom <- function(model.theta = model.theta, type = type, N = N, n=n, y = y, X = X, Z = Z,
                      link_code_mu = link_code_mu,
                      link_prior_beta = link_prior_beta,
                      hyperparam.beta = hyperparam.beta,

                      link_code_theta = link_code_theta,
                      hyper.theta.a = hyper.theta.a,
                      hyper.theta.b = hyper.theta.b,

                      link_prior_psi = link_prior_psi,
                      hyperparam.psi = hyperparam.psi,
                      n.iter, warmup.perc, n.chain, thin, verbose, ...){

  data.stan <- list(
      N = N,
      n = n,
      y = y,
      X = X,
      Z = Z,
      K = ncol(X),
      H = ncol(Z),
      link_code_mu = link_code_mu,
      link_prior_beta = link_prior_beta,
      hyperprior_beta = hyperparam.beta,

      link_code_theta = link_code_theta,
      hyper_theta_a = hyper.theta.a,
      hyper_theta_b = hyper.theta.b,

      link_prior_psi = link_prior_psi,
      hyperprior_psi = hyperparam.psi
  )

  if(model.theta) stan.model <- stanmodels[[paste(type, "_theta", sep="")]]
  else stan.model <- stanmodels[[type]]

  fit.stan = rstan::sampling(
    object = stan.model,
    data = data.stan,
    chains = n.chain,
    thin = thin,
    #init = initl,
    iter = n.iter, warmup = round(n.iter*warmup.perc),
    refresh = verbose*n.iter/10, #show an update @ each %10
    # seed=1
    ...
  )

  output <- fit.stan
  return(output)

}





















