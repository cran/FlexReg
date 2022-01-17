#' Flexible Regression Models for Binomial data
#'
#' @description The function fits some flexible regression models for binomial data via a Bayesian approach to inference based on Hamiltonian Monte Carlo algorithm.
#' Available regression models are the flexible beta-binomial (\code{type="FBB"}), the beta-binomial (\code{"type=BetaBin"}), and the binomial one (\code{"type=Bin"}).
#'
#'
#' @param formula an object of class \code{`formula`}: a symbolic description of the model to be fitted (of type \code{y ~ x} or \code{y ~ x | z}).
#' @param data an optional data frame, list, or object that is coercible to a data frame through \code{base::as.data.frame} containing the variables in the model. If not found in data, the variables in formula are taken from the environment from which the function flexreg is called.
#' @param type a character specifying the type of regression model. Current options are the flexible beta-binomial \code{"FBB"} (default), the beta-binomial \code{"BetaBin"}, and the binomial one \code{"Bin"}.
#' @param n the total number of trials.
#' @param link.mu a character specifying the link function for the mean model (mu). Currently, \code{"logit"} (default), \code{"probit"}, \code{"cloglog"}, and \code{"loglog"} are supported.
#' @param prior.beta a character specifying the prior distribution for the \code{beta} regression coefficients of the mean model. Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.beta a positive numeric (vector of length 1) specifying the hyperprior standard deviation parameter for the prior distribution of \code{beta} regression coefficients. A value of 100 is suggested if the prior is \code{"normal"}, 2.5 if \code{"cauchy"}.
#' @param link.theta a character specifying the link function for the overdispersion model (theta). Currently, \code{"identity"} (default), \code{"logit"}, \code{"probit"}, \code{"cloglog"}, and \code{"loglog"} are supported. If \code{link.theta = "identity"}, the prior distribution for \code{theta} is a beta.
#' @param hyper.theta.a a numeric (vector of length 1) specifying the first shape parameter for the beta prior distribution of \code{theta}.
#' @param hyper.theta.b a numeric (vector of length 1) specifying the second shape parameter for the beta prior distribution of \code{theta}.
#' @param prior.psi a character specifying the prior distribution for \code{psi} regression coefficients of the overdispersion model (not supported if \code{link.theta="identity"}). Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.psi a positive numeric (vector of length 1) specifying the hyperprior standard deviation parameter for the prior distribution of \code{psi} regression coefficients. A value of 100 is suggested if the prior is \code{"normal"}, 2.5 if \code{"cauchy"}.
#' @param n.iter 	a positive integer specifying the number of iterations for each chain (including warmup). The default is 5000.
#' @param burnin.perc the percentage of iterations per chain to discard.
#' @param n.chain a positive integer specifying the number of Markov chains. The default is 1.
#' @param thin a positive integer specifying the period for saving samples. The default is 1.
#' @param verbose \code{TRUE} (default) or \code{FALSE}: flag indicating whether to print intermediate output.
#' @param ... additional arguments for \code{rstan::sampling}.
#'
#' @return The \code{flexreg_binom} function returns an object of class \code{`flexreg`}, i.e. a list with the following elements:
#' \item{\code{call}}{the function call.}
#' \item{\code{formula}}{the original formula.}
#' \item{\code{link.mu}}{a character specifing the link function in the mean model.}
#' \item{\code{link.theta}}{a character specifing the link function in the overdispersion model.}
#' \item{\code{model}}{an object of class \code{`stanfit`} containing the fitted model.}
#' \item{\code{response}}{the response variable, assuming values in (0, 1).}
#' \item{\code{design.X}}{the design matrix for the mean model.}
#' \item{\code{design.Z}}{the design matrix for the overdispersion model (if defined).}
#'
#' @details Let Y be a random variable whose distribution can be specified in the \code{type} argument and \eqn{\mu} be the mean of Y/n.
#' The \code{flexreg_binom} function links the parameter \eqn{\mu} to a linear predictor through a function  \eqn{g(\cdot)} specified in \code{link.mu}:
#' \deqn{g(\mu_i) = x_i^t \bold{\beta},} where \eqn{\bold{\beta}} is the vector of regression coefficients for the mean model.
#' By default, \code{link.theta="identity"}, meaning that the overdispersion parameter \eqn{\theta} is assumed to be constant.
#' It is possible to extend the model by linking \eqn{\theta} to an additional (possibly overlapping) set of covariates through a proper link
#' function \eqn{q(\cdot)}  specified in the \code{link.theta} argument: \deqn{q(\theta_i) = z_i^t \bold{\psi},} where \eqn{\bold{\psi}} is the vector of regression coefficients for the overdispersion model.
#' In \code{flexreg_binom}, the regression model for the mean and, where appropriate, for the overdispersion parameter can be specified in the
#' \code{formula} argument with a formula of type \eqn{y \sim x_1 + x_2 | z_1 + z_2} where covariates on the left of ("|") are included in the regression model
#' for the mean and covariates on the right of ("|") are included in the regression model for the overdispersion.
#'
#' If the second part is omitted, i.e., \eqn{y \sim x_1 + x_2}, the overdispersion is assumed constant for each observation.
#'
#' @references {
#' Ascari, R., and Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @examples
#' \dontrun{
#' data(Bacteria)
#' fbb <- flexreg_binom(y~females, n=n, data=Bacteria, type="FBB")
#' }
#' @import Rcpp methods rstan faraway
#'
#' @export

flexreg_binom <- function(formula, data, type="FBB", n=NULL,
                     link.mu="logit", prior.beta = "normal", hyperparam.beta = 100,
                     hyper.theta.a=NULL, hyper.theta.b=NULL,
                     link.theta=NULL, prior.psi = NULL, hyperparam.psi = NULL,
                     n.iter=5000, burnin.perc=.5, n.chain=1, thin=1, verbose=TRUE, ...)
{
  cl <- match.call()

  if(is.na(match(type, c("FBB", "Bin", "BetaBin")))) stop("Please specify the type of model correctly.")

  if (missing(data)) data <- environment(formula) else  data <- as.data.frame(data)

  ###############################################
  ###############################################
  # Aggiungere controllo che deve cercare n in data oppure nello stesso environment di formula

  formula <- Formula::as.Formula(formula)


  if(is.character(link.theta)) link_code_theta <-
                  pmatch(link.theta, c("identity", "logit", "probit", "cloglog", "loglog"))

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
    #link_code_theta <- 1

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
  if(length(formula)[2] < 2){
    n <- model.frame(update(formula, ~. + n),
                     data)$n
  } else if(length(formula)[2] >=2){
    n <- model.frame(update(formula, . ~ . + n|.),
                     data)$n
  }


  if(N < 1)
    stop("Empty model")
  if(min(y) < 0)
    stop("Invalid dependent variable: y cannot assume negative values.")
  if(max(y/n) > 1)
    stop("Invalid dependent variable: y must be less or equal to n for every observation.")

  # Design matrix:
  X <- model.matrix(formula, data = mf, rhs = 1)
  if((model.theta == F)){
    Z <- NULL
  } else {
    Z <- model.matrix(formula, data = mf, rhs = 2)
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

                     n.iter = n.iter, burnin.perc = burnin.perc,
                     n.chain = n.chain, thin = thin, verbose=verbose, ...)


  output <- list(call=cl, formula=formula, link.mu=link.mu, link.theta=link.theta,
                 model=model, response=y, design.X=X, design.Z=Z, n=n)
  class(output)<-"flexreg"
  print(summary(output))
  invisible(output)
  #return(output)
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
                      n.iter, burnin.perc, n.chain, thin, verbose, ...){

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
    iter = n.iter, warmup = round(n.iter*burnin.perc),
    refresh = verbose*n.iter/10, #show an update @ each %10
    # seed=1
    ...
  )

  output <- fit.stan
  return(output)

}




















