#' Flexible Regression Models for Bounded Continuous Responses
#'
#' @description The function fits some flexible regression models for bounded continuous responses (e.g., proportions and rates) via a Bayesian approach to inference based on Hamiltonian Monte Carlo algorithm.
#' Available regression models are the flexible beta regression model (\code{type = "FB"}, default), the variance inflated beta (\code{type = "VIB"}),  the beta  (\code{type = "Beta"}), as well as their augmented versions.
#'
#'
#' @param formula an object of class "\code{\link{formula}}": a symbolic description of the mean model (\code{y ~ x}) or the mean and precision models (\code{y ~ x | z}) to be fitted (see Details).
#' @param zero.formula an object of class "\code{\link{formula}}": a symbolic description of the zero augmented model to be fitted (see Details).
#' @param one.formula an object of class "\code{\link{formula}}": a symbolic description of the one augmented model to be fitted (see Details).
#' @param data an optional  \code{data.frame}, list, or object that is coercible to a  \code{data.frame} through \code{\link{as.data.frame}} containing the variables in the model. If not found in \code{data}, the variables in \code{formula}, \code{zero.formula}, and \code{one.formula} are taken from the environment from which the function \code{\link{flexreg}} is called.
#' @param type a character specifying the type of regression model. Current options are \code{"FB"} (flexible beta, default), \code{"VIB"} (variance inflated beta), and \code{"Beta"}.
#' @param link.mu a character specifying the link function for the mean model (mu). Currently, \code{"logit"} (default), \code{"probit"}, \code{"cloglog"}, and \code{"loglog"} are supported.
#' @param prior.beta a character specifying the prior distribution for the  regression coefficients of the mean model, \code{beta}. Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.beta a positive numeric (vector of length 1) specifying the hyperprior scale parameter for the prior distribution of \code{beta} regression coefficients. The default is 100 if the prior is \code{"normal"}, 2.5 if it is \code{"cauchy"}.
#' @param prior.omega0 a character specifying the prior distribution for the  regression coefficients of the augmented model in zero, \code{omega0}. Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.omega0 a positive numeric (vector of length 1) specifying the hyperprior scale parameter for the prior distribution of \code{omega0} regression coefficients. The default is 100 if the prior is \code{"normal"}, 2.5 if it is \code{"cauchy"}.
#' @param prior.omega1 a character specifying the prior distribution for the  regression coefficients of the augmented model in one, \code{omega1}. Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.omega1 a positive numeric (vector of length 1) specifying the hyperprior scale parameter for the prior distribution of \code{omega1} regression coefficients. The default is 100 if the prior is \code{"normal"}, 2.5 if it is \code{"cauchy"}.
#' @param link.phi a character specifying the link function for the precision model (phi). Currently, \code{"identity"} (default), \code{"log"}, and \code{"sqrt"} are supported.
#' @param prior.phi a character specifying the prior distribution for precision parameter \code{phi} if \cr \code{link.phi = "identity"}. Currently, \code{"gamma"} (default) and \code{"unif"} are supported.
#' @param hyperparam.phi a positive numeric (vector of length 1) specifying the hyperprior parameter for the prior distribution of \code{phi}. If the prior is \code{"gamma"}, the value identifies the gamma's shape and rate parameters (the default is 0.001). If the prior is \code{"uniform"} the hyperparameter must be specified to define the upper limit of the support of \code{phi}.
#' @param prior.psi a character specifying the prior distribution for the regression coefficients of the precision model \code{psi} (not supported if \code{link.phi = "identity"}). Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.psi a positive numeric (vector of length 1) specifying the hyperprior scale parameter for the prior distribution of \code{psi} regression coefficients. The default is 100 if the prior is \code{"normal"}, 2.5 if it is \code{"cauchy"}.
#' @param n.chain a positive integer specifying the number of Markov chains. The default is 1.
#' @param n.iter 	a positive integer specifying the number of iterations for each chain (including warm-up). The default is 5000.
#' @param warmup.perc the percentage of iterations per chain to discard.
#' @param thin a positive integer specifying the period for saving samples. The default is 1.
#' @param verbose a logical (with default \code{TRUE}) indicating whether to print intermediate output.
#' @param ... additional arguments from \code{\link[rstan]{sampling}}.
#'
#' @return The \code{\link{flexreg}} function returns an object of class \code{`flexreg`}, i.e. a list with the following elements:
#' \item{\code{call}}{the function call.}
#' \item{\code{type}}{the type of regression model.}
#' \item{\code{formula}}{the overall formula.}
#' \item{\code{aug}}{a character specifing the absence of the augmentation (\code{"No"}) or the presence of augmentation in zero (\code{"0"}), one (\code{"1"}), or both (\code{"01"}).}
#' \item{\code{link.mu}}{a character specifing the link function in the mean model.}
#' \item{\code{link.phi}}{a character specifing the link function in the precision model.}
#' \item{\code{model}}{an object of class \code{`stanfit`} containing the fitted model.}
#' \item{\code{response}}{the response variable, assuming values in (0, 1).}
#' \item{\code{design.X}}{the design matrix for the mean model.}
#' \item{\code{design.Z}}{the design matrix for the precision model (if defined).}
#' \item{\code{design.X0}}{the design matrix for the augmented model in zero (if defined).}
#' \item{\code{design.X1}}{the design matrix for the augmented model in one (if defined).}
#'
#' @details Let Y be a continuous bounded random variable whose distribution can be specified in the \code{type} argument and \eqn{\mu} be the mean of Y.
#' The \code{\link{flexreg}} function links the parameter \eqn{\mu} to a linear predictor through a function  \eqn{g_1(\cdot)} specified in \code{link.mu}:
#' \deqn{g_1(\mu) = \bold{x}^t \bold{\beta},} where \eqn{\bold{\beta}} is the vector of regression coefficients for the mean model.
#' The prior distribution and the related hyperparameter of \eqn{\bold{\beta}} can be specified in \code{prior.beta} and \code{hyperparam.beta}, respectively.
#' By default, the precision parameter \eqn{\phi} is assumed to be constant.
#' The prior distribution and the related hyperparameter of \eqn{\phi} can be specified in \code{prior.phi} and \code{hyperparam.phi}.
#' It is possible to extend the model by linking \eqn{\phi} to an additional (possibly overlapping) set of covariates through a proper link
#' function \eqn{g_2(\cdot)}  specified in the \code{link.phi} argument:
#' \deqn{g_2(\phi) = \bold{z}^t \bold{\psi},} where \eqn{\bold{\psi}} is the vector of regression coefficients for the precision model.
#' The prior distribution and the related hyperparameter of \eqn{\bold{\psi}} can be specified in \code{prior.psi} and \code{hyperparam.psi}.
#' In the function \code{\link{flexreg}}, the regression model for the mean and, where appropriate, for the precision parameter can be specified in the
#' \code{formula} argument with a formula of type \code{y ~ x1 + x2 | z1 + z2} where covariates on the left of "|" are included in the regression model
#' for the mean, whereas covariates on the right of "|" are included in the regression model for the precision.
#'
#' If the second part is omitted, i.e., \code{y ~ x1 + x2}, the precision is assumed constant for each observation.
#'
#' In presence of zero values in the response, one has to link the parameter \eqn{q_0}, i.e., the probability of augmentation in zero,  to an additional (possibly overlapping) set of covariates through a logit link function:
#' \deqn{g_3(q_{0}) = \bold{x}_{0}^t \bold{\omega_0},} where \eqn{\bold{\omega_0}} is the vector of regression coefficients for the augmented model in zero.
#' The prior distribution and the related hyperparameter of \eqn{\bold{\omega_0}} can be specified in \code{prior.omega0} and \code{hyperparam.omega0}.
#' In presence of one values in the response, one has to link the parameter \eqn{q_1}, i.e., the probability of augmentation in one,  to an additional (possibly overlapping) set of covariates through a logit link function:
#' \deqn{g_4(q_{1}) = \bold{x}_{1}^t \bold{\omega_1},} where \eqn{\bold{\omega_1}} is the vector of regression coefficients for the augmented model in one.
#' The prior distribution and the related hyperparameter of \eqn{\bold{\omega_1}} can be specified in \code{prior.omega1} and \code{hyperparam.omega1}.
#' If both the augmented models in zero and one are specified, the link function is a bivariate logit.
#' In \code{\link{flexreg}} function, the augmented models in zero and/or one can be specified in the
#' \code{zero.formula} and/or \code{one.formula} arguments with a formula of type \code{ ~ x}.
#' Left hand side in \code{zero.formula} and \code{one.formula} can be omitted; if specified, they have to be the same as left hand side in \code{formula}.

#'
#' @references {
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020). Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309. doi:10.1177/1471082X18821213 \cr
#' \cr
#'  Ferrari, S.L.P., Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018) A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' }
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, type="FB")
#'
#' # Regression model with one augmentation:
#' AFB1 <- flexreg(accuracy ~ dyslexia | iq + dyslexia + iq:dyslexia,
#' one.formula = ~ iq + dyslexia, data = Reading, type="FB")
#'}
#' @import Rcpp methods rstan
#'
#' @export

flexreg <- function(formula, zero.formula=NULL, one.formula=NULL, data,
                    type="FB", link.mu="logit",
                    prior.beta = "normal", hyperparam.beta = NULL,
                    prior.omega0 = "normal", hyperparam.omega0 = NULL,
                    prior.omega1 = "normal", hyperparam.omega1 = NULL,
                    link.phi = NULL, prior.phi = NULL, hyperparam.phi = NULL,
                    prior.psi = NULL, hyperparam.psi = NULL,
                    n.chain=1, n.iter=5000, warmup.perc=.5, thin=1, verbose=TRUE, ...)
{
  cl <-   match.call()

  if(is.na(match(type, c("FB", "Beta", "VIB")))) stop("Please specify the type of model correctly.")

  if (missing(data)) data <- environment(formula) else  data <- as.data.frame(data)

  formula <- Formula::as.Formula(formula)

  if (is.character(link.phi)) link_code_phi <- pmatch(link.phi, c("identity", "log", "sqrt"))
  if(length(formula)[2] >= 2){
    model.phi <- TRUE
    if(is.null(link.phi)) link_code_phi <- 2
    if(link_code_phi <2) stop("Error: invalid link function for regression model for phi")
    if(!is.null(prior.phi) | !is.null(hyperparam.phi)) stop("Error: prior chosen for phi but regression model for phi in formula. Please define priors for coefficients psi instead.")
    #if(length(formula)[2] > 2)  warning("formula must not have more than two RHS parts")
  } else if(length(formula)[2] < 2){
    model.phi <- FALSE
    if(is.null(link.phi)) link_code_phi <- 1
    if(link_code_phi >1) stop("Error: covariates for regression model for phi must be specified")
    if(!is.null(prior.psi) | !is.null(hyperparam.psi)) stop("Error: prior for psi coefficients defined but no regression model for phi. Please define covariates for regression model for phi")
  }

  if(is.null(zero.formula)){
    #aug.zero <- F
    #formula0 <- NULL
    link_prior_omega0 <- hyperparam.omega0 <- X0 <- NULL
  } else {
    #aug.zero <- T
    formula0 <- Formula::as.Formula(zero.formula)
    if(!is.null(attr(formula0, "lhs")[[1]])){
      if(attr(formula, "lhs")[[1]]!=attr(formula0, "lhs")[[1]]) stop("Left hand side in formula and zero.formula have to be the same. Left hand side in zero.formula can be omitted ")
    }
    mf0 <- model.frame(formula0, data)
    X0 <- model.matrix(formula0, data = mf0, rhs = 1)
  }

  if(is.null(one.formula)){
    #aug.one <- F
    #formula1 <- NULL
    link_prior_omega1 <- hyperparam.omega1 <- X1 <- NULL
  } else {
    #aug.one <- T
    formula1 <- Formula::as.Formula(one.formula)
    if(!is.null(attr(formula1, "lhs")[[1]])){
      if(attr(formula, "lhs")[[1]]!=attr(formula1, "lhs")[[1]]) stop("Left hand side in formula and one.formula have to be the same. Left hand side in one.formula can be omitted ")
    }
    mf1 <- model.frame(formula1, data)
    X1 <- model.matrix(formula1, data = mf1, rhs = 1)
  }

  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  N <- length(y)

  if(N < 1) stop("Empty model")
  if(min(y) == 0 & (is.null(zero.formula))){
    stop("Presence of zeros in the response variable. Specify a model with zero augmentation.")
  } else if(min(y)> 0 & !is.null(zero.formula)){
    stop("A model with zero augmentation is specified but there are no zeros in the response variable.")
  }
  if(max(y) == 1 & (is.null(one.formula))){
    stop("Presence of ones in the response variable. Specify a model with one augmentation.")
  }else if(max(y)<1 & !is.null(one.formula)){
    stop("A model with one augmentation is specified but there are no ones in the response variable.")
  }
  if(min(y)<0 | max(y)>1){
    stop("Invalid dependent variable, all observations must be in [0, 1].")
  }

  X <- model.matrix(formula, data = mf, rhs = 1)
  if(isFALSE(model.phi)){
    Z <- NULL
  } else {
    Z <- model.matrix(formula, data = mf, rhs = 2)
  }

  if(link.mu %in%c("logit", "probit", "cloglog", "loglog")) link_code_mu <- pmatch(link.mu, c("logit", "probit", "cloglog", "loglog")) else
    stop("Invalid link function for mu.")

  if (prior.beta %in% c("normal", "cauchy")) link_prior_beta  <- pmatch(prior.beta, c("normal", "cauchy")) else
    stop("Invalid prior for beta parameter.")

  #hyperparam beta
  if(is.null(hyperparam.beta)) hyperparam.beta <- ifelse(link_prior_beta==1,100,2.5)
  if(hyperparam.beta < 0) stop("Hyperprior for beta coefficients must be positive")
  if(link_prior_beta == 1 & hyperparam.beta < 10) warning("A hyperprior of at least 10 for normal prior for beta coefficients is recommended")
  if(link_prior_beta == 2 & hyperparam.beta != 2.5) warning("A hyperprior of 2.5 for cauchy prior for beta coefficients is recommended")

  if(link_code_phi == 1) {
    prior_code_phi <- ifelse(is.null(prior.phi),1,
                             pmatch(prior.phi, c( "gamma", "unif")))
    if(prior_code_phi==1) hyper_phi <- ifelse(is.null(hyperparam.phi),0.001, hyperparam.phi) else #add warnings if g grande
      if(prior_code_phi==2) hyper_phi <-ifelse(is.null(hyperparam.phi), stop("Please specify an hyperparameter A>0 for uniform prior for phi"), hyperparam.phi)
      link_prior_psi <- NULL
  }  else if (link_code_phi != 1) {
    prior_code_phi <- NULL
    hyper_phi <- NULL
    if(is.null(prior.psi)) link_prior_psi <- 1 else
      if (is.character(prior.psi)) link_prior_psi  <- pmatch(prior.psi, c("normal", "cauchy")) else
        stop("Invalid prior for psi parameter.")

    if(is.null(hyperparam.psi)) hyperparam.psi <- ifelse(link_prior_psi==1, 100, 2.5)
    else {
      if(hyperparam.psi < 0) stop("Hyperprior for psi coefficients must be positive")
      if(link_prior_psi == 1 & hyperparam.psi < 10) warning("An hyperprior of at least 10 for normal for psi coefficients is recommended")
      if(link_prior_psi == 2 & hyperparam.psi != 2.5) warning("An hyperprior of 2.5 for cauchy prior for psi coefficients is recommended")
    }
  }


  #checks on omega1 prior
  if(prior.omega1 %in% c("normal", "cauchy")){
    link_prior_omega1  <- pmatch(prior.omega1, c("normal", "cauchy"))
  } else stop("Invalid prior for omega1 parameter.")
  if(is.null(hyperparam.omega1)) hyperparam.omega1 <- ifelse(link_prior_omega1==1, 100, 2.5)
  if(hyperparam.omega1 < 0) stop("Hyperprior for omega1 coefficients must be positive")
  if(link_prior_omega1 == 1 & hyperparam.omega1 < 10) warning("An hyperprior of at least 10 for normal prior for omega1 coefficients is recommended")
  if(link_prior_omega1 == 2 & hyperparam.omega1 != 2.5) warning("An hyperprior of 2.5 for cauchy prior for omega1 coefficients is recommended")

  #checks on omega0 prior
  if (prior.omega0 %in% c("normal", "cauchy")){
    link_prior_omega0  <- pmatch(prior.omega0, c("normal", "cauchy"))
  } else stop("Invalid prior for omega0 parameter.")
  if(is.null(hyperparam.omega0)) hyperparam.omega0 <- ifelse(link_prior_omega0==1, 100, 2.5)
  if(hyperparam.omega0 < 0) stop("Hyperprior for omega0 coefficients must be positive")
  if(link_prior_omega0 == 1 & hyperparam.omega0 < 10) warning("An hyperprior of at least 10 for normal prior for omega0 coefficients is recommended")
  if(link_prior_omega0 == 2 & hyperparam.omega0 != 2.5) warning("An hyperprior of 2.5 for cauchy prior for omega0 coefficients is recommended")


  if(is.null(one.formula) & is.null(zero.formula)){
    aug_code <- "No"
  }else if(is.null(one.formula)){
    aug_code <- "0"
  } else if (is.null(zero.formula)){
    aug_code <- "1"
  } else aug_code <- "01"


  #put all formulas together
  if(!is.null(zero.formula)|!is.null(one.formula)){
    formula.final <- paste0(deparse(attr(formula, "lhs")[[1]]),"~",
                            deparse(attr(formula, "rhs")[[1]]),"|",
                            ifelse(model.phi==F,1,deparse(attr(formula, "rhs")[[2]])),"|",
                            ifelse(is.null(zero.formula),1,deparse(attr(formula0, "rhs")[[1]])),
                            "|", ifelse(is.null(one.formula),1,deparse(attr(formula1, "rhs")[[1]])))
    formula <- Formula::as.Formula(formula.final)
  }

  ############################################

  model <- fit.model(model.phi = model.phi, type = type, N = N,  y = y,
                     X = X, X0 = X0, X1 = X1, Z = Z, link_code_mu = link_code_mu,
                     aug_code = aug_code,
                     link_prior_beta, hyperparam.beta,
                     link_prior_omega0, hyperparam.omega0,
                     link_prior_omega1, hyperparam.omega1,
                     link_code_phi = link_code_phi,
                     prior_code_phi = prior_code_phi, hyper_phi = hyper_phi,
                     link_prior_psi, hyperparam.psi,
                     n.iter, warmup.perc, n.chain, thin,
                     verbose, ...)

  #questa riga serve per summary
  link.phi <- c("identity", "log", "sqrt")[link_code_phi]
  output <- list(call=cl, type=type, formula=formula, aug = aug_code, link.mu=link.mu, link.phi=link.phi,
                 model=model, response=y, design.X=X, design.Z=Z, design.X0=X0,design.X1=X1)
  class(output)<- c("flexreg", "flexreg_bound")
  #invisible(output)
  return(output)
}

fit.model <- function(model.phi = model.phi, type = type, N = N,  y = y, X = X, X0 = X0, X1 = X1, Z = Z,
                      link_code_mu = link_code_mu, aug_code,
                      link_prior_beta, hyperparam.beta,
                      link_prior_omega0, hyperparam.omega0,
                      link_prior_omega1, hyperparam.omega1,
                      link_code_phi = link_code_phi,
                      prior_code_phi = prior_code_phi, hyper_phi = hyper_phi,
                      link_prior_psi, hyperparam.psi,
                      n.iter, warmup.perc, n.chain, thin, verbose, ...){

  data.stan <- list(
    N = N,
    y = y,
    X = X,
    X0 = X0,
    X1 = X1,
    Z = Z,
    K = ncol(X),
    K0 = ncol(X0),
    K1 = ncol(X1),
    H = ncol(Z),
    link_code_mu = link_code_mu,
    link_prior_beta = link_prior_beta,
    hyperprior_beta = hyperparam.beta,
    link_prior_omega1 = link_prior_omega1,
    hyperprior_omega1 = hyperparam.omega1,
    link_prior_omega0 = link_prior_omega0,
    hyperprior_omega0 = hyperparam.omega0,
    link_code_phi = link_code_phi,
    prior_code_phi = prior_code_phi,
    hyper_phi = hyper_phi,
    link_prior_psi = link_prior_psi,
    hyperprior_psi = hyperparam.psi
  )

  if(model.phi){
    stan.model <- stanmodels[[paste0(type,aug_code, "_phi")]]
  } else {  stan.model <- stanmodels[[paste0(type,aug_code)]]}
  fit.stan = rstan::sampling(
    object = stan.model,#
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
