#' Beta probability density function
#'
#' The function computes the probability density function of the beta distribution with a mean-precision parameterization.
#' @param x a vector of quantiles.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a positive real value.
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The beta distribution has density
#' \deqn{\frac{\Gamma{(\phi)}}{\Gamma{(\mu\phi)}\Gamma{((1-\mu)\phi)}}x^{\mu\phi-1}(1-x)^{(1-\mu)\phi-1}}
#' for \eqn{0<x<1}, where \eqn{0<\mu<1} identifies the mean and \eqn{\phi>0} is the precision parameter.
#'
#' @examples dBeta_mu(x = c(.5,.7,.8), mu = 0.3, phi = 20)
#'
#' @references{
#' Ferrari, S.L.P., and Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501
#' }
#'
#' @import stats
#'
#' @export
#'
#'
dBeta_mu <- function(x, mu, phi){
  if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi
  return(dbeta(x, shape1=alpha1, shape2=alpha2))
}


#' Flexible beta probability density function
#'
#' The function computes the probability density function of the flexible beta distribution.
#'
#' @param x a vector of quantiles.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a positive real value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among clusters. It must lie in (0, 1).
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The FB distribution is a special mixture of two beta distributions
#' \deqn{p Beta(x|\lambda_1,\phi)+(1-p)Beta(x|\lambda_2,\phi)}
#'  for \eqn{0<x<1} where \eqn{Beta(x|\cdot,\cdot)} is the beta distribution with a mean-precision parameterization.
#'  Moreover, \eqn{0<p<1} is the mixing weight, \eqn{\phi>0} is a precision parameter,
#'  \eqn{\lambda_1=\mu+(1-p)w} and \eqn{\lambda_2=\mu-pw} are the component means of the first and second component of the mixture,
#'  \eqn{0<\mu=p\lambda_1+(1-p)\lambda_2<1} is the overall mean, and \eqn{0<w<1} is the  normalized distance between clusters.
#
#' @examples dFB(x = c(.5,.7,.8), mu = 0.3, phi = 20, p = .5, w = .5)
#'
#' @references {
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079
#' }
#'
#' @import stats
#'
#' @export

dFB <- function(x, mu, phi, p, w){
  if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")

  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde
  return(p*dBeta_mu(x,lambda1,phi) + (1-p)*dBeta_mu(x,lambda2,phi))
}


#' Variance-inflated beta probability density function
#'
#' The function computes the probability density function of the variance-inflated beta distribution.
#'
#' @param x a vector of quantiles.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a positive real value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param k the extent of the variance inflation. It must lie in (0, 1).
#'
#' @return A vector with the same length as \code{x}.
#'
#'
#' @details The VIB distribution is a special mixture of two beta distributions
#' \deqn{p Beta(x|\mu,\phi k)+(1-p)Beta(x|\mu,\phi)}
#'  for \eqn{0<x<1} where \eqn{Beta(x|\cdot,\cdot)} is the beta distribution with a mean-precision parameterization.
#'  Moreover, \eqn{0<p<1} is the mixing weight, \eqn{0<\mu<1} represents the overall (as well as mixture component)
#'  mean, \eqn{\phi>0} is a precision parameter, and \eqn{0<k<1} determines the extent of the variance inflation.
#'
#' @examples dVIB(x = c(.5,.7,.8), mu = 0.3, phi = 20, p = .5, k= .5)
#'
#' @references {
#' Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020) Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309.
#' doi:10.1177/1471082X18821213
#' }
#'
#' @import stats
#'
#' @export
dVIB <- function(x, mu, phi, p, k){
  if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")

  return(p*dBeta_mu(x,mu,phi*k) + (1-p)*dBeta_mu(x,mu,phi))
}


#' Beta-binomial probability mass function
#'
#' The function computes the probability mass function of the beta-binomial distribution.
#' @param x a vector of quantiles.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion parameter. It must  lie in (0, 1).
#' @param phi the precision parameter. It is an alternative way to specify the \code{theta} parameter. It must be a positive real value.
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The beta-binomial distribution has probability mass function
#' \deqn{{n\choose y} \frac{\Gamma{(\phi)}}{\Gamma{(\mu\phi)}\Gamma{((1-\mu)\phi)}} \frac{\Gamma{(\mu\phi+y)}\Gamma{((1-\mu)\phi + n - y)}}{\Gamma{(\phi + n)}},}
#' for \eqn{x \in \lbrace 0, 1, \dots, n \rbrace}, where \eqn{0<\mu<1} identifies the mean and \eqn{\phi=(1-\theta)/\theta >0} is the precision parameter.
#'
#' @examples dBetaBin(x = 5, size = 10, mu = .3, phi = 10)
#'
#' @references{
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @import stats
#'
#' @export
#'
#'
dBetaBin <- function(x, size, mu, theta=NULL, phi=NULL){
  if (any(x < 0 | x > size) | !all.equal(x, as.integer(x))) stop("x has to be an integer between 0 and n")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(size < 0 | size != as.integer(size))) stop("size must be a non-negative integer")

  if (!is.null(theta) & !is.null(phi)) {
    if(theta != 1/(phi+1)) stop("Please specify 'theta' or 'phi' but not both") else
      warning("In dFBB() specify 'theta' or 'phi' but not both")
  } else if (is.null(theta) & is.null(phi)) {
    stop("Pleasy specify 'theta' or 'phi' (but not both)")
  } else if (is.null(theta) & !is.null(phi)) {
    if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  } else {
    if (any(theta < 0 | theta > 1)) stop("Parameter theta must lie in (0,1)")
    phi <- (1-theta)/theta
  }

  return(exp(lgamma(size+1)-lgamma(x+1)-lgamma(size-x+1) +
               lgamma(phi)+lgamma((phi*mu)+x) + lgamma((phi*(1-mu))+size-x) -
               lgamma(phi+size) - lgamma(phi*mu) - lgamma(phi*(1-mu))
             ))
}

#' Flexible beta-binomial probability mass function
#'
#' The function computes the probability mass function of the flexible beta-binomial distribution.
#'
#' @param x a vector of quantiles.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion paramete. It must  lie in (0, 1).
#' @param phi the precision parameter. It is an alternative way to specify the \code{theta} parameter. It must be a positive real value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among clusters. It must lie in (0, 1).
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The FBB distribution is a special mixture of two beta-binomial distributions
#' \deqn{p BB(x|\lambda_1,\phi)+(1-p)BB(x|\lambda_2,\phi)}
#'  for \eqn{x \in \lbrace 0, 1, \dots, n \rbrace} where \eqn{BB(x|\cdot,\cdot)} is the beta-binomial distribution with a mean-precision parameterization.
#'  Moreover, \eqn{\phi=(1-\theta)/\theta}, \eqn{0<p<1} is the mixing weight, \eqn{\phi>0} is a precision parameter,
#'  \eqn{\lambda_1=\mu+(1-p)w} and \eqn{\lambda_2=\mu-pw} are the component means of the first and second component of the mixture,
#'  \eqn{0<\mu=p\lambda_1+(1-p)\lambda_2<1} is the overall mean, and \eqn{0<w<1} is the  normalized distance between clusters.
#
#' @examples dFBB(x = c(5,7,8), size=10, mu = 0.3, phi = 20, p = .5, w = .5)
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @import stats
#'
#' @export

dFBB <- function(x, size, mu, theta=NULL, phi=NULL, p, w){
  if (any(x < 0 | x > size) | !all.equal(x, as.integer(x))) stop("x has to be an integer between 0 and n")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(size < 0 | size != as.integer(size))) stop("size must be a non-negative integer")

  if (!is.null(theta) & !is.null(phi)) {
    if(theta != 1/(phi+1)) stop("Please specify 'theta' or 'phi' but not both") else
      warning("In dFBB() specify 'theta' or 'phi' but not both")
  } else if (is.null(theta) & is.null(phi)) {
    stop("Pleasy specify 'theta' or 'phi' (but not both)")
  } else if (is.null(theta) & !is.null(phi)) {
    if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  } else {
    if (any(theta < 0 | theta > 1)) stop("Parameter theta must lie in (0,1)")
    phi <- (1-theta)/theta
  }
  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde
  return(p*dBetaBin(x=x, size=size, mu=lambda1, phi=phi) + (1-p)*dBetaBin(x=x, size=size, mu=lambda2, phi=phi))
}






#' Draw density plots
#'
#' The function draws a curve corresponding to the probability density/mass function of the specified distribution (beta, flexible beta, variance-inflated beta, binomial, beta-binomial, or flexible beta-binomial).
#' @param type a character specifying the distribution type to be plotted (\code{"Beta"}, \code{"FB"}, \code{"VIB"}, \code{"Bin"}, \code{"BetaBin"}, or \code{"FBB"}).
#' @param size the total number of trials (to be specified if \code{type} is \code{"Bin"}, \code{"BetaBin"}, or \code{"FBB"}).
#' @param mu the mean parameter of the distribution. It must lie in (0, 1).
#' @param theta the overdispersion parameter (to be specified if \code{type} is \code{"BetaBin"} or \code{"FBB"}). It must  lie in (0, 1).
#' @param phi  the precision parameter (an alternative way to specify the \code{theta} parameter if \code{type} is \code{"BetaBin"} or \code{"FBB"}). It must be a positive real value.
#' @param p  the mixing weight (to be specified if \code{type} is \code{"FB"} or \code{"VIB"}). It must lie in (0, 1).
#' @param w  the normalized distance among clusters of the FB distribution (to be specified if \code{type = "FB"}). It must lie in (0, 1).
#' @param k  the extent of the variance inflation (to be specified if \code{type = "VIB"}). It must lie in (0, 1).
#' @param ... additional arguments of \code{stat_function()}.
#'
#'
#' @examples
#' curve.density("Beta", mu=0.5, phi=20)
#' curve.density("FB", mu=0.5, phi=20, p=0.4, w=.8)
#' curve.density("VIB", mu=0.5, phi=20, p=0.9, k=.8, col=3)
#'
#' curve.density("Bin", size=10, mu=.7)
#' curve.density("BetaBin", size=10, mu=.7, phi=10)
#' curve.density("FBB", size=10, mu=.7, phi=10, p=.2,w=.7)
#'
#'
#' @references{
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018) A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' \cr
#' Ascari, R., and Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020) Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309. doi:10.1177/1471082X18821213 \cr
#' \cr
#'  Ferrari, S.L.P., and Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501
#'}
#'
#'
#' @import stats ggplot2
#'
#' @export
#'
#'
curve.density <- function(type = NULL, size = NULL, mu = NULL, theta = NULL, phi = NULL, p = NULL,
                          w = NULL, k = NULL, ...){
  if(is.null(mu) | ((is.null(phi)&(is.null(theta))) & type!="Bin")) stop("Specify parameters according to the chosen distribution")
  if(type == "FB" & (is.null(p) | is.null(w))) stop("Specify parameters according to the chosen distribution")
  if(type == "VIB" & (is.null(p) | is.null(k))) stop("Specify parameters according to the chosen distribution")

  if(type == "Beta" & !is.null(p)) stop("Specify parameters according to the chosen distribution")
  if(!(type %in% c("FB","FBB")) & !is.null(w)) stop("Specify parameters according to the chosen distribution")
  if(type != "VIB" & !is.null(k)) stop("Specify parameters according to the chosen distribution")

  if((type %in% c("FB","Beta", "VIB"))){
    if(!is.null(size)) stop("Parameter size must not be specified if type is ''Beta'', ''FB'', or ''VIB''.")
    if(type == "Beta") {
      fun <- function(x) dBeta_mu(x,mu,phi)
    } else if(type == "FB") {
      fun <- function(x) dFB(x, mu = mu, phi = phi, p = p, w = w)
      } else if(type == "VIB") {
        fun <- function(x) dVIB(x, mu = mu, phi = phi, p = p, k = k)
        } else stop("Please specify the type of distribution correctly.")


    title.default <- paste(type," density function", sep="")

    pp <- ggplot(data = data.frame(x = 0))+
      stat_function(fun = fun, ...) + xlim(0,1) + xlab("x") + ylab("Density") + ggtitle(title.default)
    } else if((type %in% c("FBB","Bin", "BetaBin"))){

      if(type=="FBB"){
        if(!is.null(theta) & is.null(phi)) phi <- (1-theta)/theta
        if(is.null(theta) & is.null(phi)) stop("Please specify a value for theta or phi")
        prob <- dFBB(0:size, size=size, mu=mu, phi=phi, w=.6, p=.5)
      } else if(type=="BetaBin"){
        if(!is.null(theta) & is.null(phi)) phi <- (1-theta)/theta
        if(is.null(theta) & is.null(phi)) stop("Please specify a value for theta or phi")
        prob <- dBetaBin(0:size, size=size, mu=mu, phi=phi)
      } else prob <- dbinom(0:size, size=size, prob=mu)

      y <- NULL
      pp <- ggplot(data=data.frame(y=0:size, prob=prob)) +
        geom_point(aes(x=y, y=prob),size=2) +
        geom_line(aes(x=y, y=prob),alpha=.3) +
        ylim(0,max(prob)+.05)

    } else stop("Error: please specify a valide type argument.")

  return(pp)
}
