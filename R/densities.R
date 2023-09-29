#' Probability density function of the beta distribution
#'
#' @description The function computes the probability density function of the beta distribution with a mean-precision parameterization.
#' It can also  compute the probability density function of the augmented beta distribution by assigning positive probabilities to zero and/or one values and a (continuous) beta density to the interval (0,1).
#'
#' @param x a vector of quantiles.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The beta distribution has density
#' \deqn{f_B(x;\mu,\phi)=\frac{\Gamma{(\phi)}}{\Gamma{(\mu\phi)}\Gamma{((1-\mu)\phi)}}x^{\mu\phi-1}(1-x)^{(1-\mu)\phi-1}}
#' for \eqn{0<x<1}, where \eqn{0<\mu<1} identifies the mean and \eqn{\phi>0} is the precision parameter.
#'
#' The augmented beta distribution has density
#' \itemize{
#' \item \eqn{q_0},  if  \eqn{x=0}
#' \item \eqn{q_1}, if  \eqn{x=1}
#' \item \eqn{(1-q_0-q_1)f_B(x;\mu,\phi)}, if \eqn{0<x<1}
#' }
#' where \eqn{0<q_0<1} identifies the augmentation in zero, \eqn{0<q_1<1} identifies the augmentation in one,
#' and \eqn{q_0+q_1<1}.
#'
#' @examples
#' dBeta(x = c(.5,.7,.8), mu = .3, phi = 20)
#' dBeta(x = c(.5,.7,.8), mu = .3, phi = 20, q0 = .2)
#' dBeta(x = c(.5,.7,.8), mu = .3, phi = 20, q0 = .2, q1= .1)
#'
#' @references{
#' Ferrari, S.L.P., Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501
#' }
#'
#' @import stats
#'
#' @export
#'
#'
dBeta <- function(x, mu, phi, q0 = NULL, q1 = NULL){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>1)) stop("The sum of q0 and q1 must be less than 1")

  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi

  fun <- (1-q0-q1)*dbeta(x, shape1=alpha1, shape2=alpha2)
  fun[which(x==0)] <- q0
  fun[which(x==1)] <- q1

  return(fun)
}


#' Probability density function of the flexible beta distribution
#'
#' @description The function computes the probability density function of the flexible beta distribution.
#' It can also  compute the probability density function of the augmented flexible beta distribution by assigning positive probabilities to zero and/or one values and a (continuous) flexible beta density to the interval (0,1).
#'
#' @param x a vector of quantiles.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among component means. It must lie in (0, 1).
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The FB distribution is a special mixture of two beta distributions with probability density function
#' \deqn{f_{FB}(x;\mu,\phi,p,w)=p f_B(x;\lambda_1,\phi)+(1-p)f_B(x;\lambda_2,\phi),}
#'  for \eqn{0<x<1}, where \eqn{f_B(x;\cdot,\cdot)} is the beta density with a mean-precision parameterization.
#'  Moreover, \eqn{0<\mu=p\lambda_1+(1-p)\lambda_2<1} is the overall mean, \eqn{\phi>0} is a precision parameter,
#'  \eqn{0<p<1} is the mixing weight,  \eqn{0<w<1} is the  normalized distance between component means, and
#'  \eqn{\lambda_1=\mu+(1-p)w} and \eqn{\lambda_2=\mu-pw} are the means of the first and second component of the mixture, respectively.
#'
#' The augmented FB distribution has density
#' \itemize{
#' \item \eqn{q_0}, if \eqn{x=0}
#' \item \eqn{q_1}, if \eqn{x=1}
#' \item \eqn{(1-q_0-q_1)f_{FB}(x;\mu,\phi,p,w)}, if \eqn{0<x<1 }
#' }
#' where \eqn{0<q_0<1} identifies the augmentation in zero, \eqn{0<q_1<1} identifies the augmentation in one,
#' and \eqn{q_0+q_1<1}.
#'
#' @examples
#' dFB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5)
#' dFB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5, q0 = .2)
#' dFB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5, q0 = .2, q1 = .1)
#' @references {
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079
#' }
#'
#' @import stats
#'
#' @export

dFB <- function(x, mu, phi, p, w, q0 = NULL, q1 = NULL){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>1)) stop("The sum of q0 and q1 must be less than 1")

  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde

  fun <- (1-q0-q1)*(p*dBeta(x,lambda1,phi) + (1-p)*dBeta(x,lambda2,phi))
  fun[which(x==0)] <- q0
  fun[which(x==1)] <- q1

  return(fun)
}


#' Probability density function of the variance-inflated beta distribution
#'
#' @description The function computes the probability density function of the variance-inflated beta distribution.
#' It can also  compute the probability density function of the augmented variance-inflated beta distribution by assigning positive probabilities to zero and/or one values and a (continuous) variance-inflated beta density to the interval (0,1).
#'
#' @param x a vector of quantiles.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param k the extent of the variance inflation. It must lie in (0, 1).
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The VIB distribution is a special mixture of two beta distributions with probability density function
#' \deqn{f_{VIB}(x;\mu,\phi,p,k)=p f_B(x;\mu,\phi k)+(1-p)f_B(x;\mu,\phi),}
#'  for \eqn{0<x<1}, where \eqn{f_B(x;\cdot,\cdot)} is the beta density with a mean-precision parameterization.
#'  Moreover, \eqn{0<p<1} is the mixing weight, \eqn{0<\mu<1} represents the overall (as well as mixture component)
#'  mean, \eqn{\phi>0} is a precision parameter, and \eqn{0<k<1} determines the extent of the variance inflation.
#' The augmented VIB distribution has density
#' \itemize{
#' \item \eqn{q_0}, if \eqn{x=0}
#' \item \eqn{q_1},  if \eqn{x=1}
#' \item \eqn{(1-q_0-q_1)f_{VIB}(x;\mu,\phi,p,k)}, if \eqn{0<x<1}
#' }
#' where \eqn{0<q_0<1} identifies the augmentation in zero, \eqn{0<q_1<1} identifies the augmentation in one,
#' and \eqn{q_0+q_1<1}.
#'
#' @examples
#' dVIB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5)
#' dVIB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q1 = .1)
#' dVIB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q0 = .2, q1 = .1)
#'
#' @references {
#' Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020). Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309.
#' doi:10.1177/1471082X18821213
#' }
#'
#' @import stats
#'
#' @export
dVIB <- function(x, mu, phi, p, k, q0 = NULL, q1 = NULL){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>1)) stop("The sum of q0 and q1 must be less than 1")

  fun <- (1-q0-q1)*(p*dBeta(x,mu,phi*k) + (1-p)*dBeta(x,mu,phi))
  fun[which(x==0)] <- q0
  fun[which(x==1)] <- q1

  return(fun)
}

#' Probability mass function of the beta-binomial distribution
#'
#' The function computes the probability mass function of the beta-binomial distribution.
#' @param x a vector of quantiles.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion parameter. It must  lie in (0, 1).
#' @param phi the precision parameter, an alternative way to specify the overdispersion parameter \code{theta}. It must be a real positive value.
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The beta-binomial distribution has probability mass function
#' \deqn{f_{BB}(x;\mu,\phi)={n\choose x} \frac{\Gamma{(\phi)}}{\Gamma{(\mu\phi)}\Gamma{((1-\mu)\phi)}} \frac{\Gamma{(\mu\phi+x)}\Gamma{((1-\mu)\phi + n - x)}}{\Gamma{(\phi + n)}},}
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

#' Probability mass function of the flexible beta-binomial distribution
#'
#' The function computes the probability mass function of the flexible beta-binomial distribution.
#'
#' @param x a vector of quantiles.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion parameter. It must  lie in (0, 1).
#' @param phi the precision parameter, an alternative way to specify the overdispersion parameter \code{theta}. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among component means. It must lie in (0, 1).
#'
#' @return A vector with the same length as \code{x}.
#'
#' @details The FBB distribution is a special mixture of two beta-binomial distributions with probability mass function
#' \deqn{f_{FBB}(x;\mu,\phi,p,w) = p BB(x;\lambda_1,\phi)+(1-p)BB(x;\lambda_2,\phi),}
#'  for \eqn{x \in \lbrace 0, 1, \dots, n \rbrace}, where \eqn{BB(x;\cdot,\cdot)} is the beta-binomial distribution with a mean-precision parameterization.
#'  Moreover, \eqn{\phi=(1-\theta)/\theta>0} is a precision parameter, \eqn{0<p<1} is the mixing weight, \eqn{0<\mu=p\lambda_1+(1-p)\lambda_2<1} is the overall mean,
#'   \eqn{0<w<1} is the  normalized distance between component means, and
#'  \eqn{\lambda_1=\mu+(1-p)w} and \eqn{\lambda_2=\mu-pw} are the scaled means of the first and second component of the mixture, respectively.
#'
#
#' @examples dFBB(x = c(5,7,8), size=10, mu = .3, phi = 20, p = .5, w = .5)
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
#' For beta, flexible beta, and variance-inflated beta, it also allows to include the representation of the probability of augmentation in zero and/or one values.
#' @param type a character specifying the distribution type to be plotted (\code{"Beta"}, \code{"FB"}, \code{"VIB"}, \code{"Bin"}, \code{"BetaBin"}, or \code{"FBB"}).
#' @param size the total number of trials (to be specified only if \code{type} is \code{"Bin"}, \code{"BetaBin"}, or \code{"FBB"}).
#' @param mu the mean parameter of the distribution. It must lie in (0, 1).
#' @param theta the overdispersion parameter (to be specified only if \code{type} is \code{"BetaBin"} or \code{"FBB"}). It must  lie in (0, 1).
#' @param phi  the precision parameter (if \code{type} is \code{"BetaBin"} or \code{"FBB"}, it represents an alternative way to specify the \code{theta} parameter). It must be a real positive value.
#' @param p  the mixing weight (to be specified only if \code{type} is \code{"FB"}, \code{"VIB"}, or \code{"FBB"}). It must lie in (0, 1).
#' @param w  the normalized distance among component means of the FB and FBB distributions (to be specified only if \code{type = "FB"}, or \code{type = "FBB"}). It must lie in (0, 1).
#' @param k  the extent of the variance inflation (to be specified only if \code{type = "VIB"}). It must lie in (0, 1).
#' @param q0 the probability of augmentation in zero (to be specified only if \code{type} is \code{"Beta"}, \code{"FB"}, or \code{"VIB"}). It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one (to be specified only if \code{type} is \code{"Beta"}, \code{"FB"}, or \code{"VIB"}). It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param ... additional arguments of \code{\link{stat_function}}.
#'
#'
#' @examples
#' curve.density("Beta", mu=.5, phi=20)
#' curve.density("Beta", mu=.5, phi=20, q1 = .3)
#' curve.density("FB", mu=.5, phi=20, p=.4, w=.8)
#' curve.density("FB", mu=.5, phi=20, p=.4, w=.8, q0= .1)
#' curve.density("VIB", mu=.5, phi=20, p=.9, k=.8, col=3)
#' curve.density("VIB", mu=.5, phi=20, p=.9, k=.8, col=3, q0=.1, q1=.3)
#'
#' curve.density("Bin", size=10, mu=.7)
#' curve.density("BetaBin", size=10, mu=.7, phi=10)
#' curve.density("FBB", size=10, mu=.7, phi=10, p=.2,w=.7)
#'
#'
#' @references{
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020). Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309. doi:10.1177/1471082X18821213 \cr
#' \cr
#' Ferrari, S.L.P., and Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501\cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#'}
#'
#'
#' @import stats ggplot2
#'
#' @export
#'
#'
curve.density <- function (type = NULL, size = NULL, mu = NULL, theta = NULL,
                           phi = NULL, p = NULL, w = NULL, k = NULL, q0 = NULL, q1 = NULL,
                           ...)
{
  if (is.null(mu) | ((is.null(phi) & (is.null(theta))) & type !=
                     "Bin"))
    stop("Specify parameters according to the chosen distribution")
  if (type == "FB" & (is.null(p) | is.null(w)))
    stop("Specify parameters according to the chosen distribution")
  if (type == "VIB" & (is.null(p) | is.null(k)))
    stop("Specify parameters according to the chosen distribution")
  if (type == "Beta" & !is.null(p))
    stop("Specify parameters according to the chosen distribution")
  if (!(type %in% c("FB", "FBB")) & !is.null(w))
    stop("Specify parameters according to the chosen distribution")
  if (type != "VIB" & !is.null(k))
    stop("Specify parameters according to the chosen distribution")
  if (!(type %in% c("Beta", "FB", "VIB")) & !is.null(c(q0,
                                                       q1)))
    stop("Parameters q0 and q1 can be defined only for Beta, FB, or VIB distributions")
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)
  if ((type %in% c("FB", "Beta", "VIB"))) {
    x <- seq(0, 1, length.out = 1000)
    if (!is.null(size))
      stop("Parameter size must not be specified if type is ''Beta'', ''FB'', or ''VIB''.")
    if (type == "Beta") {
      fun <- dBeta(x, mu, phi, q0, q1)
    }
    else if (type == "FB") {
      fun <- dFB(x, mu = mu, phi = phi, p = p, w = w, q0 = q0,
                 q1 = q1)
    }
    else if (type == "VIB") {
      fun <- dVIB(x, mu = mu, phi = phi, p = p, k = k,
                  q0 = q0, q1 = q1)
    }
    else stop("Please specify the type of distribution correctly.")
    title.default <- paste0(type, " density function ", ifelse(q0 >
                                                                 0 & q1 == 0, "with zero augmentation", ifelse(q0 ==
                                                                                                                 0 & q1 > 0, "with one augmentation", ifelse(q0 >
                                                                                                                                                               0 & q1 > 0, "with zero and one augmentation", ""))))
    pp <- ggplot(data = data.frame(x = x[-c(1,length(x))], y = fun[-c(1,length(x))])) + geom_line(aes(x = x,
                                                                                                      y = y), alpha = 1) + xlab("x") + ylab("Density") +
      ggtitle(title.default) + theme_minimal()
    if (q0 > 0) {
      pp <- pp + geom_point(aes(x = x[1], y = q0), size = 2)+
        geom_segment(aes(x = x[1], y = x[1],xend = x[1], yend = q0), linetype="dashed")
    }
    if (q1 > 0) {
      pp <- pp + geom_point(aes(x = x[length(x)], y = q1),
                            size = 2)+
        geom_segment(aes(x = x[length(x)], y = x[1], xend = x[length(x)], yend = q1), linetype="dashed")
    }
  }
  else if ((type %in% c("FBB", "Bin", "BetaBin"))) {
    if (type == "FBB") {
      if (!is.null(theta) & is.null(phi))
        phi <- (1 - theta)/theta
      if (is.null(theta) & is.null(phi))
        stop("Please specify a value for theta or phi")
      title.default <- paste0(type, " mass function ")

      prob <- dFBB(0:size, size = size, mu = mu, phi = phi,
                   w = 0.6, p = 0.5)
    }
    else if (type == "BetaBin") {
      title.default <- paste0("Beta binomial mass function ")
      if (!is.null(theta) & is.null(phi))
        phi <- (1 - theta)/theta
      if (is.null(theta) & is.null(phi))
        stop("Please specify a value for theta or phi")
      prob <- dBetaBin(0:size, size = size, mu = mu, phi = phi)
    }
    else {
      prob <- dbinom(0:size, size = size, prob = mu)
      title.default <- paste0("Binomial mass function ")

    }
    y <- NULL
    pp <- ggplot(data = data.frame(y = 0:size, prob = prob)) +
      geom_point(aes(x = y, y = prob), size = 2) + geom_line(aes(x = y,
                                                                 y = prob), alpha = 0.3) + ylim(0, max(prob) + 0.05) +
      xlab("x") + ylab("Probability") + ggtitle(title.default)+ theme_minimal()
  }
  else stop("Error: please specify a valide type argument.")
  return(pp)
}


