#' Random generator from the beta distribution
#'
#' @description The function  generates random values from the beta distribution with a mean-precision parameterization, or from the augmented beta distribution.
#'
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#'
#' @return A vector of length \code{n}.
#'
#' @examples
#' rBeta(n = 100, mu = .5, phi = 30)
#' rBeta(n = 100, mu = .5, phi = 30, q0 = .2, q1 = .1)

#'
#' @references{
#' Ferrari, S.L.P.,  Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501
#' }
#'
#' @import stats
#'
#' @export
#'

rBeta <- function(n, mu, phi, q0 = NULL, q1 = NULL){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")

  if (any(q0+q1>1)) stop("The sum of q0 and q1 must be less than 1")

  n <- floor(n)

  v <- sample(c(0,1,2), n, replace=T, prob=c(1-q0-q1,q0,q1))
  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi

  x <- vector(mode="numeric", length = n)
  x[v==0] <- rbeta(length(which(v==0)),alpha1,alpha2)
  x[v==1] <- 0
  x[v==2] <- 1
  return(x)
}

#' Random generator from the flexible beta distribution
#'
#' @description The function generates random values from the flexible beta distribution, or from the augmented flexible beta distribution.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among clusters. It must lie in (0, 1).
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#'
#' @return A vector of length  \code{n}.
#'
#' @examples
#' rFB(n = 100, mu = .5, phi = 30,p = .3, w = .6)
#' rFB(n = 100, mu = .5, phi = 30,p = .3, w = .6, q0 = .2, q1 = .1)
#' @references {
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079
#' }
#'
#' @import stats
#'
#' @export

rFB <- function(n, mu, phi, p, w,q0 = NULL, q1 = NULL){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")

  if (any(q0+q1>1)) stop("The sum of q0 and q1 must be less than 1")

  n <- floor(n)
  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde
  x <- vector(mode="numeric", length = n)
  v.aug <- sample(c(0,1,2), n, replace=T, prob=c(1-q0-q1,q0,q1))
  x[v.aug==1] <- 0
  x[v.aug==2] <- 1
  v <- rbinom(length(which(v.aug==0)),1,prob=p)
  x[v.aug==0][v==1] <- rBeta(length(which(v==1)),lambda1,phi)
  x[v.aug==0][v==0] <- rBeta(length(which(v==0)),lambda2,phi)
  return(x)
}

#' Random generation from the variance-inflated beta distribution
#'
#' @description The function generates random values from the variance-inflated beta distribution, or from the augmented variance-inflated beta distribution.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param k the extent of the variance inflation. It must lie in (0, 1).
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#'
#' @return A vector of length  \code{n}.
#'
#' @examples
#' rVIB(n = 100, mu = .5, phi = 30, p = .3, k = .6)
#' rVIB(n = 100, mu = .5, phi = 30, p = .3, k = .6, q0 = .2, q1 = .1)
#'
#' @references{
#'  Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020). Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309. doi:10.1177/1471082X18821213
#' }
#'
#' @import stats
#'
#' @export

rVIB <- function(n, mu, phi, p, k, q0 = NULL, q1 = NULL){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")

  if (any(q0+q1>1)) stop("The sum of q0 and q1 must be less than 1")

  n <- floor(n)
  x <- vector(mode="numeric", length = n)

  v.aug <- sample(c(0,1,2), n, replace=T, prob=c(1-q0-q1,q0,q1))
  x[v.aug==1] <- 0
  x[v.aug==2] <- 1
  v <- rbinom(length(which(v.aug==0)),1,prob=p)
  x[v.aug==0][v==1] <- rBeta(length(which(v==1)),mu,phi*k)
  x[v.aug==0][v==0] <- rBeta(length(which(v==0)),mu,phi)

  return(x)
}


#' Random generator from the beta-binomial distribution
#'
#' The function  generates random values from the beta-binomial distribution.
#' @param n the number of values to generate. If \code{length(n)} > 1, the length is taken to be the number required.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion parameter. It must  lie in (0, 1).
#' @param phi the precision parameter, an alternative way to specify the overdispersion parameter \code{theta}.  It must be a real positive value.
#'
#' @return A vector of length \code{n}.
#'
#' @examples
#' rBetaBin(n = 100, size = 40, mu = .5, theta = .4)
#' rBetaBin(n = 100, size = 40, mu = .5, phi = 1.5)
#'
#' @references{
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @import stats
#'
#' @export
#'

rBetaBin <- function(n, size=NULL, mu=NULL, theta=NULL, phi=NULL){
  if (length(n)>1) n <- length(n)
  if (any(is.null(mu) | mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(is.null(size) | size < 0 | size != as.integer(size))) stop("size must be a non-negative integer")

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

  probs <- rBeta(n, mu = mu, phi = phi)
  return(rbinom(n, size=size, prob = probs))
}

#' Random generator from the flexible beta-binomial distribution
#'
#' The function generates random values from the flexible beta-binomial distribution.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion parameter. It must lie in (0, 1).
#' @param phi the precision parameter, an alternative way to specify the overdispersion parameter \code{theta}. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among clusters. It must lie in (0, 1).
#'
#' @return A vector of length  \code{n}.
#'
#' @examples
#' rFBB(n = 100, size = 40, mu = .5, theta = .4, p = .3, w = .6)
#' rFBB(n = 100, size = 40, mu = .5, phi = 1.5, p = .3, w = .6)
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @import stats
#'
#' @export

rFBB <- function(n, size=NULL, mu, theta=NULL, phi=NULL, p, w){
  if (length(n)>1) n <- length(n)
  if (any(is.null(mu) | mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(is.null(p) | p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(is.null(w) | w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(is.null(size) | size < 0 | size != as.integer(size))) stop("size must be a non-negative integer")

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

  n <- floor(n)
  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde
  v <- rbinom(n,1,prob=p)
  x <- vector(mode="numeric", length = n)
  x[v==1] <- rBetaBin(n=length(which(v==1)), size=size, mu = lambda1, phi=phi)
  x[v==0] <- rBetaBin(n=length(which(v==0)), size=size, mu = lambda2, phi=phi)
  return(x)
}
