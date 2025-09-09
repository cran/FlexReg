#' The Flexible Beta-Binomial Distribution
#'
#' @description Mass function, distribution function, quantile function, and random generation
#' for the flexible beta-binomial distribution.
#'
#' @param x,q a vector of quantiles.
#' @param prob a vector of probabilities.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion parameter. It must  lie in (0, 1).
#' @param phi the precision parameter, an alternative way to specify the overdispersion parameter \code{theta}. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among component means. It must lie in (0, 1).
#' @param log logical; if TRUE, probabilities are returned on log-scale.
#' @param log.prob logical; if TRUE, probabilities \code{prob} are given as log(prob).
#'
#' @return The function \code{dFBB} returns a vector with the same length as \code{x} containing the probability mass values.
#' The function \code{pFBB} returns a vector with the same length as \code{q} containing the values of the distribution function.
#' The function \code{qFBB} returns a vector with the same length as \code{prob} containing the quantiles.
#' The function \code{rFBB} returns a vector of length \code{n} containing the generated random values.
#'
#' @details The FBB distribution is a special mixture of two beta-binomial distributions with probability mass function
#' \deqn{f_{FBB}(x;\mu,\phi,p,w) = p BB(x;\lambda_1,\phi)+(1-p)BB(x;\lambda_2,\phi),}
#'  for \eqn{x \in \lbrace 0, 1, \dots, n \rbrace}, where \eqn{BB(x;\cdot,\cdot)} is the beta-binomial distribution with a mean-precision parameterization.
#'  Moreover, \eqn{\phi=(1-\theta)/\theta>0} is a precision parameter, \eqn{0<p<1} is the mixing weight, \eqn{0<\mu=p\lambda_1+(1-p)\lambda_2<1} is the overall mean,
#'   \eqn{0<w<1} is the  normalized distance between component means, and
#'  \eqn{\lambda_1=\mu+(1-p)w} and \eqn{\lambda_2=\mu-pw} are the scaled means of the first and second component of the mixture, respectively.
#'
#
#' @examples
#' dFBB(x = c(5,7,8), size=10, mu = .3, phi = 20, p = .5, w = .5)
#' dFBB(x = c(5,7,8), size=10, mu = .3, theta = 1/(20+1), p = .5, w = .5)
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @import stats
#'
#' @rdname FBB
#'
#' @export

dFBB <- Vectorize(function(x, size, mu, theta = NULL, phi = NULL, p, w, log = FALSE){
  #if (any(x < 0 | x > size) | !all.equal(x, as.integer(x))) stop("x has to be an integer between 0 and n")
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

  fun <-
    ifelse((x < 0 | x > size), 0,
           ifelse(!(x == as.integer(x)), 0,
                  p*dBetaBin(x=x, size=size, mu=lambda1, phi=phi) + (1-p)*dBetaBin(x=x, size=size, mu=lambda2, phi=phi)))

  if(log) fun <- log(fun)
  return(fun)
}, vectorize.args = c("x", "size", "mu", "theta", "phi", "p", "w"))




#' @examples
#' qFBB(prob = .5, size=10, mu = .3, phi = 20, p = .5, w = .5)
#' qFBB(prob = .5, size=10, mu = .3, theta = 1/(20+1), p = .5, w = .5)
#'
#' @import stats
#'
#' @rdname FBB
#'
#' @export
#'

qFBB <- Vectorize(function(prob, size, mu, theta = NULL, phi = NULL, p, w, log.prob = FALSE){

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

  if(log.prob) prob <- exp(prob)

  for(x in 0:size){
    if(pFBB(x, size=size, mu=mu, phi=phi, p=p, w=w) >= prob){
      fun <- x
      break()
    }
  }

  return(fun)
}, vectorize.args = c("prob", "size", "mu", "theta", "phi", "p", "w"))




#' @examples
#' pFBB(q = c(5,7,8), size=10, mu = .3, phi = 20, p = .5, w = .5)
#' pFBB(q = c(5,7,8), size=10, mu = .3, theta = 1/(20+1), p = .5, w = .5)
#'
#' @import stats
#'
#' @rdname FBB
#'
#' @export
pFBB <- Vectorize(function(q, size, mu, theta = NULL, phi = NULL, p, w, log.prob = FALSE){

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

  if(length(q) == 1){
    p_temp <- sum(dFBB(0:floor(q), size=size, mu=mu, phi=phi, p=p, w=w))
  } else {
    p_temp <- vector(mode="numeric", length = length(q))
    for(x in 1:length(q)){
      p_temp[x] <- sum(dFBB(0:floor(q[x]), size=size, mu=mu, phi=phi, p=p, w=w))
    }
  }

  fun <- p_temp
  if(log.prob) fun <- log(fun)


  return(fun)
}, vectorize.args = c("q", "size", "mu", "phi", "p", "w"))


#' @examples
#' rFBB(n = 100, size = 40, mu = .5, phi = 5, p = .3, w = .6)
#' rFBB(n = 100, size = 40, mu = .5, theta = 1/(5+1), p = .3, w = .6)
#'
#' @import stats
#'
#' @rdname FBB
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

  probs <- rFB(n, mu = mu, phi = phi, p=p, w=w)
  return(as.numeric(rbinom(n, size=size, prob = probs)))

}
