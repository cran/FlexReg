#' The Beta-Binomial Distribution
#'
#' @description Mass function, distribution function, quantile function, and random generation
#' for the beta-binomial distribution.
#'
#' @param x,q a vector of quantiles.
#' @param prob a vector of probabilities.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param size the total number of trials.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param theta the overdispersion parameter. It must  lie in (0, 1).
#' @param phi the precision parameter, an alternative way to specify the overdispersion parameter \code{theta}. It must be a real positive value.
#' @param log logical; if TRUE, probabilities are returned on log-scale.
#' @param log.prob logical; if TRUE, probabilities \code{prob} are given as log(prob).
#'
#' @return The function \code{dBetaBin} returns a vector with the same length as \code{x} containing the probability mass values.
#' The function \code{pBetaBin} returns a vector with the same length as \code{q} containing the values of the distribution function.
#' The function \code{qBetaBin} returns a vector with the same length as \code{prob} containing the quantiles.
#' The function \code{rBetaBin} returns a vector of length \code{n} containing the generated random values.
#'
#' @details The beta-binomial distribution has probability mass function
#' \deqn{f_{BB}(x;\mu,\phi)={n\choose x} \frac{\Gamma{(\phi)}}{\Gamma{(\mu\phi)}\Gamma{((1-\mu)\phi)}} \frac{\Gamma{(\mu\phi+x)}\Gamma{((1-\mu)\phi + n - x)}}{\Gamma{(\phi + n)}},}
#' for \eqn{x \in \lbrace 0, 1, \dots, n \rbrace}, where \eqn{0<\mu<1} identifies the mean and \eqn{\phi=(1-\theta)/\theta >0} is the precision parameter.
#'
#' @examples
#' dBetaBin(x = 5, size = 10, mu = .3, phi = 10)
#' dBetaBin(x = 5, size = 10, mu = .3, theta = 1/(10+1))
#'
#' @references{
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005
#' }
#'
#' @import stats
#'
#' @rdname BetaBin
#'
#' @export
#'

dBetaBin <- Vectorize(function(x, size, mu, theta = NULL, phi = NULL, log = FALSE){
  #if (any(x < 0 | x > size) | !all.equal(x, as.integer(x))) stop("x has to be an integer between 0 and n")
  #if (!all.equal(x, as.integer(x))) warning("x has to be an integer between 0 and n")

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

  fun <-
    ifelse((x < 0 | x > size), 0,
           ifelse(!(x == as.integer(x)), 0, exp(lgamma(size+1)-lgamma(x+1)-lgamma(size-x+1) +
                                                  lgamma(phi)+lgamma((phi*mu)+x) + lgamma((phi*(1-mu))+size-x) -
                                                  lgamma(phi+size) - lgamma(phi*mu) - lgamma(phi*(1-mu)))))

  if(log) fun <- log(fun)
  return(fun)
}, vectorize.args = c("x", "size", "mu", "theta", "phi"))


#' @examples
#' qBetaBin(prob = .5, size = 10, mu = .3, phi = 10)
#' qBetaBin(prob = .5, size = 10, mu = .3, theta = 1/(10+1))
#'
#' @import stats
#'
#' @rdname BetaBin
#'
#' @export
#'

qBetaBin <- Vectorize(function(prob, size, mu, theta = NULL, phi = NULL, log.prob = FALSE){

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
    if(pBetaBin(x, size=size, mu=mu, phi=phi) >= prob){
      fun <- x
      break()
    }
  }

  return(fun)
}, vectorize.args = c("prob", "size", "mu", "theta", "phi"))


#' @examples
#' pBetaBin(q = 5, size = 10, mu = .3, phi = 10)
#' pBetaBin(q = 5, size = 10, mu = .3, theta = 1/(10+1))
#'
#' @import stats
#'
#' @rdname BetaBin
#'
#' @export
#'

pBetaBin <- Vectorize(function(q, size, mu, theta = NULL, phi = NULL, log.prob = FALSE){

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

  if(length(q) == 1){
    p_temp <- sum(dBetaBin(0:floor(q), size=size, mu=mu, phi=phi))
  } else {
    p_temp <- vector(mode="numeric", length = length(q))
    for(x in 1:length(q)){
      p_temp[x] <- sum(dBetaBin(0:floor(q[x]), size=size, mu=mu, phi=phi))
    }
  }

  fun <- p_temp
  if(log.prob) fun <- log(fun)

  return(fun)
}, vectorize.args = c("q", "size", "mu", "phi"))



#' @examples
#' rBetaBin(n = 100, size = 40, mu = .5, theta = .4)
#' rBetaBin(n = 100, size = 40, mu = .5, phi = 1.5)
#'
#' @import stats
#'
#' @rdname BetaBin
#'
#' @export

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
  return(as.numeric(rbinom(n, size=size, prob = probs)))
}
