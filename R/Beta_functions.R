#' The Mean-Precision Parameterized Beta Distribution
#'
#' @description Density function, distribution function, quantile function, and random generation
#' for the (augmented) beta distribution with the mean-precision parameterization.
#'
#' @param x,q a vector of quantiles.
#' @param prob a vector of probabilities.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param log logical; if TRUE, densities are returned on log-scale.
#' @param log.prob logical; if TRUE, probabilities \code{prob} are given as log(prob).
#'
#' @return The function \code{dBeta} returns a vector with the same length as \code{x} containing the density values.
#' The function \code{pBeta} returns a vector with the same length as \code{q} containing the values of the distribution function.
#' The function \code{qBeta} returns a vector with the same length as \code{prob} containing the quantiles.
#' The function \code{rBeta} returns a vector of length \code{n} containing the generated random values.
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
#' @rdname Beta
#'
#'
dBeta <- Vectorize(function(x, mu, phi, q0 = NULL, q1 = NULL, log = FALSE){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  #if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi

  fun <- (1-q0-q1)*dbeta(x, shape1=alpha1, shape2=alpha2)
  fun[which(x==0)] <- q0
  fun[which(x==1)] <- q1

  if (any(x < 0 | x > 1)) fun <- 0

  if(log) fun <- log(fun)
  return(fun)
}, vectorize.args = c("x", "mu", "phi", "q0", "q1"))



#'
#' @examples
#' qBeta(prob = c(.5,.7,.8), mu = .3, phi = 20)
#' qBeta(prob = c(.5,.7,.8), mu = .3, phi = 20, q0 = .2)
#' qBeta(prob = c(.5,.7,.8), mu = .3, phi = 20, q0 = .2, q1= .1)
#'
#' @import stats
#'
#' @export
#'
#' @rdname Beta
#'


qBeta <- Vectorize(function(prob, mu, phi, q0 = NULL, q1 = NULL, log.prob = FALSE){
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)

  #if (any(q < 0 | q > 1)) stop("q has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi

  if(log.prob) prob <- exp(prob)

  fun <- uniroot(function(x) pBeta(x, mu = mu, phi = phi, q0=q0, q1=q1) - prob,
                 lower = 0, upper = 1, tol = 1e-9)$root

  return(fun)
}, vectorize.args = c("prob", "mu", "phi", "q0", "q1"))



#'
#' @examples
#' pBeta(q = c(.5,.7,.8), mu = .3, phi = 20)
#' pBeta(q = c(.5,.7,.8), mu = .3, phi = 20, q0 = .2)
#' pBeta(q = c(.5,.7,.8), mu = .3, phi = 20, q0 = .2, q1= .1)
#'
#' @import stats
#'
#' @rdname Beta
#'
#' @export
#'

pBeta <- Vectorize(function(q, mu, phi, q0 = NULL, q1 = NULL, log.prob = FALSE){
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)

  #if (any(q < 0 | q > 1)) stop("q has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi

  fun <- q0*(q >= 0) + (1-q0-q1)*pbeta(q, alpha1, alpha2) + q1*(q>=1)

  if(log.prob) fun <- log(fun)

  return(fun)
}, vectorize.args = c("q", "mu", "phi", "q0", "q1"))





#'
#' @examples
#' rBeta(n = 100, mu = .5, phi = 30)
#' rBeta(n = 100, mu = .5, phi = 30, q0 = .2, q1 = .1)
#'
#' @import stats
#'
#' @rdname Beta
#'
#' @export
#'

rBeta <- function(n, mu, phi, q0 = NULL, q1 = NULL){

  if(is.null(q0)){
    q0 <- 0
  }
  if(is.null(q1)){
    q1 <- 0
  }

  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")

  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  n <- floor(n)

  l.mu <- length(mu)
  l.phi <- length(phi)
  l.q0 <- length(q0)
  l.q1 <- length(q1)

  L.max <- max(l.mu, l.phi, l.q0, l.q1)
  ## recycling check
  if (L.max %% l.mu != 0 | L.max %% l.phi != 0 |
      L.max %% l.q0 != 0 | L.max %% l.q1 != 0)
    warning("longer object length is not a multiple of shorter object length")

  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  q0 <- rep(q0, length.out = n)
  q1 <- rep(q1, length.out = n)

  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi

  v.aug <- unlist(lapply(1:n, function(l) sample(c(0,1,2),1,replace = T, prob = c(1-q0[l]-q1[l],q0[l], q1[l]))))

  out <- vector(mode="numeric", length = n)

  if(sum(v.aug == 0)>0){
    out[v.aug==0] <- rbeta(sum(v.aug==0),
                           alpha1[which(v.aug==0)],
                           alpha2[which(v.aug==0)])
  }

  out[v.aug==1] <- 0
  out[v.aug==2] <- 1
  return(out)
}

