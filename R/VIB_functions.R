#' The Variance-Inflated Beta Distribution
#'
#' @description Density function, distribution function, quantile function, and random generation
#' for the (augmented) variance-inflated beta distribution.
#'
#' @param x,q a vector of quantiles.
#' @param prob a vector of probabilities.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param k the extent of the variance inflation. It must lie in (0, 1).
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param log logical; if TRUE, densities are returned on log-scale.
#' @param log.prob logical; if TRUE, probabilities \code{prob} are given as log(prob).
#'
#' @return The function \code{dVIB} returns a vector with the same length as \code{x} containing the density values.
#' The function \code{pVIB} returns a vector with the same length as \code{q} containing the values of the distribution function.
#' The function \code{qVIB} returns a vector with the same length as \code{prob} containing the quantiles.
#' The function \code{rVIB} returns a vector of length \code{n} containing the generated random values.
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
#' @rdname VIB
#'
#' @export
dVIB <- Vectorize(function(x, mu, phi, p, k, q0 = NULL, q1 = NULL, log = FALSE){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  #if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  fun <- (1-q0-q1)*(exp(log(p) + dBeta(x, mu, phi*k, log = T)) +
                      exp(log(1-p) + dBeta(x, mu, phi, log = T)))
  fun[which(x==0)] <- q0
  fun[which(x==1)] <- q1

  if (any(x < 0 | x > 1)) fun <- 0

  if(log) fun <- log(fun)
  return(fun)
}, vectorize.args = c("x", "mu", "phi", "p", "k", "q0", "q1"))



#'
#' @examples
#' qVIB(prob = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5)
#' qVIB(prob = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q1 = .1)
#' qVIB(prob = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q0 = .2, q1 = .1)
#'
#' @import stats
#'
#' @rdname VIB
#'
#' @export

qVIB <- Vectorize(function(prob, mu, phi, p, k, q0 = NULL, q1 = NULL, log.prob = FALSE){
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)

  #if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  if(log.prob) prob <- exp(prob)

  fun <- uniroot(function(x) pVIB(x, mu = mu, phi = phi, p = p, k = k, q0=q0, q1=q1) - prob,
                 lower = 0, upper = 1, tol = 1e-9)$root

  return(fun)
}, vectorize.args = c("prob", "mu", "phi", "p", "k", "q0", "q1"))






#'
#' @examples
#' pVIB(q = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5)
#' pVIB(q = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q1 = .1)
#' pVIB(q = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q0 = .2, q1 = .1)
#'
#' @import stats
#'
#' @rdname VIB
#'
#' @export
pVIB <- Vectorize(function(q, mu, phi, p, k, q0 = NULL, q1 = NULL, log.prob = FALSE){
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)

  #if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  fun <- q0*(q >= 0) +
    (1-q0-q1)*(p*pBeta(q,mu,phi*k) +
                 (1-p)*pBeta(q,mu,phi)) +
    q1*(q >= 1)

  if(log.prob) fun <- log(fun)

  return(fun)
}, vectorize.args = c("q", "mu", "phi", "p", "k", "q0", "q1"))





#'
#' @examples
#' rVIB(n = 100, mu = .5, phi = 30, p = .3, k = .6)
#' rVIB(n = 100, mu = .5, phi = 30, p = .3, k = .6, q0 = .2, q1 = .1)
#'
#' @import stats
#'
#' @rdname VIB
#'
#' @export

rVIB <- function(n, mu, phi, p, k, q0 = NULL, q1 = NULL){
  if(is.null(q0)){
    q0 <- 0
  }
  if(is.null(q1)){
    q1 <- 0
  }

  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")

  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  n <- floor(n)

  l.mu <- length(mu)
  l.phi <- length(phi)
  l.p <- length(p)
  l.k <- length(k)
  l.q0 <- length(q0)
  l.q1 <- length(q1)

  L.max <- max(l.mu, l.phi, l.p, l.k, l.q0, l.q1)
  ## recycling check
  if (L.max %% l.mu != 0 | L.max %% l.phi != 0 |
      L.max %% l.p != 0 | L.max %% l.k != 0 |
      L.max %% l.q0 != 0 | L.max %% l.q1 != 0)
    warning("longer object length is not a multiple of shorter object length")

  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  p <- rep(p, length.out = n)
  k <- rep(k, length.out = n)
  q0 <- rep(q0, length.out = n)
  q1 <- rep(q1, length.out = n)


  out <- vector(mode="numeric", length = n)

  v.aug <- unlist(lapply(1:n, function(l) sample(c(0,1,2),1,replace = T, prob = c(1-q0[l]-q1[l],q0[l], q1[l]))))

  out[v.aug==1] <- 0
  out[v.aug==2] <- 1

  if(sum(v.aug == 0) > 0){
    mu_v0 <- mu[v.aug==0]
    k_v0 <- k[v.aug==0]
    phi_v0 <- phi[v.aug==0]

    for(i in 1:sum(v.aug == 0)){
      v.mixt <- rbinom(1, 1, prob=p)

      if(v.mixt == 1){
        out[v.aug==0][i] <- rBeta(1, mu_v0[i], k_v0[i]*phi_v0[i])
      } else out[v.aug==0][i] <- rBeta(1, mu_v0[i], phi_v0[i])
    }
  }

  return(out)
}
