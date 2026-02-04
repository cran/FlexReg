#' The Flexible Beta Distribution
#'
#' @description  Density function, distribution function, quantile function, and random generation
#' for the (augmented) flexible beta distribution.
#'
#' @param x,q a vector of quantiles.
#' @param prob a vector of probabilities.
#' @param n the number of values to generate. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu the mean parameter. It must lie in (0, 1).
#' @param phi the precision parameter. It must be a real positive value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among component means. It must lie in (0, 1).
#' @param q0 the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param q1 the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).
#' @param log logical; if TRUE, densities are returned on log-scale.
#' @param log.prob logical; if TRUE, probabilities \code{prob} are given as log(prob).
#'
#' @return The function \code{dFB} returns a vector with the same length as \code{x} containing the density values.
#' The function \code{pFB} returns a vector with the same length as \code{q} containing the values of the distribution function.
#' The function \code{qFB} returns a vector with the same length as \code{prob} containing the quantiles.
#' The function \code{rFB} returns a vector of length \code{n} containing the generated random values.
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
#'
#' @references {
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079
#' }
#'
#' @import stats
#'
#' @rdname FB
#'
#' @export

dFB <- Vectorize(function(x, mu, phi, p, w, q0 = NULL, q1 = NULL, log = FALSE){
  q0 <- ifelse(is.null(q0),0,q0)
  q1 <- ifelse(is.null(q1),0,q1)
  #if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde

  fun <- (1-q0-q1)*(exp(log(p) + dBeta(x, lambda1, phi, log = T)) +
    exp(log(1-p) + dBeta(x, lambda2, phi, log = T)))
  fun[which(x==0)] <- q0
  fun[which(x==1)] <- q1

  if (any(x < 0 | x > 1)) fun <- 0

  if(log) fun <- log(fun)
  return(fun)
}, vectorize.args = c("x", "mu", "phi", "p", "w", "q0", "q1"))







#'
#' @examples
#' qFB(prob = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5)
#' qFB(prob = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5, q0 = .2)
#' qFB(prob = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5, q0 = .2, q1 = .1)
#'
#' @import stats
#'
#' @rdname FB
#'
#' @export

qFB <- Vectorize(function(prob, mu, phi, p, w, q0 = NULL, q1 = NULL, log.prob = FALSE){
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)

  #if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde

  if(log.prob) prob <- exp(prob)

  fun <- uniroot(function(x) pFB(x, mu = mu, phi = phi, p = p, w = w, q0=q0, q1=q1) - prob,
                 lower = 0, upper = 1, tol = 1e-9)$root

  return(fun)
}, vectorize.args = c("prob", "mu", "phi", "p", "w", "q0", "q1"))








#' @examples
#' pFB(q = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5)
#' pFB(q = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5, q0 = .2)
#' pFB(q = c(.5,.7,.8), mu = .3, phi = 20, p = .5, w = .5, q0 = .2, q1 = .1)
#'
#' @import stats
#'
#' @rdname FB
#'
#' @export

pFB <- Vectorize(function(q, mu, phi, p, w, q0 = NULL, q1 = NULL, log.prob = FALSE){
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)

  #if (any(x < 0 | x > 1)) stop("x has to be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")
  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde

  fun <- q0*(q >= 0) +
    (1-q0-q1)*(p*pBeta(q,lambda1,phi) +
                 (1-p)*pBeta(q,lambda2,phi)) +
    q1*(q >= 1)

  if(log.prob) fun <- log(fun)

  return(fun)
}, vectorize.args = c("q", "mu", "phi", "p", "w", "q0", "q1"))


#' @examples
#' rFB(n = 100, mu = .5, phi = 30,p = .3, w = .6)
#' rFB(n = 100, mu = .5, phi = 30,p = .3, w = .6, q0 = .2, q1 = .1)
#'
#' @import stats
#'
#' @rdname FB
#'
#' @export

rFB <- function(n, mu, phi, p, w, q0 = NULL, q1 = NULL){
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
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  if (any(q0 < 0 | q0 > 1)) stop("Parameter q0 has to be between 0 and 1")
  if (any(q1 < 0 | q1 > 1)) stop("Parameter q1 has to be between 0 and 1")

  if (any(q0+q1>=1)) stop("The sum of q0 and q1 must be less than 1")

  n <- floor(n)

  l.mu <- length(mu)
  l.phi <- length(phi)
  l.p <- length(p)
  l.w <- length(w)
  l.q0 <- length(q0)
  l.q1 <- length(q1)

  L.max <- max(l.mu, l.phi, l.p, l.w, l.q0, l.q1)
  ## recycling check
  if (L.max %% l.mu != 0 | L.max %% l.phi != 0 |
      L.max %% l.p != 0 | L.max %% l.w != 0 |
      L.max %% l.q0 != 0 | L.max %% l.q1 != 0)
    warning("longer object length is not a multiple of shorter object length")

  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  p <- rep(p, length.out = n)
  w <- rep(w, length.out = n)
  q0 <- rep(q0, length.out = n)
  q1 <- rep(q1, length.out = n)


  wtilde <- w*pmin(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde

  out <- vector(mode="numeric", length = n)

  v.aug <- unlist(lapply(1:n, function(l) sample(c(0,1,2),1,replace = T, prob = c(1-q0[l]-q1[l],q0[l], q1[l]))))

  out[v.aug==1] <- 0
  out[v.aug==2] <- 1

  if(sum(v.aug == 0) > 0){
    lambda1_v0 <- lambda1[v.aug==0]
    lambda2_v0 <- lambda2[v.aug==0]
    phi_v0 <- phi[v.aug==0]

    for(i in 1:sum(v.aug == 0)){
      v.mixt <- rbinom(1, 1, prob=p)

      if(v.mixt == 1){
        out[v.aug==0][i] <- rBeta(1, lambda1_v0[i], phi_v0[i])
      } else out[v.aug==0][i] <- rBeta(1, lambda2_v0[i], phi_v0[i])
    }
  }

  return(out)
}
