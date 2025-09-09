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
#' @param ... additional arguments of \code{\link[ggplot2]{stat_function}}.
#'
#'
#' @examples
#' \dontrun{
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
#' }
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
  if (!(type %in% c("Beta", "FB", "VIB")) & !is.null(c(q0,q1)))
    stop("Parameters q0 and q1 can be defined only for Beta, FB, or VIB distributions")
  q0 <- ifelse(is.null(q0), 0, q0)
  q1 <- ifelse(is.null(q1), 0, q1)
  if ((type %in% c("FB", "Beta", "VIB"))) {
    x <- seq(0, 1, length.out = 1000)
    if (!is.null(size))
      stop("Parameter size must not be specified if type is ''Beta'', ''FB'', or ''VIB''.")
    if (type == "Beta") {
      fun <- dBeta(x, mu, phi, q0, q1)
    } else if (type == "FB") {
      fun <- dFB(x, mu = mu, phi = phi, p = p, w = w, q0 = q0,
                 q1 = q1)
    } else if (type == "VIB") {
      fun <- dVIB(x, mu = mu, phi = phi, p = p, k = k,
                  q0 = q0, q1 = q1)
    } else stop("Please specify the type of distribution correctly.")
    title.default <- paste0(type, " density function ",
                            ifelse(q0 > 0 & q1 == 0, "with zero augmentation",
                                   ifelse(q0 == 0 & q1 > 0, "with one augmentation",
                                          ifelse(q0 > 0 & q1 > 0, "with zero and one augmentation", ""))))

    pp <- ggplot(data = data.frame(x = x[-c(1,length(x))], y = fun[-c(1,length(x))])) + geom_line(aes(x = x,
                                                                                                      y = y), alpha = 1) + xlab("x") + ylab("Density") +
      ggtitle(title.default) + theme_minimal()
    if (q0 > 0) {
      pp <- pp + annotate("segment", x = x[1], y = x[1],xend = x[1], yend = q0, linetype="dashed") +
        annotate("point", x = x[1], y = q0, size = 2)
    }
    if (q1 > 0) {
      pp <- pp + annotate("segment", x = x[length(x)], y = x[1],xend = x[length(x)], yend = q1, linetype="dashed") +
        annotate("point", x = x[length(x)], y = q1, size = 2)

    }
  } else if ((type %in% c("FBB", "Bin", "BetaBin"))) {
    if (type == "FBB") {
      if (!is.null(theta) & is.null(phi))
        phi <- (1 - theta)/theta
      if (is.null(theta) & is.null(phi))
        stop("Please specify a value for theta or phi")
      title.default <- paste0(type, " mass function ")

      prob <- dFBB(0:size, size = size, mu = mu, phi = phi,
                   w = w, p = p)
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
      geom_point(aes(x = y, y = prob), size = 2) +
      geom_line(aes(x = y, y = prob), alpha = 0.3) +
      ylim(0, max(prob) + 0.05) +
      xlab("x") + ylab("Probability") +
      ggtitle(title.default)+ theme_minimal()
  }
  else stop("Error: please specify a valide type argument.")
  return(pp)
}


