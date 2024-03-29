% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/densities.R
\name{dVIB}
\alias{dVIB}
\title{Probability density function of the variance-inflated beta distribution}
\usage{
dVIB(x, mu, phi, p, k, q0 = NULL, q1 = NULL)
}
\arguments{
\item{x}{a vector of quantiles.}

\item{mu}{the mean parameter. It must lie in (0, 1).}

\item{phi}{the precision parameter. It must be a real positive value.}

\item{p}{the mixing weight. It must lie in (0, 1).}

\item{k}{the extent of the variance inflation. It must lie in (0, 1).}

\item{q0}{the probability of augmentation in zero. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).}

\item{q1}{the probability of augmentation in one. It must lie in (0, 1). In case of no augmentation, it is \code{NULL} (default).}
}
\value{
A vector with the same length as \code{x}.
}
\description{
The function computes the probability density function of the variance-inflated beta distribution.
It can also  compute the probability density function of the augmented variance-inflated beta distribution by assigning positive probabilities to zero and/or one values and a (continuous) variance-inflated beta density to the interval (0,1).
}
\details{
The VIB distribution is a special mixture of two beta distributions with probability density function
\deqn{f_{VIB}(x;\mu,\phi,p,k)=p f_B(x;\mu,\phi k)+(1-p)f_B(x;\mu,\phi),}
 for \eqn{0<x<1}, where \eqn{f_B(x;\cdot,\cdot)} is the beta density with a mean-precision parameterization.
 Moreover, \eqn{0<p<1} is the mixing weight, \eqn{0<\mu<1} represents the overall (as well as mixture component)
 mean, \eqn{\phi>0} is a precision parameter, and \eqn{0<k<1} determines the extent of the variance inflation.
The augmented VIB distribution has density
\itemize{
\item \eqn{q_0}, if \eqn{x=0}
\item \eqn{q_1},  if \eqn{x=1}
\item \eqn{(1-q_0-q_1)f_{VIB}(x;\mu,\phi,p,k)}, if \eqn{0<x<1}
}
where \eqn{0<q_0<1} identifies the augmentation in zero, \eqn{0<q_1<1} identifies the augmentation in one,
and \eqn{q_0+q_1<1}.
}
\examples{
dVIB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5)
dVIB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q1 = .1)
dVIB(x = c(.5,.7,.8), mu = .3, phi = 20, p = .5, k= .5, q0 = .2, q1 = .1)

}
\references{
{
Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020). Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309.
doi:10.1177/1471082X18821213
}
}
