#' @title Predict Method for \code{`flexreg`} Objects
#'
#' @description Method that computes various types of predictions from objects of class \code{`flexreg`}.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param newdata an optional  \code{data.frame} containing variables with which to predict. If omitted, the fitted values are used.
#' @param n.new an optional vector containing the total number of trials with which to predict. It must be specified if \code{newdata} is not \code{NULL} and the
#' \code{\link{flexreg}} object is the result of the \code{\link{flexreg_binom}} function (i.e., the fitted model is binomial, beta-binomial, or flexible beta-binomial). The vector must have the same length as \code{nrow(newdata)}.
#' @param cluster a logical (with default \code{FALSE}).  The option \code{cluster = TRUE} is available only for \code{"FB"} and \code{"FBB"} models and allows to compute some component-specific predictions (see Details).
#' @param type a character indicating the type of prediction. Available options are: \code{"response"}, returning the marginal fitted mean of the response/relative response;
#' \code{"link"}, returning the linear predictor of the mean model;
#' \code{"precision"}, returning the fitted precision parameter;
#' \code{"overdispersion"}, returning the fitted overdispersion parameter;
#' \code{"variance"}, returning the fitted variance of the response.
#' @param estimate a character indicating the type of estimate. Available options are \code{"mean"} (default), \code{"median"}, and \code{"quantile"}.
#' @param q if \code{estimate = "quantile"}, a numeric value of probability in (0, 1).
#' @param ... additional arguments. Currently not used.
#'
#' @details The \code{\link{predict}} method computes various types of predictions from objects of class \code{`flexreg`}.
#' If \code{type = "response"}, the function returns the marginal mean, i.e., \eqn{\mu}.
#' In case of models for continuous bounded responses with augmentation, the function returns also the overall mean
#' \eqn{q_1+(1-q_0-q_1)\mu} and the probabilities of augmentation \eqn{q_0} and/or \eqn{q_1}.
#' If \code{type = "variance"}, the function returns \eqn{Var(Y|0<Y<1)} in case of no augmentation and
#' \eqn{(1-q_0-q_1)Var(Y|0<Y<1)+q_1^2+(1-q_0-q_1)\mu^2-(q_1+(1-q_0-q_1)\mu)^2} in case of augmentation.
#' If \code{cluster = TRUE}, for FB and FBB models, the function returns the cluster means (\eqn{\lambda_1} and \eqn{\lambda_2}) when \code{type = "response"} and the cluster variances when \code{type = "variance"}.
#'
#'The option \code{type = "overdispersion"}  is available only for beta-binomial and flexible beta-binomial models and returns the fitted overdispersion.

#'@return The function returns a \code{data.frame} of different dimensions depending on the type of prediction.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data=Reading, type="FB")
#' predict(FB, type="response", cluster=TRUE)
#' }
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' }
#'
#' @import  stats
#'
#' @export
#'

predict.flexreg <- function(object, newdata = NULL, n.new = NULL, cluster = FALSE,
                            type = "response",
                            estimate = "mean", q = NULL, ...){
  model <- object
  model.type <- model$type
  model.class <- class(model)
  formula <- model$formula
  cluster.var <- cluster

  if("flexreg_binom" %in% model.class){
    if(is.null(newdata))  {
        n <- model$n
    } else {
      n <- n.new
      if(is.null(n)) stop("Please specify n.new containing the values of n with which to predict.")
      if(length(n) != nrow(newdata)) stop("The vector in n must be of the same length as the number of row in newdata.")
    }
  } else {
    n <- NULL
  }

  #start checks
  if((model.type %in% c("Beta", "BetaBin", "Bin")) & cluster==T){
    cluster <- F
    warning("Beta, Binomial, and Beta-Binomial models are not mixture and thus clusters are not available.")
  }

  if( (model.type == "VIB") & cluster==T){
    cluster <- F
    if(type == "response")  warning("For VIB models, the component-specific means are omitted.")
  }

  if(!(estimate %in% c("mean", "median", "quantile")))
    stop("Argument `estimate` must be set equal to `mean`, `median` or `quantile`.")

  if((estimate == "quantile") & is.null(q))
    stop("Please specify the order of the quantile.")

  if(!(type %in% c("response", "link", "precision", "variance", "overdispersion")))
    stop("Argument `type` must be set equal to `response`, `link`, `precision`, `variance`, or `overdispersion`.")

  if((type == "overdispersion") & (model.type %in% c("Bin", "Beta", "FB", "VIB")))
    stop("predicted ovedispersion is available only for BetaBin and FBB models.")

  if((type == "precision") & ("flexreg_binom" %in% class(model)))
    stop("predicted precision is available only for Beta, FB, and VIB models.")


  if( cluster == T & !(type %in% c("response", "variance"))){
    cluster = F
    warning("Clusters are printed only when type is response or variance.")
  }
  #end checks

  if(!is.null(newdata)) newdata <- newdata.adjust(newdata, formula) #adjust the newdata

  #produce outcome depending on the chosen type
  if(type == "response") pred.chain <- predict_response(model, newdata, cluster, n)
  if(type == "link") pred.chain <- predict_link(model, newdata)
  if(type == "precision") pred.chain <- predict_precision(model, newdata)
  if(type == "variance") pred.chain <- predict_variance(model, newdata, cluster, cluster.var, model.type, model.class, n)
  if(type == "overdispersion") pred.chain <- predict_over(model, newdata)


  pred.chain <- pred.chain[unlist(lapply(pred.chain, is.matrix))]
  predicted <- lapply(pred.chain, function(x)
      apply(x, 2, eval(str2expression(estimate)), probs = q))
  predicted <- do.call(data.frame, predicted)

  return(as.data.frame(predicted))
}


#' @title Residuals Method for flexreg Objects
#'
#' @description Method that computes various types of residuals from objects of class \code{`flexreg`}. If the model type is  \code{FB} or \code{FBB} and \code{cluster = TRUE}, the method returns also residuals with respect to cluster means.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param type a character indicating type of residuals (\code{"raw"} or \code{"standardized"}).
#' @param cluster logical. If the model is \code{"FB"} without augmentation or \code{"FBB"}, \code{cluster = TRUE} returns the cluster means. By default \code{cluster = FALSE}.
#' @param estimate a character indicating the type of estimate: \code{"mean"} (default), \code{"median"}, or \code{"quantile"}.
#' @param q if \code{estimate = "quantile"}, a numeric value of probability in (0, 1).
#' @param ... additional arguments. Currently not used.
#'
#'@return The method returns an array with as many rows as the number of observations in the sample. If \code{cluster = FALSE}, the array has only one column containing either the raw or standardized residuals.
#' If \code{cluster = TRUE}, the array has four columns: the first column contains the raw or standardized residuals, the second and third columns contain the cluster residuals,
#' and the fourth column contains the classification labels (see Details).
#'
#' @details The \code{residuals} method  computes raw and standardized residuals from objects of class \code{`flexreg`}.
#' Raw residuals are defined as \eqn{r=y-\hat{\mu}} for bounded continuous  responses or as \eqn{r= y/n-\hat{\mu}} for bounded discrete responses.
#' Values \eqn{y} and \eqn{y/n} are  the observed
#' responses which are specified on the left-hand side of \code{formula} in the
#' \code{\link{flexreg}} and \code{\link{flexreg_binom}} functions, respectively.
#' Moreover,  \eqn{\hat{\mu}}  is the predicted value, the result of
#'   the \code{\link{predict}} function with \code{type = "response"}.
#' Standardized residuals are defined as \eqn{\frac{r}{\sqrt{\widehat{Var}(y)}}} where
#'  \eqn{\widehat{Var}(y)}
#' is the variance of the response evaluated at the posterior means
#' --by default, otherwise evaluated at the posterior quantiles of order \code{q}-- of the parameters.
#' If the model is \code{"FB"} or \code{"FBB"}, \code{type = "raw"}, and \code{cluster = TRUE}, the cluster raw residuals are computed as
#' the difference between the observed response/relative response  and the cluster means, i.e.,
#' \eqn{\hat{\lambda}_{1}} and \eqn{\hat{\lambda}_{2}}.
#' If the model is \code{"FB"} or \code{"FBB"}, \code{type = "standardized"} and \code{cluster = TRUE}, the cluster standardized residuals are computed as the
#' cluster raw residuals divided by the square root of the cluster variances.
#' Cluster residuals, either raw or standardized, can be used for classification purpose. Indeed, with \code{cluster = TRUE} the \code{residuals} method returns also a column named
#' \code{"label"} assigning values 1 or 2 to observations depending on whether they are classified in cluster 1 (if the corresponding cluster residual is smaller) or in cluster 2.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data=Reading, type="FB")
#' residuals(FB, type="raw", cluster=TRUE)
#' }
#'
#' @references {
#' Ascari, R., Migliorati, S. (2021). A new regression model for overdispersed binomial data accounting for outliers and an excess of zeros. Statistics in Medicine, \bold{40}(17), 3895--3914. doi:10.1002/sim.9005 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S. (2020). A new mixed-effects mixture model for constrained longitudinal data. Statistics in Medicine, \bold{39}(2), 129--145. doi:10.1002/sim.8406 \cr
#' \cr
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' }
#'
#' @import stats
#'
#' @export
#'
#'

residuals.flexreg <- function(object, type = "raw",
                              cluster=FALSE, estimate="mean", q=NULL, ...){
  model <- object
  posterior <- model$model
  model.type <- model$type

  #start checks
  if(!(model.type %in% c("FB", "FBB")) & cluster==T){
    cluster <- F
    warning("Cluster residuals are available only for FB and FBB models, and thus cluster is set equal to FALSE")
  }

  if(!(estimate %in% c("mean", "median", "quantile")))
    stop("Argument `estimate` must be set equal to `mean`, `median` or `quantile`")
  if(estimate=="quantile" & is.null(q))
    stop("You have to specify the order of the quantile")
  if(!(type %in% c("raw", "standardized")))
    stop("Argument `type` must be set equal to `raw` or `standardized`")

  if("flexreg_binom" %in% class(model)){
    n <- model$n
    response <- model$response/n
  } else {
    response <- model$response
  }
  #end checks

  predicted <- predict.flexreg(model, type = "response", newdata = NULL, n.new = NULL, cluster = cluster, estimate = estimate, q = q)
  residuals <- data.frame(raw = response-predicted$response)

  if(cluster == T){
    residuals <- data.frame(residuals, cluster1 = response - predicted$l1, cluster2 = response - predicted$l2)
    residuals$label <- as.numeric(abs(residuals$cluster1) > abs(residuals$cluster2))+1
  }

  if(type == "standardized"){
    variance <- predict.flexreg(object = model, type = "variance", newdata = NULL, n.new = NULL, cluster = cluster, estimate = estimate, q = q)
    residuals.old <- residuals
    residuals  <- data.frame(standardized = residuals$raw / sqrt(variance$variance))

    if( cluster == T){
      residuals <- data.frame(residuals, cluster1 = residuals.old$cluster1 / sqrt(variance$cluster1),
                                         cluster2 = residuals.old$cluster2 / sqrt(variance$cluster2))
      residuals$label <- as.numeric(abs(residuals$cluster1) > abs(residuals$cluster2))+1
    }
  }
  return(residuals)
}




