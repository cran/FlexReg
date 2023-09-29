#' @title Methods for \code{`flexreg`} Objects
#'
#' @description Methods for extracting information from fitted  regression model objects of class \code{`flexreg`}.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param digits an integer indicating the number of decimal places. Default equal to 4.
#' @param ... additional arguments.
#'
#' @details  The \code{\link{summary.flexreg}} method summarizes the results of \code{\link{flexreg}} and \code{\link{flexreg_binom}} functions, adding also information from the functions
#' \code{\link{residuals.flexreg}} and \code{\link{WAIC}}. The \code{\link{summary.flexreg}} method returns an object of class \code{`summary.flexreg`} containing the relevant summary statistics which can subsequently be
#' printed using the associated \code{\link{print.summary.flexreg}} method.
#'
#' @examples
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, n.iter = 1000)
#' summary(FB)
#'
#'
#' @method summary flexreg
#' @export
#'

summary.flexreg <- function(object, ..., digits=4){
  x <- object
  call <- x$call
  model.name <- x$model@model_name
  model.type <- x$type

  if("flexreg_bound" %in% class(x)){
    link.mu <- x$link.mu
    link.phi <- x$link.phi
    formula <- x$formula

    covariate.names.mu <- colnames(x$design.X)
    covariate.names.phi <- colnames(x$design.Z)
    covariate.names.0 <- colnames(x$design.X0)
    covariate.names.1 <- colnames(x$design.X1)

    posterior <- x$model
    pars <- extract.pars(posterior = posterior)
    n.pars <- length(pars)
    pp <- rstan::extract(posterior, pars)
    summa <- lapply(pp, function(x) c(mean(x), sd(x),quantile(x, probs = c(.025,.5,.975))))
    summa.mat <- round(matrix(unlist(summa), ncol=5, byrow = T ),digits)
    colnames(summa.mat) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
    rownames(summa.mat) <- pars

    summ.mu <- summa.mat[grep("beta",rownames(summa.mat)),]
    # if there is only the intercept term:
    if(is.null(dim(summ.mu))){
      summ.mu <- matrix(ncol=length(summ.mu),nrow=1, data=summ.mu)
      colnames(summ.mu) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
    }
    rownames(summ.mu) <- covariate.names.mu

    summ.q0 <- summa.mat[grep("omega0",rownames(summa.mat)),]
    if(is.na(summ.q0[1])){
      summ.q0  <- NULL
    } else if(is.null(dim(summ.q0))){ # if there is only the intercept term:
      summ.q0 <- matrix(ncol=length(summ.q0),nrow=1, data=summ.q0)
      colnames(summ.q0) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
    }
    rownames(summ.q0) <- covariate.names.0

    summ.q1 <- summa.mat[grep("omega1",rownames(summa.mat)),]
    if(is.na(summ.q1[1])){
      summ.q1  <- NULL
    } else if(is.null(dim(summ.q1))){ # if there is only the intercept term:
      summ.q1 <- matrix(ncol=length(summ.q1),nrow=1, data=summ.q1)
      colnames(summ.q1) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
    }
    rownames(summ.q1) <- covariate.names.1

    if(is.null(covariate.names.phi)){
        summ.phi <- summa.mat[which(rownames(summa.mat)=="phi"),]
        dim(summ.phi) <- c(1,5)
        colnames(summ.phi) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
        rownames(summ.phi) <- "phi"
      } else {
        summ.phi <- summa.mat[grep("psi",rownames(summa.mat)),]
        if(is.null(dim(summ.phi)))  dim(summ.phi) <- c(1,5) #summ.phi <- t(as.matrix(summ.phi))
        rownames(summ.phi) <- covariate.names.phi
      }

      n.parz <- nrow(summ.mu)+nrow(summ.phi)+ifelse(is.null(summ.q0),0,nrow(summ.q0))+ifelse(is.null(summ.q1),0,nrow(summ.q1))
      if( n.pars > n.parz) {
        summ.add <- summa.mat[(n.parz+1):n.pars,]
      } else summ.add <- NULL

      residuals <- residuals.flexreg(x, type = "raw", cluster = FALSE, estimate="mean")
      summ.res <- round(quantile(residuals$raw), digits)
      names(summ.res) <- c("Min", "1Q", "Median", "3Q", "Max")


    waic_out <- suppressWarnings(WAIC(x))

    output <- list(call=call,type=model.type, Model=model.name, formula=formula, link.mu=link.mu, link.phi=link.phi,
                   Summary.res=summ.res,
                   Summary.mu=summ.mu, Summary.phi=summ.phi,
                   Summary.q0=summ.q0, Summary.q1=summ.q1,
                   Summary.add=summ.add,
                   waic_out = waic_out)
  } else {
    # ELSE, if the model is for binomial data:
    link.mu <- x$link.mu
    link.theta <- x$link.theta
    formula <- x$formula

    covariate.names.mu <- colnames(x$design.X)
    covariate.names.theta <- colnames(x$design.Z)

    posterior <- x$model
    pars <- extract.pars(posterior = posterior)
    n.pars <- length(pars)
    pp <- rstan::extract(posterior, pars)
    summa <- lapply(pp, function(x) c(mean(x), sd(x),quantile(x, probs = c(.025,.5,.975))))
    summa.mat <- round(matrix(unlist(summa), ncol=5, byrow = T ),digits)
    colnames(summa.mat) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
    rownames(summa.mat) <- pars

    summ.mu <- summa.mat[grep("beta",rownames(summa.mat)),]
    if(is.null(dim(summ.mu))){ # if there is only the intercept term:
      summ.mu <- matrix(ncol=length(summ.mu),nrow=1, data=summ.mu)
      colnames(summ.mu) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
    }
    rownames(summ.mu) <- covariate.names.mu

    if(model.name != "Bin"){
      if(is.null(covariate.names.theta)){
        summ.theta <- summa.mat[which(rownames(summa.mat)=="theta"),]
        summ.theta <- rbind(summ.theta,
                            summa.mat[which(rownames(summa.mat)=="phi"),])
        dim(summ.theta) <- c(2,5)
        colnames(summ.theta) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
        rownames(summ.theta) <- c("theta", "phi")
      } else {
        summ.theta <- summa.mat[grep("psi",rownames(summa.mat)),]
        rownames(summ.theta) <- covariate.names.theta
      }

      n.parz <- nrow(summ.mu)+nrow(summ.theta)
      if( n.pars > n.parz) {
        summ.add <- summa.mat[(n.parz+1):n.pars,]
      } else summ.add <- NULL

      residuals <- residuals.flexreg(x, type = "raw", cluster = FALSE, estimate = "mean")
      summ.res <- round(quantile(residuals$raw), digits)
      names(summ.res) <- c("Min", "1Q", "Median", "3Q", "Max")
    } else {
      summ.phi <- NULL
      summ.add <- NULL
      summ.theta <- NULL

      #summ.res <- NULL

      residuals <- residuals.flexreg(x, type = "raw", cluster = FALSE, estimate = "mean")
      summ.res <- round(quantile(residuals$raw), digits)
      names(summ.res) <- c("Min", "1Q", "Median", "3Q", "Max")
    }

    waic_out <- suppressWarnings(WAIC(x))

    output <- list(call=call, type=model.type, Model=model.name, formula=formula,
                   link.mu=link.mu, link.theta=link.theta,
                   Summary.res=summ.res,
                   Summary.mu=summ.mu,
                   Summary.theta=summ.theta,
                   Summary.add=summ.add,
                   waic_out = waic_out)

  }
  class(output) <- "summary.flexreg"
  return(output)
}


#' Print Methods for summary.flexreg Objects
#'
#' @param x an object of class \code{`summary.flexreg`}.
#' @param ... additional arguments. Currently not used.
#'
#' @rdname summary.flexreg
#' @export
#'

print.summary.flexreg <- function(x, ...){
  cat("Call: ")
  print(x$call)
  cat("\nModel name: ", x$type,
      ifelse(is.null(x$call$zero.formula),"", "with zero augmentation"),
      ifelse(!is.null(x$call$zero.formula) & !is.null(x$call$one.formula), "and", ""),
      ifelse(is.null(x$call$one.formula),"", "with one augmentation"),"\n \n")

  cat("Residuals:\n")
  print(x$Summary.res)

  cat("\nCoefficients (mean model with", x$link.mu, "link):\n")
  print(x$Summary.mu)

  if(x$type %in% c("Beta", "FB",  "VIB")){
    cat("\nCoefficients (precision model with", x$link.phi, "link for phi):\n")
    print(x$Summary.phi)
  } else if(x$Model %in% c("FBB_theta", "BetaBin_theta")){
    cat("\nCoefficients (overdispersion model with", x$link.theta, "link for theta):\n")
    print(x$Summary.theta)
  } else if(x$Model %in% c("FBB", "BetaBin")){
    cat("\nOverdispersion parameter(s):\n")
    print(x$Summary.theta)
  }  # else if Bin then nothing should be printed!

  if(!is.null(x$Summary.q0)){
    cat("\nCoefficients (zero augmentation model with logit link):\n")
    print(x$Summary.q0)
  }

  if(!is.null(x$Summary.q1)){
    cat("\nCoefficients (one augmentation model with logit link):\n")
    print(x$Summary.q1)
  }

  if(!is.null(x$Summary.add)){
    cat("\nAdditional Parameters:\n")
    print(x$Summary.add)
  }

  cat("\nWaic method:")
  suppressWarnings(print(x$waic_out$waic_out))

}

#' Print Methods for flexreg Objects
#'
#' @param x an object of class \code{`flexreg`}.
#' @param ... additional arguments. Currently not used.
#'
#' @rdname print.flexreg
#' @export
#'

print.flexreg <- function(x, ...){
  summ <- summary.flexreg(x)

  cat("Call: ")
  print(summ$call)
  cat("\nModel name: ", summ$type,
      ifelse(is.null(summ$call$zero.formula),"", "with zero augmentation"),
      ifelse(!is.null(summ$call$zero.formula) & !is.null(summ$call$one.formula), "and", ""),
      ifelse(is.null(summ$call$one.formula),"", "with one augmentation"),"\n \n")


  cat("\nCoefficients (mean model with", summ$link.mu, "link):\n")
  print(summ$Summary.mu[,1])

  if(summ$type %in% c("Beta", "FB",  "VIB")){
    cat("\nCoefficients (precision model with", summ$link.phi, "link for phi):\n")
    print(summ$Summary.phi[,1])
  } else if(summ$Model %in% c("FBB_theta", "BetaBin_theta")){
    cat("\nCoefficients (overdispersion model with", summ$link.theta, "link for theta):\n")
    print(summ$Summary.theta[,1])
  } else if(summ$Model %in% c("FBB", "BetaBin")){
    cat("\nOverdispersion parameter(s):\n")
    print(summ$Summary.theta[,1])
  }  # else if Bin then nothing should be printed!

  if(!is.null(summ$Summary.q0)){
    cat("\nCoefficients (zero augmentation model with logit link):\n")
    print(summ$Summary.q0[,1])
  }

  if(!is.null(summ$Summary.q1[,1])){
    cat("\nCoefficients (one augmentation model with logit link):\n")
    print(summ$Summary.q1[,1])
  }

  if(!is.null(summ$Summary.add)){
    cat("\nAdditional Parameters:\n")
    print(summ$Summary.add[,1])
  }

}




#' @title Plot Method for \code{flexreg} Objects
#'
#' @description Method for plotting regression curves for the mean from fitted regression model objects of class \code{`flexreg`}.
#'
#' @param x an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param name.x a character containing the name of the covariate from the mean model to be plotted on the x-axis of the scatterplot.
#' @param additional.cov.default a list of additional covariates from the mean model and their value to be set as default.
#' @param smooth a logical value indicating wheater the curves should be smooth (\code{TRUE}) or piecewise linear (\code{FALSE}, default).
#' @param cluster logical. If the model is \code{"FB"} or \code{"FBB"}, \code{cluster = TRUE} plots the cluster means. By default, \code{cluster = FALSE}.
#' @param type  a vector of characters indicating the regression curves to be plotted. Available options are \code{"response"} and \code{"response.aug"} for augmented models.
#' @param ... additional arguments. Currently not used.
#'
#' @details The function produces a scatterplot of the covariate from the mean model specified in \code{name.x} and \code{y} or \code{y/n} if the response is bounded continuous or discrete, respectively. Any other variable specified in the mean model must be set to a default through the \code{additional.cov.default} argument.
#' The argument \code{type = "response"} plots the conditional mean curve (i.e., \eqn{\mu}), whereas the argument \code{type = "response.aug"}, available only for augmented models,
#' plots the augmented mean curve.
#' If the regression model is of \code{"FB"} or \code{"FBB"} type and \code{cluster = TRUE}, then the function returns two additional curves corresponding to the component means, i.e., \eqn{\lambda_1} and \eqn{\lambda_2}.
#'
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq + dyslexia, data = Reading)
#' plot(FB, name.x="iq", additional.cov.default = list("dyslexia"=1))
#'}
#'
#' @import ggplot2
#'
#' @method plot flexreg
#' @export
#'

plot.flexreg <- function(x, name.x, additional.cov.default = NA, smooth = TRUE,
                         cluster = FALSE, type = "response", ...)
  {

  if(!all(type %in% c("response", "response.aug"))) {
    stop("Argument `type` must be set equal to `response` and/or `response.aug`")
  }

  if((x$type %in% c("Beta", "BetaBin", "Bin", "VIB")) & cluster == TRUE){
    cluster <- FALSE
    warning("Clusters means are not plotted for Beta, VIB,  Binomial, and Beta-Binomal models.")
  }


  group <- Response <- NULL
  object <- x
  model.type <- object$type
  y <- object$response
  y.name <- "y"

  if("flexreg_binom" %in% class(object)){
    n <- object$n
    y <- y/n
    y.name <- "y/n"
  } else{
    n <- NULL
  }

  if(smooth == F){
    x <- object$design.X[,which(colnames(object$design.X) == name.x)]
  } else if(smooth == TRUE){
    x <- object$design.X[,which(colnames(object$design.X) == name.x)]
    N <- dim(unique(object$design.X))[1]
    x <- seq(min(x), max(x), length.out=ifelse(N<50,50,N))
  }

  newdata <- data.frame(x, additional.cov.default)
  names(newdata) <- c(name.x, names(additional.cov.default))

  newdata <- newdata.adjust(newdata, object$formula)

  mu.hat <- predict(object, newdata = newdata, n.new = rep(1, nrow(newdata)), type = "response", cluster = cluster)
  if(!all(type %in% names(mu.hat))){
    type <- "response"
    warning("The augmented response is not available for this model, the response is plotted instead.")
  }

  if(cluster == TRUE) type <-  c(type,"l1", "l2")
  mu.hat <- subset(mu.hat, select = c(type))

  cov.obs <- object$design.X[,which(colnames(object$design.X) == name.x)]

  data.plot <- data.frame(y=y, x=cov.obs)
  pp1 <- ggplot(data=data.plot, aes(x=x, y=y))+
    geom_point()+
    theme_minimal()

  dd <- data.frame(x = x, y = mu.hat)
  dd <- reshape(dd, direction = "long",
                varying = 2:ncol(dd), v.names=c('Response'))
  dd$group <- as.factor(dd$time)

  plot1 <-
    pp1 + geom_line(data=dd, aes(x=x, y=Response, group=group, colour = group)) +
    ylab(y.name)

  labels <- c(expression(mu), expression(mu[aug]), expression(lambda[1]),expression(lambda[2]))
  codice <- c("response", "response.aug", "l1", "l2")
  labels <- labels[match(names(mu.hat), codice)]
  col <- c("black", "#009E73", "#D55E00", "#0072B2")
  col <- col[match(names(mu.hat), codice)]

  pp <- plot1+ geom_point() +
    scale_color_manual(labels = labels, values = col) +
    theme(legend.position="bottom") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 15))
  return(pp)
}


