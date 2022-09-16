#' @title Methods for flexreg Objects
#'
#' @description Methods for extracting information from fitted  regression model objects of class \code{`flexreg`}.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}}.
#' @param digits an integer indicating the number of decimal places. Default equal to 4.
#' @param ... additional arguments.
#'
#' @details  The \code{summary.flexreg} method summarizes the results of \code{\link{flexreg}} and \code{\link{flexreg_binom}} functions, adding also information from the functions
#' \code{\link{residuals.flexreg}} and \code{\link{WAIC}}. The \code{summary.flexreg} method returns an object of class \code{`summary.flexreg`} containing the relevant summary statistics which can subsequently be
#' printed using the associated \code{print} method.
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
  if(model.type %in% c("Beta","FB","VIB")){
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
    if(is.null(dim(summ.mu))){ # if there is only the intercept term:
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

      residuals <- residuals.flexreg(x, type = "raw", cluster=FALSE, estimate="mean")
      summ.res <- round(quantile(residuals), digits)
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

      residuals <- residuals.flexreg(x, type = "raw", cluster=FALSE, estimate="mean")
      summ.res <- round(quantile(residuals), digits)
      names(summ.res) <- c("Min", "1Q", "Median", "3Q", "Max")
    } else {
      summ.phi <- NULL
      summ.add <- NULL
      summ.theta <- NULL

      #summ.res <- NULL

      residuals <- residuals.flexreg(x, type = "raw", cluster=FALSE, estimate="mean")
      summ.res <- round(quantile(residuals), digits)
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


#' @title Plot method for flexreg Objects
#'
#' @description Method for plotting regression curves for the mean from fitted regression model objects of class \code{`flexreg`}.
#'
#' @param x an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}} or \code{\link{flexreg_binom}}.
#' @param name.x a character containing the name of the covariate from the mean model to be plotted on the x-axis of the scatterplot.
#' @param additional.cov.default a list of additional covariates from the mean model to be set as default.
#' @param ... additional arguments. Currently not used.
#'
#' @details The function produces a scatterplot of the covariate from the mean model specified in \code{name.x} and \code{y} or \code{y/n} if the response is continuous bounded or binomial, respectively. Any other variable specified in the mean model must be set to a default through the \code{additional.cov.default} argument.
#' If the regression model is of \code{FB} without augmentation or \code{FBB} type the function returns a scatterplot with three curves, one corresponding to the overall mean and two corresponding to the component means of the FB distribution, i.e., \eqn{\lambda_1} and \eqn{\lambda_2}.
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

plot.flexreg <- function(x, name.x, additional.cov.default = NA, ...)
  {
  group <- Response <- NULL
  #assign("group", "Response", envir = .GlobalEnv)
  object <- x
  model.name <- object$model@model_name
  y <- object$response
  y.name <- "y"
  if(model.name %in% c("FBB", "FBB_theta", "BetaBin",
                                    "BetaBin_theta", "Bin")){
    n <- object$n
    y <- y/n
    y.name <- "y/n"
  }
  x <- object$design.X[,which(colnames(object$design.X) == name.x)]

 # additional.cov.default <- list("x2" = 0)#list("dyslexia"= -1, "iq:dyslexia" = 0)

  #prova con modello senza intercetta
  intercept <- ifelse("(Intercept)"  %in% colnames(object$design.X), T, F)
  if( intercept == T){
  newdata <- data.frame(1, x, additional.cov.default)
  names(newdata) <- c("(Intercept)", name.x, names(additional.cov.default))
  } else {
  newdata <- data.frame(x, additional.cov.default)
  names(newdata) <- c(name.x, names(additional.cov.default))
  }
  cluster <-   model.name %in% c("FBNo", "FBNo_phi", "FBB", "FBB_theta")
  if (cluster == T) {
  mu.hat <- predict(object, newdata = newdata, type = "response", cluster = T)
  data.plot <- data.frame(y = y, x = x, mu.hat)
  data.plot <- reshape(data.plot, direction = "long",
                       varying = 3:5, v.names=c('Response'))
  data.plot$group <- as.factor(data.plot$time)
  plot1 <- ggplot(data.plot, aes(x=x, y=y, group = group))+ geom_point()+
    ylab(y.name)+
    geom_line(aes(x=x, y=Response, colour = group))
  } else {
    mu.hat <- predict(object, newdata = newdata, type = "response", cluster = F)
    data.plot <- data.frame(y = y, x = x,  mu.hat)
    names(data.plot)[3] <- "Response"
    plot1 <- ggplot(data.plot, aes(x=x, y=y))+
      geom_line(aes(x=x, y=Response))+
      ylab(y.name)
  }
  return(plot1+ geom_point()+theme_bw()+
           scale_color_manual(labels = c(expression(mu), expression(lambda[1]),expression(lambda[2])), values = c("black", "red", "blue"))+
        theme(legend.title = element_blank()))
}





