#' Convergence plots
#'
#' The function produces some convergence plots from the Monte Carlo draws.
#' @param model an object of class \code{`flexreg`}.
#' @param file a character string giving the name of the file (including the extension .pdf) containing the convergence plots. If \code{NULL}, the convergence plots are printed in the graphics window.
#' @param plotfun  an optional character vector of diagnostics plots. The default is to compute \code{"all"} plots, otherwise one can specify a subset of plots among \code{"density"}, \code{"trace"}, \code{"intervals"}, \code{"rate"}, \code{"rhat"}, and \code{"acf"}.
#' @param pars an optional character vector of parameter names. If \code{pars} is not specified, all parameters in the regression models are evaluated.
#' @param point_est an optional character to specify the point estimate to be shown between \code{"median"} (the default), \code{"mean"}, or \code{"none"}.
#' @param prob the probability mass to be included in the inner interval (for \code{"intervals"} plot) or in the shaded region (for \code{"density"} plot). The default is 0.5.
#' @param prob_outer the probability mass to be included in the outer interval  of the \code{"intervals"} plot. The default is 0.9.
#' @param lags the number of lags to be shown in the \code{"acf"} plot. The default is 10.
#' @param width,height the width and height of the graphics region of each plot in inches. The default values are 7.
#' @param warmup a logical scalar indicating whether to include the warmup draws or not (default).
#'
#' @return A .pdf file with one plot per page.
#'
#' @import bayesplot ggplot2
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, type = "FB")
#' convergence.plot(FB, file = "Convergence_plot_Output.pdf", pars = "beta")
#'}
#'
#' @references {
#' Brooks, SP., Gelman, A. (1998). General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics, \bold{7}, 434-455. \cr
#' \cr
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#' }
#'
#' @details
#' \itemize{
#' \item \code{"density"} returns a density plot for each parameter in \code{pars} computed from the posterior draws. See \code{bayesplot::mcmc_areas} for further details.
#' \item \code{"trace"} returns a trace plot for each parameter in \code{pars} computed from the posterior draws. See \code{bayesplot::mcmc_trace} for further details.
#' \item \code{"intervals"} returns a plot of uncertainty interval for each parameter in \code{pars} computed from the posterior draws. See \code{bayesplot::mcmc_intervals} for further details.
#' \item \code{"rate"} returns a plot for each parameter in \code{pars}  with the number of iterations on the x-axis and the Monte Carlo mean until iteration i-th on the y-axis.
#' \item \code{"rhat"} returns a plot with the Rhat values for each parameter in \code{pars}. See \code{bayesplot::mcmc_rhat} for further details.
#' \item \code{"acf"} returns the autocorrelation plots (one for each parameter in \code{pars}). See \code{bayesplot::mcmc_acf} for further details.
#'
#' }
#' Moreover, the convergence plots can be further customized using the \pkg{ggplot2} package.
#'
#' @export
#'




convergence.plot <- function(model, file = "convergence-output.pdf",
                             plotfun = "all", pars = NULL, point_est = "median",
                             prob = 0.5, prob_outer = 0.9, lags = 10, warmup = F,
                             width = 7, height = 7)
{

  if(is.null(model$aug)) model$aug <- "No"

  if(model$aug %in% c("0", "1")){
    beta_index <- grep("beta", names(model$model[[2]]))
    names(model$model[[2]])[beta_index] <-
      paste0("omega", model$aug, "[", 1:length(beta_index), "]")
  }

  if(!is.null(pars)){
    omega_index <- grep("omega", pars)

    if(length(omega_index) == 0){ # if no augmentation
      pars_1 <- pars

      out1 <- convergence.plot.internal(model$model[[1]], file,
                                        plotfun, pars_1, point_est,
                                        prob, prob_outer, lags, warmup,
                                        width, height)
      out2 <- NULL
    } else {
      pars_1 <- pars[-grep("omega", pars)]

      if(model$aug == "01"){
        pars_2 <- pars[grep("omega", pars)]

        out2 <- convergence.plot.internal(model$model[[2]], file,
                                          plotfun, pars_2, point_est,
                                          prob, prob_outer, lags, warmup,
                                          width, height)
      } else {
        pars_2 <- paste0("omega", model$aug, "[", 1:length(beta_index), "]")

        out2 <- convergence.plot.internal(model$model[[2]], file,
                                          plotfun, pars_2, point_est,
                                          prob, prob_outer, lags, warmup,
                                          width, height)
      }

      if(length(pars_1) >0){
        out1 <- convergence.plot.internal(model$model[[1]], file,
                                          plotfun, pars_1, point_est,
                                          prob, prob_outer, lags, warmup,
                                          width, height)
      } else out1 <- NULL


    }
  } else { # pars == NULL
    out1 <- convergence.plot.internal(model$model[[1]], file,
                                      plotfun, pars, point_est,
                                      prob, prob_outer, lags, warmup,
                                      width, height)
    if(model$aug != "No"){
      out2 <- convergence.plot.internal(model$model[[2]], file,
                                        plotfun, pars, point_est,
                                        prob, prob_outer, lags, warmup,
                                        width, height)
    } else out2 <- NULL
  }

  output <- vector(mode = "list", 2)
  names(output) <- c("MainModel", "AugModel")
  output$MainModel <- out1
  output$AugModel <- out2
  M <- length(out1)+length(out2)


  D <- length(pars)
  file.extension.position <- regexpr("\\.([[:alnum:]]+)$",file)
  file.extension <- tolower(substr(file, file.extension.position +
                                     1, nchar(file)))
  file.name <- substr(file, 1, file.extension.position - 1)

  if(is.null(file)){

    if(!is.null(out1)){
      for(i in 1:length(out1)){
        question <- utils::menu(c("Yes", "No", "Exit"), title=paste0("Plot ", i, "/",length(out1), " of the main model", "\nDo you want to see it?"))
        #grDevices::devAskNewPage(ask = TRUE)
        if (question == 1)
          #invisible
          (utils::capture.output(print(out1[[i]])))
        else if (question == 2) next()
        else if (question ==3) break()#oppure next()
      }
    }

    if(!is.null(out2)){
      for(i in 1:length(out2)){
        question <- utils::menu(c("Yes", "No", "Exit"), title=paste0("Plot ", i, "/",length(out2)," of the augmented model", "\nDo you want to see it?"))
        #grDevices::devAskNewPage(ask = TRUE)
        if (question == 1)
          #invisible
          (utils::capture.output(print(out2[[i]])))
        else if (question == 2) next()
        else if (question ==3) break()#oppure next()
      }
    }

  } else if (file.extension == "pdf") {
    out <- c(out1, out2)
    grDevices::pdf(file, width = width, height = height)
    invisible(utils::capture.output(print(out)))#DA QUI
    garbage <- grDevices::dev.off()
  } else stop("File extension not supported")
}


#' internal function
#' @keywords internal
#'
convergence.plot.internal <- function(posterior, file,
                                      plotfun, pars, point_est,
                                      prob, prob_outer, lags, warmup,
                                      width, height)
{

  if (is.null(pars)) pars <- extract.pars(posterior = posterior)
  D <- length(pars)
  chains <- rstan::extract(posterior, pars, inc_warmup = warmup, permuted=F)#ottengo degli array
  n.iter <- dim(posterior)[1]
  n.warmup <- dim(chains)[1]-n.iter
  pars.full <- names(posterior)

  #if one writes beta, it comprises all betas
  if ("beta" %in% pars) {
    pars <- c(pars.full[grep("beta",pars.full)], pars)
    pars <- pars[-which(pars == "beta")]
  }
  if ("psi" %in% pars) {
    pars <- c(pars.full[grep("psi",pars.full)], pars)
    pars <- pars[-which(pars=="psi")]
  }
  if ("omega0" %in% pars) {
    pars <- c(pars.full[grep("omega0",pars.full)], pars)
    pars <- pars[-which(pars=="omega0")]
  }
  if ("omega1" %in% pars) {
    pars <- c(pars.full[grep("omega1",pars.full)], pars)
    pars <- pars[-which(pars=="omega1")]
  }

  ppp <- c()
  #density plot
  if (any(plotfun %in% c("all", "density"))){
    pp <- lapply(pars, function(x)  mcmc_areas(posterior, pars = x, prob = prob))
    title.plot <- lapply(pars, function(x) ggtitle(paste("Posterior distribution of ", x,
                                                         "with", point_est,"and", prob*100, "% intervals", sep=" ")))
    pp <- lapply(1:D, function(x) pp[[x]] + title.plot[[x]])
    ppp <- c(ppp, pp)
  }

  #trace plot
  if (any(plotfun %in% c("all", "trace"))){
    n_warmup <- ifelse(warmup==F, 0, n.warmup)
    pp <- lapply(pars, function(x) mcmc_trace(chains, pars = x, n_warmup = n_warmup))
    title.plot <- lapply(pars, function(x) ggtitle(paste("Traceplot of ", x, sep=" ")))
    pp <- lapply(1:D, function(x) pp[[x]] + title.plot[[x]]+theme(axis.title.y=element_text(angle=0,hjust=1)) )
    ppp <- c(ppp, pp)
  }

  #intervals plot
  if (any(plotfun %in% c("all", "intervals"))){
    pp <-  lapply(pars, function(x) mcmc_intervals(posterior, pars = x, prob = prob, prob_outer = prob_outer))
    title.plot <- lapply(pars, function(x) ggtitle(paste("Posterior interval estimates of ", x, sep=" ")))
    pp <- lapply(1:D, function(x) pp[[x]] + title.plot[[x]])
    ppp <- c(ppp, pp)
  }

  #rate plot
  if (any(plotfun %in% c("all", "rate"))){
    pp <- rate_plot(chains, pars, n.warmup = n.warmup)
    ppp <- c(ppp, pp)
  }

  if (any(plotfun %in% c("all", "rhat"))){
    pp <- mcmc_rhat(rhat = rhat(posterior, pars = pars)) + yaxis_text(hjust = .0) + ggtitle("Rhat plot")
    ppp <- c(ppp, list(pp))
  }

  #acf plot
  if (any(plotfun %in% c("all", "acf"))){
    pp <- mcmc_acf(posterior, pars = pars, lags = lags) + ggtitle("ACF plot")
    ppp <- c(ppp, list(pp))
  }

  return(ppp)

}



#' Convergence diagnostics
#'
#' The function returns some diagnostic measures to check for convergence to the equilibrium distribution of the Markov Chain(s).
#' Moreover, it prints the number (and percentage) of iterations that ended with a divergence and that saturated the max treedepth, and the E-BFMI values for each chain for which E-BFMI is less than 0.2.
#'
#' @param model an object of class \code{`flexreg`}.
#' @param diagnostics  an optional character vector of diagnostics names. The default is to compute \code{"all"} diagnostics, otherwise one can specify a selection of diagnostics among \code{"Rhat"}, \code{"geweke"}, \code{"raftery"}, \code{"heidel"}, and \code{"gelman"}.
#' @param pars an optional character vector of parameter names. If \code{pars} is not specified, all parameters in the regression models are evaluated.
#' @param additional.args a list containing additional arguments (see details)
#'
#' @return A print from \code{\link[rstan]{check_hmc_diagnostics}} function and a list of convergence diagnostics.
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, type = "FB")
#' convergence.diag(FB,  diagnostics = c("Rhat", "geweke"), pars = "beta")
#'
#'}

#'
#' @references {
#' Brooks, SP., Gelman, A. (1998). General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics, \bold{7}, 434-455. \cr
#' \cr
#' Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to calculating posterior moments. In Bayesian Statistics 4 (ed JM Bernado, JO Berger, AP Dawid and AFM Smith). Clarendon Press, Oxford, UK. \cr
#' \cr
#' Heidelberger P., Welch P.D. (1981). A spectral method for confidence interval generation and run length control in simulations. Comm. ACM. \bold{24}, 233-245.\cr
#' \cr
#' Raftery, A.E. and Lewis, S.M. (1992). One long run with diagnostics: Implementation strategies for Markov chain Monte Carlo. Statistical Science, \bold{7}, 493-497.\cr
#' \cr
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#' }
#'
#' @details \itemize{
#' \item \code{"Rhat"} returns the potential scale reduction factor on split chains. An R-hat greater than 1 is indicative of a bad mix of the chains. At convergence R-hat has to be less than 1.05. See \code{rstan::Rhat} for further details.
#' \item \code{"geweke"} returns the z-scores, one for each parameter, for a test of equality between the means of the first 10\% and last 50\% of the chain. The fraction to use from the first and last part of the chain can be edited through the additional arguments \code{frac1} and \code{frac2}. The sum of \code{frac1} and \code{frac2} has to be strictly less than 1. See \code{coda::geweke.diag} for further details.
#' \item \code{"raftery"} returns the estimate of the "dependence factor" \eqn{I}. Values of \eqn{I} greater than 5 may indicate a strong autocorrelation.
#' Additional parameters such as the quantile to be estimated (\code{q}), the desired margin of error of the estimate (\code{r}), and the probability (\code{s}) of obtaining an estimate between \eqn{q-r} and \eqn{q+r}  can be passed as a list in the \code{additional.args} argument. See \code{coda::raftery.diag} for further details.
#' \item \code{"heidel"} returns the p-values, one for each parameter, referred to a convergence test where the null hypothesis is that the sampled values come from a stationary distribution. It is possible to set the target value for ratio of halfwidth to sample mean (\code{eps}) and the significance level of the test (\code{pvalue})  into the \code{additional.args} argument. See \code{coda::heidel.diag} for further details.
#' \item \code{"gelman"} returns the estimate of the potential scale reduction factor and the upper confidence limit. At least two chains are needed to compute the Gelman and Rubin's convergence diagnostic. Additional parameters such as the confidence level (\code{confidence}), a logical flag indicating whether variables should be transformed (\code{transform}),
#' a logical flag indicating whether only the second half of the series should be used in the computation (\code{autoburnin}), and a logical flag indicating whether the multivariate potential scale reduction factor should be calculated for multivariate chains (\code{multivariate}) can be passed as a list in the \code{additional.args} argument. See \code{coda::gelman.diag} for further details.
#' }
#'
#' @import rstan
#'
#' @export
#'

convergence.diag <- function(model, diagnostics = "all", pars = NULL, additional.args=list())
{
  n.chain <- model$call$n.chain
  if(is.null(model$aug)) model$aug <- "No"

  if(model$aug %in% c("0", "1")){
    beta_index <- grep("beta", names(model$model[[2]]))
    names(model$model[[2]])[beta_index] <-
      paste0("omega", model$aug, "[", 1:length(beta_index), "]")
  }
  if(!is.null(pars)){
    omega_index <- grep("omega", pars)

    if(length(omega_index) == 0){ # if no augmentation
      pars_1 <- pars

      out1 <- convergence.diag.internal(model$model[[1]], pars = pars_1, additional.args = additional.args, aug = NA, n.chain = n.chain)
      out2 <- NULL
    } else {
      pars_1 <- pars[-grep("omega", pars)]

      if(model$aug == "01"){
        pars_2 <- pars[grep("omega", pars)]

        out2 <- convergence.diag.internal(model$model[[2]], pars = pars_2, additional.args = additional.args, aug = model$aug, n.chain = n.chain)
      } else {
        pars_2 <- paste0("omega", model$aug, "[", 1:length(beta_index), "]")

        out2 <- convergence.diag.internal(model$model[[2]], pars = pars_2, additional.args = additional.args, aug = model$aug, n.chain = n.chain)
      }

      if(length(pars_1) >0){
        out1 <- convergence.diag.internal(model$model[[1]], pars = pars_1, additional.args = additional.args, aug = NA, n.chain = n.chain)
      } else out1 <- NULL


    }
  } else { # pars == NULL
    out1 <- convergence.diag.internal(model$model[[1]], pars = pars, additional.args = additional.args,aug = NA, n.chain = n.chain)
    if(model$aug != "No"){
      out2 <- convergence.diag.internal(model$model[[2]], pars = pars, additional.args = additional.args, aug = model$aug, n.chain = n.chain)
    } else out2 = NULL
  }

  output <- vector(mode = "list", 2)
  names(output) <- c("MainModel", "AugModel")
  output$MainModel <- out1
  output$AugModel <- out2



  class(output) <- c("convergence.diag.flexreg", class(model)[2])

  return(output)
}

#' internal function
#' @keywords internal

convergence.diag.internal <- function(posterior, diagnostics = "all", pars = NULL, additional.args=list(), aug = NULL, n.chain)
{
  #posterior <- model
  hmc.diag <- rstan::check_hmc_diagnostics(posterior)

  if(is.null(pars)) pars <- extract.pars(posterior = posterior)
  mcmc <- rstan::extract(posterior, pars=pars, permuted=F)

  if(any(diagnostics == "all")) diagnostics  <- c("Rhat", "geweke", "raftery", "heidel", "gelman")
  #Rhat diag
  if (any(diagnostics %in% "Rhat")){
    rhat <- apply(mcmc, 3, function(x) rstan::Rhat(x))
    diagnostics <- diagnostics[-which(diagnostics=="Rhat")]
  } else rhat <- NULL

  mcmc <- rstan::As.mcmc.list(posterior, pars=pars)

  if(aug %in% c("0", "1")){
    beta_index <- grep("beta", colnames(mcmc[[1]]))
    colnames(mcmc[[1]])[beta_index] <-
      paste0("omega", aug, "[", 1:length(beta_index), "]")
  }

  #gestiamo gli argomenti aggiuntivi
  arg.coda <- lapply(diagnostics, function(x)
    formals(utils::getFromNamespace(paste0(x,".diag"), "coda"))[-1])
  names(arg.coda) <- diagnostics

  if(length((additional.args))>0){
    arg.coda <- unlist(arg.coda)
    arg.match <- lapply(names(additional.args), function(x) grepl(x,names(arg.coda)))
    for(i in 1:length(arg.match))  arg.coda[arg.match[[i]]] <- additional.args[i]
  }
  arg.coda <- unlist(arg.coda)
  e <- names(arg.coda)
  out <- within(as.list(arg.coda),{

    if (any(diagnostics %in%  "gelman")){
      if(is.null(n.chain)) n.chain <- 1
      if (n.chain < 2){
        diagnostics <- diagnostics[-which(diagnostics=="gelman")]
        warning("You need at least two chains to compute Gelman and Rubin's convergence diagnostic")}
    }

    other.diag <- lapply(diagnostics, function(x)
      do.call(utils::getFromNamespace(paste0(x,".diag"), "coda"),append(list(mcmc),
                                                                        lapply(e[startsWith(e,x)], function(a) get(a)))))

    names(other.diag) <- diagnostics
  })

  if(!is.null(rhat))  out$other.diag$Rhat <- rhat

  return(list(out$other.diag))

}




#' Print Methods for convergence.diag.flexreg Objects
#'
#' @param x an object of class \code{`convergence.diag.flexreg`}.
#' @param ... additional arguments. Currently not used.
#'
#' @rdname convergence.diag
#' @export
#'

print.convergence.diag.flexreg <- function(x, ...){

  if(!is.null(x$MainModel)){
    if("flexreg_bound" %in% class(x)){
      cat("Convergence Diagnostics for the continuous part of the model: \n \n")
      print(x$MainModel)
    } else {
      cat("Convergence Diagnostics: \n \n")
      print(x$MainModel)
    }
  }

  if(!is.null(x$AugModel)){
    cat("Convergence Diagnostics for the augmented part of the model: \n \n")
    print(x$AugModel)
  }
}
