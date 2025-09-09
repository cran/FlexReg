#' Construct Design Matrices for flexreg Objects
#'
#' @description Method for extracting design matrices from fitted
#' regression model objects of class \code{`flexreg`}.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of
#'  \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param ... additional arguments. Currently not used.
#'
#' @details The method returns a list containing all the design matrices involved
#'  in the regression models, namely \code{`X`} (regression on the mean), \code{`Z`}
#'  (regression on the precision or overdispersion parameters), \code{`X0`} (regression
#'  on the augmentation in zero probability), and/or \code{`X1`} (regression
#'  on the augmentation in one probability).
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, n.iter = 1000)
#' model.matrix(FB)
#' }
#'
#'
#' @method model.matrix flexreg
#'
#' @export
#'
model.matrix.flexreg <- function(object, ...){

  if("flexreg_bound" %in% class(object)){
    matrix_list <- vector(mode="list", 4)
    names(matrix_list) <- c("X", "Z", "X0", "X1")

    matrix_list$X <- object$design.X
    matrix_list$Z <- object$design.Z
    matrix_list$X0 <- object$design.X0
    matrix_list$X1 <- object$design.X1


  } else if("flexreg_binom" %in% class(object)) {
    matrix_list <- vector(mode="list", 2)
    names(matrix_list) <- c("X", "Z")

    matrix_list$X <- object$design.X
    matrix_list$Z <- object$design.Z

  }

  return(matrix_list)

}


#' Construct Design Matrices for flexreg Objects
#'
#' @description Method for extracting design matrix from fitted
#' regression model objects of class \code{`flexreg`}.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of
#'  \code{\link{flexreg}} or \code{\link{flexreg_binom}} functions.
#' @param ... additional arguments. Currently not used.
#'
#' @details The method returns a list containing all the design matrices involved
#'  in the regression models, namely \code{`X`} (regression on the mean), \code{`Z`}
#'  (regression on the precision or overdispersion parameters), \code{`X0`} (regression
#'  on the augmentation in zero probability), and/or \code{`X1`} (regression
#'  on the augmentation in one probability).
#'
#' @examples
#' \dontrun{
#' data("Reading")
#' FB <- flexreg(accuracy.adj ~ iq, data = Reading, n.iter = 1000)
#' model.frame(FB)
#' }
#'
#' @import stats
#'
#' @export
#'
model_frame <- function(object, ...){
  formula <- object$formula
  formula <- Formula::as.Formula(formula)

  # var.names <- c(as.character(attr(formula, "lhs")[[1]]),
  #                as.character(attr(formula, "rhs")[[1]]))
  #
  # if(!is.null(object$design.Z)) varnames <- c(var.names, as.character(attr(formula, "rhs")[[2]]))
  # if(!is.null(object$design.X0)) varnames <- c(var.names, as.character(attr(formula, "rhs")[[3]]))
  # if(!is.null(object$design.X1)) varnames <- c(var.names, as.character(attr(formula, "rhs")[[4]]))
  #
  # var.names <- unique(var.names)
  #
  # var.names <- setdiff(var.names, c("1", "+"))

  data <- object$call$data


  if(is.null(data)) {
    data <- environment(formula)
   # out <- model.frame(formula, data)
  } else{
  data <- eval(data)
  #out <- data[, var.names]
  }
  out <- model.frame(formula, data)
  return(out)
}

#
#
# model.frame.flexreg <- function(formula, zero.formula=NULL, one.formula=NULL, data = NULL, ...){
#
#   formula <- Formula::as.Formula(formula)
#
#   if(length(formula)[2] >= 2){
#     model.phi <- TRUE
#       } else {
#     model.phi <- FALSE
#       }
#   if(!is.null(one.formula)) formula1 <- Formula::as.Formula(one.formula)
#   if(!is.null(zero.formula)) formula0 <- Formula::as.Formula(zero.formula)
#
#     #put all formulas together
#     if(!is.null(zero.formula)|!is.null(one.formula)){
#       formula <- paste0(deparse(attr(formula, "lhs")[[1]]),"~",
#                               deparse(attr(formula, "rhs")[[1]]),"|",
#                               ifelse(model.phi==F,1,deparse(attr(formula, "rhs")[[2]])),"|",
#                               ifelse(is.null(zero.formula),1,deparse(attr(formula0, "rhs")[[1]])),
#                               "|", ifelse(is.null(one.formula),1,deparse(attr(formula1, "rhs")[[1]])))
#      # formula <- formula.final
#     }
#       formula <- Formula::as.Formula(formula)
#
#     if(is.null(data)) {
#     data <- environment(formula)
#     # out <- model.frame(formula, data)
#   } else{
#     data <- eval(data)
#     #out <- data[, var.names]
#   }
#   out <- model.frame(formula, data)
#   return(out)
# }
