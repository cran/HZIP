#' Test for Overdispersion
#'
#' Generic function for simulation-based dispersion tests.
#'
#' @param object A fitted model object.
#' @param ... Further arguments passed to methods.
#'
#' @export
testDisp <- function(object, ...) UseMethod("testDisp")

#' Test for Zero-Inflation
#'
#' Generic function for simulation-based zero-inflation tests.
#'
#' @param object A fitted model object.
#' @param ... Further arguments passed to methods.
#'
#' @export
testZI <- function(object, ...) UseMethod("testZI")
