#' Salamanders data
#'
#' This dataset is adapted from the \pkg{glmmTMB} package and contains
#' salamander counts with information on mining status and species.
#' It is intended for illustrating zero-inflated Poisson models with
#' random effects using the \code{hzip()} function.
#'
#' @format A data frame with 644 rows and 4 variables:
#' \describe{
#'   \item{Ind}{Individual identifier (factor).}
#'   \item{y}{Count response variable (integer).}
#'   \item{mined}{Mining status: \code{"yes"} or \code{"no"}.}
#'   \item{spp}{Species factor with multiple levels (e.g., GP, PR, DM, etc.).}
#' }
#'
#' @details
#' The dataset was originally included in the \pkg{glmmTMB} package
#' (Brooks et al., 2017), and has been slightly modified for testing
#' the \pkg{HZIP} package.
#'
#' @source
#' Adapted from the \pkg{glmmTMB} package.
#'
#' @examples
#' \donttest{
#' data(salamanders, package = "HZIP")
#'
#' ## Fit zero-inflated Poisson with random effects
#' fit.salamander <- hzip(y ~ mined | mined + spp + mined * spp,
#'                        data = salamanders)
#' summary(fit.salamander)
#' }
"salamanders"
