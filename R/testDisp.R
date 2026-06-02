#' Simulation-Based Dispersion Test for HZIP Models
#'
#' Compares the observed variance of the response with the distribution
#' of variances obtained from parametric bootstrap replicates simulated
#' under the fitted \code{HZIP} model.
#'
#' @param object An object of class \code{HZIP}.
#' @param nsim A positive integer giving the number of simulated
#'   replicates used to build the reference distribution.
#'   Defaults to \code{500}.
#' @param alternative Character string specifying the alternative
#'   hypothesis. One of \code{"two.sided"} (default),
#'   \code{"greater"} (overdispersion), or \code{"less"}
#'   (underdispersion).
#' @param alpha Numeric value in \code{(0, 1)} giving the significance level
#'   used for the decision. Defaults to \code{0.05}.
#' @param seed An optional integer seed for reproducibility.
#' @param plot Logical; if \code{TRUE} (default) a histogram of the
#'   simulated statistics is printed with the observed value marked.
#' @param ... Further arguments (currently ignored).
#'
#' @return An object of class \code{"hzip_test"} (invisibly), a list
#'   with components:
#'   \describe{
#'     \item{\code{statistic}}{Observed variance of the response.}
#'     \item{\code{p.value}}{Simulation-based p-value.}
#'     \item{\code{alternative}}{The alternative hypothesis used.}
#'     \item{\code{alpha}}{Significance level used for the decision.}
#'     \item{\code{simulated}}{Numeric vector of simulated variances.}
#'     \item{\code{method}}{Description of the test.}
#'     \item{\code{h0}}{Statement of the null hypothesis.}
#'   }
#'
#' @seealso \code{\link{testZI.HZIP}}, \code{\link{simulate.HZIP}}
#'
#' @examples
#' \donttest{
#' fit.salamander <- hzip(
#'   y ~ mined | mined + spp,
#'   data = salamanders
#' )
#'
#' testDisp(fit.salamander, nsim = 200, seed = 42)
#' }
#'
#' @importFrom stats var model.response model.frame
#' @importFrom Formula Formula
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_bw
#'
#' @export
testDisp.HZIP <- function(object,
                                nsim        = 500,
                                alternative = c("two.sided", "greater", "less"),
                                alpha       = 0.05,
                                seed        = NULL,
                                plot        = TRUE,
                                ...) {

  alternative <- match.arg(alternative)

  formula <- object$formula
  data    <- object$data
  Y <- model.response(model.frame(Formula::Formula(formula), data = data))

  obs_stat <- var(Y)

  sims     <- simulate(object, nsim = nsim, seed = seed)
  sim_stat <- apply(sims, 2, var)

  p.value <- switch(alternative,
    two.sided = {
      p_low  <- mean(sim_stat <= obs_stat)
      p_high <- mean(sim_stat >= obs_stat)
      2 * min(p_low, p_high)
    },
    greater = mean(sim_stat >= obs_stat),
    less    = mean(sim_stat <= obs_stat)
  )

  result <- structure(
    list(
      statistic   = obs_stat,
      p.value     = p.value,
      alternative = alternative,
      alpha       = alpha,
      simulated   = sim_stat,
      method      = "Simulation-based dispersion test for HZIP models",
      h0          = "Var(Y_obs) = E[Var(Y_sim)]",
      stat_name   = "Variance"
    ),
    class = "hzip_test"
  )

  if (plot) print(.plot_hzip_test(result))
  print(result)

  invisible(result)
}
