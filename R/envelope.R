envelopeRQR <- function(x, nsim = 100) {
  n <- length(x)
  x_sorted <- sort(x)
  x_theo <- qnorm(ppoints(n))

  sim_res <- matrix(rnorm(nsim * n), nrow = nsim, ncol = n)
  sim_sorted <- t(apply(sim_res, 1, sort))


  lower <- apply(sim_sorted, 2, min)
  upper <- apply(sim_sorted, 2, max)
  mean_env <- apply(sim_sorted, 2, mean)


  y_range <- range(x_sorted, lower, upper)

  df.enve <-data.frame(x_sorted,lower,mean_env,upper)

  ggplot(df.enve) +
    labs(
      x = "Theoretical quantiles",
      y = "Randomized quantile residuals"
    ) +
    scale_y_continuous(breaks = seq(-3, 3, by = 1)) +

    stat_qq(
      aes(sample = x_sorted),
      colour = "black",
      size = 3
    ) +

    stat_qq(
      aes(sample = lower),
      colour = "blue",
      geom = "line",
      linewidth = 1
    ) +

    stat_qq(
      aes(sample = mean_env),
      colour = "blue",
      geom = "line",
      linewidth = 1
    ) +

    stat_qq(
      aes(sample = upper),
      colour = "blue",
      geom = "line",
      linewidth = 1
    ) +

    theme_bw() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 25)
    )
}

#' Envelope Plot for HZIP Models
#'
#' Produces a normal Q-Q plot with a simulation envelope for
#' randomized quantile residuals from a Hierarchical Zero-Inflated
#' Poisson (HZIP) model fitted via \code{\link{hzip}}.
#'
#' @param object Either an object of class \code{HZIP} (typically
#'   obtained from \code{\link{hzip}}), or a numeric vector of
#'   randomized quantile residuals obtained from
#'   \code{\link{residuals.HZIP}}.
#' @param nsim Integer specifying the number of Monte Carlo simulations
#'   used to construct the simulation envelope. Default is \code{100}.
#' @param method Character string passed to \code{\link{residuals.HZIP}}
#'   when \code{object} is of class \code{HZIP}. One of \code{"ghq"}
#'   (default) or \code{"laplace"}. Ignored if \code{object} is a
#'   numeric vector.
#' @param Q Integer giving the number of Gauss-Hermite quadrature nodes
#'   when \code{method = "ghq"}. Default is \code{21}. Ignored if
#'   \code{object} is a numeric vector or \code{method = "laplace"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The envelope is constructed using Monte Carlo simulation based on
#' \code{nsim} replications of independent standard normal samples.
#' For each order statistic, the 2.5\% and 97.5\% empirical quantiles
#' define the envelope limits.
#'
#' If \code{object} is of class \code{HZIP}, residuals are computed
#' internally via \code{\link{residuals.HZIP}} before plotting.
#'
#' @return A \code{ggplot2} object containing the Q-Q plot with
#'   simulation envelope.
#'
#' @seealso \code{\link{hzip}}, \code{\link{residuals.HZIP}}
#'
#' @examples
#' \donttest{
#' fit.salamander <- hzip(
#'   y ~ mined | mined + spp,
#'   data = salamanders
#' )
#'
#' # Passing the fitted object directly
#' envelope.HZIP(fit.salamander)
#'
#' # Or passing residuals explicitly
#' res <- residuals(fit.salamander)
#' envelope.HZIP(res)
#' }
#'
#' @importFrom stats qnorm rnorm ppoints quantile residuals
#' @importFrom ggplot2 ggplot labs scale_y_continuous stat_qq
#' @importFrom ggplot2 theme_bw theme element_text aes
#'
#' @export
envelope.HZIP <- function(object, nsim = 21, method = "ghq", Q = 15, ...) {
  if (inherits(object, "HZIP")) {
    object <- residuals(object, method = method, Q = Q)
  }
  if (!is.numeric(object)) {
    stop("'object' must be an HZIP model or a numeric vector of residuals.")
  }
  envelopeRQR(object, nsim = nsim)
}
