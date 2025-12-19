envelopeRQR <- function(x, nsim = 100) {
  n <- length(x)
  x_sorted <- sort(x)
  x_theo <- qnorm(ppoints(n))

  sim_res <- matrix(rnorm(nsim * n), nrow = nsim, ncol = n)
  sim_sorted <- t(apply(sim_res, 1, sort))


  lower <- apply(sim_sorted, 2, quantile, probs = 0.025)
  upper <- apply(sim_sorted, 2, quantile, probs = 0.975)
  mean_env <- apply(sim_sorted, 2, mean)


  y_range <- range(x_sorted, lower, upper)

  df.enve <-data.frame(x_sorted,lower,mean_env,upper)

  ggplot(df.enve) +
    labs(x = "Theorical quantiles",y="Randomized quantile residuals")+
    scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
    stat_qq(aes(sample = x_sorted), colour = "black",size=3) +
    stat_qq(aes(sample = lower), colour = "blue",geom="line", size=1) +
    stat_qq(aes(sample = mean_env), colour = "blue",geom="line", size=1) +
    stat_qq(aes(sample = upper), colour = "blue",geom="line", size=1) +
    theme_bw()+
    theme(legend.position="none")+
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25))
}
#' Envelope simulation for HZIP Model
#'
#' Produces a Q-Q plot of residuals from a Hierarchical Zero-Inflated Poisson (HZIP) Model fitted via \code{\link{hzip}}.
#'
#' @param object An object of class \code{HZIP}, typically returned by \code{\link{hzip}}.
#' @param nsim Integer. Number of simulations used to construct the envelope. Default is \code{100}.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' A simulation envelope is added using Monte Carlo replications.
#'
#' @return Envelope simulation plot.
#'
#' @seealso \code{\link{hzip}}, \code{\link{residuals.HZIP}}
#'
#' @examples
#' \donttest{
#' data(salamanders)
#' fit.salamander <- hzip(y ~ mined|mined+spp,data = salamanders)
#' res <- residuals(fit.salamander)
#' envelope.HZIP(res, nsim = 21)
#' }
#'
#' @importFrom stats qnorm rnorm ppoints quantile
#' @importFrom ggplot2 ggplot labs scale_y_continuous stat_qq
#' @importFrom ggplot2 theme_bw theme element_text aes
#'
#' @export
envelope.HZIP <- function(object, nsim = 100, ...) {
  envelopeRQR(object, nsim)
}
