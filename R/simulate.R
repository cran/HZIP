#' Simulate Responses from a Fitted HZIP Model
#'
#' Generates \code{nsim} replicated response vectors from a fitted
#' \code{HZIP} model by sampling latent random effects from the
#' Generalized Log-Gamma distribution and then drawing from the
#' hurdle-Poisson component.
#'
#' @param object An object of class \code{HZIP}, typically obtained
#'   from \code{\link{hzip}}.
#' @param nsim A positive integer giving the number of response vectors
#'   to simulate. Defaults to \code{1}.
#' @param seed An optional integer seed passed to \code{\link{set.seed}}
#'   for reproducibility. If \code{NULL} (default) the current random
#'   state is used.
#' @param ... Further arguments (currently ignored).
#'
#' @return A \code{data.frame} with \code{n} rows and \code{nsim} columns,
#'   where \code{n} is the number of observations in \code{object}.
#'   Each column is one simulated response vector, named
#'   \code{sim_1}, \code{sim_2}, \ldots, \code{sim_nsim}.
#'   This follows the convention of the S3 generic
#'   \code{\link[stats]{simulate}}.
#'
#' @details
#' For each replicate and each observation \eqn{i}, the algorithm:
#' \enumerate{
#'   \item Samples \eqn{b_{1i} \sim \mathrm{LGG}(0, \sigma_1, \lambda_1)}
#'         and \eqn{b_{2i} \sim \mathrm{LGG}(0, \sigma_2, \lambda_2)}
#'         independently via \code{rgengamma()}.
#'   \item Computes
#'         \eqn{\pi_i = 1 - \exp\{-\exp(\eta_i + b_{1i})\}} and
#'         \eqn{\mu_i = \exp(\rho_i + b_{2i})}.
#'   \item Draws \eqn{Y_i}: with probability \eqn{\pi_i} sets
#'         \eqn{Y_i = 0}; otherwise draws from a zero-truncated
#'         \eqn{\mathrm{Poisson}(\mu_i)} via rejection sampling.
#' }
#'
#' @seealso \code{\link{hzip}}, \code{\link{testDisp.HZIP}},
#'   \code{\link{testZI.HZIP}}
#'
#' @examples
#' \donttest{
#' fit.salamander <- hzip(
#'   y ~ mined | mined + spp,
#'   data = salamanders
#' )
#'
#' sims <- simulate(fit.salamander, nsim = 10, seed = 42)
#' dim(sims)   # n x 10
#' head(sims)
#' }
#'
#' @importFrom stats model.frame model.matrix runif rpois simulate
#' @importFrom Formula Formula
#'
#' @method simulate HZIP
#' @export
simulate.HZIP <- function(object, nsim = 1, seed = NULL, ...) {

  if (!is.null(seed)) set.seed(seed)

  formula <- object$formula
  data    <- object$data
  theta1  <- c(object$scale_zero,  object$coefficients_zero)
  theta2  <- c(object$scale_count, object$coefficients_count)

  s1 <- theta1[1]; l1 <- theta1[1]
  s2 <- theta2[1]; l2 <- theta2[1]

  w1 <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
  w2 <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)

  eta_lin <- as.numeric(w1 %*% theta1[-1])
  rho_lin <- as.numeric(w2 %*% theta2[-1])
  n <- length(eta_lin)

  out <- replicate(nsim, {

    b1 <- rgengamma(n, mu = 0, sigma = s1, lambda = l1)
    b2 <- rgengamma(n, mu = 0, sigma = s2, lambda = l2)

    pi_i <- -expm1(-exp(eta_lin + b1))
    mu_i <- exp(rho_lin + b2)

    zero_i  <- runif(n) < pi_i
    y_count <- rpois(n, lambda = mu_i)

    retry   <- which(!zero_i & y_count == 0)
    max_try <- 100L
    iter    <- 0L
    while (length(retry) > 0 && iter < max_try) {
      y_count[retry] <- rpois(length(retry), lambda = mu_i[retry])
      retry <- retry[y_count[retry] == 0]
      iter  <- iter + 1L
    }

    ifelse(zero_i, 0L, y_count)
  })

  out <- as.data.frame(out)
  colnames(out) <- paste0("sim_", seq_len(nsim))
  out
}
