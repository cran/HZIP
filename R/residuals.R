dLGG_vec <- function(b, sigma, lambda, log = FALSE, zero = 1e-4) {
  z <- b / sigma
  if (abs(lambda) > zero) {
    lam2 <- lambda^(-2)
    res  <- log(abs(lambda)) + lam2 * log(lam2) +
      lam2 * (lambda * z - exp(lambda * z)) -
      lgamma(lam2) - log(sigma)
    if (!log) res <- exp(res)
  } else {
    res <- dnorm(z, 0, 1, log = log) / sigma
  }
  res
}

PMF_vec <- function(pi, u, Y) {
  d <- dpois(Y, lambda = u)
  ifelse(Y == 0, pi + (1 - pi) * d, (1 - pi) * d)
}

CDF_vec <- function(pi, u, Y) {
  F0 <- as.numeric(Y >= 0)
  GJ <- ppois(Y, lambda = u)
  GJ[Y < 0] <- 0
  pi * F0 + (1 - pi) * GJ
}

predictHZIP_ghq <- function(Y, w1, w2, theta1, theta2, Q = 21) {

  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Package 'statmod' is required for Gauss-Hermite quadrature. Please install it using install.packages('statmod').")
  }

  n  <- length(Y)
  gh <- statmod::gauss.quad(Q, kind = "hermite")
  nodes   <- gh$nodes
  weights <- gh$weights * exp(gh$nodes^2)

  dlgg1 <- dLGG_vec(nodes, sigma = theta1[1], lambda = theta1[1])  # comp. zero
  dlgg2 <- dLGG_vec(nodes, sigma = theta2[1], lambda = theta2[1])  # comp. count


  eta_lin <- as.numeric(w1 %*% theta1[-1])  # n
  rho_lin <- as.numeric(w2 %*% theta2[-1])  # n
  Y0 <- (Y == 0)                            # n

  num_b1 <- numeric(n)
  num_b2 <- numeric(n)
  den    <- numeric(n)

  for (k in seq_len(Q)) {
    b1k     <- nodes[k]
    wk_dlgg1 <- weights[k] * dlgg1[k]
    eta_k   <- eta_lin + b1k
    pi_k    <- -expm1(-exp(eta_k))

    for (l in seq_len(Q)) {
      b2l <- nodes[l]
      w_kl <- wk_dlgg1 * weights[l] * dlgg2[l]
      u_kl <- exp(rho_lin + b2l)

      d_pois <- dpois(Y, lambda = u_kl)

      pmf <- ifelse(Y0, pi_k + (1 - pi_k) * d_pois, (1 - pi_k) * d_pois)

      contrib <- w_kl * pmf
      den    <- den    + contrib
      num_b1 <- num_b1 + b1k * contrib
      num_b2 <- num_b2 + b2l * contrib
    }
  }

  data.frame(b1 = num_b1 / den, b2 = num_b2 / den)
}

predictHZIP_laplace <- function(Y, w1, w2, theta1, theta2) {

  n   <- length(Y)
  b1  <- numeric(n)
  b2  <- numeric(n)

  s1 <- theta1[1]; l1 <- theta1[1]
  s2 <- theta2[1]; l2 <- theta2[1]
  th_zero  <- theta1[-1]
  th_count <- theta2[-1]

  for (i in seq_len(n)) {
    Yi   <- Y[i]
    w1i  <- w1[i, , drop = TRUE]
    w2i  <- w2[i, , drop = TRUE]
    eta0 <- sum(w1i * th_zero)
    rho0 <- sum(w2i * th_count)

    neg_logpost <- function(b) {
      eta <- eta0 + b[1]
      rho <- rho0 + b[2]
      pi_h <- -expm1(-exp(eta))
      u_h  <- exp(rho)
      if (Yi == 0) {
        lpmf <- log(pi_h + (1 - pi_h) * dpois(0, u_h))
      } else {
        lpmf <- log1p(-pi_h) + dpois(Yi, u_h, log = TRUE)
      }
      lp1 <- dLGG_vec(b[1], sigma = s1, lambda = l1, log = TRUE)
      lp2 <- dLGG_vec(b[2], sigma = s2, lambda = l2, log = TRUE)
      -(lpmf + lp1 + lp2)
    }

    opt <- optim(par = c(0, 0), fn = neg_logpost, method = "BFGS",
                 control = list(reltol = 1e-8))
    b1[i] <- opt$par[1]
    b2[i] <- opt$par[2]
  }

  data.frame(b1 = b1, b2 = b2)
}

#' Residuals for HZIP Models
#'
#' Computes randomized quantile residuals for fitted
#' \code{HZIP} models.
#'
#' The residuals are obtained by integrating out the latent random
#' effects using either Gauss-Hermite quadrature or a Laplace
#' approximation.
#'
#' @param object An object of class \code{HZIP}, typically obtained
#'   from \code{\link{hzip}}.
#' @param method Character string indicating the approximation method
#'   used to estimate the posterior random effects. Possible values are:
#'   \describe{
#'     \item{\code{"ghq"}}{Gauss-Hermite quadrature approximation.}
#'     \item{\code{"laplace"}}{Laplace approximation based on posterior modes.}
#'   }
#' @param Q Integer specifying the number of Gauss-Hermite quadrature
#'   nodes when \code{method = "ghq"}. Ignored when
#'   \code{method = "laplace"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A numeric vector containing the randomized quantile residuals.
#'
#' @details
#' Let \eqn{Y_{ij}} denote the response variable for the \eqn{j}-th
#' observation of the \eqn{i}-th subject. The residuals are computed
#' using the randomized quantile residual approach proposed by
#' Dunn and Smyth (1996).
#'
#' For each observation, the conditional cumulative distribution function
#' is evaluated as
#'
#' \deqn{
#' F^*(y_{ij}) =
#' F(y_{ij} - 1) + U_{ij} P(Y_{ij} = y_{ij}),
#' }
#'
#' where \eqn{U_{ij} \sim \mathrm{Uniform}(0,1)}.
#'
#' The residuals are then obtained as
#'
#' \deqn{
#' r_{ij} = \Phi^{-1}(F^*(y_{ij})),
#' }
#'
#' where \eqn{\Phi^{-1}(\cdot)} is the standard normal quantile function.
#'
#' Random effects are approximated through either:
#' \itemize{
#'   \item Gauss-Hermite quadrature (\code{method = "ghq"});
#'   \item Laplace approximation (\code{method = "laplace"}).
#' }
#'
#' The Laplace approximation computes posterior modes using
#' \code{\link[stats]{optim}}, while the quadrature approach relies on
#' \code{\link[statmod]{gauss.quad}}.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996).
#' Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics},
#' 5(3), 236--244.
#'
#' @examples
#' \donttest{
#' fit.salamander <- hzip(
#'   y ~ mined | mined + spp,
#'   data = salamanders
#' )
#'
#' res1 <- residuals(fit.salamander, method = "ghq")
#' res2 <- residuals(fit.salamander, method = "laplace")
#'
#' head(res1)
#' head(res2)
#' }
#'
#' @importFrom stats model.response model.frame model.matrix qnorm runif dpois ppois dnorm optim
#' @importFrom Formula Formula
#'
#' @export
residuals.HZIP <- function(object, method = c("ghq", "laplace"), Q = 21, ...) {

  method <- match.arg(method)

  formula <- object$formula
  data    <- object$data
  theta1  <- c(object$scale_zero,  object$coefficients_zero)
  theta2  <- c(object$scale_count, object$coefficients_count)

  Y  <- model.response(model.frame(Formula::Formula(formula), data = data))
  w1 <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
  w2 <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)

  vB <- if (method == "ghq") {
    predictHZIP_ghq(Y, w1, w2, theta1, theta2, Q = Q)
  } else {
    predictHZIP_laplace(Y, w1, w2, theta1, theta2)
  }
  b1 <- vB[, 1]
  b2 <- vB[, 2]

  eta.hat <- as.numeric(w1 %*% theta1[-1]) + b1
  pi.hat  <- -expm1(-exp(eta.hat))

  rho.hat <- as.numeric(w2 %*% theta2[-1]) + b2
  u.hat   <- exp(rho.hat)

  n  <- length(Y)
  Fq <- CDF_vec(pi.hat, u.hat, Y - 1) + runif(n) * PMF_vec(pi.hat, u.hat, Y)
  qnorm(Fq)
}
