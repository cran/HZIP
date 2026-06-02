#' Fit a Hierarchical Zero-Inflated Poisson (HZIP) Model
#'
#' \code{hzip()} fits a longitudinal/clustered zero-inflated Poisson model with
#' subject-level random effects following a Generalized Log-Gamma (GLG)
#' distribution, by maximizing a marginal likelihood approximated via adaptive
#' Gauss-Hermite quadrature. The model uses a two-part \link[Formula]{Formula}:
#' \eqn{y \sim \text{zero part} \mid \text{count part}}, where the count
#' intensity (Poisson mean) and the zero-inflation probability are linked to
#' (possibly different) sets of covariates. Initial values are obtained from
#' \code{pscl::zeroinfl(..., dist = "poisson", link = "cloglog")}.
#'
#' @param formula A two-part \link[Formula]{Formula} of the form
#'   \code{y ~ w_zero + ... | x_count + ...}, where the right-hand side before
#'   the bar specifies covariates for the zero-inflation component and the
#'   right-hand side after the bar specifies covariates for the Poisson mean.
#' @param data A \code{data.frame} containing all variables used in
#'   \code{formula} and a subject identifier named \code{Ind} (one row per
#'   observation).
#' @param hessian Logical; if \code{TRUE} (default) the observed Hessian at the
#'   optimum is computed and used for standard-error estimates.
#' @param method Character string passed to \code{\link[stats]{optim}} for the
#'   first optimization stage (default \code{"BFGS"}).
#' @param Q Integer; number of Gauss-Hermite nodes for quadrature
#'   (default \code{15}). Larger values improve accuracy at higher computational
#'   cost.
#' @param lower Bounds on the variables for the \code{"L-BFGS-B"} method
#'   (passed to \code{\link[stats]{optim}}).
#' @param upper Bounds on the variables for the \code{"L-BFGS-B"} method
#'   (passed to \code{\link[stats]{optim}}).
#' @param control Optional \code{list} passed to \code{\link[stats]{optim}}'s
#'   \code{control=} argument (e.g., \code{list(maxit = 500)}).
#' @param start Optional numeric vector of initial values of length
#'   \eqn{p_1 + p_2 + 2} in the order
#'   \code{c(scale_zero, beta_zero, scale_count, beta_count)}. If \code{NULL}
#'   (default), initial values are obtained from \code{pscl::zeroinfl} and the
#'   scales are set to \code{0.5}.
#' @param refine Logical; if \code{TRUE}, the BFGS optimum is further refined
#'   by \code{\link[stats]{nlminb}} (a trust-region Newton method). Useful as
#'   a safeguard for difficult problems but typically does not improve the
#'   solution when BFGS already converges. Default \code{FALSE} (BFGS only).
#' @param ep_method Character; method used to compute the Hessian for the
#'   standard errors. Either \code{"optim"} (default; uses the finite-difference
#'   Hessian computed by \code{\link[stats]{optim}}, fast) or \code{"numDeriv"}
#'   (uses \code{\link[numDeriv]{hessian}} with Richardson extrapolation, more
#'   accurate but slower). For most applications the two give virtually identical
#'   standard errors; \code{"numDeriv"} is recommended only when extra accuracy
#'   is needed.
#' @param n_starts Integer; if greater than \code{1}, performs a multi-start
#'   procedure with \code{n_starts} perturbations of the initial values and
#'   keeps the best log-likelihood. Useful for difficult problems with possible
#'   local optima. Default \code{1} (no multi-start).
#' @param verbose Logical; if \code{TRUE}, prints progress information during
#'   the optimization. Default \code{FALSE}.
#' @param ... Further arguments passed to \code{\link[stats]{optim}}.
#'
#' @details
#' Let \eqn{y_{ij}} denote the count response for subject \eqn{i} at occasion
#' \eqn{j}. The HZIP model assumes
#' \deqn{P(y_{ij}=0 \mid u_i) = \pi_{ij}(u_i) + \{1-\pi_{ij}(u_i)\}\exp\{-\mu_{ij}(u_i)\},}
#' \deqn{P(y_{ij}=k \mid u_i) = \{1-\pi_{ij}(u_i)\}\frac{\mu_{ij}(u_i)^k e^{-\mu_{ij}(u_i)}}{k!},\quad k\ge 1,}
#' with linear predictors for the count and zero parts (links typically
#' \code{log} for the Poisson mean and \code{cloglog} for the zero-inflation).
#' Subject-specific random effects \eqn{u_i} follow a GLG distribution and
#' induce within-subject dependence; the marginal likelihood is approximated by
#' adaptive Gauss-Hermite quadrature with \code{Q} nodes.
#'
#' \strong{Optimization strategy.} The default workflow is:
#' \enumerate{
#'   \item Generate initial values from \code{pscl::zeroinfl} (with scales set
#'         to a neutral \code{0.5}, ensuring reproducibility);
#'   \item Run \code{\link[stats]{optim}} with the specified \code{method}
#'         (default BFGS);
#'   \item If \code{refine = TRUE}, refine the result with
#'         \code{\link[stats]{nlminb}} (off by default; useful for
#'         difficult problems);
#'   \item If \code{n_starts > 1}, repeat steps 2-3 with perturbed initial
#'         values and keep the best log-likelihood;
#'   \item Compute standard errors from the inverse Hessian using the method
#'         specified by \code{ep_method}.
#' }
#'
#' @return An object of class \code{"HZIP"}, a \code{list} with elements:
#' \item{call}{The matched call.}
#' \item{formula}{The model \code{Formula}.}
#' \item{coefficients_zero}{Estimated coefficients for the zero-inflation part.}
#' \item{coefficients_count}{Estimated coefficients for the count part.}
#' \item{scale_zero}{Estimated scale (zero part).}
#' \item{scale_count}{Estimated scale (count part).}
#' \item{loglik}{Maximized log-likelihood.}
#' \item{loglikz}{Auxiliary log-likelihood from \code{mlez_hat}.}
#' \item{AIC, BIC}{Information criteria based on \code{loglik}.}
#' \item{AICz, BICz}{Information criteria based on \code{loglikz}.}
#' \item{convergence}{Convergence code (0 = success).}
#' \item{n}{Number of observations.}
#' \item{m}{Cluster sizes per subject (vector ordered by \code{Ind}).}
#' \item{ep}{Approximate standard errors (square roots of the diagonal of the
#'   inverse Hessian).}
#' \item{iter}{Number of optimization iterations.}
#' \item{method}{Optimization method.}
#' \item{refine}{Whether \code{nlminb} refinement was applied.}
#' \item{ep_method}{Method used to compute standard errors.}
#' \item{n_starts}{Number of starting points used.}
#' \item{optim}{Raw output from the final optimization stage.}
#' \item{data}{The input \code{data}.}
#'
#' @note
#' The subject identifier must be named \code{Ind}. The internal parameter
#' vector follows the order
#' \code{c(scale_zero, beta_zero, scale_count, beta_count)}.
#'
#' For difficult problems (e.g., models with many interaction terms), if the
#' default fit produces large standard errors or the AIC suggests a worse fit
#' than a nested simpler model, try (i) \code{n_starts = 10} for multi-start
#' optimization and (ii) supplying \code{start = ...} from a fit of the simpler
#' nested model.
#'
#' @references
#' Min, Y., & Agresti, A. (2005). Random effect models for repeated measures of
#' zero-inflated count data. \emph{Statistical Modelling}, 5(1), 1-19.
#'
#' Jackman, S. (2020). \emph{pscl}: Classes and Methods for R Developed in the
#' Political Science Computational Laboratory. R package version 1.5.5.
#'
#' Zeileis, A., & Croissant, Y. (2010). Extended model formulas in R:
#' \emph{Journal of Statistical Software}, 34(1), 1-13. (\pkg{Formula})
#'
#' @examples
#' \donttest{
#' # Basic fit
#' fit.salamander <- hzip(y ~ mined | mined + spp, data = salamanders)
#' summary(fit.salamander)
#'
#' # More robust fit for difficult problems (multi-start)
#' fit.salamander2 <- hzip(y ~ mined | mined + spp + mined:spp,
#'                         data = salamanders,
#'                         n_starts = 10, verbose = TRUE)
#' summary(fit.salamander2)
#'
#' # Using the simpler model as starting point for the harder one
#' init <- with(fit.salamander, c(scale_zero, coefficients_zero,
#'                                scale_count, coefficients_count,
#'                                rep(0, 6)))  # 6 zeros for the interaction
#' fit.salamander3 <- hzip(y ~ mined | mined + spp + mined:spp,
#'                         data = salamanders, start = init)
#' summary(fit.salamander3)
#'
#' }
#'
#' @importFrom dplyr group_split group_by mutate
#' @importFrom stats model.frame model.matrix model.response optim nlminb pnorm as.formula rnorm
#' @importFrom pscl zeroinfl
#' @importFrom statmod gauss.quad
#' @importFrom numDeriv hessian
#' @import Formula
#' @export
hzip <- function(formula, data, hessian = TRUE, method = "BFGS",
                 Q = 15, lower = -Inf, upper = Inf,
                 control = NULL, start = NULL,
                 refine = FALSE, ep_method = c("optim", "numDeriv"),
                 n_starts = 1, verbose = FALSE, ...) {

  ep_method <- match.arg(ep_method)

  if (!"Ind" %in% names(data)) stop("data must contain 'Ind'")
  if (length(unique(data$Ind)) != max(as.integer(factor(data$Ind)))) {
    data <- dplyr::mutate(data, Ind = as.integer(factor(Ind)))
  }
  data_list <- dplyr::group_split(dplyr::group_by(data, Ind))

  fForm <- Formula::Formula(formula)
  xlist <- lapply(data_list, function(df) model.matrix(fForm, df, rhs = 1))
  wlist <- lapply(data_list, function(df) model.matrix(fForm, df, rhs = 2))
  ylist <- lapply(data_list, function(df) model.response(model.frame(fForm, df)))

  p1 <- ncol(xlist[[1]])
  p2 <- ncol(wlist[[1]])

  # ---- Gauss-Hermite ----
  QGauss   <- statmod::gauss.quad(Q, kind = "hermite")
  Qnodes   <- QGauss$nodes
  Qweights <- QGauss$weights
  # --------

  if (is.null(start)) {
    lhs  <- formula(fForm, lhs = 1, rhs = 0)
    rhs1 <- formula(fForm, lhs = 0, rhs = 1)
    rhs2 <- formula(fForm, lhs = 0, rhs = 2)
    fAux <- paste(deparse(lhs[[2]]), "~",
                  deparse(rhs2[[2]]), "|",
                  deparse(rhs1[[2]]))
    fit.aux <- pscl::zeroinfl(Formula::Formula(as.formula(fAux)),
                              data = data, dist = "poisson", link = "cloglog")
    initial <- c(0.5, as.numeric(fit.aux$coefficients$zero),
                 0.5, as.numeric(fit.aux$coefficients$count))
  } else {
    if (length(start) != p1 + p2 + 2) {
      stop(sprintf("'start' must have length %d (scale_zero, beta_zero[%d], scale_count, beta_count[%d])",
                   p1 + p2 + 2, p1, p2))
    }
    initial <- as.numeric(start)
  }
  if (verbose) {
    cat(sprintf("[hzip] Starting point (length=%d):\n", length(initial)))
    print(round(initial, 4))
  }

  control_user <- if (is.null(control)) list() else control

  ask_hess_inline <- (ep_method == "optim") && (!refine) && (n_starts == 1) && hessian

  fit_one <- function(init) {
    op <- tryCatch({
      optim(par = init, fn = lvero,
            ylist = ylist, xlist = xlist, wlist = wlist,
            Qnodes = Qnodes, Qweights = Qweights,
            method = method, hessian = ask_hess_inline,
            control = control_user, lower = lower, upper = upper, ...)
    }, error = function(e) NULL)
    if (is.null(op)) return(list(par = init, value = Inf, convergence = 99,
                                 counts = c(NA, NA), hessian = NULL))

    if (refine) {
      val_before <- op$value
      ref <- tryCatch({
        nlminb(start = op$par, objective = lvero,
               ylist = ylist, xlist = xlist, wlist = wlist,
               Qnodes = Qnodes, Qweights = Qweights,
               control = list(rel.tol = 1e-12, x.tol = 1e-10, iter.max = 500))
      }, error = function(e) NULL)
      if (!is.null(ref) && ref$objective < op$value) {
        op$par <- ref$par
        op$value <- ref$objective
        op$convergence <- ref$convergence
        op$counts <- c(op$counts[1] + ref$iterations, NA)
        if (verbose) cat(sprintf("[hzip] nlminb refined: %.4f -> %.4f\n",
                                 -val_before, -op$value))
      }
    }
    op
  }

  if (n_starts > 1) {
    if (verbose) cat(sprintf("[hzip] Multi-start with %d points...\n", n_starts))
    set.seed(1)
    fits <- vector("list", n_starts)
    fits[[1]] <- fit_one(initial)
    for (k in seq_len(n_starts - 1) + 1) {
      perturb <- initial + rnorm(length(initial), sd = 0.3)
      perturb[1]      <- abs(perturb[1])      + 1e-3
      perturb[p1 + 2] <- abs(perturb[p1 + 2]) + 1e-3
      fits[[k]] <- fit_one(perturb)
      if (verbose) cat(sprintf("  start %d: logLik = %.4f\n", k, -fits[[k]]$value))
    }
    best <- which.min(sapply(fits, function(f) f$value))
    op <- fits[[best]]
    if (verbose) cat(sprintf("[hzip] Best start: #%d (logLik = %.4f)\n",
                             best, -op$value))
  } else {
    op <- fit_one(initial)
  }

  if (verbose) {
    cat(sprintf("[hzip] Final logLik: %.4f (convergence=%d)\n",
                -op$value, op$convergence))
  }

  ep <- rep(NA_real_, length(op$par))
  if (hessian) {
    H <- NULL
    if (ask_hess_inline && !is.null(op$hessian)) {
      H <- op$hessian
    }
    if (is.null(H)) {
      if (ep_method == "numDeriv") {
        H <- tryCatch({
          numDeriv::hessian(func = lvero, x = op$par,
                            ylist = ylist, xlist = xlist, wlist = wlist,
                            Qnodes = Qnodes, Qweights = Qweights,
                            method = "Richardson")
        }, error = function(e) NULL)
      } else {
        op_h <- tryCatch({
          optim(par = op$par, fn = lvero,
                ylist = ylist, xlist = xlist, wlist = wlist,
                Qnodes = Qnodes, Qweights = Qweights,
                method = method, hessian = TRUE,
                control = list(maxit = 1), lower = lower, upper = upper)
        }, error = function(e) NULL)
        H <- if (!is.null(op_h)) op_h$hessian else NULL
      }
    }

    if (!is.null(H) && all(is.finite(H))) {
      ep <- tryCatch({
        Hinv <- solve(H)
        diag_v <- diag(Hinv)
        if (any(diag_v < 0)) {
          warning("Hessian not positive-definite; some standard errors will be NA.")
          diag_v[diag_v < 0] <- NA_real_
        }
        sqrt(diag_v)
      }, error = function(e) {
        warning("Hessian not invertible; standard errors set to NA.")
        rep(NA_real_, length(op$par))
      })
    } else {
      warning("Could not compute Hessian.")
    }
  }

  loglikz <- mlez_hat(op$par, xlist, wlist, ylist)
  npar    <- length(op$par)        #  p1 + p2 + 2
  npar_z  <- p1 + 1                # ell_z: beta_zero (p1) + lambda_1 (1)

  fit.hzip <- list(
    call = match.call(),
    formula = formula,
    coefficients_zero  = op$par[2:(p1 + 1)],
    coefficients_count = op$par[(p1 + 3):(p1 + p2 + 2)],
    scale_zero  = op$par[1],
    scale_count = op$par[p1 + 2],
    loglik  = -op$value,
    loglikz = loglikz,
    AIC  = -2 * (-op$value) + 2 * npar,
    BIC  = -2 * (-op$value) + log(nrow(data)) * npar,
    AICz = -2 * loglikz     + 2 * npar_z,
    BICz = -2 * loglikz     + log(nrow(data)) * npar_z,
    convergence = op$convergence,
    n   = as.numeric(max(data$Ind)), #nrow(data),
    m   = as.numeric(table(data$Ind)),
    ep  = ep,
    iter = op$counts[1],
    method = method,
    refine = refine,
    ep_method = ep_method,
    n_starts = n_starts,
    optim = op,
    data = data
  )
  class(fit.hzip) <- "HZIP"
  return(fit.hzip)
}

#' @export
print.HZIP <- function(x, ...) {
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients (Zero part):\n"); print(round(x$coefficients_zero, 4))
  cat("\nCoefficients (Count part):\n"); print(round(x$coefficients_count, 4))
  cat("\nZero scale:\n"); print(round(x$scale_zero, 4))
  cat("\nCount scale:\n"); print(round(x$scale_count, 4))
  cat("\nlogLik:", round(x$loglik, 4),
      "  (conv =", x$convergence, ", iter =", x$iter, ")\n")
  invisible(x)
}


#' @export
summary.HZIP <- function(object, ...) {

  fForm <- Formula::Formula(object$formula)
  data  <- eval(object$call$data)

  y <- model.response(model.frame(fForm, data = data))
  x <- model.matrix(fForm,  data = data, rhs = 1)
  w <- model.matrix(fForm,  data = data, rhs = 2)

  coef_zero    <- object$coefficients_zero
  coef_count   <- object$coefficients_count
  scale_zero   <- object$scale_zero
  scale_count  <- object$scale_count
  ep           <- object$ep
  iter         <- object$iter
  loglik       <- object$loglik
  loglikz      <- object$loglikz
  n            <- object$n
  convergence  <- object$convergence

  p1 <- ncol(x)
  p2 <- ncol(w)

  # Layout do vetor de parametros: (lambda1, beta_zero[p1], lambda2, beta_count[p2])
  if (!is.null(ep) && length(ep) == (p1 + p2 + 2)) {
    std_scale_zero  <- ep[1]
    std_coeff_zero  <- ep[2:(p1 + 1)]
    std_scale_count <- ep[(p1 + 2)]
    std_coeff_count <- ep[(p1 + 3):(p1 + p2 + 2)]
  } else {
    std_scale_zero  <- NA_real_
    std_coeff_zero  <- rep(NA_real_, length(coef_zero))
    std_scale_count <- NA_real_
    std_coeff_count <- rep(NA_real_, length(coef_count))
  }

  # ----Wald ----
  z_beta_zero  <- coef_zero  / std_coeff_zero
  p_beta_zero  <- 2 * (1 - pnorm(abs(z_beta_zero)))
  z_eta_count  <- coef_count / std_coeff_count
  p_eta_count  <- 2 * (1 - pnorm(abs(z_eta_count)))

  # ---- Coeficients ----
  df_coef_zero <- data.frame(
    Estimate     = coef_zero,
    `Std. Error` = std_coeff_zero,
    `z value`    = z_beta_zero,
    `Pr(>|z|)`   = p_beta_zero,
    row.names    = colnames(x),
    check.names  = FALSE
  )

  df_coef_count <- data.frame(
    Estimate     = coef_count,
    `Std. Error` = std_coeff_count,
    `z value`    = z_eta_count,
    `Pr(>|z|)`   = p_eta_count,
    row.names    = colnames(w),
    check.names  = FALSE
  )

  df_coef_scale_zero <- data.frame(
    Estimate  = scale_zero,
    Std.Error = std_scale_zero,
    row.names = "scale"
  )

  df_coef_scale_count <- data.frame(
    Estimate  = scale_count,
    Std.Error = std_scale_count,
    row.names = "scale"
  )

  # ---- Exit ----
  out <- list(
    call         = object$call,
    loglik       = loglik,
    loglikz      = loglikz,
    AIC          = object$AIC,
    BIC          = object$BIC,
    AICz         = object$AICz,
    BICz         = object$BICz,
    AIC_HZIP     = object$AIC  - object$AICz,
    BIC_HZIP     = object$BIC  - object$BICz,
    iter         = iter,
    coef_zero    = df_coef_zero,
    coef_count   = df_coef_count,
    scale_zero   = df_coef_scale_zero,
    scale_count  = df_coef_scale_count,
    convergence  = convergence
  )

  class(out) <- "summary.HZIP"
  return(out)
}

#' @importFrom stats printCoefmat
#' @export
print.summary.HZIP <- function(x, digits = 5, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\n")

  w <- 14L

  center <- function(s, width = w) {
    s <- as.character(s)
    pad <- max(0L, width - nchar(s))
    left  <- pad %/% 2L
    right <- pad - left
    paste0(strrep(" ", left), s, strrep(" ", right))
  }

  fmt_num <- function(v) center(formatC(v, digits = digits, format = "f"))
  fmt_int <- function(v) center(formatC(v, format = "d"))
  fmt_lab <- function(s) center(s)

  labs1 <- c("logLik", "AIC", "BIC", "Iter")
  vals1 <- c(fmt_num(x$loglik), fmt_num(x$AIC),
             fmt_num(x$BIC),    fmt_int(x$iter))
  cat(paste0(vapply(labs1, fmt_lab, character(1)), collapse = ""), "\n", sep = "")
  cat(paste0(vals1, collapse = ""), "\n", sep = "")

  labs2 <- c("AICz", "BICz", "AIC_HZIP", "BIC_HZIP")
  vals2 <- c(fmt_num(x$AICz), fmt_num(x$BICz),
             fmt_num(x$AIC_HZIP), fmt_num(x$BIC_HZIP))
  cat(paste0(vapply(labs2, fmt_lab, character(1)), collapse = ""), "\n", sep = "")
  cat(paste0(vals2, collapse = ""), "\n", sep = "")

  cat("Convergence:", x$convergence, "\n\n")

  # ---- Coeficients ----
  cat("Coefficients (Zero part):\n")
  printCoefmat(x$coef_zero,
               digits = digits,
               signif.stars = TRUE,
               na.print = "NA")

  cat("\nCoefficients (Count part):\n")
  printCoefmat(x$coef_count,
               digits = digits,
               signif.stars = TRUE,
               na.print = "NA")

  cat("\nZero Scale:\n")
  print(x$scale_zero,  digits = digits)

  cat("\nCount Scale:\n")
  print(x$scale_count, digits = digits)

  invisible(x)
}
