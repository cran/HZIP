#' @import utils
utils::globalVariables("Ind")
NULL


.plot_hzip_test <- function(x, ...) {
  sim_df <- data.frame(stat = x$simulated)
  ggplot2::ggplot(sim_df, ggplot2::aes(x = .data$stat)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue",
                            color = "white", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = x$statistic,
                        color = "firebrick", linewidth = 1, linetype = "dashed") +
    ggplot2::labs(
      title    = x$method,
      subtitle = paste0("p-value = ", round(x$p.value, 4),
                        "  |  alternative: ", x$alternative),
      x        = paste("Simulated", x$stat_name),
      y        = "Count"
    ) +
    ggplot2::theme_bw()
}

#' Print Method for hzip_test Objects
#'
#' Prints a formatted summary of a simulation-based test returned by
#' \code{\link{testDisp.HZIP}} or \code{\link{testZI.HZIP}},
#' including the null hypothesis, observed statistic, p-value, and decision.
#'
#' @param x An object of class \code{hzip_test}.
#' @param digits Number of significant digits to print. Default is \code{4}.
#' @param ... Further arguments (currently ignored).
#'
#' @importFrom rlang .data
#' @export
print.hzip_test <- function(x, digits = 4, ...) {

  cat("\n", x$method, "\n", sep = "")
  cat(strrep("-", 50), "\n")
  cat("H0:", x$h0, "\n\n")
  cat("Observed statistic :", round(x$statistic, digits), "\n")
  cat("p-value            :", format.pval(x$p.value, digits = digits), "\n")
  cat("Alternative        :", x$alternative, "\n")
  cat("Significance level :", x$alpha, "\n")
  cat("No. simulations    :", length(x$simulated), "\n\n")

  decision <- if (x$p.value < x$alpha) {
    paste0("Reject H0 at the ", x$alpha * 100, "% significance level.")
  } else {
    paste0("Fail to reject H0 at the ", x$alpha * 100, "% significance level.")
  }
  cat("Decision:", decision, "\n\n")

  invisible(x)
}
