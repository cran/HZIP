# rgengamma (interno)
#' @keywords internal
#' @noRd
#' @importFrom stats rgamma rnorm
rgengamma <- function(n, mu=0, sigma=1, lambda, zero=0.0001) {
  # generates n random generalized loggamma variates
  if (abs(lambda) > zero) {
    alpha <- 1/lambda^2
    ei    <- rgamma(n,shape=alpha,rate=1)
    ui    <- (log(ei)-log(alpha))/lambda
  } else {
    ui <- rnorm(n,mean=0,sd=1)
  }
  ui <- ui*sigma + mu
  return(ui)
}

#' Simulate Data from a Hierarchical Zero-Inflated Poisson (HZIP) Model
#'
#' \code{rHZIP()} generates panel/longitudinal data from a two–part
#' hierarchical zero-inflated Poisson model with subject-specific random
#' effects for both the zero-inflation and the count components. Random
#' effects are drawn from a generalized log-gamma (GLG) distribution via
#' \code{rgengamma()} (user must provide/attach a function with this name; see
#' \strong{Dependencies}).
#'
#' @param n Integer. Number of subjects.
#' @param m Integer vector of length \code{n} (or a scalar recycled to length
#'   \code{n}). Numbers of repeated measurements per subject.
#' @param para Numeric vector of parameters in the order
#'   \code{c(lambda1, beta1, lambda2, beta2)}, where:
#'   \itemize{
#'     \item \code{lambda1}: GLG scale/shape parameter for the zero part.
#'     \item \code{beta1}:  length-\code{p1} vector of coefficients for the zero part
#'           (matches \code{ncol(x1)}).
#'     \item \code{lambda2}: GLG scale/shape parameter for the count part.
#'     \item \code{beta2}:  length-\code{p2} vector of coefficients for the count part
#'           (matches \code{ncol(x2)}).
#'   }
#'   Internally, the linear predictors are
#'   \eqn{\eta^{(0)}_{ij} = x^{\top}_{1,ij}\beta_1 + b^{(0)}_i} and
#'   \eqn{\eta^{(1)}_{ij} = x^{\top}_{2,ij}\beta_2 + b^{(1)}_i}, with
#'   \eqn{p_{ij} = 1 - \exp\{-\exp(\eta^{(0)}_{ij})\}} (cloglog link for zero-inflation)
#'   and \eqn{\mu_{ij} = \exp(\eta^{(1)}_{ij})} (log link for counts).
#' @param x1 Numeric matrix of covariates for the zero-inflation part
#'   (dimension \code{sum(m)} by \code{p1}). Include an intercept column if desired.
#' @param x2 Numeric matrix of covariates for the count part
#'   (dimension \code{sum(m)} by \code{p2}). Include an intercept column if desired.
#'
#' @details
#' For subject \eqn{i=1,\dots,n} with \eqn{m_i} observations, the model draws
#' two subject-level random effects \eqn{b^{(0)}_i} and \eqn{b^{(1)}_i} independently
#' from GLG distributions parameterized by \code{lambda1} and \code{lambda2}.
#' Conditional on these effects, outcomes are generated as
#' \deqn{Y_{ij} = Z_{ij}\times C_{ij},}
#' where \eqn{Z_{ij}\sim \mathrm{Bernoulli}(p_{ij})} controls structural zeros and
#' \eqn{C_{ij}\sim \mathrm{Poisson}(\mu_{ij})} controls the count size.
#'
#' The returned data frame contains subject IDs, the response \code{y}, and the
#' supplied covariates. An attribute \code{"propZeros"} stores a small summary
#' with the percentage of structural zeros, additional Poisson zeros, and total
#' zeros.
#'
#' @return A \code{tibble} with columns:
#' \item{Ind}{Subject identifier (1..n), repeated according to \code{m}.}
#' \item{y}{Simulated response.}
#' \item{x*, w*}{The covariates from \code{x1} and \code{x2} (renamed as described below).}
#'
#' The object has an attribute \code{"propZeros"}: a 3×1 \code{data.frame} with
#' rows \code{"Zeros"} (structural zeros, \%),
#' \code{"Count"} (extra zeros from the Poisson part, \%), and
#' \code{"Total"} (overall zero percentage).
#'
#' @section Column naming:
#' Columns of \code{x1} are renamed to \code{x1, x2, ..., xp1}.
#' Columns of \code{x2} are copied except the first column is dropped and the
#' remaining are renamed \code{w1, w2, ..., w_{p2-1}}. (This mirrors the current
#' implementation that excludes the first \code{x2} column from the output.)
#'
#' @section Dependencies:
#' This function calls \code{rgengamma()} to draw GLG random effects. Ensure that
#' such a function is available on the search path (e.g., from a package that
#' provides a generalized log-gamma RNG) or provide your own implementation
#' with the signature \code{rgengamma(n, mu, sigma, lambda)}. It also uses
#' \pkg{dplyr} and \pkg{tibble}.
#'
#' @examples
#' set.seed(123)
#'
#' n <- 50
#' m <- rep(4, n)
#' N <- sum(m)
#'
#' # design matrices (with intercepts)
#' x1 <- cbind(1, rnorm(N))
#' x2 <- cbind(1, rnorm(N), rbinom(N, 1, 0.5))
#'
#' p1 <- ncol(x1); p2 <- ncol(x2)
#' lambda1 <- 0.7
#' beta1   <- c(-0.2, 0.6)
#' lambda2 <- 0.9
#' beta2   <- c( 0.3, 0.5, -0.4)
#'
#' para <- c(lambda1, beta1, lambda2, beta2)
#'
#' sim <- rHZIP(n, m, para, x1, x2)
#' head(sim)
#' attr(sim, "propZeros")
#'
#' @seealso \code{\link{hzip}} for model fitting.
#'
# rHZIP
#' @importFrom dplyr filter
#' @importFrom tibble tibble
#' @importFrom stats rbinom rpois
#' @export

rHZIP <- function(n,m,para,x1,x2){
  p1 <- ncol(x1)
  p2 <- ncol(x2)

  lambda1 <- para[1]
  beta1 <- para[2:(p1+1)]

  lambda2 <- para[(p1+2)]
  beta2 <- para[(p1+3):(p1+p2+2)]

  ind <- rep(seq_len(n), times = m)

  sample1 <- data.frame(x1,ind)
  sample2 <- data.frame(x2,ind)

  y <- vector()
  conta1 <- vector()
  conta2 <- vector()
  i <- 1
  while(i <=n){
    s.aux1 <- filter(sample1,ind==i)
    b1 <-  rgengamma(n = 1, mu = 0, sigma = lambda1, lambda =lambda1)
    eta1 <- as.matrix(s.aux1[,1:ncol(x1)])%*%beta1+b1
    p <- 1-exp(-exp(eta1))

    s.aux2 <- filter(sample2,ind==i)
    b2 <- rgengamma(n = 1, mu = 0, sigma = lambda2, lambda =lambda2)
    eta2 <- as.matrix(s.aux2[,1:ncol(x2)])%*%beta2+b2
    mu <- exp(eta2)

    y.auxi1 <- rbinom(m[i],1,p)
    y.auxi2 <- rpois(m[i],mu)

    conta1 <- c(conta1,y.auxi1)
    conta2 <- c(conta2,y.auxi2)

    y <- c(y,y.auxi1*y.auxi2) #rZIP(mi,mu=mu,sigma=p)

    i <- i+1
  }

  colnames(x1) <- paste0("x", seq_len(ncol(x1)))
  xx2 <-x2[,-1]
  colnames(xx2) <- paste0("w", seq_len(ncol(x2)-1))


  out <- tibble(
    Ind = ind,
    y = y,
    !!!as.data.frame(x1),
    !!!as.data.frame(xx2)
  )

  # Proportions

  pConta1 <- 100*(1-mean(conta1))
  pConta2 <- 100*mean(ifelse(conta2==0,1,0))
  pConta3 <- 100*mean(ifelse(y==0,1,0))

  propZeros <- data.frame(rbind(pConta1,(pConta3-pConta1),pConta3))
  colnames(propZeros) <- "Zero-proportions(%)"
  rownames(propZeros) <- c("Zeros","Count","Total")
  propZeros

  attr(out, "propZeros")  <- propZeros

  return(out)
}
