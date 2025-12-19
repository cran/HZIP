dLGG <- function(b, mu=0, sigma=1, lambda, log = FALSE, zero=0.0001) {
  # generalized loggamma density
  b <- (b-mu)/sigma
  if (abs(lambda) > zero) {
    lam2 <- lambda^(-2)
    res  <- log(abs(lambda))+lam2*log(lam2)+(lam2)*(lambda*b-exp(lambda*b))-lgamma(lam2) - log(sigma)
    if (!log)
      res <- exp(res)
  } else {
    res <- dnorm(b, mean = 0, sd = 1, log = log)/sigma
  }
  return(res)
}
PMF.HZIP <- function(vB,Y,w1,w2,theta1,theta2){

  eta.hat <- w1%*%theta1[-1]+vB[1]
  pi.hat <- clogloglink(eta.hat, bvalue = NULL, inverse = TRUE,
                        deriv = 0, short = TRUE)

  rho.hat <- w2%*%theta2[-1]+vB[2]
  u.hat <- exp(rho.hat)

  if(Y==0){
    ll <- pi.hat+(1-pi.hat)*dpois(Y,lambda=u.hat,log=FALSE)
  }else{
    ll <- (1-pi.hat)*dpois(Y,lambda=u.hat,log=FALSE)
  }
  return(ll)
}
Integrate1.b1 <- function(vB,Y,w1,w2,theta1,theta2){
  aux1 <- vB[1]*PMF.HZIP(vB,Y,w1,w2,theta1,theta2)
  aux2 <- dLGG(vB[1], mu=0, sigma=theta1[1], lambda=theta1[1],log = FALSE, zero=0.0001)
  aux3 <- dLGG(vB[2], mu=0, sigma=theta2[1], lambda=theta2[1],log = FALSE, zero=0.0001)
  ll <- aux1*aux2*aux3
  return(ll)
}
Integrate1.b2 <- function(vB,Y,w1,w2,theta1,theta2){
  aux1 <- vB[2]*PMF.HZIP(vB,Y,w1,w2,theta1,theta2)
  aux2 <- dLGG(vB[1], mu=0, sigma=theta1[1], lambda=theta1[1],log = FALSE, zero=0.0001)
  aux3 <- dLGG(vB[2], mu=0, sigma=theta2[1], lambda=theta2[1],log = FALSE, zero=0.0001)
  ll <- aux1*aux2*aux3
  return(ll)
}
Integrate2 <- function(vB,Y,w1,w2,theta1,theta2){
  aux1 <- PMF.HZIP(vB,Y,w1,w2,theta1,theta2)
  aux2 <- dLGG(vB[1], mu=0, sigma=theta1[1], lambda=theta1[1],log = FALSE, zero=0.0001)
  aux3 <- dLGG(vB[2], mu=0, sigma=theta2[1], lambda=theta2[1],log = FALSE, zero=0.0001)
  ll <- aux1*aux2*aux3
  return(ll)
}
PMF <- function(pi,u,Y){
  ll <- numeric()
  n <- length(Y)
  for(ij in 1:n){
    if(Y[ij]==0){
      ll[ij] <- pi[ij]+(1-pi[ij])*dpois(Y[ij],lambda=u[ij],log=FALSE)
    }else{
      ll[ij] <- (1-pi[ij])*dpois(Y[ij],lambda=u[ij],log=FALSE)
    }
  }
  return(ll)
}
CDF <- function(pi,u,Y){
  ll <- numeric()
  n <- length(Y)
  for(ij in 1:n){
    F0 <- ifelse(Y[ij]<0,0,1)
    GJ <- ppois(Y[ij],lambda=u[ij],lower.tail=TRUE,log.p = FALSE)
    ll[ij] <- pi[ij]*F0+(1-pi[ij])*GJ
  }
  return(ll)
}
predictHZIP <- function(Y,w1,w2,theta1,theta2,lower=c(-Inf,-Inf),
                         upper=c(Inf,Inf)){
  n <- length(Y)
  b1 <- numeric()
  b2 <- numeric()
  for(i in 1:n){
    temp1 <- hcubature(Integrate1.b1,lowerLimit=lower,upperLimit=upper,tol=1e-6,
                       Y=Y[i],w1=w1[i,],
                       w2=w2[i,],theta1,theta2)$integral

    temp2 <- hcubature(Integrate1.b2,lowerLimit=lower,upperLimit=upper,tol=1e-4,
                       Y=Y[i],w1=w1[i,],
                       w2=w2[i,],theta1,theta2)$integral
    temp3 <- hcubature(Integrate2,lowerLimit=lower,upperLimit=upper,tol=1e-4,
                       Y=Y[i],w1=w1[i,],
                       w2=w2[i,],theta1,theta2)$integral

    b1[i] <- temp1/temp3
    b2[i] <- temp2/temp3
  }
  return(data.frame(b1=b1,b2=b2))
}

#' Compute Residuals for HZIP Models
#'
#' This function calculates residuals for objects of class \code{HZIP}
#' usingrandomized quantile residuals. The computation is performed efficiently
#' using C++ functions for predicting random effects and calculating
#' residuals.
#'
#' @param object An object of class \code{HZIP}, typically returned from \code{\link{hzip}}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of residuals with length equal to the total number
#'   of observations in the dataset.
#'
#' @details
#' The function internally groups the data by individual (\code{Ind}), constructs
#' model matrices for both zero-inflation and count parts of the model, and then
#' calls the C++ functions \code{predict_HZIP_cpp_vec} and \code{r_ij_cpp_vec}
#' to efficiently compute the residuals. Random effects are integrated using
#' adaptive quadrature based on the supplied \code{nodes} and \code{weights}.
#'
#' @examples
#' \donttest{
#' fit.salamander <- hzip(y ~ mined|mined+spp,data = salamanders)
#' residuals(fit.salamander)
#' }
#'
#' @importFrom stats model.response model.frame model.matrix qnorm runif dpois ppois dnorm
#' @importFrom Formula Formula
#' @importFrom VGAM clogloglink
#' @importFrom cubature hcubature
#'
#' @export
residuals.HZIP <- function(object,...){

  formula <- object$formula
  data <-object$data
  n <- object$n
  theta1 <- c(object$scale_zero,object$coefficients_zero)
  theta2 <- c(object$scale_count,object$coefficients_count)

  Y <- model.response(model.frame(Formula(formula), data = data))
  w1 <- model.matrix(Formula(formula), data = data, rhs = 1)
  w2 <- model.matrix(Formula(formula), data = data, rhs = 2)

  vB <- predictHZIP(Y,w1,w2,theta1,theta2,lower=c(-Inf,-Inf),
                     upper=c(Inf,Inf))
  b1 <- vB[,1]
  b2 <- vB[,2]

  eta.hat <- w1%*%theta1[-1]+b1
  pi.hat <- 1-exp(-exp(eta.hat))

  rho.hat <- w2%*%theta2[-1]+b2
  u.hat <- exp(rho.hat)

  E.ZIP <- (1-pi.hat)*u.hat
  Var.ZIP <- u.hat*(1-pi.hat)*(1+u.hat*pi.hat)

  n <- length(Y)

  Fq <- CDF(pi=pi.hat,u=u.hat,Y=Y-1)+
      runif(n)*PMF(pi=pi.hat,u=u.hat,Y=Y)

  rq <- qnorm(Fq)
  return(rq)
}

