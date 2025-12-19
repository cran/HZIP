#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
NumericMatrix Kappas(IntegerVector y) {
  int n = y.size();
  IntegerVector ksizes(n);
  for (int i = 0; i < n; i++) {
    ksizes[i] = y[i] + 1;
  }

  int total = 1;
  for (int i = 0; i < n; i++) {
    total *= ksizes[i];
  }

  NumericMatrix out(total, n);
  for (int col = 0; col < n; col++) {
    int repeat_len = 1;
    for (int k = col + 1; k < n; k++) repeat_len *= ksizes[k];
    int block = 1;
    for (int k = 0; k < col; k++) block *= ksizes[k];
    for (int i = 0; i < block; i++) {
      for (int val = 0; val < ksizes[col]; val++) {
        for (int j = 0; j < repeat_len; j++) {
          int index = i * ksizes[col] * repeat_len + val * repeat_len + j;
          out(index, col) = val;
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
double dLGG(double b, double mu, double sigma, double lambda, bool log = false) {
  double lfunc = 0.0;

  if (lambda == 0.0) {
    lfunc = R::dnorm(b, mu, sigma, true); // log=TRUE
  } else {
    double phi = pow(lambda, -2.0);
    double c = std::abs(lambda) * pow(phi, phi) / R::gammafn(phi);
    double zz = (b - mu) / sigma;
    lfunc = std::log(c / sigma) + (zz / lambda) - phi * exp(lambda * zz);
  }

  if (log) {
    return (lfunc);
  } else {
    return exp(lfunc);
  }
}

// [[Rcpp::export]]
double dPoisLGG(NumericVector& theta2,
                 NumericMatrix& wi,
                 IntegerVector& yi,
                 NumericVector& k,
                 NumericVector& nodes,
                 NumericVector& weights){

  int mi = yi.size();
  int p = wi.ncol();

  double lambda = theta2[0];
  NumericVector beta(theta2.begin() + 1, theta2.end());

  double phi = pow(lambda, -2.0);

  NumericVector vMu(mi);
  for (int i = 0; i < mi; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) {
      xb += wi(i,j) * beta[j];
    }
    vMu[i] = exp(xb);
  }

  IntegerVector mZ(mi);
  int contador = 0;
  for (int i = 0; i < mi; i++) {
    mZ[i] = (yi[i] == 0) ? 1 : 0;
    contador += mZ[i];
  }

  double ld2 = 0.0;

  if (contador != 0 && contador != mi) {

    int Q2 = nodes.size();
    NumericVector sum1(Q2);
    for (int q = 0; q < Q2; q++) {
      double prod_log = 0.0;
      for (int j = 0; j < mi; j++) {
        double med_aux = vMu[j] * exp(nodes[q]);
        double t1 = (mZ[j] == 0) ? R::dpois(yi[j], med_aux, false) : 1.0;
        double t2 = exp(-med_aux * k[j]);
        double val = t1 * t2;
        prod_log += log(val);
      }

      double dLGG_val = dLGG(nodes[q], 0.0, lambda, lambda, true);

      double tempo1 = exp(prod_log + dLGG_val);
      sum1[q] = weights[q] * tempo1 / exp(-pow(nodes[q],2));
    }
    ld2 = std::accumulate(sum1.begin(), sum1.end(), 0.0);
  }

  if (contador == mi) {
    double temp1 = 0.0;
    for (int i = 0; i < mi; i++) temp1 += k[i] * vMu[i];
    double ld2_aux = phi*log(phi) - phi*log(temp1 + phi);
    ld2 = exp(ld2_aux);
  }

  if (contador == 0) {
    double ymas = std::accumulate(yi.begin(), yi.end(), 0.0);
    double l1 = R::lgammafn(phi + ymas) + phi*log(phi) -
      std::accumulate(yi.begin(), yi.end(), 0.0,
                      [](double s, int yi){ return s + R::lgammafn(yi + 1.0);})-
                        R::lgammafn(phi);

    for (int i = 0; i < mi; i++) l1 += yi[i] * log(vMu[i]);

    double sumkv = 0.0;
    for (int i = 0; i < mi; i++) sumkv += (k[i] + 1.0) * vMu[i];

    l1 -= (phi + ymas) * log(phi + sumkv);
    ld2 = exp(l1);
  }

  return ld2;
}

// [[Rcpp::export]]
double dPoisLGG_AGHQ(NumericVector& theta2,
                     NumericMatrix& wi,
                     IntegerVector& yi,
                     NumericVector& k,
                     NumericVector& nodes,
                     NumericVector& weights) {

  int mi = yi.size();
  int p = wi.ncol();

  double lambda = theta2[0];
  NumericVector beta(theta2.begin() + 1, theta2.end());

  double phi = pow(lambda, -2.0);

  // vMu = exp(X beta)
  NumericVector vMu(mi);
  for (int i = 0; i < mi; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) xb += wi(i, j) * beta[j];
    vMu[i] = std::exp(xb);
  }

  // mZ e contador
  IntegerVector mZ(mi);
  int contador = 0;
  for (int i = 0; i < mi; i++) {
    mZ[i] = (yi[i] == 0) ? 1 : 0;
    contador += mZ[i];
  }

  double ld2 = 0.0;

  // ===========================================================
  //  AGHQ
  // ===========================================================
  if (contador != 0 && contador != mi) {

    auto log_integrand = [&](double b) {
      double val = 0.0;
      for (int j = 0; j < mi; j++) {
        double med_aux = vMu[j] * std::exp(b);
        if (mZ[j] == 0) {
          val += R::dpois(yi[j], med_aux, true); // log
        }
        val += -med_aux * k[j];
      }
      val += dLGG(b, 0.0, lambda, lambda, true); // log density LGG
      return val;
    };

    auto deriv1 = [&](double b) {
      double h = 1e-5;
      return (log_integrand(b + h) - log_integrand(b - h)) / (2*h);
    };
    auto deriv2 = [&](double b) {
      double h = 1e-5;
      return (log_integrand(b + h) - 2*log_integrand(b) + log_integrand(b - h)) / (h*h);
    };

    double b_mode = 0.0; // initial chute
    for (int it = 0; it < 30; it++) {
      double f1 = deriv1(b_mode);
      double f2 = deriv2(b_mode);
      if (std::fabs(f2) < 1e-8) break;
      double step = f1 / f2;
      b_mode -= step;
      if (std::fabs(step) < 1e-8) break;
    }

    // Escala
    double f2_mode = deriv2(b_mode);
    double sigma = 1.0 / std::sqrt(-f2_mode);

    // Quadratura adaptativa
    int Q2 = nodes.size();
    NumericVector sum1(Q2);
    for (int q = 0; q < Q2; q++) {
      double bq = b_mode + sigma * nodes[q];
      double val = std::exp(log_integrand(bq));
      sum1[q] = weights[q] * val;
    }
    ld2 = sigma * std::accumulate(sum1.begin(), sum1.end(), 0.0);
  }

  // ===========================================================
  // CASO contador == mi
  // ===========================================================
  if (contador == mi) {
    double temp1 = 0.0;
    for (int i = 0; i < mi; i++) temp1 += k[i] * vMu[i];
    double ld2_aux = phi*std::log(phi) - phi*std::log(temp1 + phi);
    ld2 = std::exp(ld2_aux);
  }

  // ===========================================================
  // CASO contador == 0
  // ===========================================================
  if (contador == 0) {
    double ymas = std::accumulate(yi.begin(), yi.end(), 0.0);
    double l1 = R::lgammafn(phi + ymas) + phi*std::log(phi)
      - std::accumulate(yi.begin(), yi.end(), 0.0,
                        [](double s, int yi){ return s + R::lgammafn(yi + 1.0); })
      - R::lgammafn(phi);

    for (int i = 0; i < mi; i++) l1 += yi[i] * std::log(vMu[i]);

    double sumkv = 0.0;
    for (int i = 0; i < mi; i++) sumkv += (k[i] + 1.0) * vMu[i];

    l1 -= (phi + ymas) * std::log(phi + sumkv);
    ld2 = std::exp(l1);
  }

  return ld2;
}

// [[Rcpp::export]]
double dBerLGG(NumericVector& theta1,
                NumericMatrix& xi,
                IntegerVector& yi,
                NumericVector& k){

  int mi = yi.size();
  int p = xi.ncol();

  double lambda = theta1[0];
  NumericVector beta(theta1.begin() + 1, theta1.end());

  double phi = pow(lambda, -2.0);

  NumericVector ptil(mi);
  for (int i = 0; i < mi; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) xb += xi(i,j) * beta[j];
    ptil[i] = exp(xb);
  }

  IntegerVector mZ(mi);
  for (int i = 0; i < mi; i++) mZ[i] = (yi[i] == 0) ? 0 : 1;

  NumericMatrix mS = Kappas(mZ);
  int nc = mS.nrow();
  NumericVector s_mas(nc);
  for (int l = 0; l < nc; l++) {
    double soma = 0.0;
    for (int j = 0; j < mi; j++) soma += mS(l,j);
    s_mas[l] = soma;
  }

  double soma_ptil =0.0;
  for(int j=0;j<mi;j++){
    soma_ptil +=ptil[j];
  }

  NumericVector aux0(nc);
  for (int l = 0; l < nc; l++) {
    double soma_aux1 = 0.0;
    for (int j = 0; j < mi; j++) {
      soma_aux1 += ptil[j] * (mS(l,j) - mZ[j] + k[j]);
    }

    double aux2 = pow(soma_aux1 + phi + soma_ptil, -phi);
    double aux3 = ( ((int)s_mas[l] % 2 == 0) ? 1.0 : -1.0 ); // (-1)^(s.mas[l])
    aux0[l] = aux3 *aux2;
  }

  double ld = pow(phi, phi) * std::accumulate(aux0.begin(), aux0.end(), 0.0);
  return ld;
}


 // [[Rcpp::export]]
 double dZIP(NumericVector& theta1,
             NumericVector& theta2,
             NumericMatrix& xi,
             NumericMatrix& wi,
             IntegerVector& yi,
             NumericVector& Qnodes,
             NumericVector& Qweights,
             bool log = false) {

   int mi = yi.size();
   IntegerVector lZ(mi);
   for (int r = 0; r < mi; r++) {
     lZ[r] = (yi[r] == 0) ? 1 : 0;
   }

   NumericMatrix mS = Kappas(lZ);
   int nrowsS = mS.nrow();

   NumericVector logI12(nrowsS);
   for (int j = 0; j < nrowsS; j++) {
     NumericVector kj(mi);
     for (int r = 0; r < mi; r++) {
       kj[r] = mS(j,r);
     }

     double I1 = std::max(dBerLGG(theta1, xi, yi, kj), 1e-15);
     double I2 = std::max(dPoisLGG_AGHQ(theta2, wi, yi, kj, Qnodes, Qweights), 1e-15);

     logI12[j] = std::log(I1) + std::log(I2);
   }

   // log-sum-exp
   double maxlog = max(logI12);
   double sumexp = 0.0;
   for (int j = 0; j < nrowsS; j++) {
     sumexp += std::exp(logI12[j] - maxlog);
   }
   double ld = maxlog + std::log(sumexp);

   if (log) {
     return ld;
   } else {
     return std::exp(ld);
   }
 }

// [[Rcpp::export]]
double lvero(NumericVector theta, List xlist, List wlist, List ylist,
             NumericVector Qnodes, NumericVector Qweights){

  int n = ylist.size();

  double total = 0.0;
  for (int i = 0; i < n; i++) {
    NumericMatrix xi = xlist[i];
    NumericMatrix wi = wlist[i];
    int p1 = xi.ncol();
    int p2 = wi.ncol();

    IntegerVector yi = ylist[i];

    NumericVector theta1 = theta[Range(0, p1)];
    NumericVector theta2 = theta[Range(p1+1, p1+p2+1)];

    double logli = dZIP(theta1,theta2,xi,wi,yi,Qnodes,Qweights,true);
    total += logli;
  }
  return -total;
}

// [[Rcpp::export]]
double mlez_hat(NumericVector theta, List xlist, List wlist, List ylist){

  NumericMatrix x1_0 = xlist[0];
  NumericMatrix w1_0 = wlist[0];
  int p1 = x1_0.ncol();
  int p2 = w1_0.ncol();

  int n = ylist.size();
  NumericVector ld(n);

  for (int i = 0; i < n; i++) {
    NumericMatrix xi = xlist[i];
    NumericMatrix wi = wlist[i];
    IntegerVector yi = ylist[i];

    NumericVector theta1 = theta[Range(0, p1)];
    NumericVector theta2 = theta[Range(p1 + 1, p1 + p2 + 1)];

    double lambda2 = theta2[0];
    double phi2 = 1.0 / (lambda2 * lambda2);

    int mi = yi.size();

    NumericVector mu_i(mi);
    for (int r = 0; r < mi; r++) {
      double xb = 0.0;
      for (int c = 0; c < p2; c++) xb += wi(r, c) * theta2[c + 1];
      mu_i[r] = std::exp(xb);
    }

    IntegerVector lZ(mi);
    for (int r = 0; r < mi; r++) lZ[r] = (yi[r] == 0) ? 1 : 0;

    NumericMatrix mS = Kappas(lZ);
    int nrowsS = mS.nrow();

    NumericVector I12(nrowsS);
    double phi2pow = std::pow(phi2, phi2);

    for (int j = 0; j < nrowsS; j++) {
      NumericVector k_j(mi);
      for (int r = 0; r < mi; r++) k_j[r] = mS(j, r);

      double I1 = dBerLGG(theta1, xi, yi, k_j);

      double mu_plus = 0.0;
      for (int r = 0; r < mi; r++) mu_plus += mS(j, r) * mu_i[r];

      I12[j] = phi2pow * pow(mu_plus + phi2, -phi2) * I1;
    }

    ld[i] = std::accumulate(I12.begin(), I12.end(), 0.0);
  }

  double total = 0.0;
  for (int i = 0; i < n; i++) {
    if (ld[i] > 1e-15)
      total += std::log(ld[i]);
    else
      total += std::log(1e-15);
  }

  return total;
}
