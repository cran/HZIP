#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ============================================================================
// Kappas: gera o produto cartesiano {0,1,...,y[0]} x ... x {0,1,...,y[n-1]}
// Quando recebe o vetor z (z_j = 1 se y_j = 0, 0 c.c.), gera exatamente o
// conjunto K~ do Teorema 1 (cada k_j satisfaz 0 <= k_j <= z_j).
// ============================================================================
// [[Rcpp::export]]
NumericMatrix Kappas(IntegerVector y) {
  int n = y.size();
  IntegerVector ksizes(n);
  for (int i = 0; i < n; i++) {
    ksizes[i] = y[i] + 1;
  }

  long long total = 1;
  for (int i = 0; i < n; i++) {
    total *= ksizes[i];
  }

  NumericMatrix out(total, n);
  for (int col = 0; col < n; col++) {
    long long repeat_len = 1;
    for (int k = col + 1; k < n; k++) repeat_len *= ksizes[k];
    long long block = 1;
    for (int k = 0; k < col; k++) block *= ksizes[k];
    for (long long i = 0; i < block; i++) {
      for (int val = 0; val < ksizes[col]; val++) {
        for (long long j = 0; j < repeat_len; j++) {
          long long index = i * ksizes[col] * repeat_len + val * repeat_len + j;
          out(index, col) = val;
        }
      }
    }
  }
  return out;
}

// ============================================================================
// dLGG: densidade Generalized Log-Gamma (em log se on_log = true)
// ============================================================================
// [[Rcpp::export]]
double dLGG(double b, double mu, double sigma, double lambda, bool on_log = false) {
  double lfunc = 0.0;

  if (lambda == 0.0) {
    lfunc = R::dnorm(b, mu, sigma, true);
  } else {
    double phi = std::pow(lambda, -2.0);
    // log c = log|lambda| + phi*log(phi) - lgamma(phi)
    double log_c = std::log(std::fabs(lambda)) + phi * std::log(phi) - R::lgammafn(phi);
    double zz = (b - mu) / sigma;
    lfunc = log_c - std::log(sigma) + (zz / lambda) - phi * std::exp(lambda * zz);
  }

  if (on_log) {
    return lfunc;
  } else {
    return std::exp(lfunc);
  }
}

// ============================================================================
// dPoisLGG (versão antiga, não adaptativa) — MANTIDA por compatibilidade,
// mas a recomendada é dPoisLGG_AGHQ. Esta função usa quadratura padrão
// e só é precisa quando lambda_2 e b são pequenos.
// ============================================================================
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
  if (!R_finite(lambda) || std::fabs(lambda) < 1e-8) return 0.0;
  NumericVector beta(theta2.begin() + 1, theta2.end());
  double phi = std::pow(lambda, -2.0);

  NumericVector vMu(mi);
  for (int i = 0; i < mi; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) xb += wi(i,j) * beta[j];
    vMu[i] = std::exp(xb);
  }

  IntegerVector mZ(mi);
  int contador = 0;
  for (int i = 0; i < mi; i++) {
    mZ[i] = (yi[i] == 0) ? 1 : 0;
    contador += mZ[i];
  }

  double ld2 = 0.0;

  // Caso 2(b): quadratura (não adaptativa — chama _AGHQ se quiser melhor)
  if (contador != 0 && contador != mi) {
    int Q2 = nodes.size();
    NumericVector logterms(Q2);
    double maxlog = -std::numeric_limits<double>::infinity();

    for (int q = 0; q < Q2; q++) {
      double prod_log = 0.0;
      for (int j = 0; j < mi; j++) {
        double med_aux = vMu[j] * std::exp(nodes[q]);
        if (mZ[j] == 0) prod_log += R::dpois(yi[j], med_aux, true);
        prod_log += -med_aux * k[j];
      }
      double dLGG_log = dLGG(nodes[q], 0.0, lambda, lambda, true);
      // Para Gauss-Hermite padrão, o peso w_q já inclui exp(-nodes[q]^2),
      // então multiplicamos por exp(nodes[q]^2) para "remover" essa parte
      // implícita do peso e somar o integrando completo.
      double term_log = std::log(weights[q]) + nodes[q]*nodes[q] + prod_log + dLGG_log;
      logterms[q] = term_log;
      if (term_log > maxlog) maxlog = term_log;
    }
    double sumexp = 0.0;
    for (int q = 0; q < Q2; q++) sumexp += std::exp(logterms[q] - maxlog);
    ld2 = std::exp(maxlog + std::log(sumexp));
  }

  // Caso 2(a)-i: todos y_ij = 0
  if (contador == mi) {
    double temp1 = 0.0;
    for (int i = 0; i < mi; i++) temp1 += k[i] * vMu[i];
    double ld2_aux = phi*std::log(phi) - phi*std::log(temp1 + phi);
    ld2 = std::exp(ld2_aux);
  }

  // Caso 2(a)-ii: todos y_ij > 0 (MNB)
  if (contador == 0) {
    double ymas = 0.0;
    double sum_lgamma_y1 = 0.0;
    for (int i = 0; i < mi; i++) {
      ymas += yi[i];
      sum_lgamma_y1 += R::lgammafn(yi[i] + 1.0);
    }
    double l1 = R::lgammafn(phi + ymas) + phi*std::log(phi)
      - sum_lgamma_y1 - R::lgammafn(phi);
    for (int i = 0; i < mi; i++) l1 += yi[i] * std::log(vMu[i]);
    double sumkv = 0.0;
    for (int i = 0; i < mi; i++) sumkv += (k[i] + 1.0) * vMu[i];
    l1 -= (phi + ymas) * std::log(phi + sumkv);
    ld2 = std::exp(l1);
  }

  return ld2;
}

// ============================================================================
// dPoisLGG_AGHQ: versão correta com quadratura de Gauss-Hermite ADAPTATIVA
// Implementa os três casos do Corolário 2:
//   (a)-i  : todos y_ij = 0           -> forma fechada eq. (f08)
//   (a)-ii : todos y_ij > 0 (MNB)     -> forma fechada eq. (f8)
//   (b)    : misto                    -> AGHQ eq. (fQ)
// ============================================================================
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
  // Proteção: lambda inválido -> integrando 0 (log -inf no chamador)
  if (!R_finite(lambda) || std::fabs(lambda) < 1e-8) return 0.0;
  NumericVector beta(theta2.begin() + 1, theta2.end());
  double phi = std::pow(lambda, -2.0);

  NumericVector vMu(mi);
  for (int i = 0; i < mi; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) xb += wi(i, j) * beta[j];
    vMu[i] = std::exp(xb);
  }

  IntegerVector mZ(mi);
  int contador = 0;
  for (int i = 0; i < mi; i++) {
    mZ[i] = (yi[i] == 0) ? 1 : 0;
    contador += mZ[i];
  }

  double ld2 = 0.0;

  // -------------------------------------------------------------------------
  // Caso 2(b): existe pelo menos um y_ij = 0 e um y_ij > 0  -> AGHQ
  // -------------------------------------------------------------------------
  if (contador != 0 && contador != mi) {

    // Integrando em log (sem o fator da quadratura)
    auto log_integrand = [&](double b) {
      double val = 0.0;
      for (int j = 0; j < mi; j++) {
        double med_aux = vMu[j] * std::exp(b);
        if (mZ[j] == 0) val += R::dpois(yi[j], med_aux, true);
        val += -med_aux * k[j];
      }
      val += dLGG(b, 0.0, lambda, lambda, true);
      return val;
    };

    // Derivadas numéricas centradas
    auto deriv1 = [&](double b) {
      double h = 1e-4;
      return (log_integrand(b + h) - log_integrand(b - h)) / (2.0*h);
    };
    auto deriv2 = [&](double b) {
      double h = 1e-4;
      return (log_integrand(b + h) - 2.0*log_integrand(b) + log_integrand(b - h)) / (h*h);
    };

    // Newton-Raphson para achar a moda b_mode
    double b_mode = 0.0;
    for (int it = 0; it < 50; it++) {
      double f1 = deriv1(b_mode);
      double f2 = deriv2(b_mode);
      if (!R_finite(f1) || !R_finite(f2)) break;
      if (std::fabs(f2) < 1e-10) break;
      double step = f1 / f2;
      // amortecimento simples para passos grandes
      if (std::fabs(step) > 5.0) step = (step > 0 ? 5.0 : -5.0);
      b_mode -= step;
      if (std::fabs(step) < 1e-8) break;
    }

    double f2_mode = deriv2(b_mode);
    if (!R_finite(f2_mode) || f2_mode >= 0.0) {
      // curvatura inválida: cai para sigma neutro
      f2_mode = -1.0;
    }
    double sigma_hat = 1.0 / std::sqrt(-f2_mode);

    // Quadratura adaptativa em log-escala (Gauss-Hermite "physicist"):
    //   Mudanca de variavel b = b_mode + sqrt(2)*sigma_hat*nu  =>  db = sqrt(2)*sigma_hat dnu
    //   I = int g(b) db = sqrt(2)*sigma_hat * int g(b(nu)) dnu
    //                   = sqrt(2)*sigma_hat * int e^{-nu^2} [e^{nu^2} g(b(nu))] dnu
    //   ~~ sqrt(2)*sigma_hat * sum_q w_q * e^{nu_q^2} * g(b_mode + sqrt(2)*sigma_hat*nu_q)
    // (nao ha divisao por sqrt(pi): os pesos w_q ja somam sqrt(pi); a forma f4 do paper
    //  com sigma/sqrt(pi) corresponde a uma convencao diferente onde phi(b_q) cancela esse fator)
    int Q2 = nodes.size();
    NumericVector logterms(Q2);
    double maxlog = -std::numeric_limits<double>::infinity();
    const double sqrt2 = std::sqrt(2.0);

    for (int q = 0; q < Q2; q++) {
      double bq = b_mode + sqrt2 * sigma_hat * nodes[q];
      double li = log_integrand(bq);
      // log do termo: log(w_q) + nu_q^2 + log_integrand(b_q)
      double term_log = std::log(weights[q]) + nodes[q]*nodes[q] + li;
      logterms[q] = term_log;
      if (R_finite(term_log) && term_log > maxlog) maxlog = term_log;
    }

    if (!R_finite(maxlog)) return 0.0;

    double sumexp = 0.0;
    for (int q = 0; q < Q2; q++) {
      if (R_finite(logterms[q])) sumexp += std::exp(logterms[q] - maxlog);
    }
    // log I = log(sqrt(2)*sigma_hat) + maxlog + log(sumexp)
    double log_I = 0.5*std::log(2.0) + std::log(sigma_hat) + maxlog + std::log(sumexp);
    ld2 = std::exp(log_I);
  }

  // -------------------------------------------------------------------------
  // Caso 2(a)-i: contador == mi  (todos y_ij = 0)
  // -------------------------------------------------------------------------
  if (contador == mi) {
    double temp1 = 0.0;
    for (int i = 0; i < mi; i++) temp1 += k[i] * vMu[i];
    // I_2 = phi^phi / (sum_j k_j*mu_j + phi)^phi
    double ld2_aux = phi*std::log(phi) - phi*std::log(temp1 + phi);
    ld2 = std::exp(ld2_aux);
  }

  // -------------------------------------------------------------------------
  // Caso 2(a)-ii: contador == 0  (todos y_ij > 0)  ->  MNB
  // -------------------------------------------------------------------------
  if (contador == 0) {
    double ymas = 0.0;
    double sum_lgamma_y1 = 0.0;
    for (int i = 0; i < mi; i++) {
      ymas += yi[i];
      sum_lgamma_y1 += R::lgammafn(yi[i] + 1.0);
    }
    double l1 = R::lgammafn(phi + ymas) + phi*std::log(phi)
      - sum_lgamma_y1 - R::lgammafn(phi);
    for (int i = 0; i < mi; i++) l1 += yi[i] * std::log(vMu[i]);
    double sumkv = 0.0;
    for (int i = 0; i < mi; i++) sumkv += (k[i] + 1.0) * vMu[i];
    l1 -= (phi + ymas) * std::log(phi + sumkv);
    ld2 = std::exp(l1);
  }

  return ld2;
}

// ============================================================================
// dBerLGG: forma fechada do Corolário 1 (eq. f7)
//   I_1 = phi_1^{phi_1} * sum_{s in S~} (-1)^{s+}
//                     * { sum_j ptil_j * (1 - z_j + k_j + s_j) + phi_1 }^{-phi_1}
// Implementação em log-escala com separação de termos pares/ímpares
// para lidar com a soma alternada (cancelamento catastrófico).
// ============================================================================
// [[Rcpp::export]]
double dBerLGG(NumericVector& theta1,
               NumericMatrix& xi,
               IntegerVector& yi,
               NumericVector& k){

  int mi = yi.size();
  int p = xi.ncol();

  double lambda = theta1[0];
  if (!R_finite(lambda) || std::fabs(lambda) < 1e-8) return 0.0;
  NumericVector beta(theta1.begin() + 1, theta1.end());
  double phi = std::pow(lambda, -2.0);

  // ptil_j = exp(x_{1ij}^T beta_1)
  NumericVector ptil(mi);
  for (int i = 0; i < mi; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) xb += xi(i,j) * beta[j];
    ptil[i] = std::exp(xb);
  }

  // z_j = 1 se y_j = 0, 0 c.c.  (essa convenção é a do Teorema 1)
  IntegerVector mZ(mi);
  for (int i = 0; i < mi; i++) mZ[i] = (yi[i] == 0) ? 1 : 0;

  // Geração de S~ : s_j in {0,1}, sem restrição extra além de z_j
  // (Kappas com argumento z gera tamanho 2^{soma(z)}, com s_j=0 forçado quando z_j=0;
  //  isso bate com a hipótese natural de S~ = {s : s_j=0 se y_j>0, s_j in {0,1} se y_j=0}.)
  NumericMatrix mS = Kappas(mZ);
  int nc = mS.nrow();

  NumericVector s_mas(nc);
  for (int l = 0; l < nc; l++) {
    double soma = 0.0;
    for (int j = 0; j < mi; j++) soma += mS(l,j);
    s_mas[l] = soma;
  }

  // soma_ptil = sum_j ptil_j  (representa o "+1" de (1 - z_j + k_j + s_j))
  double soma_ptil = 0.0;
  for (int j = 0; j < mi; j++) soma_ptil += ptil[j];

  // Para cada s em S~: A_l = sum_j ptil_j*(1 - z_j + k_j + s_j) + phi
  // log(termo_l) = -phi * log(A_l)
  // Sinal = (-1)^{s+}
  // Acumular pares/ímpares separados e depois subtrair na escala original.
  double maxlog = -std::numeric_limits<double>::infinity();
  NumericVector logA(nc);

  for (int l = 0; l < nc; l++) {
    double soma_aux1 = 0.0;
    for (int j = 0; j < mi; j++) {
      // (s_j - z_j + k_j); o "+1" entra via soma_ptil somado abaixo
      soma_aux1 += ptil[j] * (mS(l,j) - mZ[j] + k[j]);
    }
    double A = soma_aux1 + soma_ptil + phi;
    if (!(A > 0.0) || !R_finite(A)) {
      // valor inválido (pode ocorrer com k inconsistente); marcar como -inf
      logA[l] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    double logterm = -phi * std::log(A);
    logA[l] = logterm;
    if (R_finite(logterm) && logterm > maxlog) maxlog = logterm;
  }

  if (!R_finite(maxlog)) return 0.0;

  double sum_pos = 0.0, sum_neg = 0.0;
  for (int l = 0; l < nc; l++) {
    if (!R_finite(logA[l])) continue;
    double e = std::exp(logA[l] - maxlog);
    if (((int)s_mas[l] % 2) == 0) sum_pos += e;
    else                          sum_neg += e;
  }
  double diff = sum_pos - sum_neg;
  if (!(diff > 0.0)) {
    // Cancelamento catastrófico ou erro de sinal: retornar valor pequeno
    // (em vez de NaN, que estraga o otimizador)
    return 0.0;
  }

  // log I_1 = phi*log(phi) + maxlog + log(diff)
  double log_I1 = phi*std::log(phi) + maxlog + std::log(diff);
  return std::exp(log_I1);
}

// ============================================================================
// dZIP: log-densidade marginal do indivíduo i (Teorema 1, eq. f4)
//   f(y_i; theta) = sum_{alpha in K~} I_1(...) * I_2(...)
// Implementado via logsumexp para estabilidade.
// ============================================================================
// [[Rcpp::export]]
double dZIP(NumericVector& theta1,
            NumericVector& theta2,
            NumericMatrix& xi,
            NumericMatrix& wi,
            IntegerVector& yi,
            NumericVector& Qnodes,
            NumericVector& Qweights,
            bool on_log = false) {

  int mi = yi.size();

  // z_j = 1 se y_j = 0  -> Kappas(z) gera K~ corretamente
  IntegerVector lZ(mi);
  for (int r = 0; r < mi; r++) lZ[r] = (yi[r] == 0) ? 1 : 0;

  NumericMatrix mS = Kappas(lZ);
  int nrowsS = mS.nrow();

  NumericVector logI12(nrowsS);
  double maxlog = -std::numeric_limits<double>::infinity();

  for (int j = 0; j < nrowsS; j++) {
    NumericVector kj(mi);
    for (int r = 0; r < mi; r++) kj[r] = mS(j,r);

    double I1 = dBerLGG(theta1, xi, yi, kj);
    double I2 = dPoisLGG_AGHQ(theta2, wi, yi, kj, Qnodes, Qweights);

    if (I1 > 0.0 && I2 > 0.0 && R_finite(I1) && R_finite(I2)) {
      logI12[j] = std::log(I1) + std::log(I2);
      if (logI12[j] > maxlog) maxlog = logI12[j];
    } else {
      logI12[j] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  if (!R_finite(maxlog)) {
    // todas as combinações falharam -> log-likelihood do indivíduo = -inf
    return on_log ? -std::numeric_limits<double>::infinity() : 0.0;
  }

  double sumexp = 0.0;
  for (int j = 0; j < nrowsS; j++) {
    if (R_finite(logI12[j])) sumexp += std::exp(logI12[j] - maxlog);
  }
  double ld = maxlog + std::log(sumexp);

  if (on_log) return ld;
  else        return std::exp(ld);
}

// ============================================================================
// lvero: -log-verossimilhança total (para minimizar com optim/nlminb)
// Convenção do vetor theta:
//   theta = ( lambda_1, beta_1[0..p1-1],  lambda_2, beta_2[0..p2-1] )
// ============================================================================
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

    // theta1 = (lambda1, beta1) tem p1+1 elementos: índices 0..p1
    NumericVector theta1 = theta[Range(0, p1)];
    // theta2 = (lambda2, beta2) tem p2+1 elementos: índices p1+1..p1+p2+1
    NumericVector theta2 = theta[Range(p1+1, p1+p2+1)];

    double logli = dZIP(theta1, theta2, xi, wi, yi, Qnodes, Qweights, true);

    // Proteção: se algum indivíduo deu -Inf/NaN, retornar penalização grande
    if (!R_finite(logli)) {
      return 1e10;
    }
    total += logli;
  }
  return -total;
}

// ============================================================================
// mlez_hat: log-verossimilhanca avaliada apenas nas posicoes com y_ij = 0,
// sob a restricao de que todos os zeros observados sao atribuidos a
// componente estrutural. Equivale a aplicar a verossimilhanca completa do
// Teorema 1 ao subconjunto y_i^(0) = 0 (vetor zerado das posicoes onde y=0).
//
// Como y_i^(0) = 0 por construcao:
//   - Todos os z_ij = 1, entao K~ tem 2^{|J_0|} elementos
//   - A parte de contagem cai no Caso 2(a)-i (forma fechada eq. f08):
//       I_2(k) = phi2^phi2 * (sum_j k_j*mu_ij + phi2)^{-phi2}
//   - A parte zero usa I_1(k) do Corolario 1 (MBerGLG)
//
// Indivíduos sem zeros (todos y_ij > 0) nao contribuem.
//
// Formula (eq. 20 do paper):
//   l_z = sum_{i: |J_0(i)|>0} log{ sum_{k in K~} I_1(k|y^(0)) * I_2^{(2a-i)}(k|y^(0)) }
//
// AIC_z = 2*(p1+1) - 2*l_z;  BIC_z = (p1+1)*log(N) - 2*l_z
// ============================================================================
// [[Rcpp::export]]
double mlez_hat(NumericVector theta, List xlist, List wlist, List ylist){

  NumericMatrix x1_0 = xlist[0];
  NumericMatrix w1_0 = wlist[0];
  int p1 = x1_0.ncol();
  int p2 = w1_0.ncol();

  int n = ylist.size();
  double total = 0.0;

  for (int i = 0; i < n; i++) {
    NumericMatrix xi_full = xlist[i];
    NumericMatrix wi_full = wlist[i];
    IntegerVector yi_full = ylist[i];
    int mi_full = yi_full.size();

    // --- Filtragem: manter apenas posicoes onde y_ij = 0 ---
    int mi = 0;
    for (int r = 0; r < mi_full; r++) if (yi_full[r] == 0) mi++;

    // Indivíduo sem zeros nao contribui
    if (mi == 0) continue;

    // Constroi xi, wi, yi filtrados
    NumericMatrix xi(mi, p1);
    NumericMatrix wi(mi, p2);
    IntegerVector yi(mi);  // todos zeros por construcao
    int idx = 0;
    for (int r = 0; r < mi_full; r++) {
      if (yi_full[r] == 0) {
        for (int c = 0; c < p1; c++) xi(idx, c) = xi_full(r, c);
        for (int c = 0; c < p2; c++) wi(idx, c) = wi_full(r, c);
        yi[idx] = 0;
        idx++;
      }
    }

    // --- Extracao dos parametros ---
    NumericVector theta1 = theta[Range(0, p1)];
    NumericVector theta2 = theta[Range(p1 + 1, p1 + p2 + 1)];

    double lambda2 = theta2[0];
    if (!R_finite(lambda2) || std::fabs(lambda2) < 1e-8) {
      total += std::log(1e-300);
      continue;
    }
    double phi2 = 1.0 / (lambda2 * lambda2);

    // mu_ij filtrado
    NumericVector mu_i(mi);
    for (int r = 0; r < mi; r++) {
      double xb = 0.0;
      for (int c = 0; c < p2; c++) xb += wi(r, c) * theta2[c + 1];
      mu_i[r] = std::exp(xb);
    }

    // Como yi = 0 para todas as mi posicoes, lZ = (1, 1, ..., 1)
    // Logo K~ tem 2^mi elementos: todas as combinacoes de {0,1}^mi
    IntegerVector lZ(mi);
    for (int r = 0; r < mi; r++) lZ[r] = 1;

    NumericMatrix mS = Kappas(lZ);
    int nrowsS = mS.nrow();

    // --- Soma sobre K~ em log-escala ---
    double log_phi2pow = phi2 * std::log(phi2);
    NumericVector logterms(nrowsS);
    double maxlog = -std::numeric_limits<double>::infinity();

    for (int j = 0; j < nrowsS; j++) {
      NumericVector k_j(mi);
      for (int r = 0; r < mi; r++) k_j[r] = mS(j, r);

      // I_1(k) do Corolario 1 (MBerGLG sobre y^(0) = 0)
      double I1 = dBerLGG(theta1, xi, yi, k_j);
      if (!(I1 > 0.0) || !R_finite(I1)) {
        logterms[j] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }

      // I_2(k) do Caso 2(a)-i: phi2^phi2 * (sum_j k_j*mu_ij + phi2)^{-phi2}
      double mu_plus = 0.0;
      for (int r = 0; r < mi; r++) mu_plus += mS(j, r) * mu_i[r];

      // log[ I_1(k) * I_2(k) ]
      double term = log_phi2pow - phi2 * std::log(mu_plus + phi2) + std::log(I1);
      logterms[j] = term;
      if (R_finite(term) && term > maxlog) maxlog = term;
    }

    // logsumexp
    if (!R_finite(maxlog)) {
      total += std::log(1e-300);
      continue;
    }
    double sumexp = 0.0;
    for (int j = 0; j < nrowsS; j++) {
      if (R_finite(logterms[j])) sumexp += std::exp(logterms[j] - maxlog);
    }
    double log_li = maxlog + std::log(sumexp);
    total += log_li;
  }

  return total;
}
