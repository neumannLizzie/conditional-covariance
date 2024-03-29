#include <covest.h>

// [[Rcpp::export(.cov_est)]]
NumericMatrix cov_est(arma::mat x, NumericVector z, double znew, double h, arma::mat cest_old, double sumK_old, arma::mat mx){
  arma::uword p = x.n_cols, n = x.n_rows, i;
  arma::mat cest(p, p), xc;
  arma::colvec K, sumK(p);

  xc = x - mx;
  K = Rcpp::dnorm(z - znew, 0, h, false);
  for (i=0; i<n; i++) cest += K[i]*xc.row(i).t()*xc.row(i);
  cest += sumK_old * cest_old;
  cest /= sum(K) + sumK_old;
  sumK.row(0) = sum(K);
  for (i = 1; i<p; i++) sumK.row(i) = NA_REAL;
return wrap(join_rows(cest, sumK));
}

// [[Rcpp::export(.mean_est)]]
NumericMatrix mean_est(arma::mat x, NumericVector z, double znew, double h, NumericVector mean_old, double sumK_old){
  arma::uword p = x.n_cols, n = x.n_rows, i;
  arma::colvec mx(p), K, sumK(p);
  arma::mat mean1(p, 2);

  K = Rcpp::dnorm(z - znew, 0, h, false);
  for (i=0; i<n; i++) mx += K[i]*x.row(i).t();
  mx += sumK_old*mean_old;
  mx /= sum(K)+sumK_old;

  sumK.row(0) = sum(K);
  for (i = 1; i<p; i++) sumK.row(i) = NA_REAL;
  return wrap(join_rows(mx, sumK));
}
