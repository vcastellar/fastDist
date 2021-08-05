#include <RcppArmadillo.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix euclidean(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), 
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  
  arma::colvec An =  sum(square(A),1);
  arma::colvec Bn =  sum(square(B),1);
  
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  return wrap(sqrt(C)); 
}

// [[Rcpp::export]]
NumericMatrix manhattan(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), 
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  arma::mat C = arma::mat(Ar.begin(), m, k, false);
  arma::mat res = arma::mat(m, m, arma::fill::zeros);
  
  for (int i = 0; i < m; i++) {
    C.each_row() -= B.row(i);
    res.col(i) = sum(abs(C), 1);
  }

  return wrap(res); 
}

// [[Rcpp::export]]
NumericMatrix minkowsky(NumericMatrix Ar, NumericMatrix Br, int p) {
  int m = Ar.nrow(), 
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  arma::mat C = arma::mat(m, n, arma::fill::zeros);
  
  double q = 1.0 / p;
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      arma::mat aux = pow(abs(A.row(i) - B.row(j)), p);
      
      arma::vec res = aux.t();
      C(i, j) = pow(sum(res), q);

    }
  }
  return wrap(C); 
}


