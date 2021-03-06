#include <RcppArmadillo.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// distancia euclidea
// [[Rcpp::export]]
NumericMatrix euclidean(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), 
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  
  arma::colvec An =  arma::sum(arma::square(A),1);
  arma::colvec Bn =  arma::sum(arma::square(B),1);
  
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  return wrap(arma::sqrt(C)); 
}

// distancia de manhattan
// [[Rcpp::export]]
NumericMatrix manhattan(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), 
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  arma::mat res = arma::mat(m, m, arma::fill::zeros);
  
  for (int i = 0; i < m; i++) {
    A.each_row() -= B.row(i);
    res.col(i) = arma::sum(arma::abs(A), 1);
  }

  return wrap(res); 
}


// distancia de minkowsky
// [[Rcpp::export]]
NumericMatrix minkowsky(NumericMatrix Ar, NumericMatrix Br, double p) {
  int m = Ar.nrow(), 
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  arma::mat res = arma::mat(m, m, arma::fill::zeros);

  double q = 1.0 / p;
  
  for (int i = 0; i < m; i++) {
    A.each_row() -= B.row(i);
    res.col(i) = arma::sum(arma::pow(arma::abs(A), p), 1);
  }
  
  return wrap(arma::pow(res, q)); 
}


