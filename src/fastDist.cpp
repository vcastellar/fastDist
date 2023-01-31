#include <RcppArmadillo.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppArmadillo;



// distancia euclidea
// [[Rcpp::export(.euclidean)]]
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
// [[Rcpp::export(.manhattan)]]
NumericMatrix manhattan(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  int x = std::min({m, n});
  if (x == n) {
    for (int j = 0; j < n; j++) {
      res.col(j) = arma::sum(arma::abs(A.each_row() - B.row(j)), 1);
    }  
  } 
  else {
    for (int j = 0; j < m; j++) {
      res.col(j) = arma::sum(arma::abs(B.each_row() - A.row(j)), 1);
    } 
  }
    
  return wrap(res);


}


// distancia de minkowsky
// [[Rcpp::export(.minkowski)]]
NumericMatrix minkowski(NumericMatrix Ar, NumericMatrix Br, double p) {
  int m = Ar.nrow(), 
      n = Br.nrow(),
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  arma::mat res = arma::mat(m, n, arma::fill::zeros);

  double q = 1.0 / p;
  
  for (int i = 0; i < n; i++) {
    res.col(i) =  arma::sum(arma::pow(arma::abs(A.each_row() - B.row(i)), p), 1);
  }
  
  res.for_each([&q](arma::mat::elem_type& val) {val = pow(val, q)28000000;});
  
  return wrap(res); 
}
