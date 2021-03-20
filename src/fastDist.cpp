#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fastPdist2(NumericMatrix Ar, NumericMatrix Br) {
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
double manhattan(NumericVector Ur, NumericVector Vr) {
  int m = Ur.length(), 
      n = Vr.length();
  if (m != n) {
    return m;
  } else {
    int l = m;
  }
  
  arma::vec U = arma::vec(Ur.begin(), m, false); 
  arma::vec V = arma::vec(Vr.begin(), n, false); 
  
  arma::vec res = abs(U - V);
  
  res = sum(res);

  return res(0); 
}

