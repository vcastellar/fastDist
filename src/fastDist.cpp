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
  C.for_each([](arma::mat::elem_type& val) {val = sqrt(val);});

  
  return wrap(C); 
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
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      res(i,j) = arma::sum(arma::abs(A.row(i) - B.row(j)));
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
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      res(i,j) =  arma::sum(arma::pow(arma::abs(A.row(i) - B.row(j)), p));
    }
  }
  
  res.for_each([&q](arma::mat::elem_type& val) {val = pow(val, q);});
  
  return wrap(res); 
}


// distancia de canberra
// [[Rcpp::export(.canberra)]]
NumericMatrix canberra(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      res(i, j) = arma::sum(arma::abs(A.row(i) - B.row(j)) / arma::abs(A.row(i) + B.row(j)));
    }
  }  
  
  return wrap(res);
}

// distancia de supremum
// [[Rcpp::export(.supremum)]]
NumericMatrix supremum(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      res(i,j) = arma::max(arma::abs(A.row(i) - B.row(j)));
    }
  }

  
  return wrap(res);
  
  
}


// distancia mahalanobis
// [[Rcpp::export(.mahalanobis)]]
NumericMatrix mahalanobis(NumericMatrix Ar) {
  int m = Ar.nrow(), 
      k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat res = arma::mat(m, m, arma::fill::zeros);
  arma::mat S(m,m);
  
  
  S = arma::inv(arma::cov(A));
  
  for (int i =0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      res(i,j) = arma::dot( (A.row(i) - A.row(j)) * S, (A.row(i) - A.row(j)).t());
    }
  }

  
  return wrap(arma::sqrt(res)); 
}

