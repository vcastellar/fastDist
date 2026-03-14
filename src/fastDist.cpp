#include <RcppArmadillo.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppArmadillo;

inline bool same_input(const NumericMatrix& Ar, const NumericMatrix& Br) {
  return (Ar.begin() == Br.begin()) && (Ar.nrow() == Br.nrow()) && (Ar.ncol() == Br.ncol());
}



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
  C.for_each([](arma::mat::elem_type& val) {val = std::sqrt(std::max(val, 0.0));});

  
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
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

#pragma omp parallel for schedule(static) if(m * n > 10000)
  for (int i = 0; i < m; i++) {
    const int jStart = symmetric ? i : 0;
    for (int j = jStart; j < n; j++) {
      double dist = 0.0;
      for (int col = 0; col < k; col++) {
        dist += std::abs(Ap[col * m + i] - Bp[col * n + j]);
      }
      res(i, j) = dist;
      if (symmetric && i != j) {
        res(j, i) = dist;
      }
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
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  double q = 1.0 / p;

#pragma omp parallel for schedule(static) if(m * n > 10000)
  for (int i = 0; i < m; i++) {
    const int jStart = symmetric ? i : 0;
    for (int j = jStart; j < n; j++) {
      double dist = 0.0;
      for (int col = 0; col < k; col++) {
        const double delta = std::abs(Ap[col * m + i] - Bp[col * n + j]);
        if (p == 1.0) {
          dist += delta;
        } else if (p == 2.0) {
          dist += delta * delta;
        } else {
          dist += std::pow(delta, p);
        }
      }
      res(i, j) = dist;
      if (symmetric && i != j) {
        res(j, i) = dist;
      }
    }
  }
  
  res.for_each([&q](arma::mat::elem_type& val) {val = std::pow(val, q);});
  
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
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

#pragma omp parallel for schedule(static) if(m * n > 10000)
  for (int i = 0; i < m; i++) {
    const int jStart = symmetric ? i : 0;
    for (int j = jStart; j < n; j++) {
      double dist = 0.0;
      for (int col = 0; col < k; col++) {
        const double a = Ap[col * m + i];
        const double b = Bp[col * n + j];
        const double den = std::abs(a) + std::abs(b);
        if (den > 0.0) {
          dist += std::abs(a - b) / den;
        }
      }
      res(i, j) = dist;
      if (symmetric && i != j) {
        res(j, i) = dist;
      }
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
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

#pragma omp parallel for schedule(static) if(m * n > 10000)
  for (int i = 0; i < m; i++) {
    const int jStart = symmetric ? i : 0;
    for (int j = jStart; j < n; j++) {
      double dist = 0.0;
      for (int col = 0; col < k; col++) {
        dist = std::max(dist, std::abs(Ap[col * m + i] - Bp[col * n + j]));
      }
      res(i, j) = dist;
      if (symmetric && i != j) {
        res(j, i) = dist;
      }
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
  arma::mat S = arma::inv_sympd(arma::cov(A));
  arma::mat AS = A * S;
  arma::vec q = arma::sum(AS % A, 1);
  arma::mat G = AS * A.t();
  arma::mat res = arma::repmat(q, 1, m) + arma::repmat(q.t(), m, 1) - 2.0 * G;
  res.transform([](double x) { return std::sqrt(std::max(x, 0.0)); });
  return wrap(res); 
}
