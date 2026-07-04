#include <RcppArmadillo.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppArmadillo;

inline bool same_input(const NumericMatrix& Ar, const NumericMatrix& Br) {
  return (Ar.begin() == Br.begin()) && (Ar.nrow() == Br.nrow()) && (Ar.ncol() == Br.ncol());
}



// euclidean distance
// [[Rcpp::export(.euclidean)]]
NumericMatrix euclidean(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();
  
  arma::colvec An = arma::sum(arma::square(A), 1);
  arma::colvec Bn = arma::sum(arma::square(B), 1);
  

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      const double ai_norm = An[i];
      for (int j = i; j < m; j++) {
        const double* a = Ap + i;
        const double* b = Bp + j;
        double dot = 0.0;
        for (int col = 0; col < k; col++) {
          dot += (*a) * (*b);
          a += m;
          b += n;
        }
        double sqDist = ai_norm + Bn[j] - 2.0 * dot;
        if (sqDist < 0.0) sqDist = 0.0;
        const double dist = std::sqrt(sqDist);
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      const double ai_norm = An[i];
      for (int j = 0; j < n; j++) {
        const double* a = Ap + i;
        const double* b = Bp + j;
        double dot = 0.0;
        for (int col = 0; col < k; col++) {
          dot += (*a) * (*b);
          a += m;
          b += n;
        }
        double sqDist = ai_norm + Bn[j] - 2.0 * dot;
        if (sqDist < 0.0) sqDist = 0.0;
        res(i, j) = std::sqrt(sqDist);
      }
    }
  }
  
  return wrap(res);
}



// manhattan distance
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
  
  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double dist = 0.0;
        for (int col = 0; col < k; col++) {
          dist += std::abs(Ap[col * m + i] - Bp[col * n + j]);
        }
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double dist = 0.0;
        for (int col = 0; col < k; col++) {
          dist += std::abs(Ap[col * m + i] - Bp[col * n + j]);
        }
        res(i, j) = dist;
      }
    }
  }
  
  
  
  return wrap(res);
  
  
}


// minkowski distance
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
  
  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
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
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
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
      }
    }
  }
  
  res.for_each([&q](arma::mat::elem_type& val) {val = std::pow(val, q);});
  
  return wrap(res); 
}


// correlation distance
// [[Rcpp::export(.correlation)]]
NumericMatrix correlation(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();
  arma::colvec meanA = arma::mean(A, 1);
  arma::colvec meanB = arma::mean(B, 1);
  arma::colvec normA(m, arma::fill::zeros);
  arma::colvec normB(n, arma::fill::zeros);

#pragma omp parallel for schedule(static) if(m > 100)
  for (int i = 0; i < m; i++) {
    double ss = 0.0;
    for (int col = 0; col < k; col++) {
      const double centered = Ap[col * m + i] - meanA[i];
      ss += centered * centered;
    }
    normA[i] = std::sqrt(ss);
  }

#pragma omp parallel for schedule(static) if(n > 100)
  for (int j = 0; j < n; j++) {
    double ss = 0.0;
    for (int col = 0; col < k; col++) {
      const double centered = Bp[col * n + j] - meanB[j];
      ss += centered * centered;
    }
    normB[j] = std::sqrt(ss);
  }

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double corr = 0.0;
        const double denom = normA[i] * normB[j];
        if (denom > 0.0) {
          double dot = 0.0;
          for (int col = 0; col < k; col++) {
            const double a = Ap[col * m + i] - meanA[i];
            const double b = Bp[col * n + j] - meanB[j];
            dot += a * b;
          }
          corr = dot / denom;
          corr = std::max(-1.0, std::min(1.0, corr));
        }
        const double dist = 1.0 - corr;
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double corr = 0.0;
        const double denom = normA[i] * normB[j];
        if (denom > 0.0) {
          double dot = 0.0;
          for (int col = 0; col < k; col++) {
            const double a = Ap[col * m + i] - meanA[i];
            const double b = Bp[col * n + j] - meanB[j];
            dot += a * b;
          }
          corr = dot / denom;
          corr = std::max(-1.0, std::min(1.0, corr));
        }
        res(i, j) = 1.0 - corr;
      }
    }
  }

  return wrap(res);
}

// cosine distance
// [[Rcpp::export(.cosine)]]
NumericMatrix cosine(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();
  arma::colvec normA = arma::sqrt(arma::sum(arma::square(A), 1));
  arma::colvec normB = arma::sqrt(arma::sum(arma::square(B), 1));

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double cos_sim = 0.0;
        const double denom = normA[i] * normB[j];
        if (denom > 0.0) {
          double dot = 0.0;
          for (int col = 0; col < k; col++) {
            dot += Ap[col * m + i] * Bp[col * n + j];
          }
          cos_sim = dot / denom;
          cos_sim = std::max(-1.0, std::min(1.0, cos_sim));
        }
        const double dist = 1.0 - cos_sim;
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double cos_sim = 0.0;
        const double denom = normA[i] * normB[j];
        if (denom > 0.0) {
          double dot = 0.0;
          for (int col = 0; col < k; col++) {
            dot += Ap[col * m + i] * Bp[col * n + j];
          }
          cos_sim = dot / denom;
          cos_sim = std::max(-1.0, std::min(1.0, cos_sim));
        }
        res(i, j) = 1.0 - cos_sim;
      }
    }
  }

  return wrap(res);
}


// canberra distance
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
  
  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
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
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
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
      }
    }
  }
  
  return wrap(res);
}

// supremum distance
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
  
  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double dist = 0.0;
        for (int col = 0; col < k; col++) {
          dist = std::max(dist, std::abs(Ap[col * m + i] - Bp[col * n + j]));
        }
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double dist = 0.0;
        for (int col = 0; col < k; col++) {
          dist = std::max(dist, std::abs(Ap[col * m + i] - Bp[col * n + j]));
        }
        res(i, j) = dist;
      }
    }
  }
  
  
  return wrap(res);
  
  
}


// squared euclidean distance
// [[Rcpp::export(.squared_euclidean)]]
NumericMatrix squared_euclidean(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  arma::colvec An = arma::sum(arma::square(A), 1);
  arma::colvec Bn = arma::sum(arma::square(B), 1);

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      const double ai_norm = An[i];
      for (int j = i; j < m; j++) {
        const double* a = Ap + i;
        const double* b = Bp + j;
        double dot = 0.0;
        for (int col = 0; col < k; col++) {
          dot += (*a) * (*b);
          a += m;
          b += n;
        }
        double sqDist = ai_norm + Bn[j] - 2.0 * dot;
        if (sqDist < 0.0) sqDist = 0.0;
        res(i, j) = sqDist;
        if (i != j) res(j, i) = sqDist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      const double ai_norm = An[i];
      for (int j = 0; j < n; j++) {
        const double* a = Ap + i;
        const double* b = Bp + j;
        double dot = 0.0;
        for (int col = 0; col < k; col++) {
          dot += (*a) * (*b);
          a += m;
          b += n;
        }
        double sqDist = ai_norm + Bn[j] - 2.0 * dot;
        if (sqDist < 0.0) sqDist = 0.0;
        res(i, j) = sqDist;
      }
    }
  }

  return wrap(res);
}


// bray-curtis distance
// [[Rcpp::export(.braycurtis)]]
NumericMatrix braycurtis(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double num = 0.0, den = 0.0;
        for (int col = 0; col < k; col++) {
          const double a = Ap[col * m + i];
          const double b = Bp[col * n + j];
          num += std::abs(a - b);
          den += std::abs(a + b);
        }
        const double dist = den > 0.0 ? num / den : 0.0;
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double num = 0.0, den = 0.0;
        for (int col = 0; col < k; col++) {
          const double a = Ap[col * m + i];
          const double b = Bp[col * n + j];
          num += std::abs(a - b);
          den += std::abs(a + b);
        }
        res(i, j) = den > 0.0 ? num / den : 0.0;
      }
    }
  }

  return wrap(res);
}


// hellinger distance
// [[Rcpp::export(.hellinger)]]
NumericMatrix hellinger(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);

  // precompute element-wise square roots (inputs are assumed non-negative)
  arma::mat sqrtA = arma::sqrt(A);
  arma::mat sqrtB = arma::sqrt(B);
  const double* Ap = sqrtA.memptr();
  const double* Bp = sqrtB.memptr();
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double delta = Ap[col * m + i] - Bp[col * n + j];
          acc += delta * delta;
        }
        const double dist = inv_sqrt2 * std::sqrt(acc);
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double delta = Ap[col * m + i] - Bp[col * n + j];
          acc += delta * delta;
        }
        res(i, j) = inv_sqrt2 * std::sqrt(acc);
      }
    }
  }

  return wrap(res);
}


// chi-squared distance
// [[Rcpp::export(.chisquared)]]
NumericMatrix chisquared(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double a = Ap[col * m + i];
          const double b = Bp[col * n + j];
          const double den = a + b;
          if (den != 0.0) {
            const double delta = a - b;
            acc += (delta * delta) / den;
          }
        }
        const double dist = 0.5 * acc;
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double a = Ap[col * m + i];
          const double b = Bp[col * n + j];
          const double den = a + b;
          if (den != 0.0) {
            const double delta = a - b;
            acc += (delta * delta) / den;
          }
        }
        res(i, j) = 0.5 * acc;
      }
    }
  }

  return wrap(res);
}


// jensen-shannon distance
// [[Rcpp::export(.jensenshannon)]]
NumericMatrix jensenshannon(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);

  // normalise each row to a probability distribution (sum to 1)
  arma::colvec sumA = arma::sum(A, 1);
  arma::colvec sumB = arma::sum(B, 1);
  arma::mat P = A;
  arma::mat Q = B;
  for (int i = 0; i < m; i++) if (sumA[i] != 0.0) P.row(i) /= sumA[i];
  for (int j = 0; j < n; j++) if (sumB[j] != 0.0) Q.row(j) /= sumB[j];
  const double* Pp = P.memptr();
  const double* Qp = Q.memptr();

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double p = Pp[col * m + i];
          const double q = Qp[col * n + j];
          const double mid = 0.5 * (p + q);
          if (p > 0.0) acc += 0.5 * p * std::log(p / mid);
          if (q > 0.0) acc += 0.5 * q * std::log(q / mid);
        }
        if (acc < 0.0) acc = 0.0;
        const double dist = std::sqrt(acc);
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double p = Pp[col * m + i];
          const double q = Qp[col * n + j];
          const double mid = 0.5 * (p + q);
          if (p > 0.0) acc += 0.5 * p * std::log(p / mid);
          if (q > 0.0) acc += 0.5 * q * std::log(q / mid);
        }
        if (acc < 0.0) acc = 0.0;
        res(i, j) = std::sqrt(acc);
      }
    }
  }

  return wrap(res);
}


// haversine (great-circle) distance, in kilometres
// [[Rcpp::export(.haversine)]]
NumericMatrix haversine(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  if (k != 2) {
    stop("haversine distance requires exactly 2 columns (latitude, longitude in degrees)");
  }
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();
  const double deg2rad = M_PI / 180.0;
  const double R = 6371.0; // mean Earth radius in km

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      const double lat1 = Ap[i] * deg2rad;
      const double lon1 = Ap[m + i] * deg2rad;
      for (int j = i; j < m; j++) {
        const double lat2 = Bp[j] * deg2rad;
        const double lon2 = Bp[n + j] * deg2rad;
        const double dlat = lat2 - lat1;
        const double dlon = lon2 - lon1;
        const double sdlat = std::sin(dlat * 0.5);
        const double sdlon = std::sin(dlon * 0.5);
        double a = sdlat * sdlat + std::cos(lat1) * std::cos(lat2) * sdlon * sdlon;
        if (a > 1.0) a = 1.0;
        const double dist = 2.0 * R * std::asin(std::sqrt(a));
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      const double lat1 = Ap[i] * deg2rad;
      const double lon1 = Ap[m + i] * deg2rad;
      for (int j = 0; j < n; j++) {
        const double lat2 = Bp[j] * deg2rad;
        const double lon2 = Bp[n + j] * deg2rad;
        const double dlat = lat2 - lat1;
        const double dlon = lon2 - lon1;
        const double sdlat = std::sin(dlat * 0.5);
        const double sdlon = std::sin(dlon * 0.5);
        double a = sdlat * sdlat + std::cos(lat1) * std::cos(lat2) * sdlon * sdlon;
        if (a > 1.0) a = 1.0;
        res(i, j) = 2.0 * R * std::asin(std::sqrt(a));
      }
    }
  }

  return wrap(res);
}


// standardized euclidean distance
// (each feature scaled by its sample variance, estimated from A)
// [[Rcpp::export(.standardized_euclidean)]]
NumericMatrix standardized_euclidean(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  // per-feature sample variance from A (n-1 denominator)
  arma::rowvec var = arma::var(A, 0, 0);
  arma::vec inv_var(k);
  for (int col = 0; col < k; col++) {
    inv_var[col] = var[col] > 0.0 ? 1.0 / var[col] : 0.0;
  }

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double delta = Ap[col * m + i] - Bp[col * n + j];
          acc += delta * delta * inv_var[col];
        }
        const double dist = std::sqrt(acc);
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double delta = Ap[col * m + i] - Bp[col * n + j];
          acc += delta * delta * inv_var[col];
        }
        res(i, j) = std::sqrt(acc);
      }
    }
  }

  return wrap(res);
}


// helper: average ranks (1-based) of every row of M, returned with the same
// column-major layout as M
inline arma::mat row_ranks(const arma::mat& M) {
  const int rows = M.n_rows;
  const int cols = M.n_cols;
  arma::mat R(rows, cols);
#pragma omp parallel for schedule(static) if(rows > 100)
  for (int i = 0; i < rows; i++) {
    std::vector<int> idx(cols);
    for (int c = 0; c < cols; c++) idx[c] = c;
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b) { return M(i, a) < M(i, b); });
    int c = 0;
    while (c < cols) {
      int j = c;
      while (j + 1 < cols && M(i, idx[j + 1]) == M(i, idx[c])) j++;
      const double avg = ((c + 1) + (j + 1)) / 2.0;
      for (int t = c; t <= j; t++) R(i, idx[t]) = avg;
      c = j + 1;
    }
  }
  return R;
}

// spearman correlation distance (1 - Spearman's rho)
// [[Rcpp::export(.spearman)]]
NumericMatrix spearman(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);

  // ranks within each row, then Pearson correlation on the ranks
  arma::mat RA = row_ranks(A);
  arma::mat RB = symmetric ? RA : row_ranks(B);
  const double* Ap = RA.memptr();
  const double* Bp = RB.memptr();
  arma::colvec meanA = arma::mean(RA, 1);
  arma::colvec meanB = arma::mean(RB, 1);
  arma::colvec normA(m, arma::fill::zeros);
  arma::colvec normB(n, arma::fill::zeros);

#pragma omp parallel for schedule(static) if(m > 100)
  for (int i = 0; i < m; i++) {
    double ss = 0.0;
    for (int col = 0; col < k; col++) {
      const double centered = Ap[col * m + i] - meanA[i];
      ss += centered * centered;
    }
    normA[i] = std::sqrt(ss);
  }

#pragma omp parallel for schedule(static) if(n > 100)
  for (int j = 0; j < n; j++) {
    double ss = 0.0;
    for (int col = 0; col < k; col++) {
      const double centered = Bp[col * n + j] - meanB[j];
      ss += centered * centered;
    }
    normB[j] = std::sqrt(ss);
  }

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double rho = 0.0;
        const double denom = normA[i] * normB[j];
        if (denom > 0.0) {
          double dot = 0.0;
          for (int col = 0; col < k; col++) {
            const double a = Ap[col * m + i] - meanA[i];
            const double b = Bp[col * n + j] - meanB[j];
            dot += a * b;
          }
          rho = dot / denom;
          rho = std::max(-1.0, std::min(1.0, rho));
        }
        const double dist = 1.0 - rho;
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double rho = 0.0;
        const double denom = normA[i] * normB[j];
        if (denom > 0.0) {
          double dot = 0.0;
          for (int col = 0; col < k; col++) {
            const double a = Ap[col * m + i] - meanA[i];
            const double b = Bp[col * n + j] - meanB[j];
            dot += a * b;
          }
          rho = dot / denom;
          rho = std::max(-1.0, std::min(1.0, rho));
        }
        res(i, j) = 1.0 - rho;
      }
    }
  }

  return wrap(res);
}


// mahalanobis distance
// [[Rcpp::export(.mahalanobis)]]
NumericMatrix mahalanobis(NumericMatrix Ar) {
  int m = Ar.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat S = arma::inv_sympd(arma::cov(A));
  arma::mat AS = A * S;
  arma::vec q = arma::sum(AS % A, 1);
  arma::mat res = arma::mat(m, m, arma::fill::zeros);
  const double* AS_p = AS.memptr();
  const double* A_p = A.memptr();

#pragma omp parallel for schedule(static) if(m * m > 10000)
  for (int i = 0; i < m; i++) {
    for (int j = i; j < m; j++) {
      double dot = 0.0;
      for (int col = 0; col < k; col++) {
        dot += AS_p[col * m + i] * A_p[col * m + j];
      }
      const double sqDist = std::max(q[i] + q[j] - 2.0 * dot, 0.0);
      const double dist = std::sqrt(sqDist);
      res(i, j) = dist;
      if (i != j) {
        res(j, i) = dist;
      }
    }
  }

  return wrap(res);
}


// hamming distance (binary/categorical)
// [[Rcpp::export(.hamming)]]
NumericMatrix hamming(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        int dist = 0;
        for (int col = 0; col < k; col++) {
          if (Ap[col * m + i] != Bp[col * n + j]) {
            dist++;
          }
        }
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        int dist = 0;
        for (int col = 0; col < k; col++) {
          if (Ap[col * m + i] != Bp[col * n + j]) {
            dist++;
          }
        }
        res(i, j) = dist;
      }
    }
  }

  return wrap(res);
}


// jaccard distance (binary/set-based)
// [[Rcpp::export(.jaccard)]]
NumericMatrix jaccard(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        int intersection = 0, union_size = 0;
        for (int col = 0; col < k; col++) {
          const bool ai = Ap[col * m + i] != 0.0;
          const bool bi = Bp[col * n + j] != 0.0;
          if (ai || bi) union_size++;
          if (ai && bi) intersection++;
        }
        const double dist = union_size > 0 ? 1.0 - (double)intersection / union_size : 0.0;
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        int intersection = 0, union_size = 0;
        for (int col = 0; col < k; col++) {
          const bool ai = Ap[col * m + i] != 0.0;
          const bool bi = Bp[col * n + j] != 0.0;
          if (ai || bi) union_size++;
          if (ai && bi) intersection++;
        }
        res(i, j) = union_size > 0 ? 1.0 - (double)intersection / union_size : 0.0;
      }
    }
  }

  return wrap(res);
}


// gower distance (mixed data types)
// [[Rcpp::export(.gower)]]
NumericMatrix gower(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(),
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  arma::mat res = arma::mat(m, n, arma::fill::zeros);
  const bool symmetric = same_input(Ar, Br);
  const double* Ap = A.memptr();
  const double* Bp = B.memptr();

  // compute range for each column (max - min from A)
  arma::vec range(k);
  for (int col = 0; col < k; col++) {
    const double minVal = A.col(col).min();
    const double maxVal = A.col(col).max();
    range[col] = maxVal - minVal;
    if (range[col] == 0.0) range[col] = 1.0; // avoid division by zero
  }

  if (symmetric) {
#pragma omp parallel for schedule(static) if(m * m > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double a = Ap[col * m + i];
          const double b = Bp[col * n + j];
          acc += std::abs(a - b) / range[col];
        }
        const double dist = acc / k;
        res(i, j) = dist;
        if (i != j) res(j, i) = dist;
      }
    }
  } else {
#pragma omp parallel for schedule(static) if(m * n > 10000)
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double acc = 0.0;
        for (int col = 0; col < k; col++) {
          const double a = Ap[col * m + i];
          const double b = Bp[col * n + j];
          acc += std::abs(a - b) / range[col];
        }
        res(i, j) = acc / k;
      }
    }
  }

  return wrap(res);
}
