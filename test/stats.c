#include"all-tests.h"
#include<math.h>

double
mean(double *xs, size_t n) {
  size_t i;
  double m = 0.0;
  
  for (i = 0; i < n; i++) {
    m += xs[i];
  }

  return m / n;
}

void
col_mean(double **xs, size_t nrow, size_t ncol, double *mu) {
  size_t i;

  for (i = 0; i < ncol; i++) {
    mu[i] = 0.0;
  }

  for (i = 0; i < nrow; i++) {
    size_t j;

    for (j = 0; j < ncol; j++) {
      mu[j] += xs[i][j];
    }
  }

  for (i = 0; i < ncol; i++) {
    mu[i] /= nrow;
  }
}

double
std(double *xs, size_t n) {
  double mu = mean(xs,n);
  double var = 0.0;
  size_t i;

  for (i = 0; i < n; i++) {
    double tmp = xs[i] - mu;
    var += tmp*tmp;
  }

  return sqrt(var / (n-1));
  
}

void
col_mean_std(double **xs, size_t nrow, size_t ncol, double *mu, double *std) {
  size_t i;

  col_mean(xs, nrow, ncol, mu);

  for (i = 0; i < ncol; i++) {
    std[i] = 0.0;
  }

  for (i = 0; i < nrow; i++) {
    size_t j;

    for (j = 0; j < ncol; j++) {
      double tmp = mu[j] - xs[i][j];

      std[j] += tmp*tmp;
    }
  }

  for (i = 0; i < ncol; i++) {
    std[i] = sqrt(std[i]/(nrow-1));
  }
}
