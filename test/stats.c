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
