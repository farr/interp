#include<math.h>
#include<gsl/gsl_rng.h>

double
log_gaussian(double x) {
  return -0.5*x*x - 0.5*log(2.0*M_PI);
}

double 
gaussian_sample(gsl_rng *rng, double x) {
  double y = x + gsl_rng_uniform(rng) - 0.5;
  double lgx = log_gaussian(x);
  double lgy = log_gaussian(y);
  double log_accept = lgy - lgx; /* Jump distribution is symmetric. */
  
  if (log(gsl_rng_uniform(rng)) < log_accept) {
    return y;
  } else {
    return x;
  }
}

