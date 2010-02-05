#include"sample.h"
#include<stdio.h>
#include<gsl/gsl_rng.h>

static double peak_fn(double x) {
  if (x < 0.5) {
    return x;
  } else {
    return 1-x;
  }
}

static double *get_sample(gsl_rng *rng) {
  double *result = malloc(2*sizeof(double));
  double x,y;

  do {
    x = gsl_rng_uniform(rng);
    y = gsl_rng_uniform(rng);
  } while (y > peak_fn(x));
  result[0] = x;

  do {
    x = gsl_rng_uniform(rng);
    y = gsl_rng_uniform(rng);
  } while (y > peak_fn(x));
  result[1] = x;

  return result;
}

int main() {
  static const int NPDFSAMP = 10000;
  static const int NSAMP = 100000;
  int i;
  double **samps = do_sample(NPDFSAMP, NSAMP, get_sample);
  
  for (i = 0; i < NSAMP; i++) {
    printf("%g %g\n", samps[i][0], samps[i][1]);
  }

  for (i = 0; i < NSAMP; i++) {
    free(samps[i]);
  }
  free(samps);

  return 0;
}
