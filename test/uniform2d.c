#include"sample.h"
#include<stdio.h>

static double *sample_uniform(gsl_rng *rng) {
  double *result = malloc(2*sizeof(double));

  result[0] = gsl_rng_uniform(rng);
  result[1] = gsl_rng_uniform(rng);

  return result;
}

int main() {
  static const int NPDFSAMP = 10000;
  static const int NSAMP = 1000;
  int i;
  double **samps = do_sample(NPDFSAMP, NSAMP, sample_uniform);
  
  for (i = 0; i < NSAMP; i++) {
    printf("%g %g\n", samps[i][0], samps[i][1]);
  }

  for (i = 0; i < NSAMP; i++) {
    free(samps[i]);
  }
  free(samps);

  return 0;
}
