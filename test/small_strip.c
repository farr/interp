#include"sample.h"
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

static double get_gaussian01(gsl_rng *rng, double sigma) {
  double x;

  do {
    x = 0.5 + gsl_ran_gaussian(rng, sigma);
  } while (x < 0.0 || x > 1.0);

  return x;
}

static double *small_strip(gsl_rng *rng) {
  static const double SIGMA0 = 0.25;
  static const double SIGMA1 = 0.025;
  double *result = malloc(2*sizeof(double));

  result[0] = get_gaussian01(rng, SIGMA0);
  result[1] = get_gaussian01(rng, SIGMA1);

  return result;
}

int main() {
  static const int NPDFSAMP = 10000;
  static const int NSAMP = 1000;
  int i;
  double **samps;

  samps = do_sample(NPDFSAMP, NSAMP, small_strip);

  for (i = 0; i < NSAMP; i++) {
    printf("%g %g\n", samps[i][0], samps[i][1]);
    free(samps[i]);
  }
  free(samps);

  return 0;
}
