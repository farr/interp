#include"../interp.h"
#include"all-tests.h"
#include<gsl/gsl_rng.h>
#include<stdio.h>
#include<math.h>#
#include<assert.h>

static double **
alloc_gaussian_samples(gsl_rng *rng, size_t n, size_t ndim) {
  double **result = alloc_matrix(n, ndim);
  size_t i;

  for (i = 0; i < ndim; i++) {
    result[0][i] = 0.0;
  }

  for (i = 0; i < n-1; i++) {
    size_t j;

    for (j = 0; j < ndim; j++) {
      result[i+1][j] = gaussian_sample(rng, result[i][j]);
    }
  }

  return result;
}

int test_gaussian () {
  FILE *out = fopen("gaussian.dat", "w");
  size_t i;
  double **samples, **interp_samples;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  const size_t nsamples = 1000000;
  const size_t ndim = 2;
  double low[ndim], high[ndim];
  tree *tr;
  double mu[ndim], std[ndim];
  int status = 1; /* Success unless proved otherwise. */
  FILE *gaussian_out = fopen("gaussian_interp.dat", "w");

  assert(gaussian_out != 0);
  assert(rng != 0);

  printf("  Testing gaussian PDF interpolation... ");

  for (i = 0; i < ndim; i++) {
    low[i] = -10.0;
    high[i] = 10.0;
  }

  randomize(rng);

  samples = alloc_gaussian_samples(rng, nsamples, ndim);

  col_mean_std(samples, nsamples, ndim, mu, std);

  tr = make_interp_tree(ndim, nsamples, samples, low, high);

  interp_samples = alloc_matrix(nsamples, ndim);

  for (i = 0; i < nsamples; i++) {
    sample_density(ndim, nsamples, samples, tr, (uniform_random) gsl_rng_uniform, rng, interp_samples[i]);
  }

  col_mean_std(interp_samples, nsamples, ndim, mu, std);

  for (i = 0; i < nsamples; i++) {
    fprintf(gaussian_out, "%g %g\n", interp_samples[i][0], interp_samples[i][1]);
  }

  fclose(gaussian_out);

  if (fabs(mu[0]) > 5e-2) status = 0;
  if (fabs(mu[1]) > 5e-2) status = 0;
  if (fabs(std[0] - 1.0) > 5e-2) status = 0;
  if (fabs(std[1] - 1.0) > 5e-2) status = 0;

  gsl_rng_free(rng);
  fclose(out);
  free_matrix(samples, nsamples, ndim);
  free_matrix(interp_samples, nsamples, ndim);

  if (status != 0) {
    printf("PASSED\n");
  } else {
    printf("FAILED\n");
  }

  return status;  
}
