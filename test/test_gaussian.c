#include"../interp.h"
#include"all-tests.h"
#include<gsl/gsl_rng.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

static double
log_gaussian(double x) {
  return -0.5*x*x - 0.5*log(2.0*M_PI);
}

static double 
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

static double **
alloc_matrix(size_t nrow, size_t ncol) {
  double **result = malloc(nrow*sizeof(double *));
  size_t i;
  
  assert(result != 0);

  for (i = 0; i < nrow; i++) {
    result[i] = malloc(ncol*sizeof(double));
    assert(result[i] != 0);
  }

  return result;
}

static void
free_matrix(double **mat, size_t nrow, size_t ncol) {
  size_t i;

  for (i = 0; i < nrow; i++) {
    free(mat[i]);
  }

  free(mat);
}

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

  for (i = 0; i < ndim; i++) {
    low[i] = -10.0;
    high[i] = 10.0;
  }

  randomize(rng);

  samples = alloc_gaussian_samples(rng, nsamples, ndim);

  col_mean_std(samples, nsamples, ndim, mu, std);

  printf("Sample Means: ");
  for (i = 0; i < ndim; i++) {
    printf("%g ", mu[i]);
  }
  printf("\nSample Stds: ");
  for (i = 0; i < ndim; i++) {
    printf("%g ", std[i]);
  }
  printf("\n");

  tr = make_interp_tree(ndim, nsamples, samples, low, high);

  interp_samples = alloc_matrix(nsamples, ndim);

  for (i = 0; i < nsamples; i++) {
    sample_density(ndim, nsamples, samples, tr, (uniform_random) gsl_rng_uniform, rng, interp_samples[i]);
  }

  col_mean_std(interp_samples, nsamples, ndim, mu, std);

  printf("Interp Means: ");
  for (i = 0; i < ndim; i++) {
    printf("%g ", mu[i]);
  }
  printf("\nInterp Stds: ");
  for (i = 0; i < ndim; i++) {
    printf("%g ", std[i]);
  }
  printf("\n");

  gsl_rng_free(rng);
  fclose(out);
  free_matrix(samples, nsamples, ndim);
  free_matrix(interp_samples, nsamples, ndim);

  return 1;  
}
