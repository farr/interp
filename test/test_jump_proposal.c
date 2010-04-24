#include"all-tests.h"
#include<gsl/gsl_rng.h>
#include<assert.h>
#include<stdlib.h>
#include"../interp.h"
#include<stdio.h>
#include<math.h>

int 
test_jump_proposal() {
  const size_t nsamples = 1000000;
  const size_t ndim = 2;
  double low[] = {-10.0, -10.0};
  double high[] = {10.0, 10.0};
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  double **samples = alloc_matrix(nsamples, ndim);
  double **interp_samples = alloc_matrix(nsamples, ndim);
  size_t i;
  tree *t;
  double mu[ndim], sigma[ndim];
  int status = 0;
  FILE *jpout = fopen("jump_proposal_gaussian.dat", "w");

  assert(rng != 0);
  assert(samples != 0);
  assert(interp_samples != 0);
  assert(jpout != 0);

  randomize(rng);

  for (i = 0; i < ndim; i++) {
    samples[0][i] = 0.0;
  }
  for (i = 0; i < nsamples-1; i++) {
    size_t j;
    for (j = 0; j < ndim; j++) {
      samples[i+1][j] = gaussian_sample(rng, samples[i][j]);
    }
  }

  t = make_interp_tree(ndim, nsamples, samples, low, high);
  assert(t != 0);

  for (i = 0; i < ndim; i++) {
    interp_samples[0][i] = 0.0;
  }
  for (i = 0; i < nsamples-1; i++) {
    double log_like = 0.0;
    double orig_log_like = 0.0;
    double log_accept;
    size_t j;
    sample_density(ndim, nsamples, samples, t, (double (*)(void *))gsl_rng_uniform, 
                   rng, interp_samples[i+1]);
    
    for (j = 0; j < ndim; j++) {
      log_like += log_gaussian(interp_samples[i+1][j]);
      orig_log_like += log_gaussian(interp_samples[i][j]);
    }

    log_accept = log_like - orig_log_like + 
      log(jump_probability(interp_samples[i], t)) - /* Jump prob backwards. */
      log(jump_probability(interp_samples[i+1], t)); /* Jump prob forwards. */

    if (log(gsl_rng_uniform(rng)) < log_accept) {
      /* Do nothing. */
    } else {
      /* Reject jump, restore original probability. */
      memcpy(interp_samples[i+1], interp_samples[i], ndim*sizeof(double));
    }
  }

  col_mean_std(interp_samples, nsamples, ndim, mu, sigma);

  if (fabs(mu[0]) > 5e-2) status = 1;
  if (fabs(mu[1]) > 5e-2) status = 2;
  if (fabs(sigma[0] - 1.0) > 5e-2) status = 3;
  if (fabs(sigma[1] - 1.0) > 5e-2) status = 4;

  for (i = 0; i < nsamples; i++) {
    fprintf(jpout, "%g %g\n", interp_samples[i][0], interp_samples[i][1]);
  }

  fclose(jpout);
  free_interp_tree(t);
  gsl_rng_free(rng);
  free_matrix(samples, nsamples, ndim);
  free_matrix(interp_samples, nsamples, ndim);

  return status;
}
