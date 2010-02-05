#include<gsl/gsl_rng.h>
#include"rng_init.h"
#include<gsl/gsl_randist.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"../denest.h"
#include<string.h>
#include<mcmc.h>

static double log_gaussian(double mu, double sigma, double x) {
  double dx = x - mu;
  double sigma2 = sigma*sigma;

  return -dx*dx/(2.0*sigma2) - 0.5*log(2.0*M_PI) - log(sigma);
}

typedef struct {
  gsl_rng *rng;
  size_t ndim;
  double *sigma;
} proposal_data;

static double *proposal(double *x, proposal_data *data) {
  size_t ndim = data->ndim;
  double *propose = malloc(ndim*sizeof(double));
  size_t i;

  memcpy(propose, x, ndim*sizeof(double));

  for (i = 0; i < ndim; i++) {
    propose[i] += gsl_ran_gaussian(data->rng, data->sigma[i]/10.0);
  }

  return propose;
}

typedef struct {
  size_t ndim;
  double *mu;
  double *sigma;
} likelihood_data;

static double log_likelihood(double *x, likelihood_data *data) {
  size_t ndim = data->ndim;
  size_t i;
  double ll = 0.0;
  
  for (i = 0; i < ndim; i++) {
    ll += log_gaussian(data->mu[i], data->sigma[i], x[i]);
  }

  return ll;
}

static double tprob(double *x, void *data) { return 0.0; }

static void bounds_of_coords(size_t ndim, size_t npts, double **pts, 
                             double *low, double *high) {
  size_t i;

  memcpy(low, pts[0], ndim*sizeof(double));
  memcpy(high, pts[0], ndim*sizeof(double));

  for (i = 1; i < npts; i++) {
    size_t j;

    for (j = 0; j < ndim; j++) {
      if (low[j] > pts[i][j]) low[j] = pts[i][j];
      if (high[j] < pts[i][j]) high[j] = pts[i][j];
    }
  }
}

int main() {
  const size_t ndim = 3;
  const size_t nsamp = 1000000;
  const double box_width = 10.0; /* In units of sigma. */
  const double epsrel = 1e-1;
  double mu[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double sigma[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  gsl_rng *rng = rng_init(gsl_rng_alloc(gsl_rng_ranlxd2));
  sample **samples;
  double **coords = malloc(nsamp*sizeof(double *));
  double *likes = malloc(nsamp*sizeof(double));
  tree *t;
  double lower_left[ndim], upper_right[ndim];
  size_t i;
  double I;
  int success;
  double *x0 = malloc(ndim*sizeof(double)); /* Will be freed by samples_free. */
  proposal_data pdata;
  likelihood_data ldata;

  pdata.rng = rng;
  pdata.sigma = sigma;
  pdata.ndim = ndim;

  ldata.ndim = ndim;
  ldata.mu = mu;
  ldata.sigma = sigma;

  memcpy(x0, mu, ndim*sizeof(double));
  samples = mcmc(rng, nsamp, x0,
                 (jump_proposal_fn *) proposal, &pdata,
                 (log_transition_probability_fn *) tprob, 0,
                 (log_posterior_fn *) log_likelihood, &ldata,
                 free);

  for (i = 0; i < nsamp; i++) {
    coords[i] = samples[i]->params;
    likes[i] = exp(samples[i]->log_posterior);
  }

  bounds_of_coords(ndim, nsamp, coords, lower_left, upper_right);

  t = make_density_tree(ndim, nsamp, coords, lower_left, upper_right);

  I = integrate_samples(ndim, nsamp, coords, likes, t);

  if (fabs(I - 1.0) < epsrel) {
    fprintf(stderr, "Success!\n");
    success = 1;
  } else {
    fprintf(stderr, "Failure: I = %g (volume = %g)\n", I, tree_volume(t));
    success = 0;
  }
  
  for (i = 0; i < nsamp; i++) {
    printf("%g %g\n", coords[i][0], coords[i][1]);
  }

  samples_free(samples, nsamp, free);
  free(coords);
  free(likes);
  free_density_tree(t);

  if (success) {
    return 0;
  } else {
    return 1;
  }
}
