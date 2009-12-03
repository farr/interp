#include"sample.h"
#include"rng_init.h"
#include"../denest.h"

double **do_sample(int npdfsamp, int nsamp, sample gen_sample) {
  int i;
  double **pdfsamps;
  double **samps;
  gsl_rng *rng;
  tree *t;
  double low[] = {0.0, 0.0}, high[] = {1.0, 1.0};

  pdfsamps = malloc(npdfsamp*sizeof(double *));
  samps = malloc(nsamp*sizeof(double*));
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  rng_init(rng);

  for (i = 0; i < npdfsamp; i++) {
    pdfsamps[i] = gen_sample(rng);
  }

  t = make_density_tree(2, npdfsamp, pdfsamps, low, high);

  for (i = 0; i < nsamp; i++) {
    samps[i] = malloc(2*sizeof(double));
    sample_density(2, npdfsamp, pdfsamps, t, (uniform_random) gsl_rng_uniform, 
                   rng, samps[i]);
  }

  for (i = 0; i < npdfsamp; i++) {
    free(pdfsamps[i]);
  }
  free(pdfsamps);

  free_density_tree(t);

  gsl_rng_free(rng);

  return samps;
}
