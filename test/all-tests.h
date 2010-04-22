#ifndef __ALL_TESTS_H__
#define __ALL_TESTS_H__

#include<gsl/gsl_rng.h>
#include<stdlib.h>

/* Memory functions */
double **
alloc_matrix(size_t nrow, size_t ncol);

void
free_matrix(double **mat, size_t nrow, size_t ncol);

/* Utility functions. */
void
randomize(gsl_rng *rng);

double
mean(double *xs, size_t n);
void
col_mean(double **xs, size_t nrow, size_t ncol, double *mu);

double
std(double *xs, size_t n);
void
col_mean_std(double **xs, size_t nrow, size_t ncol, double *mu, double *std);

/* All test functions. */
int test_gaussian();

#endif /* __ALL_TESTS_H__ */
