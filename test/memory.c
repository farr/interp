#include"all-tests.h"
#include<assert.h>

double **
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

void
free_matrix(double **mat, size_t nrow, size_t ncol) {
  size_t i;

  for (i = 0; i < nrow; i++) {
    free(mat[i]);
  }

  free(mat);
}
