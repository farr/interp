#include"all-tests.h"
#include<gsl/gsl_rng.h>
#include<assert.h>
#include<stdio.h>

void
randomize(gsl_rng *rng) {
  FILE *devrandom = fopen("/dev/random", "r");
  unsigned long seed;
  size_t status;

  assert(devrandom != 0);
  
  status = fread(&seed, sizeof(unsigned long), 1, devrandom);
  assert(status == 1); /* One object read. */
  
  gsl_rng_set(rng, seed);

  fclose(devrandom);
}
