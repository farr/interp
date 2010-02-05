#include"rng_init.h"
#include<stdio.h>

gsl_rng *rng_init(gsl_rng *rng) {
  FILE *devrandom = fopen("/dev/random", "r");
  unsigned long int seed;

  fread(&seed, sizeof(unsigned long int), 1, devrandom);

  gsl_rng_set(rng, seed);

  fclose(devrandom);

  return rng;
}
