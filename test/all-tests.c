#include"../interp.h"
#include"all-tests.h"
#include<stdio.h>

typedef int (*test)();

typedef struct {
  test t;
  char *name;
} test_info;

int main() {
  const size_t ntests = 1;
  test_info tests[ntests];
  size_t i = 0;
  size_t nfail = 0;

  tests[0].t = test_gaussian;
  tests[0].name = "interpolate 2-D Gaussian PDF";

  printf("Running tests...\n");

  for (i = 0; i < ntests; i++) {
    int status;
    printf("  Running test: %s...", tests[i].name);
    fflush(stdout);
    
    status = tests[i].t();

    if (status) {
      printf("FAILED (with code %d)\n", status);
      fflush(stdout);
    } else {
      printf("PASSED\n");
      fflush(stdout);
    }
  }

  printf("%ld failed, %ld succeeded out of %ld run\n", nfail, i-nfail, i);

  return nfail;
}
