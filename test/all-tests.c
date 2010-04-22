#include"../interp.h"
#include"all-tests.h"
#include<stdio.h>

typedef int (*test)();

typedef struct {
  test t;
  char *name;
} test_info;

int main() {
  const size_t ntests = 2;
  test_info tests[ntests];
  size_t i = 0;
  size_t nfail = 0;

  tests[0].t = test_gaussian;
  tests[0].name = "interpolate 2-D Gaussian PDF";

  tests[1].t = test_jump_proposal;
  tests[1].name = "jump proposal test";

  printf("Running tests...\n");

  for (i = 0; i < ntests; i++) {
    int status;
    printf("  Running test: %s...", tests[i].name);
    fflush(stdout);
    
    status = tests[i].t();

    if (status) {
      printf("FAILED (with code %d)\n", status);
      fflush(stdout);
      nfail++;
    } else {
      printf("PASSED\n");
      fflush(stdout);
    }
  }

  printf("%ld failed, %ld succeeded out of %ld run\n", nfail, ntests-nfail, ntests);

  return nfail;
}
