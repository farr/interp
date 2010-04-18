#include"../interp.h"
#include"all-tests.h"
#include<stdio.h>

typedef int (*test)();

int main() {
  test tests[] = 
    { test_gaussian, 
      0 /* Terminated by NULL pointer. */};
  size_t i = 0;
  size_t nfail = 0;

  printf("Running tests...\n");

  while(tests[i] != 0) {
    if (!((tests[i])())) {
      nfail++;
    }
    i++;
  }

  printf("%ld failed out of %ld run\n", nfail, i);

  return nfail;
}
