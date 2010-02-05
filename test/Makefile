GSL_PREFIX = /sw
MCMC_PREFIX = /Users/farr/Documents/code/mcmc

CFLAGS += -I$(GSL_PREFIX)/include -I$(MCMC_PREFIX) -L$(GSL_PREFIX)/lib -lgslcblas -lgsl
LFLAGS += -L$(GSL_PREFIX)/lib -lgslcblas -lgsl -lm

OBJS = sample.o rng_init.o ../denest.o $(MCMC_PREFIX)/mcmc.o

CHECK_PROGS = uniform2d linear_peak small_strip mcmc_gaussian_integral

$(CHECK_PROGS): sample.h $(OBJS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f ../denest.o
	rm -f *.dat
	rm -f $(CHECK_PROGS)