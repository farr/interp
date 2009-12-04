GSL_PREFIX = /sw

CFLAGS += -I$(GSL_PREFIX)/include -L$(GSL_PREFIX)/lib -lgslcblas -lgsl
LFLAGS += -L$(GSL_PREFIX)/lib -lgslcblas -lgsl

OBJS = sample.o rng_init.o ../denest.o

CHECK_PROGS = uniform2d linear_peak small_strip

$(CHECK_PROGS): sample.h $(OBJS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f ../denest.o
	rm -f *.dat
	rm -f $(CHECK_PROGS)