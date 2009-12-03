GSL_PREFIX = /sw

CFLAGS += -I$(GSL_PREFIX)/include -L$(GSL_PREFIX)/lib -lgslcblas -lgsl
LFLAGS += -L$(GSL_PREFIX)/lib -lgslcblas -lgsl

OBJS = sample.o rng_init.o ../denest.o

uniform2d: sample.h $(OBJS)
linear_peak: sample.h $(OBJS)

.PHONY: clean
clean:
	rm -rf *.o
	rm uniform2d