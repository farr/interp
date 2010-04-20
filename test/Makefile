# Configuration options:
GSL_DIR = /opt/local
ARCH_OPT = -march=core2 # comment out if not on core2 processor. 

GSL_LIB_DIR = $(GSL_DIR)/lib
GSL_INC_DIR = $(GSL_DIR)/include

GSL_LIBS = -lgslcblas -lgsl

COPTFLAGS = -O3 $(ARCH_OPT)
CDEBUGFLAGS = -ggdb

#######################################
# End config---stuff below shouldn't need to be changed.

TEST_OBJS = test_gaussian.o

UTILITY_OBJS = randomize.o stats.o

CFLAGS = $(COPTFLAGS) $(CDEBUGFLAGS) -I$(GSL_INC_DIR)

LDFLAGS = -L$(GSL_LIB_DIR)

LDLIBS = $(GSL_LIBS)

all-tests: ../interp.o $(UTILITY_OBJS) $(TEST_OBJS)

%.o: ../interp.h utility.h

.PHONY: clean
clean:
	rm -f ../interp.o
	rm -f $(TEST_OBJS)
	rm -f all-tests.o all-tests