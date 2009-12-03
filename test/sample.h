#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#include<gsl/gsl_rng.h>

typedef double *(*sample)(gsl_rng *);

double **do_sample(int npdfsample, int nsample, sample samp);

#endif /* __SAMPLE_H__ */
