#ifndef __SAMPLING_H_
#define __SAMPLING_H_

#include <stdlib.h>

void build_multinomial_cdf(double *cdf, double *pdf, size_t pdf_size);
size_t sample_multinomial_cdf(double *cdf, size_t cdf_size);

int sample_from_discrete_log_pd(double *log_pd, int log_pd_size);

void kitagawa_resampling(double *w, size_t n_weights,
                         size_t *s_res, size_t n_res);

size_t sample_uniform(size_t pdf_size);

#endif // __SAMPLING_H_

