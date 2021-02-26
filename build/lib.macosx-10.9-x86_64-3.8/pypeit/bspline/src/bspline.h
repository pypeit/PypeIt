
#ifndef _BSPLINE_H_
#define _BSPLINE_H_

#include <stdbool.h>
#include <stdint.h>

int32_t column_to_row_major_index(int32_t k, int32_t nr, int32_t nc);
void flat_row_major_indices(int32_t k, int32_t nr, int32_t nc, int32_t *i, int32_t *j);
int32_t* upper_triangle(int32_t kn, bool upper_left);
void bspline_model(double *action, int64_t *lower, int64_t *upper, double *coeff,
                   int32_t n, int32_t nord, int32_t npoly, int32_t nd, double *yfit);
void intrv(int32_t nord, double *breakpoints, int32_t nb, double *x, int32_t nx, int64_t *indx);
void solution_arrays(int32_t nn, int32_t npoly, int32_t nord, int32_t nd, double *ydata,
                     double *ivar, double *action, int64_t *upper, int64_t *lower,
                     double *alpha, int32_t ar, double *beta, int32_t bn);
void cholesky_solve(double *a, int32_t ar, int32_t ac, double *b, int32_t bn);
int cholesky_band(double *lower, int32_t lr, int32_t lc);

#endif // _BSPLINE_H_

