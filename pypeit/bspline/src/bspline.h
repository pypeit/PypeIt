
#ifndef _BSPLINE_H_
#define _BSPLINE_H_

#include <stdbool.h>

int column_to_row_major_index(int k, int nr, int nc);
void flat_row_major_indices(int k, int nr, int nc, int *i, int *j);
int* upper_triangle(int kn, bool upper_left);
void bspline_model(double *action, long *lower, long *upper, double *coeff, int n, int nord,
                   int npoly, int nd, double *yfit);
void intrv(int nord, double *breakpoints, int nb, double *x, int nx, long *indx);
void solution_arrays(int nn, int npoly, int nord, int nd, double *ydata, double *ivar,
                     double *action, long *upper, long *lower, double *alpha, int ar,
                     double *beta, int bn);
void cholesky_solve(double *a, int ar, int ac, double *b, int bn);
int cholesky_band(double *lower, int lr, int lc);

#endif // _BSPLINE_H_

