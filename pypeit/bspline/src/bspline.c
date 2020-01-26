/*
   Support algorithms for bspline.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "bspline.h"

void cholesky_solve(double *a, int ar, int ac, double *b, int bn) {
    /*
       Solve the equation Ax=b where A is a Cholesky-banded matrix.

    Args:
        a: 
            The flattened array A used in the equation A x = b.  The
            number of columns (2nd axis) in a is given by ``ac``.

        ar:
            The number of rows in the A matrix.

        ac:
            The number of columns in the A matrix.

        b:
            The vector b in the equation A x = b.  The values are
            replaced by the solution.

        bn:
            The number of elements in the b vector.
    */
    int n = bn - ar;
    int kn = ar - 1;
    int i, j;
    double s;
    for (j = 0; j < n; ++j) {
        b[j] /= a[j];
        for (i = 1; i < kn+1; ++i)
            b[j+i] -= b[j]*a[i*ac + j];
    }
    for (j = n-1; j > -1; --j) {
        s = 0;
        for (i = 1; i < kn+1; ++i)
            s += a[i*ac+j] * b[j+i];
        b[j] = (b[j] - s)/a[j];
    }
}

int cholesky_band(double *lower, int lr, int lc) {
    /*
       Compute the Cholesky decomposition of banded matrix.

    Args:
        lower: 
            The flattened matrix on which to perform the Cholesky
            decomposition.  The input matrix is replaced.

        lr:
            Number of rows (1st axis) in the matrix.

        lc:
            Number of columns (2st axis) in the matrix.

    Returns:
        int: Returns an integer.  If the decomposition was successful,
        the integer is -1.  Otherwise, the integer is the index of the
        column that contains a problem for the decomposition.
    */
    int i, j, k;
    int kn = lr - 1;
    int n = lc - lr;

    int nbi = kn*(kn+1)/2;

    int *bi = (int*) malloc (nbi * sizeof(int));
    int s = 0;
    for (i = kn; i > 0; --i)
        for (j = 0; j < i; ++j) {
            bi[s] = j+(kn-i)*(kn+1);
            s += 1;
        }

    int *here = (int*) malloc (nbi*n * sizeof(int));
    for (i = 0; i < nbi; ++i)
        for (j = 0; j < n; ++j) {
            k = i*n+j;
            here[k] = bi[i] + (j+1)*lr;
            // Because I've given up on trying to figure out how to more
            // simply write the indices in the untransposed matrix! This
            // converts the indices from a column-major order to a
            // row-major order.
            here[k] = (here[k] - (here[k]/lr)*lr)*lc + here[k]/lr;
        }
            
    double *hmm = (double*) malloc (nbi * sizeof(double));
    for (j = 0; j < n; ++j) {
        if (lower[j] <= 0) {
            // Would result in a NaN or Inf calculation
            // Free the memory
            free(bi);
            free(here);
            free(hmm);
            // Return the column that causes the problem
            return j;
        }
        lower[j] = sqrt(lower[j]);
        for (i = 1; i < lr; ++i)
            lower[i*lc + j] /= lower[j];
        s = 0;
        for (i = 0; i < kn; ++i)
            for (k = i; k < kn; ++k)
                hmm[s++] = lower[(k+1)*lc+j] * lower[(i+1)*lc+j];
        for (i = 0; i < nbi; ++i)
            lower[here[i*n+j]] -= hmm[i];
    }

    // Free the memory
    free(bi);
    free(here);
    free(hmm);
    // Return success
    return -1;
}



