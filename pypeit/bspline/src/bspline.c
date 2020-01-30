/*
Support algorithms for bspline.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "bspline.h"


int column_to_row_major_index(int k, int nr, int nc) {
    /*
    Convert a flattened index in a column-major stored array into the
    flattened index in row-major stored array.
    
    Args:
        k:
            Index in the flattened column-major array.
        nr:
            Number of array rows.
        nc:
            Number of array columns.

    Returns:
        Index in the flattened row-major array.
    */
    return (k - (k/nr)*nr)*nc + k/nr;
}


void flat_row_major_indices(int k, int nr, int nc, int *i, int *j) {
    /*
    Convert the flattened index in a row-major stored array into the
    indices along the two axes.

    Args:
        k:
            Index in the flattened row-major array.
        nr:
            Number of array rows.
        nc:
            Number of array columns.
        i:
            Index along the first axis of the row-major stored array.
            Value pointed to is replaced.
        j:
            Index along the second axis of the row-major stored
            array. Value pointed to is replaced.
    */
    *i = k/nc;
    *j = k - *i*nc;
}


int* upper_triangle(int n, bool upper_left) {
    /*
    Compute the indices in a flattened 2D array that contain the
    upper triangle.

    To get the indices for the *upper-left* triangle, set
    upper_left=true.

    For example, the upper-right triangle of a 3x3 matrix has
    flattened indices of: 0, 1, 2, 4, 5, 8. The upper-left indices
    are: 0, 1, 2, 3, 4, 6.

    Args:
        n:
            Size of the square array.
        upper_left:
            Return the upper-left triangle instead of the
            upper-right. The upper triangle is typically defined as
            the upper-right triangle.

    Returns:
        Returns the pointer to the integer array with the indices.
    */
    int i, j;
    int nbi = n*(n+1)/2;
    int* bi = (int*) malloc (nbi * sizeof(int));
    int s = 0;
    for (i = n; i > 0; --i)
        for (j = 0; j < i; ++j) {
            bi[s] = upper_left ? j+(n-i)*n : j+(n-i)*(n+1);
            s += 1;
        }
    return bi;
}


void bspline_model(double *action, long *lower, long *upper, double *coeff, int n, int nord,
                   int npoly, int nd, double *yfit) {
    /*
    Calculate the bspline model.

    Args:
        action:
            Action matrix. See pypeit.bspline.bspline.bspline.action.
            The shape of the array is expected to be ``nd`` by
            ``npoly*nord``.
        lower:
            Vector with the starting indices along the second axis of
            action used to construct the model.
        upper:
            Vector with the (inclusive) ending indices along the
            second axis of action used to construct the model.
        coeff:
            The model coefficients used for each action.
        n:
            Number of unmasked measurements included in the fit.
        nord:
            Fit order.
        npoly:
            Polynomial per fit order.
        nd:
            Total number of data points.
        yfit:
            Pointer to the memory location for the bspline model.
            Memory must have already been allocated.
    */
    int nn = npoly*nord;    // This is the same as the number of columns in action
    int mm = n - nord+1;
    int i, j, k;
    for (i = 0; i < mm; ++i) {
        if (!(upper[i]+1 > lower[i]))
            continue;
        for (j = lower[i]; j <= upper[i]; ++j) {
            yfit[j] = 0;
            for (k = 0; k < nn; ++k)
                yfit[j] += action[k*nd + j] * coeff[i*npoly + k];
        }
    }
}

void intrv(int nord, double *breakpoints, int nb, double *x, int nx, long *indx) {
    /*
    Find the segment between breakpoints which contain each value in
    the array x.

    The minimum breakpoint is nbkptord -1, and the maximum is nbkpt -
    nbkptord - 1.

    Args:
        nord:
            Order of the fit.
        breakpoints:
            Locations of good breakpoints
        x:
            Data values, assumed to be monotonically increasing.
        indx:
            Replaced on output: the break-point segments.
    */
    int n = nb - nord;
    int ileft = nord - 1;
    int i;
    for (i = 0; i < nx; ++i) {
        while ((x[i] > breakpoints[ileft+1]) & (ileft < n - 1))
            ileft += 1;
        indx[i] = ileft;
    }
}


void solution_arrays(int nn, int npoly, int nord, int nd, double *ydata, double *ivar,
                     double *action, long *upper, long *lower, double *alpha, int ar,
                     double *beta, int bn) {
    /*
    Support function that builds the arrays for Cholesky
    decomposition.

    Array size requirements:

        - action is 2D, first axis is the same length as ydata and
          ivar, the second axis is npoly*nord.

        - The second axis of alpha and the length of beta are the
          same length as the second axis of action (npoly*nord).

        - BEWARE that the type of upper and lower must match the
          input: np.int32 for int and np.int64 for long. Current
          input must be long.

    Args:
        nn:
            Number of good break points.
        npoly:
            Polynomial per fit order.
        nord:
            Fit order.
        nd:
            Total number of data points.
        ydata:
            Data to fit
        ivar:
            Inverse variance in the data to fit.
        action:
            Action matrix. See pypeit.bspline.bspline.bspline.action.
            The shape of the array is expected to be ``nd`` by
            ``npoly*nord``.
        upper:
            Vector with the (inclusive) ending indices along the
            second axis of action used to construct the model.
        lower:
            Vector with the starting indices along the second axis of
            action used to construct the model.
        alpha:
            Solution matrix for Cholesky decomposition. Memory must
            have already been allocated.
        ar:
            Number of rows (first axis) in alpha.
        beta:
            Solution vector for Cholesky decomposition. Memory must
            have already been allocated.
        bn:
            Number of elements in beta. Same as the number of columns
            in alpha.
    */
    // Get the upper triangle indices
    int bw = npoly * nord;      // This is the length of the second axis of action
    int nbi = bw*(bw+1)/2;
    int *bi = upper_triangle(bw, false);
    int *bo = upper_triangle(bw, true);

    int i, j, k;
    int ii, jj, kk;
    int itop;

    // Convenience data
    // TODO: These are big allocations.  Can we avoid them?
    double *ierr = (double*) malloc (nd * sizeof(double));
    double *a2 = (double*) malloc (nd*bw * sizeof(double));
    for (i = 0; i < nd; ++i) {
        ierr[i] = sqrt(ivar[i]);
        for (j = 0; j < bw; ++j)
            a2[i*bw + j] = action[j*nd + i] * ierr[i];
    }

    // Zero input arrays
    // TODO: Just assume this?
    for (i = 0; i < ar; ++i)
        for (j = 0; j < bn; ++j) {
            beta[j] = 0;
            alpha[i*bn+j] = 0;
        }

    // Construct alpha and beta
    for (k = 0; k < nn-nord+1; ++k) {
        if (!(upper[k]+1 > lower[k]))
            continue;

        itop = k*npoly;
        for (i = 0; i < nbi; ++i) {
            kk = column_to_row_major_index(bo[i]+itop*bw, ar, bn);
            flat_row_major_indices(bi[i], nd, bw, &ii, &jj);
            for (j = lower[k]; j <= upper[k]; ++j)
                alpha[kk] += a2[j*bw+ii] * a2[j*bw+jj];
        }
        for (i = 0; i < bw; ++i)
            for (j = lower[k]; j <= upper[k]; ++j)
                beta[itop+i] += ydata[j] * ierr[j] * a2[j*bw + i];
    }
    // Free memory
    free(a2);
    free(ierr);
    free(bo);
    free(bi);
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
    int i, j, k, s;
    int kn = lr - 1;
    int n = lc - lr;

    int nbi = kn*(kn+1)/2;
    int *bi = upper_triangle(kn, false);

    int *here = (int*) malloc (nbi*n * sizeof(int));
    for (i = 0; i < nbi; ++i)
        for (j = 0; j < n; ++j) {
            k = i*n+j;
            here[k] = bi[i] + (j+1)*lr;
            // Because I've given up on trying to figure out how to more
            // simply write the indices in the untransposed matrix! This
            // converts the indices from a column-major order to a
            // row-major order.
            here[k] = column_to_row_major_index(here[k], lr, lc);
        }

    double *hmm = (double*) malloc (nbi * sizeof(double));
    for (j = 0; j < n; ++j) {
        if (lower[j] <= 0) {
            // Would result in a NaN or Inf calculation
            // Free the memory
            free(hmm);
            free(here);
            free(bi);
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
    free(hmm);
    free(here);
    free(bi);
    // Return success
    return -1;
}


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

