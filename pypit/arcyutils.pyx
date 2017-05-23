# To get this running, you must do the following at the command line:
# python arcycomb_setup.py build_ext --inplace
# although I'm not really sure what the --inplace does, I think it means "only valid for this directory"

import numpy as np
cimport numpy as np
cimport cython
DTYPE = np.float64
ctypedef np.float_t DTYPE_t
ITYPE = np.int64
ctypedef np.int_t ITYPE_t

cdef extern from "math.h":
    double csqrt "sqrt" (double)
    double cexp "exp" (double)
    double ccos "cos" (double)
    double csin "sin" (double)
    double cpow "pow" (double, double)
    double cmin "fmin" (double, double)

cdef extern from "gsl/gsl_multifit.h":
    ctypedef struct gsl_matrix:
        pass
    ctypedef struct gsl_vector:
        pass
    ctypedef struct gsl_multifit_linear_workspace:
        pass
    gsl_matrix *gsl_matrix_alloc(int,int)
    gsl_vector *gsl_vector_alloc(int)
    gsl_multifit_linear_workspace *gsl_multifit_linear_alloc(int,int)
    void gsl_matrix_set(gsl_matrix *A,int,int,double)
    void gsl_vector_set(gsl_vector *v,int,double)
    double gsl_vector_get(gsl_vector *v,int)
    void gsl_multifit_linear(gsl_matrix *A, gsl_vector *u, gsl_vector *v, gsl_matrix *B, double *d, gsl_multifit_linear_workspace *w)
    void gsl_multifit_linear_est(gsl_vector *u, gsl_vector *c, gsl_matrix *A, double *d, double *e)
    void gsl_multifit_wlinear(gsl_matrix *A, gsl_vector *e, gsl_vector *u, gsl_vector *v, gsl_matrix *B, double *d, gsl_multifit_linear_workspace *w)
    void gsl_multifit_linear_free(gsl_multifit_linear_workspace *w)
    void gsl_matrix_free(gsl_matrix *A)
    void gsl_vector_free(gsl_vector *v)

cdef extern from "gsl/gsl_bspline.h":
    ctypedef struct gsl_bspline_workspace

    gsl_bspline_workspace *gsl_bspline_alloc(int,int)
    void gsl_bspline_knots_uniform(double, double, gsl_bspline_workspace *bw)
    void gsl_bspline_eval(double, gsl_vector *v, gsl_bspline_workspace *bw)
    void gsl_bspline_free(gsl_bspline_workspace *bw)

cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng_type
    ctypedef struct gsl_rng

    cdef gsl_rng_type *gsl_rng_default
    gsl_rng *gsl_rng_alloc ( gsl_rng_type * T) nogil
    void gsl_rng_free (gsl_rng * r) nogil
    gsl_rng_type * gsl_rng_env_setup () nogil


#@cython.boundscheck(False)
def func2d_fit_val(np.ndarray[DTYPE_t, ndim=1] yarr not None,
                    np.ndarray[DTYPE_t, ndim=2] frame not None,
                    np.ndarray[DTYPE_t, ndim=2] weights not None,
                    int polyorder):

    cdef int x, sz_x, y, sz_y, j
    cdef int degree = polyorder+1
    cdef double mval, chisq

    sz_x = frame.shape[0]
    sz_y = frame.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] image = np.zeros((sz_x,sz_y), dtype=DTYPE)

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat

    cdef gsl_vector *yvec
    cdef gsl_vector *wgt
    cdef gsl_vector *c

    xmat = gsl_matrix_alloc(sz_y, degree)
    wgt  = gsl_vector_alloc(sz_y)
    yvec = gsl_vector_alloc(sz_y)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)
    ws = gsl_multifit_linear_alloc(sz_y, degree)

    for y in range(sz_y):
        gsl_matrix_set(xmat, y, 0, 1.0)
        for j in range(1,degree):
            gsl_matrix_set(xmat, y, j, cpow(yarr[y], j))

    for x in range(sz_x):
        # Set the arrays
        for y in range(sz_y):
            gsl_vector_set(yvec, y, frame[x,y])
            gsl_vector_set(wgt, y, weights[x,y])
        # Fit with a polynomial
        gsl_multifit_wlinear(xmat, wgt, yvec, c, cov, &chisq, ws)
        # Obtain the model values for this row
        for y in range(sz_y):
            mval = 0.0
            for j in range(0,degree): mval += gsl_vector_get(c, j) * cpow(yarr[y], j)
            image[x,y] = mval

    # Free some of these arrays from memory
    gsl_multifit_linear_free(ws)
    gsl_matrix_free(xmat)
    gsl_matrix_free(cov)
    gsl_vector_free(yvec)
    gsl_vector_free(c)

    # Return an image
    return image


#@cython.boundscheck(False)
def grow_masked(np.ndarray[DTYPE_t, ndim=2] img not None,
              double grow, double growval):

    cdef int x, y, sz_x, sz_y
    cdef int i, j, mnx, mxx, mny, mxy
    sz_x = img.shape[0]
    sz_y = img.shape[1]
    cdef int gval = <int>(1.0+grow)

    cdef np.ndarray[DTYPE_t, ndim=2] imgout = img.copy()

    # Grow any masked values by the specified amount
    for x in range(0,sz_x):
        for y in range(0,sz_y):
            if img[x,y] != growval: continue
            # Set x limits
            mnx = x-gval
            mxx = x+gval+1
            if mnx < 0: mnx = 0
            if mxx > sz_x: mxx = sz_x
            # Set y limits
            mny = y-gval
            mxy = y+gval+1
            if mny < 0: mny = 0
            if mxy > sz_y: mxy = sz_y
            for i in range(mnx,mxx):
                for j in range(mny, mxy):
                    if csqrt(<double>((i-x)*(i-x)+(j-y)*(j-y))) <= grow:
                        imgout[i,j] = growval
    return imgout


#@cython.boundscheck(False)
def mean_weight(np.ndarray[DTYPE_t, ndim=1] array not None,
                np.ndarray[DTYPE_t, ndim=1] weight not None,
                int rejhilo, double maskval):
    cdef int sz_x
    cdef int x, y
    cdef double temp
    cdef double wght = 0.0

    sz_x = array.shape[0]

    if rejhilo != 0:
        # Sort the array
        for x in range(sz_x-1):
            for y in range(x,sz_x):
                # Sort the array
                if array[y] < array[x]:
                    temp = array[x]
                    array[x] = array[y]
                    array[y] = temp
    temp = 0.0
    for x in range(rejhilo,sz_x-rejhilo):
        temp += (weight[x]*array[x])
        wght += weight[x]
    if wght != 0.0: return temp/wght
    return maskval


#@cython.boundscheck(False)
def order_pixels(np.ndarray[DTYPE_t, ndim=3] pixlocn not None,
                np.ndarray[DTYPE_t, ndim=2] lord not None,
                np.ndarray[DTYPE_t, ndim=2] rord not None):
    """
    Based on physical pixel locations, determine which pixels are within the orders
    """
    cdef int x, y, o
    cdef int sz_x, sz_y, sz_o

    sz_x = pixlocn.shape[0]
    sz_y = pixlocn.shape[1]
    sz_o = lord.shape[1]

    cdef np.ndarray[ITYPE_t, ndim=2] outfr = np.zeros((sz_x, sz_y), dtype=ITYPE)

    # First determine the size of the final arrays
    for x in range(sz_x):
        for y in range(sz_y):
            for o in range(sz_o):
                if (lord[x,o] < rord[x,o]) and (pixlocn[x,y,1]>lord[x,o]) and (pixlocn[x,y,1]<rord[x,o]):
                    outfr[x,y] = o+1
                    break # Speed up the routine a little by only assigning a single order to a given pixel
                elif (lord[x,o] > rord[x,o]) and (pixlocn[x,y,1]<lord[x,o]) and (pixlocn[x,y,1]>rord[x,o]):
                    outfr[x,y] = o+1
                    break # Speed up the routine a little by only assigning a single order to a given pixel
    return outfr


#@cython.boundscheck(False)
def smooth_x(np.ndarray[DTYPE_t, ndim=2] array not None,
            np.ndarray[DTYPE_t, ndim=2] weight not None,
            int fact, int rejhilo, double maskval):
    cdef int sz_x, sz_y
    cdef int x, y, b
    fact /= 2

    sz_x = array.shape[0]
    sz_y = array.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] smtarr = np.zeros((sz_x,sz_y), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] medarr = np.zeros((2*fact+1), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] wgtarr = np.zeros((2*fact+1), dtype=DTYPE)

    for y in range(sz_y):
        for x in range(sz_x):
            for b in range(2*fact+1):
                if (x+b-fact < 0) or (x+b-fact >= sz_x): wgtarr[b] = 0.0
                else:
                    medarr[b] = array[x+b-fact,y]
                    wgtarr[b] = weight[x+b-fact,y]
            smtarr[x,y] = mean_weight(medarr, wgtarr, rejhilo, maskval)
    return smtarr
