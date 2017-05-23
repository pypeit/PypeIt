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
    void gsl_multifit_wlinear(gsl_matrix *A, gsl_vector *e, gsl_vector *u, gsl_vector *v, gsl_matrix *B, double *d, gsl_multifit_linear_workspace *w)
    void gsl_multifit_linear_free(gsl_multifit_linear_workspace *w)
    void gsl_matrix_free(gsl_matrix *A)
    void gsl_vector_free(gsl_vector *v)


#@cython.boundscheck(False)
def cr_screen(np.ndarray[DTYPE_t, ndim=2] frame not None, double maskval):
    """
    Perofrm a masked median filter along the spatial
    direction to identify candidate cosmic rays.
    """
    cdef int sz_x, sz_y
    cdef int i, j, cnt
    cdef int x, y
    cdef double temp

    sz_x = frame.shape[0]
    sz_y = frame.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] medarr = np.zeros((sz_x), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] madarr = np.zeros((sz_x), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] arrt = np.zeros((sz_y),dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] sigimg = np.zeros((sz_x,sz_y),dtype=DTYPE)

    # Get the median value along the spatial direction
    for x in range(sz_x):
        # Fill in the array
        cnt = 0
        for y in range(sz_y):
            if frame[x,y] != maskval:
                arrt[cnt] = frame[x,y]
                cnt += 1
        if cnt == 0:
            medarr[x] = maskval
            continue
        # Sort the array
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if arrt[j] < arrt[i]:
                    temp = arrt[i]
                    arrt[i] = arrt[j]
                    arrt[j] = temp
        # Find the median value
        if cnt%2==0:
            medarr[x] = 0.5*(arrt[cnt/2] + arrt[cnt/2 - 1])
        else:
            medarr[x] = arrt[(cnt-1)/2]

    # Get the median absolute deviation along the spatial direction
    for x in range(sz_x):
        # Fill in the array
        cnt = 0
        for y in range(sz_y):
            if frame[x,y] != maskval:
                temp = frame[x,y]-medarr[x]
                if temp<0.0: temp *= -1.0
                arrt[cnt] = temp
                cnt += 1
        if cnt == 0:
            madarr[x] = maskval
            continue
        # Sort the array
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if arrt[j] < arrt[i]:
                    temp = arrt[i]
                    arrt[i] = arrt[j]
                    arrt[j] = temp
        # Find the median value
        if cnt%2==0:
            madarr[x] = 1.4826 * 0.5*(arrt[cnt/2] + arrt[cnt/2 - 1])
        else:
            madarr[x] = 1.4826 * arrt[(cnt-1)/2]

    # Return an image of significance detections
    for x in range(sz_x):
        for y in range(sz_y):
            if madarr[x] == maskval: sigimg[x,y] = maskval
            else:
                temp = (frame[x,y]-medarr[x])/madarr[x]
                if temp < 0.0: temp *= -1.0
                sigimg[x,y] = temp
    return sigimg
