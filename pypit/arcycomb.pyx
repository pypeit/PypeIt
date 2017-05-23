# To get this running, you must do the following at the command line:
# python arcycomb_setup.py build_ext --inplace
# although I'm not really sure what the --inplace does, I think it means "only valid for this directory"

import numpy as np
cimport numpy as np
cimport cython
DTYPE = np.float64
ctypedef np.float_t DTYPE_t

cdef extern from "math.h":
    double csqrt "sqrt" (double)

#@cython.boundscheck(False)
def masked_limitget(np.ndarray[DTYPE_t, ndim=3] array not None,
                  double limvalue,
                  int limtype):

    if limtype < -2 or limtype > 2:
        raise ValueError("Arg 3 of masked_limit can only be only of: -2 (<), -1 (<=), 0 (==), 1 (>=), 2 (>)")

    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double sumv, intv

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    for x in range(sz_x):
        for y in range(sz_y):
            for n in range(nfr):
                if limtype == -2:
                    if array[x,y,n] < limvalue:
                        array[x,y,n] = 1
                    else:
                        array[x,y,n] = 0
                elif limtype == -1:
                    if array[x,y,n] <= limvalue:
                        array[x,y,n] = 1
                    else:
                        array[x,y,n] = 0
                elif limtype == 0:
                    if array[x,y,n] == limvalue:
                        array[x,y,n] = 1
                    else:
                        array[x,y,n] = 0
                elif limtype == 1:
                    if array[x,y,n] >= limvalue:
                        array[x,y,n] = 1
                    else:
                        array[x,y,n] = 0
                elif limtype == 2:
                    if array[x,y,n] > limvalue:
                        array[x,y,n] = 1
                    else:
                        array[x,y,n] = 0
    return array


#@cython.boundscheck(False)
def masked_limitset(np.ndarray[DTYPE_t, ndim=3] array not None,
                  double limvalue,
                  int limtype,
                  double maskvalue):

    if limtype < -2 or limtype > 2:
        raise ValueError("Arg 3 of masked_limit can only be only of: -2 (<), -1 (<=), 0 (==), 1 (>=), 2 (>)")

    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double sumv, intv

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    for x in range(sz_x):
        for y in range(sz_y):
            for n in range(nfr):
                if limtype == -2:
                    if array[x,y,n] < limvalue:
                        array[x,y,n] = maskvalue
                elif limtype == -1:
                    if array[x,y,n] <= limvalue:
                        array[x,y,n] = maskvalue
                elif limtype == 0:
                    if array[x,y,n] == limvalue:
                        array[x,y,n] = maskvalue
                elif limtype == 1:
                    if array[x,y,n] >= limvalue:
                        array[x,y,n] = maskvalue
                elif limtype == 2:
                    if array[x,y,n] > limvalue:
                        array[x,y,n] = maskvalue
    return array


#@cython.boundscheck(False)
def masked_limitsetarr(np.ndarray[DTYPE_t, ndim=3] array not None,
                  np.ndarray[DTYPE_t, ndim=2] limvalue not None,
                  int limtype,
                  double maskvalue):

    if limtype < -2 or limtype > 2:
        raise ValueError("Arg 3 of masked_limit can only be only of: -2 (<), -1 (<=), 0 (==), 1 (>=), 2 (>)")

    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double sumv, intv

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    for x in range(sz_x):
        for y in range(sz_y):
            for n in range(nfr):
                if limtype == -2:
                    if array[x,y,n] < limvalue[x,y]:
                        array[x,y,n] = maskvalue
                elif limtype == -1:
                    if array[x,y,n] <= limvalue[x,y]:
                        array[x,y,n] = maskvalue
                elif limtype == 0:
                    if array[x,y,n] == limvalue[x,y]:
                        array[x,y,n] = maskvalue
                elif limtype == 1:
                    if array[x,y,n] >= limvalue[x,y]:
                        array[x,y,n] = maskvalue
                elif limtype == 2:
                    if array[x,y,n] > limvalue[x,y]:
                        array[x,y,n] = maskvalue
    return array


#@cython.boundscheck(False)
def masked_mean(np.ndarray[DTYPE_t, ndim=3] array not None,
                  double maskvalue):
    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double sumv, intv

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] meanarr = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(sz_x):
        for y in range(sz_y):
            # Fill in the array
            intv = 0.0
            sumv = 0.0
            for n in range(nfr):
                if array[x,y,n] != maskvalue:
                    sumv += array[x,y,n]
                    intv += 1.0
            if intv == 0.0:
                meanarr[x,y] = maskvalue
            else:
                meanarr[x,y] = sumv/intv
    return meanarr


#@cython.boundscheck(False)
def masked_median(np.ndarray[DTYPE_t, ndim=3] array not None,
                  double maskvalue):
    cdef int sz_x, sz_y, nfr
    cdef int i, j, cnt
    cdef int x, y, n
    cdef double temp

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] medarr = np.zeros((sz_x,sz_y), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] arrt = np.zeros(nfr,dtype=DTYPE)

    for x in range(sz_x):
        for y in range(sz_y):
            # Fill in the array
            cnt = 0
            for n in range(nfr):
                if array[x,y,n] != maskvalue:
                    arrt[cnt] = array[x,y,n]
                    cnt += 1
            if cnt == 0:
                medarr[x,y] = maskvalue
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
                medarr[x,y] = 0.5*(arrt[cnt/2] + arrt[cnt/2 - 1])
            else:
                medarr[x,y] = arrt[(cnt-1)/2]
    return medarr


#@cython.boundscheck(False)
def masked_replace(np.ndarray[DTYPE_t, ndim=2] array not None,
                  np.ndarray[DTYPE_t, ndim=2] repvalue not None,
                  double maskvalue):

    cdef int sz_x, sz_y
    cdef int x, y

    sz_x = array.shape[0]
    sz_y = array.shape[1]

    for x in range(sz_x):
        for y in range(sz_y):
            if array[x,y] == maskvalue:
                array[x,y] = repvalue[x,y]
    return array


#@cython.boundscheck(False)
def masked_weightmean(np.ndarray[DTYPE_t, ndim=3] array not None,
                  double maskvalue):
    # This routine uses the weighted mean formula described
    # by Mighell (1999), ApJ, 518, 380
    # and is given by their Eq. (18).

    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double sumv, intv

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] meanarr = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(sz_x):
        for y in range(sz_y):
            # Fill in the array
            intv = 0.0
            sumv = 0.0
            for n in range(nfr):
                if array[x,y,n] != maskvalue:
                    if array[x,y,n] <= 1.0: # Deal with spurious cases
                        sumv += 0.0
                        intv += 1.0
                    else:
                        sumv += csqrt(array[x,y,n])*array[x,y,n]
                        intv += csqrt(array[x,y,n])
#					else:
#						sumv += 1.0
#						intv += 1.0/(1.0+array[x,y,n])
            if intv == 0.0:
                meanarr[x,y] = maskvalue
            else:
                meanarr[x,y] = sumv/intv
    return meanarr


#@cython.boundscheck(False)
def maxnonsat(np.ndarray[DTYPE_t, ndim=3] array not None,
           double satval):

    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double temp, minv

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] mmarr = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(sz_x):
        for y in range(sz_y):
            # Sort the array
            temp = 0.0
            minv = satval
            for n in range(nfr):
                if array[x,y,n] > temp and temp < satval:
                    temp = array[x,y,n]
                if array[x,y,n] < minv:
                    minv = array[x,y,n]
            if temp == 0.0:
                mmarr[x,y] = minv
            else:
                mmarr[x,y] = temp
    return mmarr


#@cython.boundscheck(False)
def mean(np.ndarray[DTYPE_t, ndim=3] array not None):
    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double sumv, intv

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] meanarr = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(sz_x):
        for y in range(sz_y):
            # Fill in the array
            sumv = 0.0
            intv = 0.0
            for n in range(nfr):
                sumv += array[x,y,n]
                intv += 1.0
            meanarr[x,y] = sumv/intv
    return meanarr


#@cython.boundscheck(False)
def median(np.ndarray[DTYPE_t, ndim=3] array not None):
    cdef int sz_x, sz_y, nfr
    cdef int x, y, n, j
    cdef double temp

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] medarr = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(sz_x):
        for y in range(sz_y):
            # Sort the array
            for n in range(nfr-1):
                for j in range(n+1,nfr):
                    if array[x,y,j] < array[x,y,n]:
                        temp = array[x,y,n]
                        array[x,y,n] = array[x,y,j]
                        array[x,y,j] = temp
            # Find the median value
            if nfr%2==0:
                medarr[x,y] = 0.5*(array[x,y,nfr/2] + array[x,y,nfr/2 - 1])
            else:
                medarr[x,y] = array[x,y,(nfr-1)/2]
    return medarr

#@cython.boundscheck(False)
def minmax(np.ndarray[DTYPE_t, ndim=3] array not None,
           int mm):

    if mm != 0 and mm != 1:
        raise ValueError("Arg 2 of minmax can only be 0 (for min) or 1 (for max)")

    cdef int sz_x, sz_y, nfr
    cdef int x, y, n
    cdef double temp

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    nfr  = array.shape[2]

    cdef np.ndarray[DTYPE_t, ndim=2] mmarr = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(sz_x):
        for y in range(sz_y):
            # Sort the array
            temp = array[x,y,0]
            for n in range(1,nfr):
                if mm == 0:
                    if array[x,y,n] < temp:
                        temp = array[x,y,n]
                else:
                    if array[x,y,n] > temp:
                        temp = array[x,y,n]
            mmarr[x,y] = temp
    return mmarr
