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
    double clog "log" (double)
    double cpow "pow" (double, double)


#@cython.boundscheck(False)
def order_saturation(np.ndarray[ITYPE_t, ndim=2] satmask not None,
                    np.ndarray[ITYPE_t, ndim=2] ordcen not None,
                    np.ndarray[ITYPE_t, ndim=2] ordwid not None):

    cdef int sz_x, sz_o, sz_y, x, y, o
    cdef int flag, xmin, xmax
    sz_x = satmask.shape[1]
    sz_y = satmask.shape[0]
    sz_o = ordcen.shape[1]

    cdef np.ndarray[ITYPE_t, ndim=2] ordsat = np.zeros((sz_y, sz_o), dtype=ITYPE)

    for o in range(sz_o):
        for y in range(sz_y):
            xmin = ordcen[y,o] - ordwid[y,o]
            xmax = ordcen[y,o] + ordwid[y,o] + 1
            if xmin < 0:
                xmin = 0
            if xmax >= sz_x:
                xmax=sz_x
            for x in range(xmin, xmax):
                if satmask[y,x] == 1:
                    ordsat[y,o] = 1
                    break
    return ordsat


#@cython.boundscheck(False)
def saturation_mask(np.ndarray[DTYPE_t, ndim=2] array not None,
                    double satlevel):
   
    cdef int x, y, sx, sy
    cdef int sz_x, sz_y
    cdef double localx, localy

    cdef double satdown = 1.001
    sz_x = array.shape[0]
    sz_y = array.shape[1]

    cdef np.ndarray[ITYPE_t, ndim=2] satmask = np.zeros((array.shape[0],array.shape[1]), dtype=ITYPE)

    for y in range (0,sz_y):
        for x in range(0,sz_x):
            if array[x,y] >= satlevel and satmask[x,y] == 0:
                satmask[x,y] = 1
                sy = 0
                localy = array[x,y+sy]
                while True:
                    satmask[x,y+sy] = 1
                    sx = 1
                    localx = array[x+sx,y+sy]
                    while True:
                        satmask[x+sx,y+sy] = 1
                        sx += 1
                        if x+sx > sz_x-1: break
                        if array[x+sx,y+sy] >= localx/satdown and array[x+sx,y+sy]<satlevel: break
                        localx = array[x+sx,y+sy]
                    sx = -1
                    localx = array[x+sx,y+sy]
                    while True:
                        satmask[x+sx,y+sy] = 1
                        sx -= 1
                        if x+sx < 0: break
                        if array[x+sx,y+sy] >= localx/satdown and array[x+sx,y+sy]<satlevel: break
                        localx = array[x+sx,y+sy]
                    sy += 1
                    if array[x,y+sy] >= localy/satdown and array[x,y+sy] < satlevel: break
                    localy = array[x,y+sy]
                sy = -1
                localy = array[x,y+sy]
                while True:
                    satmask[x,y+sy] = 1
                    sx = 1
                    localx = array[x+sx,y+sy]
                    while True:
                        satmask[x+sx,y+sy] = 1
                        sx += 1
                        if x+sx > sz_x-1: break
                        if array[x+sx,y+sy] >= localx/satdown and array[x+sx,y+sy]<satlevel: break
                        localx = array[x+sx,y+sy]
                    sx = -1
                    localx = array[x+sx,y+sy]
                    while True:
                        satmask[x+sx,y+sy] = 1
                        sx -= 1
                        if x+sx < 0: break
                        if array[x+sx,y+sy] >= localx/satdown and array[x+sx,y+sy]<satlevel: break
                        localx = array[x+sx,y+sy]
                    sy -= 1
                    if array[x,y+sy] >= localy/satdown and array[x,y+sy] < satlevel: break
                    localy = array[x,y+sy]
    return satmask
