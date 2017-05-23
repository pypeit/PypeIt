import numpy as np
cimport numpy as np
cimport cython
DTYPE = np.float64
ctypedef np.float_t DTYPE_t
ITYPE = np.int64
ctypedef np.int_t ITYPE_t


cdef extern from "math.h":
    double csqrt "sqrt" (double)
    double catan "atan" (double)
    double ccos "cos" (double)
    double cpow "pow" (double, double)
    double cerf "erf" (double)

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
def point_in_poly(np.ndarray[DTYPE_t, ndim=2] poly not None,
                    double x, double y):

    cdef int p, sz_p, inside
    cdef double tst, p1x, p2x, p1y, p2y, xints

    sz_p = poly.shape[0]
    inside = 0

    p1x = poly[0,0]
    p1y = poly[0,1]
    for p in range(sz_p+1):
        p2x = poly[p%sz_p,0]
        p2y = poly[p%sz_p,1]
        if p2y < p1y:
            tst = p2y
        else:
            tst = p1y
        if y > tst:
            if p2y < p1y:
                tst = p1y
            else:
                tst = p2y
            if y <= tst:
                if p1x < p2x:
                    tst = p2x
                else:
                    tst = p1x
                if x <= tst:
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        if inside == 0: inside = 1
                        else: inside = 0
        p1x,p1y = p2x,p2y
    return inside


#@cython.boundscheck(False)
def poly_area(np.ndarray[DTYPE_t, ndim=2] p not None):
    cdef int x, sz_x
    sz_x = p.shape[0]

    cdef double area = p[sz_x-1,0]*p[0,1]
    for x in range(0,sz_x-1):
        area += p[x,0]*p[x+1,1]
    for x in range(0,sz_x-1):
        area -= p[x+1,0]*p[x,1]
    area -= p[0,0]*p[sz_x-1,1]
    if area < 0.0: area *= -0.5
    else: area *= 0.5
    return area


#@cython.boundscheck(False)
def polyfit(np.ndarray[DTYPE_t, ndim=1] x not None,
        np.ndarray[DTYPE_t, ndim=1] y not None,
        int pmin, int pmax, int degree, double offset,
        np.ndarray[DTYPE_t, ndim=1] coeffs not None):

    cdef int sz_x
    cdef int i, j
    cdef double chisq

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat
    cdef gsl_vector *yvec
    cdef gsl_vector *c

    sz_x = pmax-pmin

    #cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(degree, dtype=DTYPE)

    xmat = gsl_matrix_alloc(sz_x, degree)
    yvec = gsl_vector_alloc(sz_x)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)

    for i in range(0,sz_x):
        gsl_matrix_set(xmat, i, 0, 1.0)
        for j in range(0,degree):
            gsl_matrix_set(xmat, i, j, cpow(x[pmin+i]-offset, j))
        gsl_vector_set(yvec, i, y[pmin+i])

    # Fit with a polynomial
    ws = gsl_multifit_linear_alloc(sz_x, degree)
    gsl_multifit_linear(xmat, yvec, c, cov, &chisq, ws)

    # Store the result in the coeffs array
    for i in range(0,degree):
        coeffs[i] = gsl_vector_get(c, i)

    # Free some of these arrays from memory
    gsl_multifit_linear_free(ws)
    gsl_matrix_free(xmat)
    gsl_matrix_free(cov)
    gsl_vector_free(yvec)
    gsl_vector_free(c)

    # Return the best-fitting chi-squared value
    return chisq


#@cython.boundscheck(False)
def rectify(np.ndarray[DTYPE_t, ndim=2] frame not None,
            np.ndarray[ITYPE_t, ndim=2] mask not None,
            np.ndarray[ITYPE_t, ndim=1] pixcen not None,
            np.ndarray[ITYPE_t, ndim=1] pixledg not None,
            np.ndarray[ITYPE_t, ndim=1] pixredg not None,
            int pixwid, double maskval):

    cdef int x, y, sz_x, sz_y
    cdef int yidx, ymin, ymax, sz_p
    cdef int shift = pixwid/2

    sz_x  = frame.shape[0]
    sz_y  = frame.shape[1]
    sz_p = 2*shift+1

    cdef np.ndarray[DTYPE_t, ndim=2] rectify = np.zeros((sz_x,sz_p), dtype=DTYPE)

    # Rectify the image
    for x in range(0,sz_x):
        for y in range(0,sz_p):
            if pixcen[x] <= 0:
                yidx = y+(pixredg[x]-shift)-shift
            elif pixcen[x] >= sz_y-1:
                yidx = y+(pixledg[x]+shift)-shift
            else:
                yidx = y+pixcen[x]-shift
            if pixredg[x] <= 0 or pixledg[x] >= sz_y-1 or yidx < 0 or yidx >= sz_y:
                rectify[x,y] = maskval
            else:
                if mask[x,yidx] == 0: rectify[x,y] = maskval
                else: rectify[x,y] = frame[x,yidx]
    return rectify


#@cython.boundscheck(False)
def rectify_undo(np.ndarray[DTYPE_t, ndim=2] recframe not None,
                np.ndarray[ITYPE_t, ndim=1] pixcen not None,
                np.ndarray[ITYPE_t, ndim=1] pixledg not None,
                np.ndarray[ITYPE_t, ndim=1] pixredg not None,
                int pixwid, double maskval, int sz_x, int sz_y):
    """
    To be used in conjunction with the rectify function directly above.
    sz_x and sz_y correspond to the frame.shape[0] and frame.shape[1]
    respectively.
    """
    cdef int x, y, ymin, ymax
    cdef int shift = pixwid/2

    cdef np.ndarray[DTYPE_t, ndim=2] unrec = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(0,sz_x):
        if pixcen[x] <= 0:
            ymin = (pixredg[x]-shift)-shift
            ymax = (pixredg[x]-shift)+shift+1
        elif pixcen[x] >= sz_y-1:
            ymin = (pixledg[x]+shift)-shift
            ymax = (pixledg[x]+shift)+shift+1
        else:
            ymin = pixcen[x]-shift
            ymax = pixcen[x]+shift+1
        if ymin < 0: ymin = 0
        elif pixledg[x] >= sz_y-1: continue
        if pixredg[x] <= 0: continue
        elif ymax > sz_y: ymax = sz_y
        for y in range(ymin,ymax):
            unrec[x,y] = recframe[x,y-ymin]
    return unrec


##################################
#  Sutherland-Hodgman Algorithm  #
##################################

#@cython.boundscheck(False)
def SH_poly_clip_area(np.ndarray[DTYPE_t, ndim=2] poly not None,
              np.ndarray[DTYPE_t, ndim=2] pixl not None,
              np.ndarray[DTYPE_t, ndim=2] p1 not None,
              np.ndarray[DTYPE_t, ndim=2] p2 not None):

    cdef int i, j, ii
    cdef int sz_p, sz_e, sz_a, dir
    cdef double sumv, intv, dd
    cdef double tx, ty, x

    cdef int ls, sz_r, side0, side1
    cdef double v1x, v1y, v0x, v0y

    sz_p = pixl.shape[0]
    sz_e = poly.shape[0]
    sz_a = p1.shape[0]

    for i in range(0,sz_a):
        if p1[i,0] == -999999.9 and p2[i,0] == -999999.9: break
        p1[i,0] = -999999.9
        p1[i,1] = -999999.9
        p2[i,0] = -999999.9
        p2[i,1] = -999999.9

    x = (pixl[1,0]-pixl[0,0])*(pixl[2,1]-pixl[1,1])-(pixl[1,1]-pixl[0,1])*(pixl[2,0]-pixl[1,0])
    if x < 0.0: dir = -1
    elif x > 0.0: dir = 1
    else: dir = 0
    v0x = poly[sz_e-1,0]
    v0y = poly[sz_e-1,1]
    x = (pixl[0,0]-pixl[sz_p-1,0])*(v0y-pixl[0,1])-(pixl[0,1]-pixl[sz_p-1,1])*(v0x-pixl[0,0])
    if x < 0.0: side0 = -1
    elif x > 0.0: side0 = 1
    else: side0 = 0
    sz_r = 0
    if side0 != -dir:
        p2[sz_r,0]=v0x
        p2[sz_r,1]=v0y
        sz_r += 1
    for ii in range(0,sz_e):
        v1x = poly[ii,0]
        v1y = poly[ii,1]
        x = (pixl[0,0]-pixl[sz_p-1,0])*(v1y-pixl[0,1])-(pixl[0,1]-pixl[sz_p-1,1])*(v1x-pixl[0,0])
        if x < 0.0: side1 = -1
        elif x > 0.0: side1 = 1
        else: side1 = 0
        if side0+side1 == 0 and side0:
            dd = (v1x-v0x)*(pixl[0,1]-pixl[sz_p-1,1])-(v1y-v0y)*(pixl[0,0]-pixl[sz_p-1,0])
            if not dd:
                ls = 0
                tx = 0.0
                ty = 0.0
            else:
                dd = ((pixl[sz_p-1,0]-v0x)*(pixl[0,1]-pixl[sz_p-1,1])-(pixl[sz_p-1,1]-v0y)*(pixl[0,0]-pixl[sz_p-1,0]))/dd
                if dd <= 0.0 or dd >= 1.0:
                    ls = 0
                    tx = 0.0
                    ty = 0.0
                else:
                    ls = 1
                    tx = v0x + dd*(v1x-v0x)
                    ty = v0y + dd*(v1y-v0y)
            if ls == 1:
                p2[sz_r,0] = tx
                p2[sz_r,1] = ty
                sz_r += 1
        if ii == sz_e-1:
            for j in range(sz_r,sz_a):
                if p2[j,0]==-999999.9: break
                p2[j,0]=-999999.9
                p2[j,1]=-999999.9
            break
        if side1 != -dir:
            p2[sz_r,0] = v1x
            p2[sz_r,1] = v1y
            sz_r += 1
        v0x = v1x
        v0y = v1y
        side0 = side1
    for i in range(0, sz_p-1):
        sz_e=-1
        for j in range(0,sz_a):
            if p2[j,0]==-999999.9 and p2[j,1]==-999999.9 and sz_e==-1:
                sz_e = j
            tx = p2[j,0]
            ty = p2[j,1]
            p2[j,0] = p1[j,0]
            p2[j,1] = p1[j,1]
            p1[j,0] = tx
            p1[j,1] = ty
        v0x = p1[sz_e-1,0]
        v0y = p1[sz_e-1,1]
        x = (pixl[i+1,0]-pixl[i,0])*(v0y-pixl[i+1,1])-(pixl[i+1,1]-pixl[i,1])*(v0x-pixl[i+1,0])
        if x < 0.0: side0 = -1
        elif x > 0.0: side0 = 1
        else: side0 = 0
        sz_r = 0
        if side0 != -dir:
            p2[sz_r,0]=v0x
            p2[sz_r,1]=v0y
            sz_r += 1
        for ii in range(0,sz_e):
            v1x = p1[ii,0]
            v1y = p1[ii,1]
            x = (pixl[i+1,0]-pixl[i,0])*(v1y-pixl[i+1,1])-(pixl[i+1,1]-pixl[i,1])*(v1x-pixl[i+1,0])
            if x < 0.0: side1 = -1
            elif x > 0.0: side1 = 1
            else: side1 = 0
            if side0+side1 == 0 and side0:
                dd = (v1x-v0x)*(pixl[i+1,1]-pixl[i,1])-(v1y-v0y)*(pixl[i+1,0]-pixl[i,0])
                if not dd:
                    ls = 0
                    tx = 0.0
                    ty = 0.0
                else:
                    dd = ((pixl[i,0]-v0x)*(pixl[i+1,1]-pixl[i,1])-(pixl[i,1]-v0y)*(pixl[i+1,0]-pixl[i,0]))/dd
                    if dd <= 0.0 or dd >= 1.0:
                        ls = 0
                        tx = 0.0
                        ty = 0.0
                    else:
                        ls = 1
                        tx = v0x + dd*(v1x-v0x)
                        ty = v0y + dd*(v1y-v0y)
                if ls == 1:
                    p2[sz_r,0] = tx
                    p2[sz_r,1] = ty
                    sz_r += 1
            if ii == sz_e-1:
                for j in range(sz_r,sz_a):
                    if p2[j,0]==-999999.9: break
                    p2[j,0]=-999999.9
                    p2[j,1]=-999999.9
                break
            if side1 != -dir:
                p2[sz_r,0] = v1x
                p2[sz_r,1] = v1y
                sz_r += 1
            v0x = v1x
            v0y = v1y
            side0 = side1

    # Calculate the size of the array
    for i in range(sz_a):
        if p2[i,0] == -999999.9:
            sz_a = i
            break
    # Calculate the area of the polygon
    dd = p2[sz_a-1,0]*p2[0,1]
    for i in range(0,sz_a-1):
        dd += p2[i,0]*p2[i+1,1]
    for i in range(0,sz_a-1):
        dd -= p2[i+1,0]*p2[i,1]
    dd -= p2[0,0]*p2[sz_a-1,1]
    if dd < 0.0: dd *= -0.5
    else: dd *= 0.5
    return dd
