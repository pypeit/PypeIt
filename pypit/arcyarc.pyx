# To get this running, you must do the following at the command line:
# python arcytrace_setup.py build_ext --inplace
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
    double clog "log" (double)
    double cpow "pow" (double, double)

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
def detections_sigma(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                    np.ndarray[DTYPE_t, ndim=1] errors not None,
                    np.ndarray[ITYPE_t, ndim=1] satmask not None,
                    double thresh, double detect):
    cdef int sz_p, flag
    cdef int p, c, d, npc, inum
    cdef double psum, esum, mnum

    sz_p = pixels.shape[0]
    cdef np.ndarray[ITYPE_t, ndim=1] pixcen = np.zeros(sz_p, dtype=ITYPE)

    c = 0
    npc = 0
    p = 1
    pixcen[0] = -1
    while p < sz_p-1:
        pixcen[p] = -1
        # If you have a saturated pixel, ignore it
        if satmask[p] == 1:
            c = 0
            while satmask[p+c] == 1:
                pixcen[p+c] = -1
                c += 1
                if p+c >= sz_p:
                    break
            p += c
            if p >= sz_p: break
            else: continue
        if pixels[p]/errors[p] >= thresh and pixels[p+1] > pixels[p]:
            psum = pixels[p]
            esum = errors[p]*errors[p]
            mnum = pixels[p]
            inum = p
            c=1
            while pixels[p-c] < pixels[p-c+1]:
                psum += pixels[p-c]
                esum += errors[p-c]*errors[p-c]
                c += 1
                if p-c < 0: break
                if pixels[p-c] <= 0.0: break
            d=1
            flag = 0
            while True:
                # Test if the peak of the emission line has been reached, and we are now decreasing in flux.
                if p+d >= (sz_p-1):
                    break
                if pixels[p+d] > pixels[p+d-1] and pixels[p+d] > pixels[p+d+1]:
                    flag = 1
                pixcen[p+d] = -1
                psum += pixels[p+d]
                esum += errors[p+d]*errors[p+d]
                if pixels[p+d]>mnum:
                    mnum = pixels[p+d]
                    inum = p+d
                d += 1
                if p+d >= sz_p:
                    break
                if pixels[p+d] > pixels[p+d-1] and flag == 1: break
                if pixels[p+d] <= 0.0: break
            p += d
            #print p, npc, psum, csqrt(esum), c, d, psum/csqrt(esum)
            if psum/csqrt(esum) >= detect and d+c >= 6:
                pixcen[npc] = inum
                npc += 1
                psum = 0.0
            if p >= sz_p: break
        else:
            p += 1
            if p >= sz_p: break
    return pixcen, npc


#@cython.boundscheck(False)
def fit_arcorder(np.ndarray[DTYPE_t, ndim=1] xarray not None,
                np.ndarray[DTYPE_t, ndim=1] yarray not None,
                np.ndarray[ITYPE_t, ndim=1] pixt not None,
                int fitp):
    cdef int p, pp, i, sz_p, sz_a
    cdef int pmin, pmax, flag
    cdef double tamp, tcen, twid

    sz_p = pixt.shape[0]
    sz_a = yarray.shape[0]
    pp = 0

    cdef np.ndarray[DTYPE_t, ndim=1] coeff = np.zeros(3, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] ampl = np.zeros(sz_a, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] cent = np.zeros(sz_a, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] widt = np.zeros(sz_a, dtype=DTYPE)

    for p in range(sz_p):
        pmin = pixt[p]-(fitp-1)/2
        pmax = pixt[p]-(fitp-1)/2 + fitp
        if pmin < 0: pmin=0
        if pmax > sz_a: pmax = sz_a
        if pmin == pmax: continue
        # Check that all pixel values are positive and reset pmin and pmax as required
        for i in range(pmin,pixt[p]):
            if yarray[pixt[p]-1-i+pmin] < 0.0:
                pmin = 1+pixt[p]-1-i+pmin
                break
        for i in range(pixt[p],pmax):
            if yarray[i] < 0.0:
                pmax = i-1
                break
        if pixt[p]-pmin <= 1 or pmax-pixt[p]<=1: continue # Probably won't be a good solution
        # Fit the gaussian
        tamp, tcen, twid, flag = fit_gauss(xarray, yarray, coeff, pmin, pmax, <double>(pixt[p]))
        if flag == 1: # The fit is acceptable
            ampl[pp] = tamp
            cent[pp] = tcen
            widt[pp] = twid
            pp += 1
    for p in range(pp,sz_p):
        ampl[p] = -1.0
        cent[p] = -1.0
        widt[p] = -1.0
    return ampl, cent, widt, pp


#@cython.boundscheck(False)
def fit_gauss(np.ndarray[DTYPE_t, ndim=1] xarray not None,
                np.ndarray[DTYPE_t, ndim=1] yarray not None,
                np.ndarray[DTYPE_t, ndim=1] coeff not None,
                int pmin, int pmax, double cent):
    cdef int i, sz_x
    cdef double da, db, dc
    cdef double l, r, dx
    cdef double Maa, Mab, Mac, Mba, Mbb, Mbc, Mca, Mcb, Mcc
    cdef double amp, cen, wid
    cdef double chisq

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat
    cdef gsl_vector *yvec
    cdef gsl_vector *c

    sz_x = xarray.shape[0]
    da = 0.0
    db = 0.0
    dc = 0.0
    Maa = 0.0
    Mab = 0.0
    Mac = 0.0
    Mba = 0.0
    Mbb = 0.0
    Mbc = 0.0
    Mca = 0.0
    Mcb = 0.0
    Mcc = 0.0

    for i in range(pmin,pmax):
        if i == 0:
            l = 0.5*(3.0*(xarray[i]-cent)-(xarray[i+1]-cent))
            r = 0.5*((xarray[i+1]-cent)+(xarray[i]-cent))
        elif i == sz_x-1:
            l = 0.5*((xarray[i]-cent)+(xarray[i-1]-cent))
            r = 0.5*(3.0*(xarray[i]-cent)-(xarray[i-1]-cent))
        else:
            l = 0.5*((xarray[i]-cent)+(xarray[i-1]-cent))
            r = 0.5*((xarray[i+1]-cent)+(xarray[i]-cent))
        dx = r-l
        da += clog(yarray[i])*(r-l)*dx
        db += clog(yarray[i])*(r**2-l**2)*dx
        dc += clog(yarray[i])*(r**3-l**3)*dx
        Maa += (r-l) * (r**3-l**3)
        Mab += (r-l) * (r**2-l**2)
        Mac += (r-l) * (r-l)
        Mba += (r**2-l**2) * (r**3-l**3)
        Mbb += (r**2-l**2) * (r**2-l**2)
        Mbc += (r**2-l**2) * (r-l)
        Mca += (r**3-l**3) * (r**3-l**3)
        Mcb += (r**3-l**3) * (r**2-l**2)
        Mcc += (r**3-l**3) * (r-l)

    Maa /= 3.0
    Mba /= 3.0
    Mca /= 3.0
    Mab /= 2.0
    Mbb /= 2.0
    Mcb /= 2.0

    xmat = gsl_matrix_alloc(3, 3)
    yvec = gsl_vector_alloc(3)
    c    = gsl_vector_alloc(3)
    cov  = gsl_matrix_alloc(3, 3)

    gsl_matrix_set(xmat, 0, 0, Maa)
    gsl_matrix_set(xmat, 0, 1, Mab)
    gsl_matrix_set(xmat, 0, 2, Mac)
    gsl_matrix_set(xmat, 1, 0, Mba)
    gsl_matrix_set(xmat, 1, 1, Mbb)
    gsl_matrix_set(xmat, 1, 2, Mbc)
    gsl_matrix_set(xmat, 2, 0, Mca)
    gsl_matrix_set(xmat, 2, 1, Mcb)
    gsl_matrix_set(xmat, 2, 2, Mcc)
    gsl_vector_set(yvec, 0, da)
    gsl_vector_set(yvec, 1, db)
    gsl_vector_set(yvec, 2, dc)

    # Fit with a polynomial
    ws = gsl_multifit_linear_alloc(3, 3)
    gsl_multifit_linear(xmat, yvec, c, cov, &chisq, ws)

    # Store the result in the coeffs array
    for i in range(0,3):
        coeff[i] = gsl_vector_get(c, i)

    # Free some of these arrays from memory
    gsl_multifit_linear_free(ws)
    gsl_matrix_free(xmat)
    gsl_matrix_free(cov)
    gsl_vector_free(yvec)
    gsl_vector_free(c)

    if coeff[0] >= 0.0:
        return 0.0, 0.0, 0.0, 0
    else:
        wid = csqrt(-0.5/coeff[0])
        cen = cent-0.5*coeff[1]/coeff[0]
        amp = cexp(coeff[2] - 0.25*coeff[1]**2/coeff[0])
        return amp, cen, wid, 1

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
def remove_similar(np.ndarray[ITYPE_t, ndim=1] array not None,
                    int numl):

    cdef int i, ii, sz_i, skipone

    ii=0
    sz_i = array.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=1] outarr = np.zeros(sz_i, dtype=ITYPE)

    skipone = 0
    for i in range(sz_i-1):
        if skipone == 1:
            skipone = 0
            continue
        if array[i+1]-array[i] > 5:
            outarr[ii] = array[i]
            ii += 1
        else:
            outarr[ii] = array[i]
            ii += 1
            skipone=1
    if skipone == 0:
        outarr[ii] = array[sz_i-1]
        ii += 1
    for i in range(ii,sz_i):
        outarr[i] = -1
    return outarr


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
