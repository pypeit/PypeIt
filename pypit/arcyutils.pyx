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

#######
#  A  #
#######



#######
#  B  #
#######

#@cython.boundscheck(False)
def bin_x(np.ndarray[DTYPE_t, ndim=2] array not None,
            int fact, int tcomb):
    # tcomb = 0: median
    # tcomb = 1: mean
    cdef int sz_x, sz_y, sz_b
    cdef int x, y, b

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    sz_b = sz_x/fact

    cdef np.ndarray[DTYPE_t, ndim=2] binarr = np.zeros((sz_b,sz_y), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] medarr = np.zeros((fact), dtype=DTYPE)

    for y in range(sz_y):
        for b in range(sz_b):
            for x in range(fact):
                medarr[x] = array[b*fact+x,y]
            if tcomb == 0:
                binarr[b,y] = median(medarr)
            elif tcomb == 1:
                binarr[b,y] = mean(medarr)
    return binarr


#@cython.boundscheck(False)
def bin_x_weight(np.ndarray[DTYPE_t, ndim=2] array not None,
                np.ndarray[DTYPE_t, ndim=2] weight not None,
                int fact, double maskval):
    cdef int sz_x, sz_y, sz_b
    cdef int x, y, b

    sz_x = array.shape[0]
    sz_y = array.shape[1]
    sz_b = sz_x/fact

    cdef np.ndarray[DTYPE_t, ndim=2] binarr = np.zeros((sz_b,sz_y), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] medarr = np.zeros((fact), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] wgtarr = np.zeros((fact), dtype=DTYPE)

    for y in range(sz_y):
        for b in range(sz_b):
            for x in range(fact):
                medarr[x] = array[b*fact+x,y]
                wgtarr[x] = weight[b*fact+x,y]
            binarr[b,y] = mean_weight(medarr, wgtarr, maskval)
    return binarr


#@cython.boundscheck(False)
def bin_y(np.ndarray[DTYPE_t, ndim=2] array not None,
            int fact, int tcomb):
    # tcomb = 0: median
    # tcomb = 1: mean
    cdef int sz_x, sz_y, sz_b
    cdef int x, y, b

    sz_x = array.shape[1]
    sz_y = array.shape[0]
    sz_b = sz_x/fact

    cdef np.ndarray[DTYPE_t, ndim=2] binarr = np.zeros((sz_y,sz_b), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] medarr = np.zeros((fact), dtype=DTYPE)

    for y in range(sz_y):
        for b in range(sz_b):
            for x in range(fact):
                medarr[x] = array[y,b*fact+x]
            if tcomb == 0:
                binarr[y,b] = median(medarr)
            elif tcomb == 1:
                binarr[y,b] = mean(medarr)
    return binarr


#@cython.boundscheck(False)
def bspline_fit(np.ndarray[DTYPE_t, ndim=1] marray not None,
                np.ndarray[DTYPE_t, ndim=1] xarray not None,
                np.ndarray[DTYPE_t, ndim=1] yarray not None,
                np.ndarray[DTYPE_t, ndim=1] earray not None,
                double minval, double maxval, int ncoeffs, int k):
   
    cdef int nbreak = ncoeffs + 2 - k
    cdef int sz_x, sz_m, i, j
    cdef double xi, yi, yerr, Bj
    cdef double dy, chisq, Rsq, dof, tss

    sz_x = xarray.shape[0]
    sz_m = marray.shape[0]

    cdef gsl_bspline_workspace *bw
    cdef gsl_vector *B
    cdef gsl_rng *r
    cdef gsl_vector *c
    cdef gsl_vector *w
    cdef gsl_vector *x
    cdef gsl_vector *y
    cdef gsl_matrix *xm
    cdef gsl_matrix *cov
    cdef gsl_multifit_linear_workspace *mw

    # Create the model array
    cdef np.ndarray[DTYPE_t, ndim=1] model = np.zeros(sz_m, dtype=DTYPE)

    gsl_rng_env_setup()
    r = gsl_rng_alloc(gsl_rng_default)

    # allocate a cubic bspline workspace
    bw = gsl_bspline_alloc(k, nbreak)
    B = gsl_vector_alloc(ncoeffs)

    x = gsl_vector_alloc(sz_x)
    y = gsl_vector_alloc(sz_x)
    w = gsl_vector_alloc(sz_x)
    xm = gsl_matrix_alloc(sz_x, ncoeffs)
    c = gsl_vector_alloc(ncoeffs)
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs)
    mw = gsl_multifit_linear_alloc(sz_x, ncoeffs)

    for i in range(0,sz_x):
        gsl_vector_set(x, i, xarray[i])
        gsl_vector_set(y, i, yarray[i])
        gsl_vector_set(w, i, 1.0 / (earray[i] * earray[i]))

    # Use uniform break points on the interval
    gsl_bspline_knots_uniform(minval, maxval, bw)

    for i in range(0,sz_x):
        xi = gsl_vector_get(x, i)
        # Compute B_j(xi) for all j
        # May want to subpixellate the bspline and sum (to mimic integration)
        gsl_bspline_eval(xi, B, bw);
        for j in range(0, ncoeffs):
            Bj = gsl_vector_get(B, j)
            gsl_matrix_set(xm, i, j, Bj)

    # Perform the fit
    gsl_multifit_wlinear(xm, w, y, c, cov, &chisq, mw)

    dof = sz_x - ncoeffs
    #tss = gsl_stats_wtss((*w).data, 1, (*y).data, 1, (*y).size)
    tss = 1.0
    Rsq = 1.0 - chisq / tss

    # Evaluate the b-spline
    for i in range(sz_m):
        gsl_bspline_eval(marray[i], B, bw)
        gsl_multifit_linear_est(B, c, cov, &yi, &yerr)
        model[i] = yi

    gsl_rng_free(r)
    gsl_bspline_free(bw)
    gsl_vector_free(B)
    gsl_vector_free(x)
    gsl_vector_free(y)
    gsl_matrix_free(xm)
    gsl_vector_free(c)
    gsl_vector_free(w)
    gsl_matrix_free(cov)
    gsl_multifit_linear_free(mw)

    # Return the best-fitting model
    return model


#@cython.boundscheck(False)
def bspline_fitmod(np.ndarray[DTYPE_t, ndim=1] xarray not None,
                    np.ndarray[DTYPE_t, ndim=1] yarray not None,
                    np.ndarray[DTYPE_t, ndim=1] earray not None,
                    double minval, double maxval, int ncoeffs, int k,
                    np.ndarray[DTYPE_t, ndim=3] pixmap not None,
                    np.ndarray[DTYPE_t, ndim=1] tilts not None,
                    np.ndarray[DTYPE_t, ndim=1] dtrc not None,
                    np.ndarray[DTYPE_t, ndim=1] scen not None,
                    np.ndarray[ITYPE_t, ndim=1] dpix not None,
                    np.ndarray[ITYPE_t, ndim=1] spix not None):
   
    cdef int nbreak = ncoeffs + 2 - k
    cdef int sz_x, sz_m, i, j
    cdef double xi, yi, yerr, Bj
    cdef double dy, chisq, Rsq, dof, tss
    cdef int p, b, donebl, donetl
    cdef int sz_p, sz_b
    cdef double fa, fb, la, lb, vb # Note, va=0 since the origin is defined as the bottom-left corner of the pixel
    cdef double dtl, stl, dbl, sbl
    cdef double cv, sva, dva, svb, dvb

    sz_x = xarray.shape[0]
    sz_m = pixmap.shape[0]
    sz_n = pixmap.shape[1]
    sz_b = dtrc.shape[0]
    sz_p = dpix.shape[0]

    cdef gsl_bspline_workspace *bw
    cdef gsl_vector *B
    cdef gsl_rng *r
    cdef gsl_vector *c
    cdef gsl_vector *w
    cdef gsl_vector *x
    cdef gsl_vector *y
    cdef gsl_matrix *xm
    cdef gsl_matrix *cov
    cdef gsl_multifit_linear_workspace *mw

    # Create the model array
    cdef np.ndarray[DTYPE_t, ndim=2] model = np.zeros((sz_m,sz_n), dtype=DTYPE)

    gsl_rng_env_setup()
    r = gsl_rng_alloc(gsl_rng_default)

    # allocate a cubic bspline workspace
    bw = gsl_bspline_alloc(k, nbreak)
    B = gsl_vector_alloc(ncoeffs)

    x = gsl_vector_alloc(sz_x)
    y = gsl_vector_alloc(sz_x)
    w = gsl_vector_alloc(sz_x)
    xm = gsl_matrix_alloc(sz_x, ncoeffs)
    c = gsl_vector_alloc(ncoeffs)
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs)
    mw = gsl_multifit_linear_alloc(sz_x, ncoeffs)

    for i in range(0,sz_x):
        gsl_vector_set(x, i, xarray[i])
        gsl_vector_set(y, i, yarray[i])
        gsl_vector_set(w, i, 1.0 / (earray[i] * earray[i]))

    # Use uniform break points on the interval
    gsl_bspline_knots_uniform(minval, maxval, bw);

    print "DONE!"
    for i in range(0,sz_x):
        xi = gsl_vector_get(x, i)
        # Compute B_j(xi) for all j
        # May want to subpixellate the bspline and sum (to mimic integration)
        gsl_bspline_eval(xi, B, bw);
        for j in range(0, ncoeffs):
            Bj = gsl_vector_get(B, j)
            gsl_matrix_set(xm, i, j, Bj)
    print "SLOW!"

    # Perform the fit
    gsl_multifit_wlinear(xm, w, y, c, cov, &chisq, mw);
    #gsl_multifit_linear(xm, y, c, cov, &chisq, mw);
    print "FIT!"

    dof = sz_x - ncoeffs
    #tss = gsl_stats_wtss((*w).data, 1, (*y).data, 1, (*y).size)
    tss = 1.0
    Rsq = 1.0 - chisq / tss

    # Construct the model image
    for p in range(sz_p):
        # get the bottom-left (i.e. the origin) and top-left corner of the pixel
        dbl = pixmap[dpix[p],spix[p],0] - 0.5*pixmap[dpix[p],spix[p],2]
        sbl = pixmap[dpix[p],spix[p],1] - 0.5*pixmap[dpix[p],spix[p],3]
        dtl = pixmap[dpix[p],spix[p],0] + 0.5*pixmap[dpix[p],spix[p],2]
        stl = pixmap[dpix[p],spix[p],1] + 0.5*pixmap[dpix[p],spix[p],3]
        donebl = 0
        donetl = 0
        for b in range(1,sz_b):
            # Bottom-left corner of pixel
            cv = dtrc[b] - tilts[b]*scen[b]
            sva = (sbl + tilts[b]*(dbl-cv))/(1.0 + tilts[b]*tilts[b])
            dva = tilts[b]*sva + cv
            if (dva > dbl) and (donebl == 0):
                # We have found the location above the pixel corner, now linearly interpolate to the pixel corner
                cv = dtrc[b-1] - tilts[b-1]*scen[b-1]
                svb = (sbl + tilts[b-1]*(dbl-cv))/(1.0 + tilts[b-1]*tilts[b-1])
                dvb = tilts[b-1]*svb + cv
                la = dtrc[b] - (dtrc[b]-dtrc[b-1])*(dva-dbl)/(dva-dvb)
                donebl = 1
            # Top-left corner of pixel
            cv = dtrc[b] - tilts[b]*scen[b]
            sva = (stl + tilts[b]*(dtl-cv))/(1.0 + tilts[b]*tilts[b])
            dva = tilts[b]*sva + cv
            if (dva > dtl) and (donetl == 0):
                cv = dtrc[b-1] - tilts[b-1]*scen[b-1]
                svb = (stl + tilts[b-1]*(dtl-cv))/(1.0 + tilts[b-1]*tilts[b-1])
                dvb = tilts[b-1]*svb + cv
                lb = dtrc[b] - (dtrc[b]-dtrc[b-1])*(dva-dtl)/(dva-dvb)
                vb = (dtl-dbl)/csqrt(1.0 + tilts[b]*tilts[b-1]) # This should really by tilts^2 at the tilt value corresponding to lb, but tilts[b]*tilts[b-1] should be close enough
                donetl = 1
            if (donebl==1) and (donetl==1):
                # Get the tilts value
                sva = tilts[b]
                break
        if la==lb:
            # Pixel cannot be given the same wavelength --> Fails at the ends of the chip
            continue
        # Use the la, lb and vb to describe how dtrc varies over the pixel
        # Convert la, lb to fa, fb (the model fluxes)
        gsl_bspline_eval(la, B, bw)
        gsl_multifit_linear_est(B, c, cov, &fa, &yerr)
        gsl_bspline_eval(lb, B, bw)
        gsl_multifit_linear_est(B, c, cov, &fb, &yerr)
        # Assuming fa = a*la + b, then
        # a=(fa-fb)/(la-lb)  (assign this to the variable dva) and
        # b=fa-a*la  (assign this to the variable dvb)
        dva = (fa-fb)/(la-lb)
        dvb = fa - dva*la
        # and then replace the following m and c by m->dva*m  and c->dva*c + dvb
        # Assuming dtrc = m*v + c, then c=la (i.e. v=0 for la), and m = (lb-la)/vb (assign this to be the variable cv)
        cv = (lb-la)/vb
        svb = csqrt(sva*sva + 1.0) / (sva*sva - 1.0)
        # Finally, calculate the appropriate value for xfit and yfit
        model[dpix[p],spix[p]] = (la*dva + dvb) + 0.5*(cv*dva)*svb*(pixmap[dpix[p],spix[p],2] + sva*pixmap[dpix[p],spix[p],3])
# 	for i in range(sz_m):
# 		gsl_bspline_eval(marray[i], B, bw)
# 		gsl_multifit_linear_est(B, c, cov, &yi, &yerr)
# 		model[i] = yi

    gsl_rng_free(r)
    gsl_bspline_free(bw)
    gsl_vector_free(B)
    gsl_vector_free(x)
    gsl_vector_free(y)
    gsl_matrix_free(xm)
    gsl_vector_free(c)
    gsl_vector_free(w)
    gsl_matrix_free(cov)
    gsl_multifit_linear_free(mw)

    # Return the best-fitting model
    return model


#######
#  C  #
#######


#@cython.boundscheck(False)
def checkmatch(np.ndarray[DTYPE_t, ndim=2] arrayA not None,
            np.ndarray[DTYPE_t, ndim=2] arrayB not None,
            double maskvalue):
    cdef int sz_x, sz_y
    cdef int x, y
    cdef double midv, ra, ea

    sz_x = arrayA.shape[0]
    sz_y = arrayA.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] ratarr = np.zeros((sz_x*sz_y), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] errarr = np.zeros((sz_x*sz_y), dtype=DTYPE)

    cdef int n = 0
    for x in range(sz_x):
        for y in range(sz_y):
            if arrayA[x,y] > 0.0 and arrayB[x,y] > 0.0:
                ratarr[x*sz_y+y] = arrayA[x,y]/arrayB[x,y]
                errarr[x*sz_y+y] = ratarr[x*sz_y+y] * csqrt(1.0/arrayA[x,y] + 1.0/arrayB[x,y])
                n += 1
            else:
                ratarr[x*sz_y+y] = maskvalue
                errarr[x*sz_y+y] = maskvalue

    # Do an efficient median calculation
    mediansort(ratarr,errarr,n/2)
    if n%2==0:
        midv = 0.5*(ratarr[n/2] + ratarr[n/2 - 1])
    else:
        midv = ratarr[(n-1)/2]
    # Calculate the chi-squared value
    cdef double chisq = 0.0
    for x in range(n):
        if ratarr[x] != maskvalue:
            chisq += ((ratarr[x]-midv)/errarr[x])**2
    return chisq

#@cython.boundscheck(False)
def crreject(np.ndarray[DTYPE_t, ndim=2] frame not None):
    """
    Identifies cosmic rays as sudden jumps in neighbouring pixel fluxes
    """
    cdef int sz_x, sz_y, x, y
    cdef double maxdiff, tst

    sz_x = frame.shape[0]
    sz_y = frame.shape[1]

    # Create the model array
    cdef np.ndarray[DTYPE_t, ndim=2] crmask = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for x in range(0,sz_x):
        for y in range(0, sz_y):
            if y==0:
                maxdiff = frame[x,y+1]-frame[x,y]
                if maxdiff<0.0: maxdiff *= -1.0
            elif y==sz_y-1:
                maxdiff = frame[x,y]-frame[x,y-1]
                if maxdiff<0.0: maxdiff *= -1.0
            else:
                maxdiff = frame[x,y+1]-frame[x,y]
                tst = frame[x,y]-frame[x,y-1]
                if maxdiff<0.0: maxdiff *= -1.0
                if tst<0.0: tst *= -1.0
                if tst>maxdiff: maxdiff = tst
            crmask[x,y] = maxdiff*frame[x,y]
    # Return the mask of cosmic rays
    return crmask


#######
#  F  #
#######

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


#######
#  G  #
#######

#@cython.boundscheck(False)
def get_splknots(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                    np.ndarray[ITYPE_t, ndim=1] knotarr not None,
                    double minv, double maxv, int num):

    cdef int x, i, sz_x, cnt
    cdef double diff, sdiff
    cdef np.ndarray[DTYPE_t, ndim=1] knots = np.zeros((num), dtype=DTYPE)

    sz_x = xarr.shape[0]
    diff = (maxv-minv)/<double>(sz_x-1)
    cnt = 0
    for x in range(0,sz_x):
        sdiff = diff/<double>(2*knotarr[x])
        for i in range(0,knotarr[x]):
            if x==0: knots[cnt+i] = minv + diff*<double>(2*i)/<double>(2*knotarr[x]+1)
            elif x==sz_x-1: knots[cnt+i] = maxv - 2.0*diff*<double>(knotarr[x]-i-1)/<double>(2*knotarr[x]+1)
            else: knots[cnt+i] = minv + diff*<double>(x) + sdiff*<double>(2*i + 1)
        cnt += knotarr[x]
    return knots


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


#######
#  H  #
#######


def hough(np.ndarray[DTYPE_t, ndim=2] frame not None,
            np.ndarray[DTYPE_t, ndim=2] mask not None,
            int ntx, int mry):
    # Calculate the Hough transform
    cdef int nimx, mimy, jx, iy, jtx, iry
    cdef double rmax, r, dr, th, dth

    nimx = frame.shape[0]
    mimy = frame.shape[1]

    mry = (mry/2)*2   # Make sure that this is even

    cdef np.ndarray[DTYPE_t, ndim=2] hought = np.zeros((ntx, mry), dtype=DTYPE)

    rmax = csqrt(<double>((nimx-1)*(nimx-1) + (mimy-1)*(mimy-1)))
    dr = rmax / <double>(mry/2)
    dth = 3.14159265358979 / <double>(ntx)

    for jx in xrange(nimx):
        for iy in xrange(mimy):
            if (mask[jx,iy]==1.0): continue
            for jtx in xrange(ntx):
                th = dth * <double>(jtx)
                r = <double>(jx)*ccos(th) + <double>(iy)*csin(th)
                if r < 0.0: iry = mry/2 + <int>(r/dr-0.5)
                else: iry = mry/2 + <int>(r/dr+0.5)
                hought[jtx, iry-1] += 1.0
    return hought


#######
#  M  #
#######

#@cython.boundscheck(False)
def mean(np.ndarray[DTYPE_t, ndim=1] array not None):
    cdef int sz_x
    cdef int x
    cdef double temp = 0.0

    sz_x = array.shape[0]

    for x in range(sz_x):
        temp += array[x]
    return temp/<double>(sz_x)


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
def median(np.ndarray[DTYPE_t, ndim=1] array not None):
    cdef int sz_x
    cdef int x, y
    cdef double temp

    sz_x = array.shape[0]

    for x in range(sz_x-1):
        for y in range(x,sz_x):
            # Sort the array
            if array[y] < array[x]:
                temp = array[x]
                array[x] = array[y]
                array[y] = temp
    # Find the median value
    if sz_x%2 == 0:
        return 0.5*(array[sz_x/2] + array[sz_x/2 - 1])
    else:
        return array[(sz_x-1)/2]


#@cython.boundscheck(False)
def mediansort(np.ndarray[DTYPE_t, ndim=1] arrayA not None,
                np.ndarray[DTYPE_t, ndim=1] arrayB not None,
                int k):
   
    cdef int i, j, l, m
    cdef DTYPE_t x, tmpA, tmpB
    l = 0
    m = <int> arrayA.shape[0] - 1
    with nogil:
        while l < m:
            x = arrayA[k]
            i = l
            j = m
            while 1:
                while arrayA[i] < x: i += 1
                while x < arrayA[j]: j -= 1
                if i <= j:
                    tmpA = arrayA[i]
                    tmpB = arrayB[i]
                    arrayA[i] = arrayA[j]
                    arrayB[i] = arrayB[j]
                    arrayA[j] = tmpA
                    arrayB[j] = tmpB
                    i += 1
                    j -= 1
                if i > j: break
            if j < k: l = i
            if k < i: m = j
    return

#######
#  N  #
#######

#######
#  O  #
#######

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


#######
#  P  #
#######

#@cython.boundscheck(False)
def polydiff(np.ndarray[DTYPE_t, ndim=1] x not None,
        np.ndarray[DTYPE_t, ndim=1] y not None,
        np.ndarray[DTYPE_t, ndim=1] c not None):

    cdef int sz_x, sz_c, i, j
    cdef double diff
    sz_x = x.shape[0]
    sz_c = c.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] md = np.zeros(sz_x, dtype=DTYPE)

    for i in range(sz_x):
        diff = 0.0
        for j in range(sz_c):
            diff += c[j] * cpow(x[i], j)
        md[i] = diff-y[i]
        # Now take the absolute value
        if md[i] < 0.0: md[i] *= -1.0
    return md


#@cython.boundscheck(False)
def poly_filter(np.ndarray[DTYPE_t, ndim=1] yarray not None,
                int order, int npix, double maskval):
    """
    Performs a polynomial filtering.
    """
    cdef int sz_x, x, i, nmask
    cdef double value, chisq

    sz_x = yarray.shape[0]

    if order + 1 >= 2*npix+1:
        print "WARNING: Polynomial order exceeds number of pixels"
        return yarray

    cdef np.ndarray[DTYPE_t, ndim=1] yfilt = np.zeros(sz_x, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coeff = np.zeros(order+1, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] pixfitx = np.arange(-npix,npix+1, dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] pixfitm = np.zeros(2*npix+1, dtype=ITYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] pixfity = np.zeros(2*npix+1, dtype=DTYPE)

    for x in range(sz_x):
        # Deal with endpoints
        if x < npix or x > sz_x-npix:
            yfilt[x] = yarray[x]
            continue
        nmask = 0
        for i in range(-npix,npix+1):
            pixfity[i+npix] = yarray[i+x]
            pixfitm[i+npix] = 0
            if yarray[i+x] == maskval:
                nmask += 1
                pixfitm[i+npix] = 1
        if order + 1 >= 2*npix+1 - nmask:
            yfilt[x] = yarray[x]
            continue
        # Fit the data
        chisq = polyfit_mask(pixfitx,pixfity,pixfitm,(2*npix+1 - nmask),order+1,coeff)
        # Calculate the central value for this pixel
        yfilt[x] = coeff[0]
    return yfilt


#@cython.boundscheck(False)
def polyfit_mask(np.ndarray[DTYPE_t, ndim=1] x not None,
                np.ndarray[DTYPE_t, ndim=1] y not None,
                np.ndarray[ITYPE_t, ndim=1] mask not None,
                int nval, int degree, np.ndarray[DTYPE_t, ndim=1] coeffs not None):

    cdef int sz_x
    cdef int i, j, k
    cdef double chisq

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat
    cdef gsl_vector *yvec
    cdef gsl_vector *c

    sz_x = x.shape[0]

    #cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(degree, dtype=DTYPE)

    xmat = gsl_matrix_alloc(nval, degree)
    yvec = gsl_vector_alloc(nval)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)

    k = 0
    for i in range(0,sz_x):
        if mask[i] == 1: continue
        gsl_matrix_set(xmat, k, 0, 1.0)
        for j in range(0,degree):
            gsl_matrix_set(xmat, k, j, cpow(x[i], j))
        gsl_vector_set(yvec, k, y[i])
        k += 1

    # Fit with a polynomial
    ws = gsl_multifit_linear_alloc(nval, degree)
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

    # Return the best-fitting coefficients
    return chisq


#@cython.boundscheck(False)
def polyfit(np.ndarray[DTYPE_t, ndim=1] x not None,
        np.ndarray[DTYPE_t, ndim=1] y not None,
        int degree, np.ndarray[DTYPE_t, ndim=1] coeffs not None):

    cdef int sz_x
    cdef int i, j
    cdef double chisq

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat
    cdef gsl_vector *yvec
    cdef gsl_vector *c

    sz_x = x.shape[0]

    #cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(degree, dtype=DTYPE)

    xmat = gsl_matrix_alloc(sz_x, degree)
    yvec = gsl_vector_alloc(sz_x)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)

    for i in range(0,sz_x):
        gsl_matrix_set(xmat, i, 0, 1.0)
        for j in range(0,degree):
            gsl_matrix_set(xmat, i, j, cpow(x[i], j))
        gsl_vector_set(yvec, i, y[i])

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

    # Return the best-fitting coefficients
    return chisq
#	return coeffs


#@cython.boundscheck(False)
def polyfit_scan(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                    np.ndarray[DTYPE_t, ndim=1] yarr not None,
                    np.ndarray[DTYPE_t, ndim=1] warr not None,
                    double maskval, int polyorder, int polypoints, int repeat):

    cdef int x, sz_x, i, j, r, ers, cnt, upto
    cdef int degree = polyorder+1
    cdef double mval, chisq
    cdef int sz_p = polypoints/2
    cdef int sz_pp = 2*sz_p + 1

    sz_x = xarr.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] model = np.zeros(sz_x, dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] mask  = np.zeros(sz_pp, dtype=ITYPE)

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat

    cdef gsl_vector *yvec
    cdef gsl_vector *wgt
    cdef gsl_vector *c

    xmat = gsl_matrix_alloc(sz_pp, degree)
    wgt  = gsl_vector_alloc(sz_pp)
    yvec = gsl_vector_alloc(sz_pp)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)
    ws = gsl_multifit_linear_alloc(sz_pp, degree)

    for r in range(repeat):
        # Set the initial arrays
        cnt = sz_pp
        for i in range(sz_pp):
            if (i-sz_p < 0):
                for j in range(0,degree):
                    gsl_matrix_set(xmat, i, j, 0.0)
                gsl_vector_set(yvec, i, 0.0)
                gsl_vector_set(wgt, i, 0.0)
                mask[i] = 1
                cnt -= 1
            elif (yarr[i-sz_p] == maskval):
                for j in range(0,degree):
                    gsl_matrix_set(xmat, i, j, 0.0)
                gsl_vector_set(yvec, i, 0.0)
                gsl_vector_set(wgt, i, 0.0)
                mask[i] = 1
                cnt -= 1
            else:
                gsl_matrix_set(xmat, i, 0, 1.0)
                for j in range(1,degree):
                    gsl_matrix_set(xmat, i, j, cpow(xarr[i-sz_p], j))
                gsl_vector_set(yvec, i, yarr[i-sz_p])
                gsl_vector_set(wgt, i, warr[i-sz_p])
                mask[i] = 0
        ers  = sz_pp-1
        upto = sz_pp-1-sz_p
        for x in range(0,sz_x):
            if (cnt < degree*3):
                model[x] = maskval
            else:
                # Fit with a polynomial
                #gsl_multifit_linear(xmat, yvec, c, cov, &chisq, ws)
                gsl_multifit_wlinear(xmat, wgt, yvec, c, cov, &chisq, ws)
                # Obtain the model value at this location
                mval = 0.0
                for i in range(0,degree): mval += gsl_vector_get(c, i) * cpow(xarr[x], i)
                model[x] = mval

            # Update the arrays to be fit with the next iteration
            upto += 1
            ers  += 1
            if ers == sz_pp: ers = 0
            if (upto >= sz_x):
                if mask[ers] == 0:
                    #for j in range(0,degree):
                    #    gsl_matrix_set(xmat, ers, j, 0.0)
                    #gsl_vector_set(yvec, ers, 0.0)
                    gsl_vector_set(wgt, ers, 0.0)
                    mask[ers] = 1
                    cnt -= 1
            elif (yarr[upto] == maskval):
                if mask[ers] == 0:
                    #for j in range(0,degree):
                    #    gsl_matrix_set(xmat, ers, j, 0.0)
                    #gsl_vector_set(yvec, ers, 0.0)
                    gsl_vector_set(wgt, ers, 0.0)
                    mask[ers] = 1
                    cnt -= 1
            else:
                gsl_matrix_set(xmat, ers, 0, 1.0)
                for j in range(1,degree):
                    gsl_matrix_set(xmat, ers, j, cpow(xarr[upto], j))
                gsl_vector_set(yvec, ers, yarr[upto])
                gsl_vector_set(wgt, ers, warr[upto])
                if mask[ers] == 1:
                    mask[ers] = 0
                    cnt += 1

        if r != repeat-1:
            # Update yarr for the next loop
            for x in range(sz_x): yarr[x] = model[x]

    # Free some of these arrays from memory
    gsl_multifit_linear_free(ws)
    gsl_matrix_free(xmat)
    gsl_matrix_free(cov)
    gsl_vector_free(yvec)
    gsl_vector_free(c)

    return model


#@cython.boundscheck(False)
def polyfit_scan_lim(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                    np.ndarray[DTYPE_t, ndim=1] yarr not None,
                    np.ndarray[DTYPE_t, ndim=1] warr not None,
                    double maskval, int polyorder, int polypoints,
                    int repeat, double limval):

    cdef int x, sz_x, i, j, r, cnt, lmv, flag
    cdef int degree = polyorder+1
    cdef double mval, chisq
    cdef int sz_p = polypoints/2
    cdef int sz_pp = 2*sz_p + 1

    sz_x = xarr.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] model = np.zeros(sz_x, dtype=DTYPE)

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat

    cdef gsl_vector *yvec
    cdef gsl_vector *wgt
    cdef gsl_vector *c

    xmat = gsl_matrix_alloc(sz_pp, degree)
    wgt  = gsl_vector_alloc(sz_pp)
    yvec = gsl_vector_alloc(sz_pp)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)
    ws = gsl_multifit_linear_alloc(sz_pp, degree)

    for r in range(repeat):
        for x in range(0,sz_x):
            cnt = 0
            # Set the central value
            if yarr[x] == maskval:
                gsl_vector_set(wgt, sz_p, 0.0)
            else:
                gsl_matrix_set(xmat, sz_p, 0, 1.0)
                for j in range(1,degree):
                    gsl_matrix_set(xmat, sz_p, j, cpow(xarr[x], j))
                gsl_vector_set(yvec, sz_p, yarr[x])
                gsl_vector_set(wgt, sz_p, warr[x])
                cnt += 1
            # Step backwards until limval is reached
            lmv = x-sz_p
            if lmv<0: lmv = 0
            flag = 0
            for i in range(0,x-lmv):
                if (xarr[x-i]-xarr[x-i-1]>limval) or (flag==1):
                    gsl_vector_set(wgt, sz_p-i-1, 0.0)
                    flag = 1
                else:
                    gsl_matrix_set(xmat, sz_p-i-1, 0, 1.0)
                    for j in range(1,degree):
                        gsl_matrix_set(xmat, sz_p-i-1, j, cpow(xarr[x-i-1], j))
                    gsl_vector_set(yvec, sz_p-i-1, yarr[x-i-1])
                    gsl_vector_set(wgt, sz_p-i-1, warr[x-i-1])
                    cnt += 1
            # Step forwards until limval is reached
            lmv = x+sz_p
            if lmv>sz_x-1: lmv = sz_x-1
            flag = 0
            for i in range(0,lmv-x):
                if (xarr[x+i+1]-xarr[x+i]>limval) or (flag==1):
                    gsl_vector_set(wgt, sz_p+i+1, 0.0)
                    flag = 1
                else:
                    gsl_matrix_set(xmat, sz_p+i+1, 0, 1.0)
                    for j in range(1,degree):
                        gsl_matrix_set(xmat, sz_p+i+1, j, cpow(xarr[x+i+1], j))
                    gsl_vector_set(yvec, sz_p+i+1, yarr[x+i+1])
                    gsl_vector_set(wgt, sz_p+i+1, warr[x+i+1])
                    cnt += 1
            # Fit the data
            if (cnt < degree*3):
                model[x] = maskval
            else:
                # Fit with a polynomial
                gsl_multifit_wlinear(xmat, wgt, yvec, c, cov, &chisq, ws)
                # Obtain the model value at this location
                mval = 0.0
                for i in range(0,degree): mval += gsl_vector_get(c, i) * cpow(xarr[x], i)
                model[x] = mval
        if r != repeat-1:
            # Update yarr for the next loop
            for x in range(sz_x): yarr[x] = model[x]

    # Free some of these arrays from memory
    gsl_multifit_linear_free(ws)
    gsl_matrix_free(xmat)
    gsl_matrix_free(cov)
    gsl_vector_free(yvec)
    gsl_vector_free(c)

    return model


#@cython.boundscheck(False)
def polyfit_scan_intext(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                        np.ndarray[DTYPE_t, ndim=1] yarr not None,
                        np.ndarray[DTYPE_t, ndim=1] warr not None,
                        np.ndarray[ITYPE_t, ndim=1] mask not None,
                        int polyorder, int polypoints, int repeat,
                        double maskval):
    """
    A scanning polynomial that will interpolate/extrapolate over
    values where:
    mask[i] = 1
    """

    cdef int x, sz_x, i, j, r, cnt, num
    cdef int degree = polyorder + 1
    cdef double mval, chisq
    cdef int sz_p = polypoints/2
    cdef int sz_pp = 2*sz_p + 1
    cdef int stend

    sz_x = xarr.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] model = np.zeros(sz_x, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] ycopy = np.zeros(sz_x, dtype=DTYPE)

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat

    cdef gsl_vector *yvec
    cdef gsl_vector *wgt
    cdef gsl_vector *c

    xmat = gsl_matrix_alloc(sz_pp, degree)
    wgt  = gsl_vector_alloc(sz_pp)
    yvec = gsl_vector_alloc(sz_pp)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)
    ws = gsl_multifit_linear_alloc(sz_pp, degree)

    # Make a copy of yarr
    for x in range(0, sz_x):
        ycopy[x] = yarr[x]

    for r in range(repeat):
        stend = 0
        for x in range(0,sz_x):
            cnt = 0
            if stend == 0:
                if mask[x] == 0:
                    # We've reached the first non-masked value
                    stend = 1
                else:
                    model[x] = maskval
                    continue
            # Set the central value
            if mask[x] == 1:
                gsl_vector_set(wgt, sz_p, 0.0)
            else:
                gsl_matrix_set(xmat, sz_p, 0, 1.0)
                for j in range(1,degree):
                    gsl_matrix_set(xmat, sz_p, j, cpow(xarr[x], j))
                gsl_vector_set(yvec, sz_p, ycopy[x])
                gsl_vector_set(wgt, sz_p, warr[x])
                cnt += 1
            # Step backwards until sz_p unmasked pixels have been included in the fit
            num = sz_p-1
            for i in range(1, x+1):
                if (mask[x-i] == 1):
                    continue
                else:
                    gsl_matrix_set(xmat, num, 0, 1.0)
                    for j in range(1, degree):
                        gsl_matrix_set(xmat, num, j, cpow(xarr[x-i], j))
                    gsl_vector_set(yvec, num, ycopy[x-i])
                    gsl_vector_set(wgt, num, warr[x-i])
                    cnt += 1
                    num -= 1
                if num < 0:
                    break
            if num > 0:
                # Not enough pixels backwards to include in the fit -> mask those values
                for i in range(0, num):
                    gsl_vector_set(wgt, sz_p, 0.0)
            # Step forwards until sz_p unmasked pixels have been included in the fit
            num = sz_p+1
            for i in range(1, sz_x-x):
                if (mask[x+i] == 1):
                    continue
                else:
                    gsl_matrix_set(xmat, num, 0, 1.0)
                    for j in range(1, degree):
                        gsl_matrix_set(xmat, num, j, cpow(xarr[x+i], j))
                    gsl_vector_set(yvec, num, ycopy[x+i])
                    gsl_vector_set(wgt, num, warr[x+i])
                    cnt += 1
                    num += 1
                if num >= sz_pp:
                    break
            if num > 0:
                # Not enough pixels backwards to include in the fit -> mask those values
                for i in range(0, num):
                    gsl_vector_set(wgt, sz_p, 0.0)
            # Fit the data
            if (cnt < degree*3) and (stend == 1):
                stend = 0
                model[x] = maskval
                continue
            elif (cnt < degree*3):
                # Obtain the model value using the previous set of coefficients
                mval = 0.0
                for i in range(0,degree): mval += gsl_vector_get(c, i) * cpow(xarr[x], i)
                model[x] = mval
            else:
                # Fit with a polynomial
                gsl_multifit_wlinear(xmat, wgt, yvec, c, cov, &chisq, ws)
                # Obtain the model value at this location
                mval = 0.0
                for i in range(0,degree): mval += gsl_vector_get(c, i) * cpow(xarr[x], i)
                model[x] = mval
            # Step backwards if we we've reached the first unmasked value
            if stend == 1:
                for i in range(0, x):
                    # Obtain the model value at this location
                    mval = 0.0
                    for j in range(0,degree): mval += gsl_vector_get(c, j) * cpow(xarr[i], j)
                    model[i] = mval
                stend = 2
        if r != repeat-1:
            # Update ycopy for the next loop
            for x in range(sz_x): ycopy[x] = model[x]

    # Free some of these arrays from memory
    gsl_multifit_linear_free(ws)
    gsl_matrix_free(xmat)
    gsl_matrix_free(cov)
    gsl_vector_free(yvec)
    gsl_vector_free(c)

    return model


#@cython.boundscheck(False)
def prepare_bsplfit(np.ndarray[DTYPE_t, ndim=2] arcfr not None,
                    np.ndarray[DTYPE_t, ndim=3] pixmap not None,
                    np.ndarray[DTYPE_t, ndim=1] tilts not None,
                    np.ndarray[DTYPE_t, ndim=1] dtrc not None,
                    np.ndarray[DTYPE_t, ndim=1] scen not None,
                    np.ndarray[ITYPE_t, ndim=1] dpix not None,
                    np.ndarray[ITYPE_t, ndim=1] spix not None):

    cdef int p, b, donebl, donetl
    cdef int sz_p, sz_b
    cdef double la, lb, vb # Note, va=0 since the origin is defined as the bottom-left corner of the pixel
    cdef double dtl, stl, dbl, sbl
    cdef double cv, sva, dva, svb, dvb

    sz_b = dtrc.shape[0]
    sz_p = dpix.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] xfit = np.zeros(sz_p, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] yfit = np.zeros(sz_p, dtype=DTYPE)

    for p in range(sz_p):
        # get the bottom-left (i.e. the origin) and top-left corner of the pixel
        dbl = pixmap[dpix[p],spix[p],0] - 0.5*pixmap[dpix[p],spix[p],2]
        sbl = pixmap[dpix[p],spix[p],1] - 0.5*pixmap[dpix[p],spix[p],3]
        dtl = pixmap[dpix[p],spix[p],0] + 0.5*pixmap[dpix[p],spix[p],2]
        stl = pixmap[dpix[p],spix[p],1] + 0.5*pixmap[dpix[p],spix[p],3]
        donebl = 0
        donetl = 0
        for b in range(1,sz_b):
            # Bottom-left corner of pixel
            cv = dtrc[b] - tilts[b]*scen[b]
            sva = (sbl + tilts[b]*(dbl-cv))/(1.0 + tilts[b]*tilts[b])
            dva = tilts[b]*sva + cv
            if (dva > dbl) and (donebl == 0):
                # We have found the location above the pixel corner, now linearly interpolate to the pixel corner
                cv = dtrc[b-1] - tilts[b-1]*scen[b-1]
                svb = (sbl + tilts[b-1]*(dbl-cv))/(1.0 + tilts[b-1]*tilts[b-1])
                dvb = tilts[b-1]*svb + cv
                la = dtrc[b] - (dtrc[b]-dtrc[b-1])*(dva-dbl)/(dva-dvb)
                donebl = 1
            # Top-left corner of pixel
            cv = dtrc[b] - tilts[b]*scen[b]
            sva = (stl + tilts[b]*(dtl-cv))/(1.0 + tilts[b]*tilts[b])
            dva = tilts[b]*sva + cv
            if (dva > dtl) and (donetl == 0):
                cv = dtrc[b-1] - tilts[b-1]*scen[b-1]
                svb = (stl + tilts[b-1]*(dtl-cv))/(1.0 + tilts[b-1]*tilts[b-1])
                dvb = tilts[b-1]*svb + cv
                lb = dtrc[b] - (dtrc[b]-dtrc[b-1])*(dva-dtl)/(dva-dvb)
                vb = (dtl-dbl)/csqrt(1.0 + tilts[b]*tilts[b-1]) # This should really by tilts^2 at the tilt value corresponding to lb, but tilts[b]*tilts[b-1] should be close enough
                donetl = 1
            if (donebl==1) and (donetl==1):
                # Get the tilts value
                sva = tilts[b]
                break
        # Use the la, lb and vb to describe how dtrc varies over the pixel
        # Assuming dtrc = m*v + c, then c=la (i.e. v=0 for la), and m = (lb-la)/vb (assign this to be the variable cv)
        cv = (lb-la)/vb
        svb = csqrt(sva*sva + 1.0) / (sva*sva - 1.0)
        # Finally, calculate the appropriate value for xfit and yfit
        xfit[p] = la + 0.5*cv*svb*(pixmap[dpix[p],spix[p],2] + sva*pixmap[dpix[p],spix[p],3])
        yfit[p] = arcfr[dpix[p],spix[p]]
    return xfit, yfit


#######
#  Q  #
#######

#@cython.boundscheck(False)
def quicksort(np.ndarray[DTYPE_t, ndim=1] array not None,
                int k):

    cdef int i, j, l, m
    cdef DTYPE_t x, tmp
    l = 0
    m = <int> array.shape[0] - 1
    with nogil:
        while l < m:
            x = array[k]
            i = l
            j = m
            while 1:
                while array[i] < x: i += 1
                while x < array[j]: j -= 1
                if i <= j:
                    tmp = array[i]
                    array[i] = array[j]
                    array[j] = tmp
                    i += 1
                    j -= 1
                if i > j: break
            if j < k: l = i
            if k < i: m = j
    return


#######
#  R  #
#######

#@cython.boundscheck(False)
def resample_taper(np.ndarray[DTYPE_t, ndim=1] xnewedg not None,
                    np.ndarray[DTYPE_t, ndim=1] xold not None,
                    np.ndarray[DTYPE_t, ndim=1] yold not None):

    cdef int n, sz_n, o, sz_o
    cdef double repval, cntr

    sz_n = xnewedg.shape[0]
    sz_o = xold.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] ynew = np.zeros(sz_n-1, dtype=DTYPE)

    o = 0
    for n in range(sz_n-1):
        if xnewedg[n] < xold[0] or xnewedg[n] > xold[sz_o-1]:
            # Could set some value for ynew[n] here...
            continue
        repval = 0.0
        cntr = 0.0
        while xold[o] < xnewedg[n+1]:
            repval += yold[o]
            cntr += 1.0
            o += 1
        ynew[n] = repval/cntr
    return ynew


#@cython.boundscheck(False)
def robust_regression_full(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                np.ndarray[DTYPE_t, ndim=1] yarr not None,
                np.ndarray[ITYPE_t, ndim=2] rarr not None,
                int order):
    cdef int sz_x, sz_r, nrun
    cdef int n, r, i
    cdef double minmad, madt, chisq

    sz_x = xarr.shape[0]
    nrun = rarr.shape[0]
    sz_r = rarr.shape[1]

    # Setup the coefficients array
    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(order+1, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coefft = np.zeros(order+1, dtype=DTYPE)

    # Create some temporary arrays
    cdef np.ndarray[DTYPE_t, ndim=1] xtmp = np.zeros(sz_r, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] ytmp = np.zeros(sz_r, dtype=DTYPE)

    # Perform the first fit
    for r in range(sz_r):
        xtmp[r] = xarr[rarr[0,r]]
        ytmp[r] = yarr[rarr[0,r]]
    chisq = polyfit(xtmp,ytmp,order+1,coeffs)
    minmad = median(polydiff(xtmp,ytmp,coeffs))

    # Fit portions of the data during each realization
    for n in range(1,nrun):
        for r in range(sz_r):
            xtmp[r] = xarr[rarr[n,r]]
            ytmp[r] = yarr[rarr[n,r]]
        chisq = polyfit(xtmp,ytmp,order+1,coefft)
        madt = median(polydiff(xarr,yarr,coefft))
        if madt < minmad:
            minmad = madt
            for i in range(order+1):
                coeffs[i] = coefft[i]
    return coeffs


#######
#  S  #
#######

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


#@cython.boundscheck(False)
def spline_interpolate(np.ndarray[DTYPE_t, ndim=1] xnew not None,
                np.ndarray[DTYPE_t, ndim=1] coeffs not None,
                double a, double b):
    cdef int sz_x, x, i, l, m, gl, gu
    cdef double s, grd, icp
    cdef int n = coeffs.shape[0] - 3

    sz_x = xnew.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] ynew = np.zeros(sz_x, dtype=DTYPE)

    cdef double h = (b-a)/n
    gl=-1
    gu=-1
    for x in range(0,sz_x):
        if xnew[x] < a:
            gl = x
        elif xnew[x] > b:
            if gu == -1:
                gu = x
                grd = (ynew[x-1]-ynew[x-2])/(xnew[x-1]-xnew[x-2])
            s = ynew[gu-1] + grd*(xnew[x]-xnew[gu-1])
        else:
            s = 0.0
            l = <int>((xnew[x]-a)//h) + 1
            m = <int>(cmin(l+3,n+3))
            for i in xrange(l,m+1):
                s += coeffs[i-1] * spline_u(xnew[x], i, a, h)
        ynew[x] = s
    for x in range(0,gl+1):
        grd = (ynew[gl+1]-ynew[gl+2])/(xnew[gl+1]-xnew[gl+2])
        ynew[x] = ynew[gl+1] + grd*(xnew[x]-xnew[gl+1])
    return ynew

#@cython.boundscheck(False)
def spline_u(double x, int k, double a, double h):
    cdef double t = (x-a)/h - <double>(k-2)
    if t < 0.0: t *= -1.0
    if t <= 1.0:
        return (4.0 - 6.0*cpow(t,2) + 3.0*cpow(t,3))
    elif t <= 2.0:
        return cpow((2.0-t),3)
    else:
        return 0.0


#######
#  T  #
#######

#######
#  U  #
#######

#@cython.boundscheck(False)
def unband(np.ndarray[DTYPE_t, ndim=2] ab not None,
            int l, int u):
    cdef int i, j, sz_x, sz_y, im, ix
    sz_x = ab.shape[0]
    sz_y = ab.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] amat = np.zeros((sz_y,sz_y), dtype=DTYPE)
    for j in range(sz_y):
        im = u-j
        if im < 0: im = 0
        ix = im+sz_x
        if ix >= sz_x: ix=sz_x
        if j >= sz_y-l: ix -= 1
        for i in range(im,ix):
            amat[i+j-u,j]= ab[i,j]
    return amat

#######
#  V  #
#######

#######
#  W  #
#######

#######
#  X  #
#######

#######
#  Y  #
#######

#######
#  Z  #
#######


