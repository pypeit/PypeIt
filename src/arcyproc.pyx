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


#######
#  A  #
#######



#######
#  B  #
#######

@cython.boundscheck(False)
def background_model(np.ndarray[DTYPE_t, ndim=1] fitfunc not None,
                    np.ndarray[ITYPE_t, ndim=1] xpix not None,
                    np.ndarray[ITYPE_t, ndim=1] ypix not None,
                    int shA, int shB):
    cdef int p, sz_p
    sz_p = fitfunc.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=2] model = np.zeros((shA,shB), dtype=DTYPE)
    for p in range(sz_p): model[xpix[p],ypix[p]] = fitfunc[p]
    return model


@cython.boundscheck(False)
def blaze_fitends(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                    np.ndarray[DTYPE_t, ndim=1] yarr not None,
                    np.ndarray[DTYPE_t, ndim=1] marr not None,
                    np.ndarray[DTYPE_t, ndim=1] secderv not None,
                    double sigma, int polyord, int nfit):
    """
    Check the ends of the fit to the blaze function.
    Continue to mask end points until an acceptable
    N pixels are found. Extrapolate a low-order polynomial
    fit to the end pixels based on the good pixels
    """

    cdef int x, y, sz_x, nlo, nhi, tst, ngood
#	cdef int nfit = 32
    cdef double mval
    sz_x = xarr.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] modelfit = np.arange((sz_x), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] xfit = np.arange((nfit), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] yfit = np.arange((nfit), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.arange((polyord+1), dtype=DTYPE)

    # First fill in the modelfit array with the best-fitting model
    for x in range(0,sz_x): modelfit[x] = marr[x]

    # Fix the beginning of the order
    ngood = 0
    for x in range(0,sz_x-5):
        if (secderv[x] < 5.0*sigma) and (secderv[x] > -5.0*sigma):
            xfit[ngood] = xarr[x]
            yfit[ngood] = yarr[x]
            ngood += 1
        else:
            tst = 0
            for y in range(5):
                if (secderv[x+y] < 5.0*sigma) and (secderv[x+y] > -5.0*sigma):
                    tst += 1
            if tst <= 2:
                ngood = 0
        if ngood >= nfit:
            nlo = x
            break
    if nlo > nfit+5:
        # Perform a fit to the good pixels, otherwise, the blaze is probably good enough
        polyfit(xfit,yfit,polyord+1,coeffs)
        # Using these coefficients, fill in the array at bad pixel locations
        ngood = 0
        for x in range(0,nlo):
            mval = coeffs[0]
            for y in range(1,polyord+1):
                mval += coeffs[y] * cpow(xarr[x],y)
            modelfit[x] = mval

    # Fix the end of the order
    ngood = 0
    for x in range(1,sz_x-4):
        if (secderv[sz_x-x] < 5.0*sigma) and (secderv[sz_x-x] > -5.0*sigma):
            xfit[ngood] = xarr[sz_x-x]
            yfit[ngood] = yarr[sz_x-x]
            ngood += 1
        else:
            tst = 0
            for y in range(5):
                if (secderv[sz_x-x-y] < 5.0*sigma) and (secderv[sz_x-x-y] > -5.0*sigma):
                    tst += 1
            if tst <= 2:
                ngood = 0
        if ngood >= nfit:
            nhi = x
            break
    if sz_x-nhi > nfit+5:
        # Perform a fit to the good pixels, otherwise, the blaze is probably good enough
        polyfit(xfit,yfit,polyord+1,coeffs)
        # Using these coefficients, fill in the array at bad pixel locations
        ngood = 0
        for x in range(sz_x-nhi,sz_x):
            mval = coeffs[0]
            for y in range(1,polyord+1):
                mval += coeffs[y] * cpow(xarr[x],y)
            modelfit[x] = mval
    return modelfit


#######
#  C  #
#######

@cython.boundscheck(False)
def combine_nrmflat(np.ndarray[DTYPE_t, ndim=2] msframe not None,
                    np.ndarray[DTYPE_t, ndim=2] frame not None,
                    np.ndarray[ITYPE_t, ndim=1] pixcen not None,
                    np.ndarray[ITYPE_t, ndim=1] pixledg not None,
                    np.ndarray[ITYPE_t, ndim=1] pixredg not None,
                    int pixwid, double maskval, int dispaxis):
    """
    Replace the master normalized flat frame values with the values for this
    order. If a value already exists in the master frame take the average
    of the two values.
    """

    cdef int x, y, sz_x, sz_y
    cdef int yidx, ymin, ymax, sz_p
    cdef int shift = pixwid/2

    sz_x  = frame.shape[dispaxis]
    sz_y  = frame.shape[1-dispaxis]
    sz_p = 2*shift+1

    for x in range(0,sz_x):
        for y in range(0,sz_p):
            if pixcen[x] <= 0:
                yidx = y+(pixredg[x]-shift)-shift
            elif pixcen[x] >= sz_y-1:
                yidx = y+(pixledg[x]+shift)-shift
            else:
                yidx = y+pixcen[x]-shift
            if pixredg[x] <= 0 or pixledg[x] >= sz_y-1 or yidx < 0 or yidx >= sz_y:
                continue
            else:
                if dispaxis == 0:
                    if msframe[x,yidx] != maskval:
                        msframe[x,yidx] += frame[x,yidx]
                        msframe[x,yidx] /= 2.0
                    else: msframe[x,yidx] = frame[x,yidx]
                else:
                    if msframe[yidx,x] != maskval:
                        msframe[yidx,x] += frame[yidx,x]
                        msframe[yidx,x] /= 2.0
                    else: msframe[yidx,x] = frame[yidx,x]
    return msframe


#######
#  M  #
#######

#######
#  N  #
#######

#######
#  O  #
#######



#######
#  P  #
#######

@cython.boundscheck(False)
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
    return


#@cython.boundscheck(False)
def polyscan_fitsky(np.ndarray[DTYPE_t, ndim=2] xarr not None,
                    np.ndarray[DTYPE_t, ndim=2] yarr not None,
                    np.ndarray[DTYPE_t, ndim=2] warr not None,
                    double maskval, int polyorder, int polypoints,
                    int nsmth, int repeat):

    cdef int x, sz_x, y, sz_y, i, j, k, r, cnt, ngood
    cdef int degree = polyorder+1
    cdef double mval, chisq
    cdef int sz_p = polypoints/2
    cdef int sz_pp = 2*sz_p + 1

    # Make sure nsmth is odd
    nsmth /= 2   # later, (2*nsmth + 1)  is used

    sz_x = xarr.shape[0]
    sz_y = xarr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] model = np.zeros((sz_x,sz_y), dtype=DTYPE)

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat

    cdef gsl_vector *yvec
    cdef gsl_vector *wgt
    cdef gsl_vector *c

    xmat = gsl_matrix_alloc(sz_pp*(2*nsmth+1), degree)
    wgt  = gsl_vector_alloc(sz_pp*(2*nsmth+1))
    yvec = gsl_vector_alloc(sz_pp*(2*nsmth+1))
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)
    ws   = gsl_multifit_linear_alloc(sz_pp*(2*nsmth+1), degree)

    for r in range(repeat):
        for x in range(0,sz_x):
            for y in range(0,sz_y):
                if yarr[x,y] == maskval:
                    model[x,y] = maskval
                    continue
                # Fill in the arrays
                cnt = 0
                ngood = 0
                for i in range(-sz_p,sz_p+1):
                    if (x+i < 0) or (x+i >=sz_x):
                        for j in range(0,2*nsmth+1):
                            gsl_vector_set(wgt, cnt, 0.0)
                            cnt += 1
                        continue
                    for j in range(-nsmth,nsmth+1):
                        if (y+j >= 0) and (y+j < sz_y):
                            if (yarr[x+i, y+j] == maskval):
                                gsl_vector_set(wgt, cnt, 0.0)
                            else:
                                gsl_matrix_set(xmat, cnt, 0, 1.0)
                                for k in range(1,degree):
                                    gsl_matrix_set(xmat, cnt, k, cpow(xarr[x+i,y+j], k))
                                gsl_vector_set(yvec, cnt, yarr[x+i,y+j])
                                gsl_vector_set(wgt, cnt, warr[x+i,y+j])
                                ngood += 1
                        else: gsl_vector_set(wgt, cnt, 0.0)
                        cnt += 1
                # Fit the data
                if (ngood < degree*nsmth):
                    model[x,y] = maskval
                else:
                    # Fit with a polynomial
                    gsl_multifit_wlinear(xmat, wgt, yvec, c, cov, &chisq, ws)
                    # Obtain the model value at this location
                    mval = 0.0
                    for i in range(0,degree): mval += gsl_vector_get(c, i) * cpow(xarr[x,y], i)
                    model[x,y] = mval
            # Update the arrays to be fit with the next iteration
        if r != repeat-1:
            # Update yarr for the next loop
            for x in range(sz_x):
                for y in range(sz_y):
                    yarr[x,y] = model[x,y]

    # Free some of these arrays from memory
    gsl_multifit_linear_free(ws)
    gsl_matrix_free(xmat)
    gsl_matrix_free(cov)
    gsl_vector_free(yvec)
    gsl_vector_free(c)

    return model


#@cython.boundscheck(False)
def prepare_bgmodel(np.ndarray[DTYPE_t, ndim=2] arcfr not None,
                    np.ndarray[DTYPE_t, ndim=3] pixmap not None,
                    np.ndarray[DTYPE_t, ndim=1] tilts not None,
                    np.ndarray[DTYPE_t, ndim=1] dtrc not None,
                    np.ndarray[DTYPE_t, ndim=1] scen not None,
                    np.ndarray[ITYPE_t, ndim=1] dpix not None,
                    np.ndarray[ITYPE_t, ndim=1] spix not None,
                    int dispdir):

    cdef int p, b, donebl, donetl
    cdef int sz_p, sz_b
    cdef double la, lb, vb # Note, va=0 since the origin is defined as the bottom-left corner of the pixel
    cdef double dtl, stl, dbl, sbl
    cdef double cv, sva, dva, svb, dvb

    sz_b = dtrc.shape[0]
    sz_p = dpix.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] xfit = np.zeros(sz_p, dtype=DTYPE)

    for p in range(sz_p):
        # get the bottom-left (i.e. the origin) and top-left corner of the pixel
        if dispdir == 0:
            dbl = pixmap[dpix[p],spix[p],0] - 0.5*pixmap[dpix[p],spix[p],2]
            sbl = pixmap[dpix[p],spix[p],1] - 0.5*pixmap[dpix[p],spix[p],3]
            dtl = pixmap[dpix[p],spix[p],0] + 0.5*pixmap[dpix[p],spix[p],2]
            stl = pixmap[dpix[p],spix[p],1] + 0.5*pixmap[dpix[p],spix[p],3]
        else:
            dbl = pixmap[spix[p],dpix[p],1] - 0.5*pixmap[spix[p],dpix[p],3]
            sbl = pixmap[spix[p],dpix[p],0] - 0.5*pixmap[spix[p],dpix[p],2]
            dtl = pixmap[spix[p],dpix[p],1] + 0.5*pixmap[spix[p],dpix[p],3]
            stl = pixmap[spix[p],dpix[p],0] + 0.5*pixmap[spix[p],dpix[p],2]
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
        if dispdir == 0:
            xfit[p] = la + 0.5*cv*svb*(pixmap[dpix[p],spix[p],2] + sva*pixmap[dpix[p],spix[p],3])
        else:
            xfit[p] = la + 0.5*cv*svb*(pixmap[spix[p],dpix[p],3] + sva*pixmap[spix[p],dpix[p],2])
    return xfit

#######
#  Q  #
#######


#######
#  R  #
#######

@cython.boundscheck(False)
def remove_maskedends(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                    np.ndarray[DTYPE_t, ndim=1] yarr not None,
                    np.ndarray[DTYPE_t, ndim=1] marr not None,
                    double maskval):
    cdef int i, nl, nr, sz
    sz = xarr.shape[0]

    # Return the input if there is nothing to strip
    if yarr[0] != maskval and yarr[sz-1] != maskval: return xarr, yarr, marr, 0, sz

    nl = 0
    i = 0
    while yarr[i] == maskval:
        nl += 1
        i += 1
    nr = 0
    i=sz-1
    while yarr[i] == maskval:
        nr += 1
        i -= 1

    # Fill in the array
    cdef np.ndarray[DTYPE_t, ndim=1] outxarr = np.zeros(sz-nl-nr, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] outyarr = np.zeros(sz-nl-nr, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] outmarr = np.zeros(sz-nl-nr, dtype=DTYPE)
    for i in range(nl,sz-nr):
        outxarr[i-nl] = xarr[i]
        outyarr[i-nl] = yarr[i]
        outmarr[i-nl] = marr[i]

    return outxarr, outyarr, outmarr, nl, sz-nr


#######
#  S  #
#######

@cython.boundscheck(False)
def scale_blaze(np.ndarray[DTYPE_t, ndim=2] recsort not None,
                    double maskval):
    cdef int d, sz_d, s, sz_s
    cdef int num
    cdef double mval, cntr

    sz_d = recsort.shape[0]
    sz_s = recsort.shape[1]

    # mean along the dispersion direction
    cdef np.ndarray[DTYPE_t, ndim=1] recmean = np.zeros(sz_s, dtype=DTYPE)

    for s in range(sz_s):
        # Determine the number of non-masked values
        num = 0
        str = -1
        for d in range(sz_d):
            if recsort[d,s] != maskval:
                if str == -1: str = d
                num += 1
        mval = 0.0
        cntr = 0.0
        for d in range(str+num/4,str+(3*num)/4):
            mval += recsort[d,s]
            cntr += 1.0
        recmean[s] = mval/cntr
    return recmean


@cython.boundscheck(False)
def smooth_gaussmask(np.ndarray[DTYPE_t, ndim=2] flux not None,
                    double maskval, double nsig):
    """
    Smooths a spectrum along the spectral direction

    nsig = number of pixels in 1 standard deviation
    """
    cdef int d, sz_d, s, sz_s, dmin, dmax, m
    cdef double mval, wght, temp
    cdef int smwid = <int>(6.0*nsig) # Go to +/-6 sigma

    sz_d = flux.shape[0]
    sz_s = flux.shape[1]

    # smoothed version (along the dispersion direction)
    cdef np.ndarray[DTYPE_t, ndim=2] fluxsmooth = np.zeros((sz_d,sz_s), dtype=DTYPE)

    for s in range(sz_s):
        # Determine the number of non-masked values
        for d in range(sz_d):
            dmin = d-smwid
            dmax = d+smwid+1
            if dmin < 0: dmin = 0
            if dmax > sz_d: dmax = sz_d
            mval = 0.0
            wght = 0.0
            for m in range(dmin,dmax):
                if flux[m,s] != maskval:
                    temp = cexp( -<double>((m-d)*(m-d))/(2.0*nsig*nsig) )
                    mval += flux[m,s]*temp
                    wght += temp
            if wght == 0.0: fluxsmooth[d,s] = maskval
            else: fluxsmooth[d,s] = mval/wght
    return fluxsmooth

#######
#  T  #
#######

#######
#  U  #
#######

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


