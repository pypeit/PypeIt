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

#######
#  A  #
#######


#######
#  B  #
#######


#######
#  C  #
#######


@cython.boundscheck(False)
def calc_angperpix(np.ndarray[DTYPE_t, ndim=1] order not None,
                    np.ndarray[DTYPE_t, ndim=1] app not None,
                    double stepsize, double startval,
                    int nord, int nsteps):
    cdef int o, n, x, sz_x, xstr, u, no
    cdef int found, minx, minn
    cdef double curval, minv, tstv, tstmv, minm

    cdef np.ndarray[DTYPE_t, ndim=1] closest = np.zeros(nord, dtype=DTYPE)

    sz_x = order.shape[0]

    for n in range(0,nsteps):
        curval = startval + stepsize * <double>(n)
        xstr = 0
        for o in range(nord):
            found = 0
            for x in range(xstr,sz_x):
                if order[x] != o:
                    xstr = x
                    break
                # Find the closest value to curval in this order
                tstv = curval-app[x]
                if tstv < 0.0: tstv *= -1.0
                if found == 0:
                    minv = tstv
                    minx = x
                else:
                    if tstv < minv:
                        minv = tstv
                        minx = x
                found = 1
            if found == 0:
                closest[o] = -999999.9 # Mask this order out, it contains no values.
            else:
                closest[o] = minv
        if n == 0:
            tstmv = median_mask(closest,-999999.9)
            if tstmv < 0.0: tstmv *= -1.0
            minn = n
            minm = tstmv
        else:
            tstmv = median_mask(closest,-999999.9)
            if tstmv < 0.0: tstmv *= -1.0
            if tstmv < minm:
                minm = tstmv
                minn = n

    # For each order, find the 'b' closest values to the best curval
    cdef int b = 2
    cdef np.ndarray[ITYPE_t, ndim=1] retarr = np.zeros(nord*b, dtype=ITYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] clsarr = np.zeros(b, dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] clsirr = np.zeros(b, dtype=ITYPE)

    xstr = 0
    no = 0
    curval = startval + stepsize * <double>(minn)
    for o in range(nord):
        found = 0
        # Mask the clsarr
        for n in range(b):
            clsarr[n] = -999999.9
        for x in range(xstr,sz_x):
            if order[x] != o:
                xstr = x
                break
            # Find the closest value to curval in this order
            tstv = curval-app[x]
            if tstv < 0.0: tstv *= -1.0
            for n in range(b):
                if clsarr[n] == -999999.9:
                    clsarr[n] = tstv
                    clsirr[n] = x
                    break
                elif tstv < clsarr[n]: # Shift everything down one
                    for u in range(0,b-n-1):
                        clsarr[b-1-u] = clsarr[b-2-u]
                        clsirr[b-1-u] = clsirr[b-2-u]
                    clsarr[n] = tstv
                    clsirr[n] = x
                    break
        # Fill in the returned array
        for n in range(b):
            if clsarr[n] != -999999.9:
                retarr[no] = clsirr[n]
                no += 1
    # Fill in the unused sections in the array with a masked value
    for n in range(no,nord*b):
        retarr[n] = -999999
    return retarr


@cython.boundscheck(False)
def calculate_lineprob_bffit(np.ndarray[DTYPE_t, ndim=1] pixls not None,
                        np.ndarray[DTYPE_t, ndim=1] carr not None,
                        np.ndarray[DTYPE_t, ndim=1] waves not None,
                        double aval, double bval, double pmean, double crit, double maxpix):

    cdef int sz_w, sz_p, sz_i
    cdef int w, p, i, j, ibst
    cdef double wbst, tst, medv, madv
    cdef int ipar, jpar
    cdef double minmad, pval, pprev

    sz_w  = waves.shape[0]
    sz_p  = pixls.shape[0]
    sz_i  = carr.shape[0]

    # Store the most likely wavelengths
    cdef np.ndarray[DTYPE_t, ndim=1] wvcls   = np.zeros((sz_p), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] prbpixl = np.zeros((sz_w), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] prbmtrx = np.zeros((sz_w), dtype=DTYPE)

    for i in range(sz_i):
        for j in range(sz_i):
            # Calculate wvcls
            for p in range(sz_p):
                pval = (pixls[p]-pmean)/maxpix
                wvcls[p] = aval + bval*pval + carr[j]*cpow(pval,2.0) + carr[i]*cpow(pval,3.0)
            for w in range(0,sz_w):
                wbst = wvcls[0]-waves[w]
                if wbst < 0.0: wbst *= -1.0
                ibst = 0
                pprev = wbst
                for p in range(1,sz_p):
                    tst = wvcls[p]-waves[w]
                    if tst < 0.0: tst *= -1.0
                    #if tst > pprev: break
                    if tst < wbst:
                        wbst = tst
                        ibst = p
                    pprev = tst
                prbpixl[w] = wbst
            # Calculate the MAD
            medv, madv = medianmad_qs(prbpixl)
            if i == 0 and j == 0:
                minmad = madv
                ipar = i
                jpar = j
            else:
                if madv < minmad:
                    minmad = madv
                    ipar = i
                    jpar = j

    # Construct the best-fitting values
    for p in range(sz_p):
        pval = (pixls[p]-pmean)/maxpix
        wvcls[p] = aval + bval*pval + carr[jpar]*cpow(pval,2.0) + carr[ipar]*cpow(pval,3.0)
    for w in range(0,sz_w):
        wbst = wvcls[0]-waves[w]
        if wbst < 0.0: wbst *= -1.0
        ibst = 0
        pprev = wbst
        for p in range(1,sz_p):
            tst = wvcls[p]-waves[w]
            if tst < 0.0: tst *= -1.0
            #if tst > pprev: break
            if tst < wbst:
                wbst = tst
                ibst = p
            pprev = tst
        if wbst < crit:
            prbpixl[w] = pixls[ibst]
            prbmtrx[w] = 1.0
        else:
            prbpixl[w] = pixls[ibst]
            prbmtrx[w] = -1.0

    # return the best values
    return prbpixl, prbmtrx, ipar, jpar


@cython.boundscheck(False)
def calculate_lineprob_iter(np.ndarray[DTYPE_t, ndim=1] pixls not None,
                        np.ndarray[DTYPE_t, ndim=2] pararr not None,
                        np.ndarray[DTYPE_t, ndim=1] waves not None,
                        np.ndarray[ITYPE_t, ndim=1] ford not None,
                        np.ndarray[ITYPE_t, ndim=1] oidx not None,
                        int order, double crit, double maxpix):

    cdef int sz_w, sz_p, sz_r, sz_i
    cdef int w, p, i, r, ibst
    cdef double wbst, tst, medv, madv
    cdef int ipar
    cdef double minmad

    sz_w  = waves.shape[0]
    sz_p  = pixls.shape[0]
    sz_r  = pararr.shape[0]
    sz_i  = pararr.shape[1]

    # Store the most likely wavelengths
    cdef np.ndarray[DTYPE_t, ndim=1] wvcls   = np.zeros((sz_p), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] pxsnd   = np.zeros((sz_p), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] ordarr  = np.zeros((sz_p), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] prbpixl = np.zeros((sz_w), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] prbmtrx = np.zeros((sz_w), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] params  = np.zeros((sz_r), dtype=DTYPE)

    for i in range(sz_p):
        ordarr[i] = <double>(order)
        pxsnd[i]  = pixls[i]/maxpix

    for i in range(sz_i):
        for r in range(sz_r):
            params[r] = pararr[r,i]
        # Reset and then calculate wvcls
        for p in range(sz_p): wvcls[p] = 0.0
        func_pcsurf_ret(pxsnd, ordarr, params, ford, oidx, wvcls)
        for w in range(0,sz_w):
            wbst = wvcls[0]-waves[w]
            if wbst < 0.0: wbst *= -1.0
            ibst = 0
            for p in range(1,sz_p):
                tst = wvcls[p]-waves[w]
                if tst < 0.0: tst *= -1.0
                if tst < wbst:
                    wbst = tst
                    ibst = p
            prbpixl[w] = wbst
        # Calculate the MAD
        medv, madv = medianmad(prbpixl)
        if i != 0:
            if madv < minmad:
                minmad = madv
                ipar = i
        else:
            minmad = madv
            ipar = i

    # Construct the best-fitting values
    for r in range(sz_r):
        params[r] = pararr[r,ipar]
    # Reset and then calculate wvcls
    for p in range(sz_p): wvcls[p] = 0.0
    func_pcsurf_ret(pxsnd, ordarr, params, ford, oidx, wvcls)
    for w in range(0,sz_w):
        wbst = wvcls[0]-waves[w]
        if wbst < 0.0: wbst *= -1.0
        ibst = 0
        for p in range(1,sz_p):
            tst = wvcls[p]-waves[w]
            if tst < 0.0: tst *= -1.0
            if tst < wbst:
                wbst = tst
                ibst = p
        if wbst < crit:
            prbpixl[w] = pixls[ibst]
            prbmtrx[w] = 1.0
        else:
            prbpixl[w] = pixls[ibst]
            prbmtrx[w] = -1.0
    # return the best values
    return prbpixl, prbmtrx, ipar


@cython.boundscheck(False)
def calibrate_order(np.ndarray[DTYPE_t, ndim=1] wcls not None,
                    np.ndarray[DTYPE_t, ndim=1] waves not None):
    cdef int w, sz_w,  c, sz_c
    cdef int bidx
    cdef double bdiff, tval

    sz_w = waves.shape[0]
    sz_c = wcls.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] bwave = np.zeros(sz_c, dtype=DTYPE)

    for c in range(0,sz_c):
        bdiff = waves[0]-wcls[c]
        if bdiff < 0.0: bdiff *= -1.0
        bidx = 0
        for w in range(1,sz_w):
            tval = waves[w]-wcls[c]
            if tval < 0.0: tval *= -1.0
            if tval < bdiff:
                bdiff = tval
                bidx = w
        bwave[c] = waves[bidx]
    return bwave


#@cython.boundscheck(False)
def centre_wave(np.ndarray[DTYPE_t, ndim=1] pixarr not None,
                np.ndarray[DTYPE_t, ndim=1] wavarr not None,
                double maxpix, int pixcnt):

    cdef int i
    cdef double chisq
    # Setup the arrays for fitting
    cdef np.ndarray[DTYPE_t, ndim=1] pixfit = np.zeros(pixcnt, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] wavfit = np.zeros(pixcnt, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(2, dtype=DTYPE)
    for i in range(pixcnt):
        pixfit[i] = pixarr[i]
        wavfit[i] = wavarr[i]
    # Perform linear regression
    chisq = polyfit(pixfit,wavfit,2,coeffs)  # i.e. y=mx+c
    cdef double wavecen = coeffs[0] + maxpix*coeffs[1]/2.0
    return wavecen


#######
#  D  #
#######

@cython.boundscheck(False)
def detections(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                np.ndarray[ITYPE_t, ndim=1] satmask not None):
    cdef int sz_p
    cdef int p, c, d, npc, pc

    sz_p = pixels.shape[0]
    cdef np.ndarray[ITYPE_t, ndim=1] pixcen = np.zeros(sz_p, dtype=ITYPE)

    c = 0
    d = 0
    npc = 0
    pixcen[sz_p-1] = -1
    for p in range(1,sz_p-1):
        if pixels[p] == 0.0 or pixels[p-1] == 0.0 or pixels[p+1] == 0.0: continue
        pixcen[p] = -1
        if c != 0 or d != 0:
            # If you have a saturated pixel, ignore it
            if satmask[p] == 1:
                c = 0
                d = 0
        if pixels[p]-pixels[p-1] >= 0.0:
            c += 1
        if (d == 0) and (c >= 5) and (pixels[p+1]-pixels[p] < 0.0): # Candidate detection for up
            c = 0
            d += 1
            pc = p
        elif d > 0:
            if (pixels[p+1]-pixels[p] < 0.0): d += 1
            else: d = 0
        if d >= 5:
            pixcen[npc] = 1+pc
            npc += 1
            d = 0
    return pixcen

@cython.boundscheck(False)
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
                if pixels[p+d] > pixels[p+d-1] and pixels[p+d] > pixels[p+d+1]:
                    flag = 1
                pixcen[p+d] = -1
                psum += pixels[p+d]
                esum += errors[p+d]*errors[p+d]
                if pixels[p+d]>mnum:
                    mnum = pixels[p+d]
                    inum = p+d
                d += 1
                if pixels[p+d] > pixels[p+d-1] and flag == 1: break
                if pixels[p+d] <= 0.0: break
                if p+d >= sz_p: break
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



@cython.boundscheck(False)
def detections_allorders(np.ndarray[DTYPE_t, ndim=2] pixels not None,
                        np.ndarray[ITYPE_t, ndim=2] satmask not None):
    cdef int sz_o, sz_p
    cdef int o, p, c, d, npc, pc

    sz_p = pixels.shape[0]
    sz_o = pixels.shape[1]
    cdef np.ndarray[ITYPE_t, ndim=2] pixcen = np.zeros((sz_p,sz_o), dtype=ITYPE)

    for o in range(0,sz_o):
        c = 0
        d = 0
        npc = 0
        pixcen[0,o] = -1
        pixcen[sz_p-1,o] = -1
        for p in range(1,sz_p-1):
            if pixels[p,o] == 0.0 or pixels[p-1,o] == 0.0 or pixels[p+1,o] == 0.0: continue
            pixcen[p,o] = -1
            if c != 0 or d != 0:
                # If you have a saturated pixel, ignore it
                if satmask[p,o] == 1:
                    c = 0
                    d = 0
            if pixels[p,o]-pixels[p-1,o] >= 0.0:
                c += 1
            else:
                c = 0
            if (d == 0) and (c >= 4) and (pixels[p+1,o]-pixels[p,o] < 0.0): # Candidate detection for up
                c = 0
                d += 1
                pc = p
            elif d > 0:
                if (pixels[p+1,o]-pixels[p,o] < 0.0): d += 1
                else: d = 0
            if d >= 4:
                pixcen[npc,o] = 1+pc
                npc += 1
                d = 0
    return pixcen

#######
#  E  #
#######


#######
#  F  #
#######

@cython.boundscheck(False)
def find_closest(np.ndarray[DTYPE_t, ndim=1] arclist not None,
                double waveval, int p):
    cdef int sz_w
    cdef int w
    cdef double tempv

    cdef double clswav = arclist[p]-waveval
    if clswav == 0.0: return arclist[p]
    sz_w = arclist.shape[0]
    cdef int tp = p
    cdef int dir
    # Set the direction to be searching (1 is increase, -1 is decrease)
    if clswav > 0.0:
        dir = -1
    else:
        clswav *= -1.0
        dir = 1
    cdef int wa = 0
    while True:
        tp += dir
        if tp == -1 or tp == sz_w: return arclist[wa]
        tempv = arclist[tp]-waveval
        if tempv < 0.0: tempv *= -1.0
        if tempv < clswav:
            clswav = tempv
            wa = tp
        else:
            return arclist[wa]

#@cython.boundscheck(False)
def find_nearest(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                np.ndarray[DTYPE_t, ndim=1] ppix not None,
                int p):
    cdef int sz_p, sz_pp, pp
    cdef double clswav, tempv

    sz_pp = ppix.shape[0]
    sz_p  = pixels.shape[0]
    cdef np.ndarray[ITYPE_t, ndim=1] apix = np.zeros(sz_pp, dtype=ITYPE)

    cdef int tp, dir, wa
    for pp in range(sz_pp):
        if ppix[pp] == 0.0:
            apix[pp] = -1
            continue
        clswav = pixels[p]-ppix[pp]
        if clswav == 0.0:
            apix[pp] = p
            continue
        tp = p
        # Set the direction to be searching (1 is increase, -1 is decrease)
        if clswav > 0.0:
            dir = -1
        else:
            clswav *= -1.0
            dir = 1
        wa = 0
        while True:
            tp += dir
            if tp == -1 or tp == sz_p:
                apix[pp] = wa
                break
            tempv = pixels[tp]-ppix[pp]
            if tempv < 0.0: tempv *= -1.0
            if tempv < clswav:
                clswav = tempv
                wa = tp
            else:
                apix[pp] = wa
                break
    return apix

@cython.boundscheck(False)
def find_patterns(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                np.ndarray[DTYPE_t, ndim=2] patterns not None):
    cdef int sz_m, sz_p
    cdef int m, pa, pb
    cdef int i, j, k
    cdef double chisq, cnum

    sz_m = patterns.shape[0]
    sz_p = pixels.shape[0]

    # the arrays that are used to estimate the pixels
    cdef np.ndarray[DTYPE_t, ndim=1] ppix = np.zeros(4, dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] apix

    # During this calculation, store a probability matrix of values that
    #represent the most likely wavelengths
    cdef int strn = 5
    cdef np.ndarray[DTYPE_t, ndim=2] prbmtrx = np.zeros((sz_p,strn), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] prbwave = np.zeros((sz_p,strn), dtype=DTYPE)
    for p in range(sz_p):
        prbmtrx[p,0] = -1.0
        prbmtrx[p,1] = -1.0
        prbmtrx[p,2] = -1.0
        prbmtrx[p,3] = -1.0
        prbmtrx[p,4] = -1.0

    for m in range(sz_m):
        print m, "/", sz_m-1
        for pa in range(sz_p-2):
            for pb in range(pa+1,sz_p-1):
                # First search assuming increasing pixels corresponds to increasing wavelength
                diff = (pixels[pb]-pixels[pa])
                for i in range(1,5):
                    if patterns[m,i] == 0.0: ppix[i-1] = 0.0
                    else: ppix[i-1] = pixels[pa] + patterns[m,i]*diff
                apix = find_nearest(pixels, ppix, pb)
                chisq = 0.0
                cnum = 0.0
                for i in range(4):
                    if apix[i] == -1: continue
                    chisq += (pixels[apix[i]]-ppix[i])**2
                    cnum += 1.0
                if cnum == 0.0: continue
                chisq /= cnum
                for i in range(4):
                    if apix[i] == -1: continue
                    if (prbmtrx[apix[i],strn-1] > chisq) or (prbmtrx[apix[i],strn-1] == -1.0):
                        for j in range(strn):
                            if prbmtrx[apix[i],j] == -1.0:
                                prbmtrx[apix[i],j] = chisq
                                prbwave[apix[i],j] = patterns[m,7+i]
                                break
                            elif chisq < prbmtrx[apix[i],j]:
                                # Shift everything down one
                                for k in range(0,strn-1-j):
                                    prbmtrx[apix[i],strn-1-k] = prbmtrx[apix[i],strn-k-2]
                                    prbwave[apix[i],strn-1-k] = prbwave[apix[i],strn-k-2]
                                prbmtrx[apix[i],j] = chisq
                                prbwave[apix[i],j] = patterns[m,7+i]
                                break
                if (prbmtrx[pa,strn-1] > chisq) or (prbmtrx[pa,strn-1] == -1.0):
                    for j in range(strn):
                        if prbmtrx[pa,j] == -1.0:
                            prbmtrx[pa,j] = chisq
                            prbwave[pa,j] = patterns[m,5]
                            break
                        elif chisq < prbmtrx[pa,j]:
                            # Shift everything down one
                            for k in range(0,strn-1-j):
                                prbmtrx[pa,strn-1-k] = prbmtrx[pa,strn-k-2]
                                prbwave[pa,strn-1-k] = prbwave[pa,strn-k-2]
                            prbmtrx[pa,j] = chisq
                            prbwave[pa,j] = patterns[m,5]
                            break
                if (prbmtrx[pb,strn-1] > chisq) or (prbmtrx[pb,strn-1] == -1.0):
                    for j in range(strn):
                        if prbmtrx[pb,j] == -1.0:
                            prbmtrx[pb,j] = chisq
                            prbwave[pb,j] = patterns[m,6]
                            break
                        elif chisq < prbmtrx[pb,j]:
                            # Shift everything down one
                            for k in range(0,strn-1-j):
                                prbmtrx[pb,strn-1-k] = prbmtrx[pb,strn-1-k-1]
                                prbwave[pb,strn-1-k] = prbwave[pb,strn-1-k-1]
                            prbmtrx[pb,j] = chisq
                            prbwave[pb,j] = patterns[m,6]
                            break
#			# Now assume decreasing pixels corresponds to increasing wavelength

    return prbmtrx, prbwave


#@cython.boundscheck(False)
def find_patterns_iter(np.ndarray[DTYPE_t, ndim=2] pixels not None,
                    np.ndarray[ITYPE_t, ndim=1] psort not None,
                    np.ndarray[DTYPE_t, ndim=2] patterns not None,
                    int searchcnt, double maxpix):
    # vamp,eamp,vpix,epix,vwid,ewid
    cdef int sz_m, sz_p
    cdef int m, p, q, pa, pb
    cdef int i, j, k, b, mincnt
    cdef double chisq, cnum, diff
    cdef double minval, wavecen, minpx, maxpx

    sz_m = patterns.shape[0]
    sz_p = pixels.shape[0]

    # the arrays that are used to estimate the pixels
    cdef np.ndarray[DTYPE_t, ndim=1] ppix = np.zeros(4, dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] apix

    # Arrays used to store the brightest pixels
    cdef np.ndarray[DTYPE_t, ndim=1] pixnew

    # During this calculation, store a probability matrix of values that
    #represent the most likely wavelengths
    cdef int strn = 5
    cdef np.ndarray[DTYPE_t, ndim=2] prbmtrx = np.zeros((sz_p,strn), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] prbwave = np.zeros((sz_p,strn), dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] prbdirc = np.zeros((sz_p,strn), dtype=ITYPE)
    for p in range(sz_p):
        for q in range(strn):
            prbmtrx[p,q] = -1.0

    # The arrays used to estimate the approximate central wavelength of the order
    cdef np.ndarray[DTYPE_t, ndim=1] pixfit = np.zeros((6), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] wavfit = np.zeros((6), dtype=DTYPE)

    cdef int searchnum = searchcnt
    if sz_p < searchcnt: searchnum = sz_p
    cdef int complete = 0
    while complete == 0:
        if searchnum == sz_p+5: break
        if sz_p < searchnum: searchnum = sz_p
        # Reset the probability matrix
        for p in range(searchnum):
            for q in range(strn):
                prbmtrx[p,q] = -1.0
        for b in range(2): # b=0 assumes increasing pixels = increasing wavelength (b=1 assumes the opposite)
            if b == 0:
                pixnew = nbrightest(pixels,psort,searchnum)
            else:
                reverse_array(pixnew, maxpix)
            # Go through all patterns
            for m in range(sz_m):
                for pa in range(searchnum-4):
                    for pb in range(pa+1,searchnum-3):
                        # First search assuming increasing pixels corresponds to increasing wavelength
                        diff = (pixnew[pb]-pixnew[pa])
                        for i in range(1,5):
                            if patterns[m,i] == 0.0: ppix[i-1] = 0.0
                            else: ppix[i-1] = pixnew[pa] + patterns[m,i]*diff
                        apix = find_nearest(pixnew, ppix, pb)
                        chisq = 0.0
                        cnum = 0.0
                        for i in range(4):
                            if apix[i] == -1: continue
                            chisq += (pixnew[apix[i]]-ppix[i])**2
                            cnum += 1.0
                        if cnum == 0.0: continue
                        chisq /= cnum
                        for i in range(4):
                            if apix[i] == -1: continue
                            if (prbmtrx[apix[i],strn-1] > chisq) or (prbmtrx[apix[i],strn-1] == -1.0):
                                for j in range(strn):
                                    if prbmtrx[apix[i],j] == -1.0:
                                        prbmtrx[apix[i],j] = chisq
                                        prbwave[apix[i],j] = patterns[m,7+i]
                                        prbdirc[apix[i],j] = b
                                        break
                                    elif chisq < prbmtrx[apix[i],j]:
                                        # Shift everything down one
                                        for k in range(0,strn-1-j):
                                            prbmtrx[apix[i],strn-1-k] = prbmtrx[apix[i],strn-k-2]
                                            prbwave[apix[i],strn-1-k] = prbwave[apix[i],strn-k-2]
                                            prbdirc[apix[i],strn-1-k] = prbdirc[apix[i],strn-k-2]
                                        prbmtrx[apix[i],j] = chisq
                                        prbwave[apix[i],j] = patterns[m,7+i]
                                        prbdirc[apix[i],j] = b
                                        break
                        if (prbmtrx[pa,strn-1] > chisq) or (prbmtrx[pa,strn-1] == -1.0):
                            for j in range(strn):
                                if prbmtrx[pa,j] == -1.0:
                                    prbmtrx[pa,j] = chisq
                                    prbwave[pa,j] = patterns[m,5]
                                    prbdirc[pa,j] = b
                                    break
                                elif chisq < prbmtrx[pa,j]:
                                    # Shift everything down one
                                    for k in range(0,strn-1-j):
                                        prbmtrx[pa,strn-1-k] = prbmtrx[pa,strn-k-2]
                                        prbwave[pa,strn-1-k] = prbwave[pa,strn-k-2]
                                        prbdirc[pa,strn-1-k] = prbdirc[pa,strn-k-2]
                                    prbmtrx[pa,j] = chisq
                                    prbwave[pa,j] = patterns[m,5]
                                    prbdirc[pa,j] = b
                                    break
                        if (prbmtrx[pb,strn-1] > chisq) or (prbmtrx[pb,strn-1] == -1.0):
                            for j in range(strn):
                                if prbmtrx[pb,j] == -1.0:
                                    prbmtrx[pb,j] = chisq
                                    prbwave[pb,j] = patterns[m,6]
                                    prbdirc[pb,j] = b
                                    break
                                elif chisq < prbmtrx[pb,j]:
                                    # Shift everything down one
                                    for k in range(0,strn-1-j):
                                        prbmtrx[pb,strn-1-k] = prbmtrx[pb,strn-k-2]
                                        prbwave[pb,strn-1-k] = prbwave[pb,strn-k-2]
                                        prbdirc[pb,strn-1-k] = prbdirc[pb,strn-k-2]
                                    prbmtrx[pb,j] = chisq
                                    prbwave[pb,j] = patterns[m,6]
                                    prbdirc[pb,j] = b
                                    break
            # Check if we've found an acceptable match (check it for both pixnew)
            minval = 3.0
            mincnt = 0
            minpx = 0.0
            maxpx = 0.0
            for pa in range(searchnum-3):
                if prbmtrx[pa,0] == minval and minval != 3.0:
                    pixfit[mincnt] = pixnew[pa]
                    wavfit[mincnt] = prbwave[pa,0]
                    mincnt += 1
                elif prbmtrx[pa,0] < minval:
                    minval = prbmtrx[pa,0]
                    pixfit[0] = pixnew[pa]
                    wavfit[0] = prbwave[pa,0]
                    mincnt = 1
                    minpx = pixnew[pa]
                    maxpx = pixnew[pa]
            if mincnt >= 5 and minval < 1.0:
                complete = 1+b
                wavecen = centre_wave(pixfit,wavfit,maxpix,mincnt)
                break
        if complete == 0:
            # Add a few new wavelengths to the next search iteration
            searchnum += 5
    # Return the central wavelength and pixel identifications
#	return pixfit, wavfit, wavecen, complete
    return prbwave, prbmtrx, prbdirc, pixnew

@cython.boundscheck(False)
def find_cont(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                np.ndarray[DTYPE_t, ndim=1] yarr not None,
                np.ndarray[ITYPE_t, ndim=1] marr not None,
                int nn):

    cdef int sz_x
    cdef int x, i, n, nmnx
    cdef double mnx

    sz_x = xarr.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] medarr = np.zeros(nn, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] xfit   = np.zeros(1+sz_x/nn, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] yfit   = np.zeros(1+sz_x/nn, dtype=DTYPE)

    n = 0
    i = 0
    nmnx = 0
    for x in range(0,sz_x):
        if marr[x] == 0:
            mnx += xarr[x]
            nmnx += 1
            medarr[i] = yarr[x]
        i += 1
        if i == nn:
            if nmnx != 0:
                xfit[n] = mnx/<double>(nmnx)
                yfit[n] = np.median(medarr[:nmnx])
                n += 1
                mnx = 0.0
                nmnx = 0
            i = 0
    if i != 0: # Append the last set of values
        if nmnx != 0:
            xfit[n] = mnx/<double>(nmnx)
            yfit[n] = np.median(medarr[:nmnx])
            n += 1
    return xfit[:n], yfit[:n]


@cython.boundscheck(False)
def find_wavecen(np.ndarray[DTYPE_t, ndim=1] waves not None,
                np.ndarray[DTYPE_t, ndim=1] orders not None):

    cdef int sz_o
    cdef int o, m, n, tnum
    cdef int i, j, k
    cdef double chisq, ewave

    sz_o = orders.shape[0]

    # During this calculation, store a probability matrix of values that
    #represent the most likely wavelengths
    cdef int strn = 5
    cdef np.ndarray[DTYPE_t, ndim=2] prbmtrx = np.zeros((sz_o,strn), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] prbwave = np.zeros((sz_o,strn), dtype=DTYPE)
    for o in range(sz_o):
        for n in range(strn):
            prbmtrx[o,n] = -1.0

    # Search from the leftmost order to rightmost order
    for o in range(0,sz_o-2):
        for m in range(o+1,sz_o-1):
            if orders[m] == orders[o]: continue
            elif orders[m] > orders[o]+1.0: break
            ewave = 2.0*waves[m]-waves[o]   # waves[o] + 2.0*(waves[m]-waves[o])
            # Find the nearest element to this prediction
            for n in range(m+1,sz_o):
                if orders[n] == orders[m]: continue
                elif orders[n] > orders[m]+1.0: break
                chisq = ewave-waves[n]
                if chisq < 0.0: chisq = chisq*-1.0
                for i in range(3):
                    if i == 0: tnum = o
                    elif i == 1: tnum = m
                    elif i == 2: tnum = n
                    if (prbmtrx[tnum,strn-1] > chisq) or (prbmtrx[tnum,strn-1] == -1.0):
                        for j in range(strn):
                            if prbmtrx[tnum,j] == -1.0:
                                prbmtrx[tnum,j] = chisq
                                prbwave[tnum,j] = waves[tnum]
                                break
                            elif chisq < prbmtrx[tnum,j]:
                                # Shift everything down one
                                for k in range(0,strn-1-j):
                                    prbmtrx[tnum,strn-1-k] = prbmtrx[tnum,strn-k-2]
                                    prbwave[tnum,strn-1-k] = prbwave[tnum,strn-k-2]
                                prbmtrx[tnum,j] = chisq
                                prbwave[tnum,j] = waves[tnum]
                                break
    # Search from the rightmost order to lefttmost order
    for o in range(1,sz_o-2):
        for m in range(o+1,sz_o-1):
            if orders[sz_o-m] == orders[sz_o-o]: continue
            elif orders[sz_o-m] < orders[sz_o-o]-1.0: break
            ewave = 2.0*waves[sz_o-m] - waves[sz_o-o]
            # Find the nearest element to this prediction
            chisq = -1.0
            for n in range(m+1,sz_o):
                if orders[sz_o-n] == orders[sz_o-m]: continue
                elif orders[sz_o-n] < orders[sz_o-m]-1.0: break
                chisq = ewave-waves[sz_o-n]
                if chisq < 0.0: chisq = chisq*-1.0
                for i in range(3):
                    if i == 0: tnum = sz_o-o
                    elif i == 1: tnum = sz_o-m
                    elif i == 2: tnum = sz_o-n
                    if (prbmtrx[tnum,strn-1] > chisq) or (prbmtrx[tnum,strn-1] == -1.0):
                        for j in range(strn):
                            if prbmtrx[tnum,j] == -1.0:
                                prbmtrx[tnum,j] = chisq
                                prbwave[tnum,j] = waves[tnum]
                                break
                            elif chisq < prbmtrx[tnum,j]:
                                # Shift everything down one
                                for k in range(0,strn-1-j):
                                    prbmtrx[tnum,strn-1-k] = prbmtrx[tnum,strn-k-2]
                                    prbwave[tnum,strn-1-k] = prbwave[tnum,strn-k-2]
                                prbmtrx[tnum,j] = chisq
                                prbwave[tnum,j] = waves[tnum]
                                break
    return prbwave, prbmtrx


@cython.boundscheck(False)
def fit_arcorder(np.ndarray[DTYPE_t, ndim=1] xarray not None,
                np.ndarray[DTYPE_t, ndim=1] yarray not None,
                np.ndarray[ITYPE_t, ndim=1] pixt not None,
                int fitp):
    cdef int p, pp, i, sz_p, sz_a
    cdef int pmin, pmax, flag
    cdef double tamp, tcen, twid

    sz_p = pixt.shape[0]
    sz_a = yarray.shape[0]

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


@cython.boundscheck(False)
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

    if coeff[0] == 0.0:
        return 0.0, 0.0, 0.0, 0
    else:
        wid = csqrt(-0.5/coeff[0])
        cen = cent-0.5*coeff[1]/coeff[0]
        amp = cexp(coeff[2] - 0.25*coeff[1]**2/coeff[0])
        return amp, cen, wid, 1


@cython.boundscheck(False)
def fit_gauss_old(np.ndarray[DTYPE_t, ndim=1] xarray not None,
                np.ndarray[DTYPE_t, ndim=1] yarray not None,
                np.ndarray[ITYPE_t, ndim=1] pixt not None,
                int fitp):
    """
    This routine uses a maximum likelihood method, but doesn't work very well.
    """
    cdef int p, pp, y, sz_p, sz_a
    cdef int pmin, pmax
    cdef double mval, wght

    sz_p = pixt.shape[0]
    sz_a = yarray.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] cent = np.zeros(sz_a, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] widt = np.zeros(sz_a, dtype=DTYPE)

    for p in range(sz_p):
        pmin = pixt[p]-(fitp-1)/2
        pmax = pixt[p]-(fitp-1)/2 + fitp
        if pmin < 0: pmin=0
        if pmax > sz_a: pmax = sz_a
        if pmin == pmax: continue
        # Find centre of this arc line
        mval = 0.0
        wght = 0.0
        # Calculate the MLE mean value
        for y in range(pmin,pmax):
            mval += <double>(y)*yarray[y]
            wght += yarray[y]
        cent[pp] = mval/wght
        # Calculate the MLE sigma value
        mval = 0.0
        for y in range(pmin,pmax):
            mval += yarray[y]*(xarray[y]-cent[pp])**2
        widt[pp] = csqrt(mval/wght)
        pp += 1
    return cent, widt, pp


@cython.boundscheck(False)
def func_pcsurf(np.ndarray[DTYPE_t, ndim=1] x not None,
                np.ndarray[DTYPE_t, ndim=1] ord not None,
                np.ndarray[DTYPE_t, ndim=1] p not None,
                np.ndarray[ITYPE_t, ndim=1] ford not None,
                np.ndarray[ITYPE_t, ndim=1] idx not None):
    """
    Returns a surface fit to polynomial coefficients
    """
    cdef int i, j, k, sz_x, sz_f
    cdef double pc

    sz_x = x.shape[0]
    sz_f = ford.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] yfit = np.zeros(sz_x, dtype=DTYPE)

    for i in range(0,sz_x):
        for j in range(0,sz_f):
            pc = 0.0
            for k in range(0,ford[j]+1):
                pc += p[idx[j]+k] * cpow(ord[i],<double>(k))
            yfit[i] +=  pc * cpow(x[i],<double>(j))
    return yfit

@cython.boundscheck(False)
def func_pcsurf_ret(np.ndarray[DTYPE_t, ndim=1] x not None,
                    np.ndarray[DTYPE_t, ndim=1] ord not None,
                    np.ndarray[DTYPE_t, ndim=1] p not None,
                    np.ndarray[ITYPE_t, ndim=1] ford not None,
                    np.ndarray[ITYPE_t, ndim=1] idx not None,
                    np.ndarray[DTYPE_t, ndim=1] yfit not None):
    """
    Returns a surface fit to polynomial coefficients
    """
    cdef int i, j, k, sz_x, sz_f
    cdef double pc

    sz_x = x.shape[0]
    sz_f = ford.shape[0]

    for i in range(0,sz_x):
        for j in range(0,sz_f):
            pc = 0.0
            for k in range(0,ford[j]+1):
                pc += p[idx[j]+k] * cpow(ord[i],<double>(k))
            yfit[i] +=  pc * cpow(x[i],<double>(j))
    return yfit


#######
#  G  #
#######

@cython.boundscheck(False)
def get_model(np.ndarray[DTYPE_t, ndim=1] arclines not None,
                double sigma):
    cdef int sz_a, sz_p
    cdef int i, j
    cdef double sz_w, sz_s

    sz_a = arclines.shape[0]
    sz_w = arclines[sz_a-1]-arclines[0]
    sz_p = 3*<int>(1.02*sz_w/sigma)
    sz_s = (<double>(sz_p)-(3.0*sz_w/sigma))/2.0
    cdef np.ndarray[DTYPE_t, ndim=1] xfit = np.zeros(sz_p, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] yfit = np.zeros(sz_p, dtype=DTYPE)

    cdef double wstart = arclines[0] - sz_s*sz_w/<double>(sz_a)
    cdef double wdiff  = (2.0*(arclines[0]-wstart)+sz_w)/<double>(sz_p-1)
    for i in range(sz_p):
        print i, "/", sz_p-1
        xfit[i] = wstart + <double>(i)*wdiff
        for j in range(sz_a):
            yfit[i] += cexp( -((xfit[i] - arclines[j])/sigma)**2 )
    return yfit

@cython.boundscheck(False)
def get_signal(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                double sigma):
    cdef int sz_p, sz_w
    cdef int i, j
    cdef double siz, sz_s

    sz_w = pixels.shape[0]
    siz = pixels[sz_w-1]-pixels[0]
    if siz < 0.0: siz *= -1.0
    sz_p = 10*<int>(1.02*siz/sigma)
    sz_s = (<double>(sz_p)-(10.0*siz/sigma))/2.0
    cdef np.ndarray[DTYPE_t, ndim=1] xfit = np.zeros(sz_p, dtype=DTYPE)
    for i in range(sz_p):
        print i, "/", sz_p-1
        for j in range(sz_w):
            xfit[i] += cexp( -((<double>(i)/10.0 - pixels[j] - sz_s)/sigma)**2 )
    return xfit

#######
#  H  #
#######

#######
#  I  #
#######

@cython.boundscheck(False)
def identify(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                np.ndarray[DTYPE_t, ndim=1] arclist not None,
                double srcrng):
    print "I think this is a deprecated function..."
    cdef int sz_p, sz_w
    cdef int p, w, wa, wb, wc
    cdef int i, j, k
    cdef double chisq, grad, wavcut
    cdef int last = 0

    sz_p = pixels.shape[0]
    sz_w = arclist.shape[0]

    # idlist is the list returned with the wavelength ids
    cdef np.ndarray[DTYPE_t, ndim=1] idlist = np.zeros(sz_p, dtype=DTYPE)

    # the arrays that are used in polyfit
    cdef np.ndarray[DTYPE_t, ndim=1] xfit = np.zeros(4, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] yfit = np.zeros(4, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(2, dtype=DTYPE)

    # During this calculation, store a probability matrix of values that
    #represent the most likely wavelengths
    cdef np.ndarray[DTYPE_t, ndim=2] prbmtrx = np.zeros((sz_p,5), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] prbwave = np.zeros((sz_p,5), dtype=DTYPE)
    for p in range(sz_p):
        prbmtrx[p,0] = -1.0
        prbmtrx[p,1] = -1.0
        prbmtrx[p,2] = -1.0
        prbmtrx[p,3] = -1.0
        prbmtrx[p,4] = -1.0

    # Generate a mask of the arc lines currently being searched
    cdef np.ndarray[ITYPE_t, ndim=1] mask = np.zeros(sz_w, dtype=ITYPE)
    for p in range(sz_p-3):
        print p, "/", sz_p-4
        idlist[p] = -1.0
        for w in range(sz_w-3):
            # Update the mask of arc lines currently being inspected
            if w == 0: mask_generate(arclist,mask,srcrng,last)
            else: mask_update(arclist,mask,w,last,srcrng)
            for wa in range(w+1,sz_w-2):
                if mask[wa+1] == 0: break
                grad = (arclist[wa]-arclist[w])/(pixels[p+1]-pixels[p])
                for wb in range(wa+1,sz_w-1):
                    if mask[wb] == 0: break
                    if 1.02*(arclist[w] + grad*(pixels[p+2]-pixels[p])) < arclist[wb]: break
                    for wc in range(wb+1,sz_w):
                        if mask[wb] == 0: break
                        if 1.02*(arclist[w] + grad*(pixels[p+3]-pixels[p])) < arclist[wc]: break
                        xfit[0] = pixels[p]
                        xfit[1] = pixels[p+1]
                        xfit[2] = pixels[p+2]
                        xfit[3] = pixels[p+3]
                        yfit[0] = arclist[w]
                        yfit[1] = arclist[wa]
                        yfit[2] = arclist[wb]
                        yfit[3] = arclist[wc]
                        chisq = polyfit(xfit,yfit,2,coeffs)  # i.e. y=mx+c
                        for k in range(4):
                            if (prbmtrx[p+k,4] > chisq) or (prbmtrx[p+k,4] == -1.0):
                                for i in range(5):
#									print p, k, i
#									print prbmtrx[p+k,:]
                                    if prbmtrx[p+k,i] == -1.0:
                                        prbmtrx[p+k,i] = chisq
                                        if k == 0: prbwave[p+k,i] = arclist[w]
                                        elif k == 1: prbwave[p+k,i] = arclist[wa]
                                        elif k == 2: prbwave[p+k,i] = arclist[wb]
                                        elif k == 3: prbwave[p+k,i] = arclist[wc]
                                        break
                                    elif chisq < prbmtrx[p+k,i]:
                                        # Shift everything down one
                                        for j in range(0,4-i):
                                            prbmtrx[p+k,4-j] = prbmtrx[p+k,4-j-1]
                                            prbwave[p+k,4-j] = prbwave[p+k,4-j-1]
                                        prbmtrx[p+k,i] = chisq
                                        if k == 0: prbwave[p+k,i] = arclist[w]
                                        elif k == 1: prbwave[p+k,i] = arclist[wa]
                                        elif k == 2: prbwave[p+k,i] = arclist[wb]
                                        elif k == 3: prbwave[p+k,i] = arclist[wc]
#										print prbmtrx[p+k,:]
#										print "-------------------"
                                        break
    return prbmtrx, prbwave
    # Flatten the vectors and set up the arrays for fitting
#	cdef np.ndarray[ITYPE_t, ndim=1] finmask = np.ones(sz_p*5, dtype=ITYPE)
#	cdef np.ndarray[DTYPE_t, ndim=1] wgtflat = np.zeros(sz_p*5, dtype=DTYPE)
#	cdef np.ndarray[DTYPE_t, ndim=1] wavflat = np.zeros(sz_p*5, dtype=DTYPE)
#	cdef np.ndarray[DTYPE_t, ndim=1] pixflat = np.zeros(sz_p*5, dtype=DTYPE)
#	for p in range(sz_p):
#		for k in range(5):
#			wgtflat[5*p+k] = 1.0/prbmtrx[p,k]
#			wavflat[5*p+k] = prbwave[p,k]
#			pixflat[5*p+k] = pixels[p]
#
#	# Make a first cut on what the likely true identifications are:
#	# More than likely, there should be a large number of the "true"
#	# identifications, and a fairly uniform random set. Try to remove
#	# the most uncertain of these.
#	cdef double medval, sigval
#	medval, sigval = medianmad(wavflat)
#	for p in range(5*sz_p):
#		if pixflat[p] < medval - 7*sigval or pixflat[p] > medval + 7.0*sigval:
#			finmask[p] = 0
#
#	You need to find a way to change the size of the arrays that are sent to polyfit_weighted
#
#	polyfit_weighted
#	return prbmtrx, prbwave


@cython.boundscheck(False)
def identifytwo(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                np.ndarray[DTYPE_t, ndim=1] arclist not None,
                double srcrng):
    print "I think this is a deprecated function..."
    cdef int sz_p, sz_w
    cdef int p, pa, w, wa
    cdef int i, j
    cdef double chisq, grad

    sz_p = pixels.shape[0]
    sz_w = arclist.shape[0]

    # idlist is the list returned with the wavelength ids
    cdef np.ndarray[DTYPE_t, ndim=1] idlist = np.zeros(sz_p, dtype=DTYPE)

    # the arrays that are used to estimate wavelengths
    cdef np.ndarray[DTYPE_t, ndim=1] wpix = np.zeros(sz_p, dtype=DTYPE)
    cdef double epix

    # During this calculation, store a probability matrix of values that
    #represent the most likely wavelengths
    cdef np.ndarray[DTYPE_t, ndim=2] prbmtrx = np.zeros((sz_p,5), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] prbwave = np.zeros((sz_p,5), dtype=DTYPE)
    for p in range(sz_p):
        prbmtrx[p,0] = -1.0
        prbmtrx[p,1] = -1.0
        prbmtrx[p,2] = -1.0
        prbmtrx[p,3] = -1.0
        prbmtrx[p,4] = -1.0

    cdef int inthere
    cdef double unq
    for p in range(sz_p-1):
        print p, "/", sz_p-2
        idlist[p] = -1.0
        if pixels[p+1] - pixels[p] < 40.0: continue
        for w in range(sz_w-1):
            # First search assuming increasing pixels corresponds to increasing wavelength
            for wa in range(w+1,sz_w):
                if arclist[wa]-arclist[w] > srcrng: break
                if arclist[wa]-arclist[w] < 1.0: break
                grad = (arclist[wa]-arclist[w])/(pixels[p+1]-pixels[p])
                # Assuming these two pixels are correctly identified,
                # calculate the estimated wavelengths for all pixels,
                # and also find the closest wavelength in arclist.
                chisq = 0.0
                for pa in range(sz_p):
                    epix     = arclist[w] + grad*(pixels[pa]-pixels[p])
                    wpix[pa] = find_closest(arclist,epix,w)
                    chisq += ((epix-wpix[pa])/wpix[pa])**2
                    inthere = 0
                    if pa == 0: unq = 1.0
                    else:
                        for i in range(pa):
                            if wpix[i] == wpix[pa]: inthere=1
                        if inthere == 0: unq += 1.0
                if unq < sz_p/2: continue
                for pa in range(sz_p):
                    if (prbmtrx[pa,4] > chisq) or (prbmtrx[pa,4] == -1.0):
                        for i in range(5):
                            if prbmtrx[pa,i] == -1.0:
                                prbmtrx[pa,i] = chisq
                                prbwave[pa,i] = wpix[pa]
                                break
                            elif chisq < prbmtrx[pa,i]:
                                # Shift everything down one
                                for j in range(0,4-i):
                                    prbmtrx[pa,4-j] = prbmtrx[pa,4-j-1]
                                    prbwave[pa,4-j] = prbwave[pa,4-j-1]
                                prbmtrx[pa,i] = chisq
                                prbwave[pa,i] = wpix[pa]
                                break
#			# Now assume decreasing pixels corresponds to increasing wavelength

    return prbmtrx, prbwave



@cython.boundscheck(False)
def identifythree(np.ndarray[DTYPE_t, ndim=1] pixels not None,
                np.ndarray[DTYPE_t, ndim=1] arclist not None,
                double srcrng):
    print "I think this is a deprecated function..."
    cdef int sz_p, sz_w
    cdef int p, pa, w, wa
    cdef int i, j, unq, inthere
    cdef double chisq, grad

    sz_p = pixels.shape[0]
    sz_w = arclist.shape[0]

    # idlist is the list returned with the wavelength ids
    cdef np.ndarray[DTYPE_t, ndim=1] idlist = np.zeros(sz_p, dtype=DTYPE)

    # the arrays that are used to estimate wavelengths
    cdef np.ndarray[DTYPE_t, ndim=1] wpix = np.zeros(sz_p, dtype=DTYPE)
    cdef double epix

    # During this calculation, store a probability matrix of values that
    #represent the most likely wavelengths
    cdef np.ndarray[DTYPE_t, ndim=2] prbmtrx = np.zeros((sz_p,5), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] prbwave = np.zeros((sz_p,5), dtype=DTYPE)
    for p in range(sz_p):
        prbmtrx[p,0] = -1.0
        prbmtrx[p,1] = -1.0
        prbmtrx[p,2] = -1.0
        prbmtrx[p,3] = -1.0
        prbmtrx[p,4] = -1.0

    cdef int gpts = 100
    cdef double gdiff = srcrng/(<double>(gpts)*(pixels[sz_p-1]-pixels[0]))
    for p in range(sz_p):
        print p, "/", sz_p-2
        idlist[p] = -1.0
        if pixels[p+1] - pixels[p] < 40.0: continue
        for w in range(sz_w-1):
            for wa in range(gpts):
                grad = <double>(wa+1)*gdiff
                chisq = 0.0
                unq = 0
                for pa in range(sz_p):
                    epix     = arclist[w] + grad*(pixels[pa]-pixels[p])
                    wpix[pa] = find_closest(arclist,epix,w)
                    chisq += ((epix-wpix[pa])/wpix[pa])**2
                    inthere = 0
                    if pa == 0: unq += 1
                    else:
                        for i in range(pa):
                            if wpix[i] == wpix[pa]: inthere=1
                        if inthere == 0: unq += 1
                if unq < sz_p/2: continue
#					print p, w, wa, arclist[w], arclist[wa], epix, wpix[pa]
                for pa in range(sz_p):
                    if (prbmtrx[pa,4] > chisq) or (prbmtrx[pa,4] == -1.0):
                        for i in range(5):
                            if prbmtrx[pa,i] == -1.0:
                                prbmtrx[pa,i] = chisq
                                prbwave[pa,i] = wpix[pa]
                                break
                            elif chisq < prbmtrx[pa,i]:
                                # Shift everything down one
                                for j in range(0,4-i):
                                    prbmtrx[pa,4-j] = prbmtrx[pa,4-j-1]
                                    prbwave[pa,4-j] = prbwave[pa,4-j-1]
                                prbmtrx[pa,i] = chisq
                                prbwave[pa,i] = wpix[pa]
                                break
#			# Now assume decreasing pixels corresponds to increasing wavelength
    return prbmtrx, prbwave


@cython.boundscheck(False)
def identify_lines(np.ndarray[DTYPE_t, ndim=1] pixls not None,
                    np.ndarray[DTYPE_t, ndim=1] wvguess not None,
                    np.ndarray[DTYPE_t, ndim=1] waves not None):

    cdef int sz_g, sz_w
    cdef int w, g, ibst
    cdef double wbst, tst

    sz_w = waves.shape[0]
    sz_g = wvguess.shape[0]

    # Store the most likely wavelengths
    cdef np.ndarray[DTYPE_t, ndim=1] prbpixl = np.zeros((sz_w), dtype=DTYPE)

    for w in range(0,sz_w):
        # Associate this arc line with a nearby one from the list
        wbst = wvguess[0]-waves[w]
        if wbst < 0.0: wbst *= -1.0
        ibst = 0
        for g in range(1,sz_g):
            tst = wvguess[g]-waves[w]
            if tst < 0.0: tst *= -1.0
            if tst < wbst:
                wbst = tst
                ibst = g
        prbpixl[w] = pixls[ibst]
    # Return the best pixel identifications
    return prbpixl


#######
#  J  #
#######

#######
#  K  #
#######


#######
#  L  #
#######

@cython.boundscheck(False)
def linreg_test(np.ndarray[DTYPE_t, ndim=2] array not None,
                int edx, int edy, double threshold):
    cdef int sz_x
    cdef int x
    cdef double ax, ay
    cdef double eddist, tdist, std

    sz_x = array.shape[0]

    # Fit a polynomial of order 1 (i.e. mx+c) to the sz_x closest points to (edx,edy)
    cdef np.ndarray[DTYPE_t, ndim=2] coeffs = polyfit_xy(array,2)
    cdef np.ndarray[DTYPE_t, ndim=1] pxdist = np.zeros(sz_x, dtype=DTYPE)

    # Derive the perpendicular distance the test point is from the best-fit line
    ax = (edy - coeffs[0] + edx/coeffs[1])/(coeffs[1] + 1.0/coeffs[1])
    ay = coeffs[1]*ax + coeffs[0]
    eddist = csqrt( (edx-ax)**2 + (edy-ay)**2 )

    # Derive the perpendicular distance each acceptable point is from the best-fit line
    for x in range(sz_x):
        ax = (array[x,1] - coeffs[0] + array[x,0]/coeffs[1])/(coeffs[1] + 1.0/coeffs[1])
        ay = coeffs[1]*ax + coeffs[0]
        tdist = csqrt( (array[x,0]-ax)**2 + (array[x,1]-ay)**2 )
        if array[x,1] < coeffs[1]*array[x,0]+coeffs[0]:
            pxdist[x] = -tdist
        else:
            pxdist[x] = tdist

    # Find the median absolute deviation of the points about the best-fit line.
    std = mad(pxdist)

    # Determine if the point is deviant and return
    if eddist/std>=threshold:
        return 1 # It's a deviant point
    else:
        return 0


#######
#  M  #
#######

@cython.boundscheck(False)
def mask_update(np.ndarray[DTYPE_t, ndim=1] arclist not None,
                np.ndarray[ITYPE_t, ndim=1] mask not None,
                int start, int last, double srcrng):

    cdef int x
    cdef int sz_w = arclist.shape[0]
    mask[start-1] = 0
    if mask[sz_w-1] == 1:
        return
    else:
        for x in range(start+1,sz_w):
            if arclist[x]-arclist[start] > srcrng:
                last = x
                break
            mask[x] = 1
        return


@cython.boundscheck(False)
def mask_generate(np.ndarray[DTYPE_t, ndim=1] arclist not None,
                    np.ndarray[ITYPE_t, ndim=1] mask not None,
                    double srcrng, int last):
    cdef int x
    cdef int sz_w = arclist.shape[0]

    last = sz_w-1
    mask[sz_w-1] = 0
    mask[sz_w-2] = 0
    mask[sz_w-3] = 0
    for x in range(sz_w):
        if arclist[x]-arclist[0] > srcrng:
            last = x
            break
        mask[x] = 1
    return mask



@cython.boundscheck(False)
def mad(np.ndarray[DTYPE_t, ndim=1] madarr not None):
    cdef int sz_x
    cdef int x, j
    cdef double temp, madval

    sz_x = madarr.shape[0]

    # Now get the median of madarr
    for x in range(sz_x-1):
        for j in range(x+1,sz_x):
            if madarr[j] < madarr[x]:
                temp = madarr[x]
                madarr[x] = madarr[j]
                madarr[j] = temp
    # Find the median value (and multiply by 1.4826 to get an estimation of the standard deviation)
    if sz_x%2==0:
        madval = 1.4826*0.5*(madarr[sz_x/2] + madarr[sz_x/2 - 1])
    else:
        madval = 1.4826*madarr[(sz_x-1)/2]
    # Return the median absolute deviation
    return madval

@cython.boundscheck(False)
def median(np.ndarray[DTYPE_t, ndim=1] array not None):
    cdef int sz_x
    cdef int x, j
    cdef double temp

    sz_x = array.shape[0]

    for x in range(sz_x-1):
        for j in range(x+1,sz_x):
            if array[j] < array[x]:
                temp = array[x]
                array[x] = array[j]
                array[j] = temp
    # Find and return the median value
    if sz_x%2==0:
        return 0.5*(array[sz_x/2] + array[sz_x/2 - 1])
    else:
        return array[(sz_x-1)/2]


@cython.boundscheck(False)
def median_mask(np.ndarray[DTYPE_t, ndim=1] array not None,
                double maskval):
    cdef int sz_x
    cdef int x, j, nnm
    cdef double temp

    sz_x = array.shape[0]

    nnm = 0
    for x in range(sz_x):
        if array[x] != maskval: nnm += 1

    cdef np.ndarray[DTYPE_t, ndim=1] arrmod = np.zeros(nnm, dtype=DTYPE)
    nnm = 0
    for x in range(sz_x):
        if array[x] != maskval:
            arrmod[nnm] = array[x]
            nnm += 1

    for x in range(nnm-1):
        for j in range(x+1,nnm):
            if arrmod[j] < arrmod[x]:
                temp = arrmod[x]
                arrmod[x] = arrmod[j]
                arrmod[j] = temp
    # Find and return the median value
    if nnm%2==0:
        return 0.5*(arrmod[nnm/2] + array[nnm/2 - 1])
    else:
        return arrmod[(nnm-1)/2]


@cython.boundscheck(False)
def medianmad(np.ndarray[DTYPE_t, ndim=1] array not None):
    cdef int sz_x
    cdef int x, j
    cdef double temp
    cdef double medval, madval

    sz_x = array.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] madarr = np.zeros(sz_x, dtype=DTYPE)

    for x in range(sz_x-1):
        for j in range(x+1,sz_x):
            if array[j] < array[x]:
                temp = array[x]
                array[x] = array[j]
                array[j] = temp
    # Find the median value
    if sz_x%2==0:
        medval = 0.5*(array[sz_x/2] + array[sz_x/2 - 1])
    else:
        medval = array[(sz_x-1)/2]
    # Calculate the Median absolute deviation
    for x in range(0,sz_x):
        temp = array[x]-medval
        if temp < 0.0:
            madarr[x] = -temp
        else:
            madarr[x] = temp
    # Now get the median of madarr
    for x in range(sz_x-1):
        for j in range(x+1,sz_x):
            if madarr[j] < madarr[x]:
                temp = madarr[x]
                madarr[x] = madarr[j]
                madarr[j] = temp
    # Find the median value (and multiply by 1.4826 to get an estimate of the standard deviation)
    if sz_x%2==0:
        madval = 1.4826*0.5*(madarr[sz_x/2] + madarr[sz_x/2 - 1])
    else:
        madval = 1.4826*madarr[(sz_x-1)/2]
    # Return the median and median absolute deviation
    return medval, madval


@cython.boundscheck(False)
def medianmad_qs(np.ndarray[DTYPE_t, ndim=1] array not None):
    cdef int sz_x
    cdef int x, j
    cdef double temp
    cdef double medval, madval

    sz_x = array.shape[0]

    # Find the median value
    if sz_x%2==0:
        quicksort(array,sz_x/2)
        medval = 0.5*(array[sz_x/2] + array[sz_x/2 - 1])
    else:
        quicksort(array,(sz_x-1)/2)
        medval = array[(sz_x-1)/2]
    # Calculate the Median absolute deviation
    for x in range(0,sz_x):
        temp = array[x]-medval
        if temp < 0.0:
            array[x] = -temp
        else:
            array[x] = temp
    # Find the median value
    if sz_x%2==0:
        quicksort(array,sz_x/2)
        madval = 0.5*(array[sz_x/2] + array[sz_x/2 - 1])
    else:
        quicksort(array,(sz_x-1)/2)
        madval = array[(sz_x-1)/2]
    # Return the median and median absolute deviation
    return medval, madval

@cython.boundscheck(False)
def medianmad_dimen(np.ndarray[DTYPE_t, ndim=1] array not None, int di):
    cdef int sz_x
    cdef int x, j
    cdef double temp
    cdef double medval, madval

    sz_x = array.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] madarr = np.zeros(sz_x, dtype=DTYPE)

    for x in range(sz_x-1):
        for j in range(x+1,sz_x):
            if array[j,di] < array[x,di]:
                temp = array[x,di]
                array[x,di] = array[j,di]
                array[j,di] = temp
    # Find the median value
    if sz_x%2==0:
        medval = 0.5*(array[sz_x/2,di] + array[sz_x/2 - 1,di])
    else:
        medval = array[(sz_x-1)/2,di]
    # Calculate the Median absolute deviation
    for x in range(0,sz_x):
        temp = array[x,di]-medval
        if temp < 0.0:
            madarr[x] = -temp
        else:
            madarr[x] = temp
    # Now get the median of madarr
    for x in range(sz_x-1):
        for j in range(x+1,sz_x):
            if madarr[j] < madarr[x]:
                temp = madarr[x]
                madarr[x] = madarr[j]
                madarr[j] = temp
    # Find the median value (and multiply by 1.4826 to get an estimate of the standard deviation)
    if sz_x%2==0:
        madval = 1.4826*0.5*(madarr[sz_x/2] + madarr[sz_x/2 - 1])
    else:
        madval = 1.4826*madarr[(sz_x-1)/2]
    # Return the median and median absolute deviation
    return medval, madval

#######
#  N  #
#######

#@cython.boundscheck(False)
def nbrightest(np.ndarray[DTYPE_t, ndim=2] pixels not None,
                np.ndarray[ITYPE_t, ndim=1] psort not None,
                int searchcnt):
    # vamp,eamp,vpix,epix,vwid,ewid
    cdef int p, j
    cdef double temp

    # An array that contains the brightest pixels in order of increasing pixel location
    cdef np.ndarray[DTYPE_t, ndim=1] pixnbrt = np.zeros(searchcnt, dtype=DTYPE)

    for p in range(searchcnt):
        pixnbrt[p] = pixels[psort[p],2]

    # Now sort the array by pixel location
    for p in range(searchcnt-1):
        for j in range(p+1,searchcnt):
            if pixnbrt[j] < pixnbrt[p]:
                temp = pixnbrt[p]
                pixnbrt[p] = pixnbrt[j]
                pixnbrt[j] = temp
    # Return the sorted N brightest pixel array
    return pixnbrt

#######
#  O  #
#######

@cython.boundscheck(False)
def order_saturation(np.ndarray[ITYPE_t, ndim=2] satmask not None,
                    np.ndarray[ITYPE_t, ndim=2] ordcen not None,
                    np.ndarray[ITYPE_t, ndim=2] ordwid not None,
                    int dispdir):

    cdef int sz_x, sz_o, sz_y, x, y, o
    cdef int flag, xmin, xmax
    sz_x = satmask.shape[1-dispdir]
    sz_y = satmask.shape[dispdir]
    sz_o = ordcen.shape[1]

    cdef np.ndarray[ITYPE_t, ndim=2] ordsat = np.zeros((sz_y,sz_o), dtype=ITYPE)

    for o in range(sz_o):
        for y in range(sz_y):
            xmin = ordcen[y,o]-ordwid[y,o]
            xmax = ordcen[y,o]+ordwid[y,o]+1
            if xmin < 0: xmin = 0
            if xmax >= sz_x: xmax=sz_x
            for x in range(xmin,xmax):
                if dispdir == 0:
                    if satmask[y,x] == 1:
                        ordsat[y,o] = 1
                        break
                else:
                    if satmask[x,y] == 1:
                        ordsat[y,o] = 1
                        break
    return ordsat

#######
#  P  #
#######

@cython.boundscheck(False)
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
        for j in range(1,degree):
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

@cython.boundscheck(False)
def polyfit_integral(np.ndarray[DTYPE_t, ndim=1] x not None,
                    np.ndarray[DTYPE_t, ndim=1] y not None,
                    np.ndarray[DTYPE_t, ndim=1] dx not None,
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
        for j in range(0,degree):
            gsl_matrix_set(xmat, i, j, (cpow(x[i]+0.5*dx[i], j+1)-cpow(x[i]-0.5*dx[i], j+1))/(dx[i]*<double>(j+1)))
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

@cython.boundscheck(False)
def polyfit_weighted(np.ndarray[DTYPE_t, ndim=1] x not None,
        np.ndarray[DTYPE_t, ndim=1] y not None,
        np.ndarray[DTYPE_t, ndim=1] w not None,
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
    wvec = gsl_vector_alloc(sz_x)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)

    for i in range(0,sz_x):
#		gsl_matrix_set(xmat, i, 0, 1.0)
        for j in range(0,degree):
            gsl_matrix_set(xmat, i, j, cpow(x[i], j))
        gsl_vector_set(yvec, i, y[i])
        gsl_vector_set(wvec, i, w[i])

    # Fit with a polynomial
    ws = gsl_multifit_linear_alloc(sz_x, degree)
    gsl_multifit_wlinear(xmat, wvec, yvec, c, cov, &chisq, ws)

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


@cython.boundscheck(False)
def polyfit_xy(np.ndarray[DTYPE_t, ndim=2] xy not None,
        int degree):

    cdef int sz_x
    cdef int i, j
    cdef double chisq

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat
    cdef gsl_vector *yvec
    cdef gsl_vector *c

    sz_x = xy.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(degree, dtype=DTYPE)

    xmat = gsl_matrix_alloc(sz_x, degree)
    yvec = gsl_vector_alloc(sz_x)
    c    = gsl_vector_alloc(degree)
    cov  = gsl_matrix_alloc(degree, degree)

    for i in range(0,sz_x):
#		gsl_matrix_set(xmat, i, 0, 1.0)
        for j in range(0,degree):
            gsl_matrix_set(xmat, i, j, cpow(xy[i,0], j))
        gsl_vector_set(yvec, i, xy[i,0])

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
    return coeffs


#######
#  Q  #
#######

@cython.boundscheck(False)
def quicksort(np.ndarray[DTYPE_t, ndim=1] array not None,
                int k):

    cdef int i, j, l, m
    cdef double x, tmp
    l = 0
    m = array.shape[0] - 1
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

@cython.boundscheck(False)
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
def reverse_array(np.ndarray[DTYPE_t, ndim=1] array not None,
                double maxpix):
    cdef int p, r
    cdef double temp

    cdef int sz_p = array.shape[0]
    if sz_p%2==0:
        r = sz_p/2
    else:
        r = (sz_p-1)/2
        array[r] = maxpix-array[r]

    for p in range(r):
        temp = maxpix-array[p]
        array[p] = maxpix-array[sz_p-p-1]
        array[sz_p-p-1] = temp
    return



@cython.boundscheck(False)
def rofunc(np.ndarray[DTYPE_t, ndim=1] xarr not None,
            np.ndarray[DTYPE_t, ndim=1] yarr not None,
            double b):
    cdef int j
    cdef double a, abdev, d, dv
    cdef double sum = 0.0

    sz_a = xarr.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] array = np.zeros(sz_a, dtype=DTYPE)

    for j in range(sz_a):
        array[j] = yarr[j]-b*xarr[j]
    if ((sz_a & 1) == 1):
        a = select(array,(sz_a-1)>>1)
    else:
        j=(sz_a>>1)
        a = 0.5*(select(array,j-1)+select(array,j))
    abdev = 0.0
    for j in range(sz_a):
        d=yarr[j]-(b*xarr[j]+a)
        if d > 0.0:
            abdev += d
        else:
            abdev -= d
        if (yarr[j] != 0.0):
            if yarr[j] < 0.0:
                dv = -yarr[j]
            else:
                dv = yarr[j]
            d /= dv
        else:
            if d > 1.0E-15: # Machine Precision
                sum += xarr[j]
            elif d < -1.0E-15: # Machine Precision
                sum -= xarr[j]
    return sum, a, abdev


@cython.boundscheck(False)
def robust_regression(np.ndarray[DTYPE_t, ndim=1] xarr not None,
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



@cython.boundscheck(False)
def robust_gridreg(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                np.ndarray[DTYPE_t, ndim=1] yarr not None,
                int order):

    cdef int i, shft, flag, fin
    cdef double madn, mado
    cdef double chisq, tval, cm

    # Set the perturbation level
    cdef double perturb = 1.0E-3

    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(order+1, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coefft = np.zeros(order+1, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coeffn = np.zeros(order+1, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] scale  = np.ones(order+1, dtype=DTYPE)
#	cdef np.ndarray[DTYPE_t, ndim=1] pertv  = np.array([-1200.0,1000.0,-667.0,-333.0,-100.0,-80.0,-60.0,-40.0,-20.0,-10.0,-8.7,-6.3,-5.0,-4.0,-3.0,-2.0,-1.0,1.0,2.0,3.0,4.0,5.0,6.7,8.3,10.0,20.0,40.0,60.0,80.0,100.0,333.0,667.0,1000.0,1200.0], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] pertv  = np.array([-2.0,-1.0,1.0,2.0], dtype=DTYPE)

    cdef int sz_ck = pertv.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=1] indx   = np.zeros(order+1, dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] bindx  = np.zeros(order+1, dtype=ITYPE)

    # First perform an ordinary least-squares regression to get reasonable starting parameters
    chisq = polyfit(xarr,yarr,order+1,coeffs)

    # Setup the model, and calculate the absolute value of the difference between the model and data:
    cdef np.ndarray[DTYPE_t, ndim=1] mdiff
    mdiff = polydiff(xarr,yarr,coeffs)
    # Calculate the median of this (i.e. the MAD)
    mado = median(mdiff)

    # Now minimize the mad by an iterative grid search
    cdef int conv = 0
    cdef int updt = 0
    cdef int nfail = 0
    while conv == 0:
        print mado, scale, coeffs
        if mado == 0.0: break
        # Set the coefficients array
        if updt == 1:
            for i in range(order+1):
                coeffs[i] = coeffn[i]
            updt = 0
        for i in range(order+1):
            indx[i] = 0
        fin  = 0
        cm   = mado
        while fin == 0:
            flag=0
            for i in range(order+1):
                if coeffs[i] == 0.0:
                    if pertv[indx[i]] < 0.0:
                        coefft[i] = -1.0E-14
                    else:
                        coefft[i] = +1.0E-14
                else:
                    if scale[i]*pertv[indx[i]]*perturb == -1.0:
                        flag = 1
                        scale[i] *= 1.56754
                    coefft[i] = coeffs[i]*(1.0+scale[i]*pertv[indx[i]]*perturb)
            if flag == 1: continue
            mdiff = polydiff(xarr,yarr,coefft)
            madn = median(mdiff)
            if madn < cm:
                cm = madn
                for i in range(order+1):
                    bindx[i] = indx[i]
            # Shift index
            shft = 0
            for i in range(order+1):
                if indx[i] == sz_ck-1 and shft == 1:
                    indx[i] = 0
                elif indx[i] == sz_ck-1:
                    indx[i] = 0
                    shft = 1
                else:
                    indx[i] += 1
                    break
            # Test if the loop is finished
            shft = 0
            for i in range(order+1):
                if indx[i] != 0: shft = 1
            if shft == 0: break
        # Set the best coefficients
        if cm == mado:
            for i in range(order+1):
                scale[i] /= 3.0
            nfail += 1
            continue
        for i in range(order+1):
            if coeffs[i] == 0.0:
                if pertv[bindx[i]] < 0.0:
                    coeffn[i] = -1.0E-14
                else:
                    coeffn[i] = +1.0E-14
            else:
                tval = scale[i]*pertv[bindx[i]]*perturb
                if tval == 1.0:
                    coeffn[i] = coeffs[i]*(1.0+scale[i]*pertv[bindx[i]-1]*perturb)
                elif tval == -1.0:
                    coeffn[i] = coeffs[i]*(1.0+scale[i]*pertv[bindx[i]+1]*perturb)
                else:
                    coeffn[i] = coeffs[i]*(1.0+tval)
        # Test for convergence
        mdiff=polydiff(xarr,yarr,coeffn)
        cm = median(mdiff)
        if cm > mado:
            for i in range(order+1):
                scale[i] /= 3.0
            nfail += 1
        elif mado-cm<1.0E-4 and nfail >= 5: conv = 1
        elif mado-cm<1.0E-4:
            updt = 1
            nfail += 1
        elif cm < mado:
            updt = 1
            nfail = 0
        mado = cm
    return coeffn




@cython.boundscheck(False)
def robust_linreg(np.ndarray[DTYPE_t, ndim=1] xarr not None,
                np.ndarray[DTYPE_t, ndim=1] yarr not None):
    cdef int sz_a = xarr.shape[0]
    cdef double a, b
    cdef double abdev = 0.0

    cdef int i
    cdef double b1, b2, dec, f, f1, f2, sigb, temp, sv, bv
    cdef double sx = 0.0
    cdef double sy = 0.0
    cdef double sxy = 0.0
    cdef double sxx = 0.0
    cdef double chisq = 0.0

    # First perform an ordinary least-squares regression to get reasonable starting parameters
    for i in range(sz_a):
        sx += xarr[i]
        sy += yarr[i]
        sxy += xarr[i]*yarr[i]
        sxx += xarr[i]*xarr[i]
    dec = <double>(sz_a)*sxx - sx*sx
    a = (sxx*sy - sx*sxy)/dec
    b = (<double>(sz_a)*sxy-sx*sy)/dec
    for i in range(sz_a):
        temp = yarr[i]-(a+b*xarr[i])
        chisq += temp*temp
    sigb = csqrt(chisq/dec)
    b1 = b
    f1, a, abdev = rofunc(xarr,yarr,b1)
    if sigb > 0.0:
        if f1 > 0.0:
            sv = 3.0*sigb
        else:
            sv = -3.0*sigb
        b2 = b + sv
        f2, a, abdev = rofunc(xarr,yarr,b2)
        if b1 == b2:
            abdev /= <double>(sz_a)
        # Bracketing
        while f1*f2 > 0.0:
            b = b2 + 1.6*(b2-b1)
            b1 = b2
            f1 = f2
            b2 = b
            f2, a, abdev = rofunc(xarr,yarr,b2)
        # Bisection
        sigb=0.1*sigb
        bv = b2-b1
        if bv < 0.0: bv *= -1.0
        while bv > sigb:
            b = b1 + 0.5*(b2-b1)
            if (b == b1) or (b == b2): break
            f, a, abdev = rofunc(xarr,yarr,b)
            if f*f1 >= 0.0:
                f1 = f
                b1 = b
            else:
                f2 = f
                b2 = b
    abdev /= <double>(sz_a)
    return a, b, abdev


#######
#  S  #
#######

@cython.boundscheck(False)
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



@cython.boundscheck(False)
def select(np.ndarray[DTYPE_t, ndim=1] array not None,
            int k):
   
    cdef int i, ir, j, l, mid, it
    cdef int n = array.shape[0]
    cdef double x, tval
    l = 0
    ir = n - 1

#	cdef np.ndarray[DTYPE_t, ndim=1] array = np.zeros(n, dtype=DTYPE)
#	for i in range(n): array[i] = arr[i]

    while 1:
        if (ir <= l+1):
            if ir == l+1 and array[ir] < array[l]:
                tval = array[ir]
                array[ir] = array[l]
                array[l] = tval
            return array[k]
        else:
            mid = (l+ir)>>1
            tval = array[mid]
            array[mid] = array[l+1]
            array[l+1] = tval
            if array[l] > array[ir]:
                tval = array[ir]
                array[ir] = array[l]
                array[l] = tval
            if array[l+1] > array[ir]:
                tval = array[ir]
                array[ir] = array[l+1]
                array[l+1] = tval
            if array[l] > array[l+1]:
                tval = array[l+1]
                array[l+1] = array[l]
                array[l] = tval
            i = l+1
            j = ir
            x = array[l+1]
            while 1:
                it = 0
                while (it==0):
                    i += 1
                    if array[i] > x: it = 1
                it = 0
                while (it==0):
                    j -= 1
                    if array[j] < x: it = 1
                if j < i: break
                tval = array[j]
                array[j] = array[i]
                array[i] = tval
            array[l+1] = array[j]
            array[j]   = x
            if j >= k: ir = j-1
            if j <= k: l = i

@cython.boundscheck(False)
def strip(np.ndarray[DTYPE_t, ndim=1] array not None,
            double maskval):
   
    cdef int i, nl, nr, sz
    sz = array.shape[0]

    # Return the input if there is nothing to strip
    if array[0] != maskval and array[sz-1] != maskval: return array

    nl = 0
    i = 0
    while array[i] == maskval:
        nl += 1
        i += 1
    nr = 0
    i=sz-1
    while array[i] == maskval:
        nr += 1
        i -= 1

    # Fill in the array
    cdef np.ndarray[DTYPE_t, ndim=1] outarr = np.zeros(sz-nl-nr, dtype=DTYPE)
    for i in range(nl,sz-nr):
        outarr[i-nl] = array[i]
    return outarr


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


