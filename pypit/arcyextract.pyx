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


@cython.boundscheck(False)
def bkgrd_polyfit(np.ndarray[DTYPE_t, ndim=1] sky not None,
                    int npix, int nord):
    cdef int p, sz_p, n, c
    cdef int pmin, pmax
    cdef double chisq, smval

    sz_p = sky.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] xarr = np.arange((sz_p), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coeff = np.zeros((nord+1), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] skysmth = np.zeros((sz_p), dtype=DTYPE)
    n = (npix-1)/2

    for p in range(sz_p):
        pmin = p-n
        pmax = p+n+1
        if pmin < 0:
            pmin = 0
        if pmax >= sz_p:
            pmax = sz_p
        chisq = polyfit(xarr,sky,pmin,pmax,nord,xarr[p],coeff)
        skysmth[p] = coeff[0]
    return skysmth


@cython.boundscheck(False)
def clip_arr(np.ndarray[DTYPE_t, ndim=2] parr not None,
            double cval):
    cdef int p, pp, sz_p
    sz_p = parr.shape[0]
    pp=0
    for p in range(sz_p):
        if parr[p,0] == cval: break
        else: pp+=1
    cdef np.ndarray[DTYPE_t, ndim=2] pout = np.zeros((pp,2), dtype=DTYPE)
    for p in range(0,pp):
        pout[p,0] = parr[p,0]
        pout[p,1] = parr[p,1]
    return pout

@cython.boundscheck(False)
def cr_maskmedian(np.ndarray[DTYPE_t, ndim=2] frame not None,
                    double maskval, double sigmadet, double sigmarep,
                    np.ndarray[ITYPE_t, ndim=2] maskcr not None):
    cdef int r, sz_r, c, sz_c
    cdef int cnt, i, j
    cdef double temp

    sz_r = frame.shape[0]
    sz_c = frame.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] rarr = np.zeros((sz_r), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] carr = np.zeros((sz_c), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] rmed = np.zeros((sz_r,2), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] cmed = np.zeros((sz_c,2), dtype=DTYPE)
#	cdef np.ndarray[DTYPE_t, ndim=2] maskcr = np.ones((sz_r,sz_c), dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] stack = np.ones((sz_r*sz_c/10,2), dtype=ITYPE)
#	cdef np.ndarray[DTYPE_t, ndim=1] masksky = np.ones((sz_r), dtype=DTYPE)

    # Calculate the median value and MAD along a column
    for c in range(sz_c):
        cnt = 0
        for r in range(sz_r):
            temp = frame[r,c] - maskval
            if temp < 0.0: temp *= -1.0
            if temp > 0.1:
                rarr[cnt] = frame[r,c]
                cnt += 1
        if cnt == 0:
            cmed[c,0] = maskval
            cmed[c,1] = maskval
            continue
        # Sort the array
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if rarr[j] < rarr[i]:
                    temp = rarr[i]
                    rarr[i] = rarr[j]
                    rarr[j] = temp
        # Find the median value
        if cnt%2==0:
            cmed[c,0] = 0.5*(rarr[cnt/2] + rarr[cnt/2 - 1])
        else:
            cmed[c,0] = rarr[(cnt-1)/2]
        # Calculate the Median absolute deviation
        for i in range(0,cnt):
            temp = rarr[i]-cmed[c,0]
            if temp < 0.0:
                rarr[i] = -temp
            else:
                rarr[i] = temp
        # Now get the median of madarr
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if rarr[j] < rarr[i]:
                    temp = rarr[i]
                    rarr[i] = rarr[j]
                    rarr[j] = temp
        # Find the median value (and multiply by 1.4826 to get an estimate of the standard deviation)
        if cnt%2==0:
            cmed[c,1] = 1.4826*0.5*(rarr[cnt/2] + rarr[cnt/2 - 1])
        else:
            cmed[c,1] = 1.4826*rarr[(cnt-1)/2]

    # Calculate the median value along a row
    for r in range(sz_r):
        cnt = 0
        for c in range(sz_c):
            temp = frame[r,c] - maskval
            if temp < 0.0: temp *= -1.0
            if temp > 0.1:
                carr[cnt] = frame[r,c]
                cnt += 1
        if cnt == 0:
            rmed[r,0] = maskval
            rmed[r,1] = maskval
            continue
        # Sort the array
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if carr[j] < carr[i]:
                    temp = carr[i]
                    carr[i] = carr[j]
                    carr[j] = temp
        # Find the median value
        if cnt%2==0:
            rmed[r,0] = 0.5*(carr[cnt/2] + carr[cnt/2 - 1])
        else:
            rmed[r,0] = carr[(cnt-1)/2]
        # Calculate the Median absolute deviation
        for i in range(0,cnt):
            temp = carr[i]-rmed[r,0]
            if temp < 0.0:
                carr[i] = -temp
            else:
                carr[i] = temp
        # Now get the median of madarr
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if carr[j] < carr[i]:
                    temp = carr[i]
                    carr[i] = carr[j]
                    carr[j] = temp
        # Find the median value (and multiply by 1.4826 to get an estimate of the standard deviation)
        if cnt%2==0:
            rmed[r,1] = 1.4826*0.5*(carr[cnt/2] + carr[cnt/2 - 1])
        else:
            rmed[r,1] = 1.4826*carr[(cnt-1)/2]
    #return rmed, cmed
    # Loop through the entire array and identify the cosmic rays.
    # Replace the cosmic rays with a zero value and generate a mask.
    for r in range(0,sz_r):
        for c in range(0,sz_c):
            temp = frame[r,c]-maskval
            if temp < 0.0: temp *= -1.0
            if temp < 0.1: continue
            if frame[r,c] > (rmed[r,0] + sigmadet*rmed[r,1]) and frame[r,c] > (cmed[c,0] + sigmadet*cmed[c,1]):
                if (rmed[r,0] + sigmarep*rmed[r,1]) > (cmed[c,0] + sigmarep*cmed[c,1]):
                    temp = rmed[r,0] + sigmarep*rmed[r,1]
                else:
                    temp = cmed[c,0] + sigmarep*cmed[c,1]
                floodfill_fast(frame, maskcr, stack, r, c, sz_r, sz_c, temp, 0.0)
#				floodfill(frame, maskcr, r, c, sz_r, sz_c, temp, 0.0)
#			elif frame[r,c] > (cmed[c,0] + 0.5*(sigmarep+sigmadet)*cmed[c,1]):
#				masksky[r] = 0.0
    return# maskcr#, masksky


@cython.boundscheck(False)
def floodfill(np.ndarray[DTYPE_t, ndim=2] frame not None,
                np.ndarray[DTYPE_t, ndim=2] mask not None,
                int r, int c, int sz_r, int sz_c,
                double cond, double replace):
    # Invoke the same routine on all neighbouring cells recursively:
    if frame[r,c] >= cond:
        frame[r,c] = replace
        mask[r,c] = 0.0
        if r > 0:
            floodfill(frame,mask,r-1,c,sz_r,sz_c,cond,replace)
        if r < sz_r - 1:
            floodfill(frame,mask,r+1,c,sz_r,sz_c,cond,replace)
        if c > 0:
            floodfill(frame,mask,r,c-1,sz_r,sz_c,cond,replace)
        if c < sz_c - 1:
            floodfill(frame,mask,r,c+1,sz_r,sz_c,cond,replace)
    return

@cython.boundscheck(False)
def floodfill_fast(np.ndarray[DTYPE_t, ndim=2] frame not None,
                np.ndarray[ITYPE_t, ndim=2] mask not None,
                np.ndarray[ITYPE_t, ndim=2] stack not None,
                int r, int c, int sz_r, int sz_c,
                double cond, double replace):
    # assume surface is a 2D image and surface[x][y] is the color at x, y.
    cdef int n = 0
    cdef int sz_n = stack.shape[0]
    stack[n,0] = r
    stack[n,1] = c
    n += 1
    while n > 0:
        n -= 1
        r = stack[n,0]
        c = stack[n,1]
        if frame[r,c] < cond:
            continue
        if n >= sz_n:
            print "WARNING: Not all cosmic rays were masked - a stricter replacement"
            print "criteria is recommended so false positives are not identified."
            break
        frame[r,c] = replace
        mask[r,c] = 0
        # right
        stack[n,0] = r + 1
        stack[n,1] = c
        n += 1
        # left
        stack[n,0] = r - 1
        stack[n,1] = c
        n += 1
        # down
        stack[n,0] = r
        stack[n,1] = c + 1
        n += 1
        # up
        stack[n,0] = r
        stack[n,1] = c - 1
        n += 1
    return


@cython.boundscheck(False)
def extract_2d(np.ndarray[DTYPE_t, ndim=2] frame not None,
                np.ndarray[DTYPE_t, ndim=2] error not None,
                np.ndarray[DTYPE_t, ndim=2] bpix not None,
                np.ndarray[ITYPE_t, ndim=1] pixcen not None,
                np.ndarray[ITYPE_t, ndim=1] piycen not None,
                np.ndarray[DTYPE_t, ndim=3] pixloc not None,
                np.ndarray[DTYPE_t, ndim=1] tilts not None,
                np.ndarray[DTYPE_t, ndim=1] ordxcen not None,
                np.ndarray[DTYPE_t, ndim=1] ordycen not None,
                np.ndarray[DTYPE_t, ndim=1] ordwid not None,
                int ordwnum, int ordlen, int ordnum, int dispaxis):

    cdef int x, sz_x, sz_y, sz_o, w
    cdef double area, uarea, darea, totarea, totflux, toterror
    cdef double grad, intc, xstp, ystp, down
    cdef int pixs, pixll, pixl, pixu, pixd, i, j, donesome
#	cdef int dme = 0

    down = <double>(ordwnum)
    sz_o = ordxcen.shape[0]
    sz_x = frame.shape[dispaxis]
    sz_y = frame.shape[1-dispaxis]

    cdef np.ndarray[DTYPE_t, ndim=2] extspec = np.zeros((ordnum,ordwnum), dtype=DTYPE) # The extracted spectrum
    cdef np.ndarray[DTYPE_t, ndim=2] errspec = np.zeros((ordnum,ordwnum), dtype=DTYPE) # The extracted spectrum
    cdef np.ndarray[DTYPE_t, ndim=2] pixlarr = np.zeros((4,2), dtype=DTYPE) # The pixel coordinates
#	cdef np.ndarray[DTYPE_t, ndim=3] pixlarr = np.zeros((4,2,2000), dtype=DTYPE) # The pixel coordinates
    cdef np.ndarray[DTYPE_t, ndim=2] polyarr = np.zeros((4,2), dtype=DTYPE) # The polygon coordinates
#	cdef np.ndarray[DTYPE_t, ndim=2] polyout
    cdef int nsiz = 11
    cdef np.ndarray[DTYPE_t, ndim=2] outlist = np.zeros((nsiz,2), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] inlist = np.zeros((nsiz,2), dtype=DTYPE)

    cdef int rend = 0
    pixll = 0
    pixl = 0
    for x in range(sz_o): # For each row in the dispersion direction
        pixs = pixll
        for w in range(ordwnum): # and each pixel in the spatial direction
            # Set the polygon extraction array
            if rend == 1: break # Reached the end of the order, the remaining values are masked
            if x == 0:
                # Calculate the equation of a line to the negative of this pixel
                grad = tilts[x]
                xstp = (ordwid[x])/((1.0+grad**2)*down)
                ystp = (ordwid[x])*grad/((1.0+grad**2)*down)
                polyarr[0,0] = ordxcen[x] + xstp*<double>(w+0) - 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[0,1] = ordycen[x] + ystp*<double>(w+0) - 0.5*ordwid[x]*grad/(1.0+grad**2)
                polyarr[1,0] = ordxcen[x] + xstp*<double>(w+1) - 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[1,1] = ordycen[x] + ystp*<double>(w+1) - 0.5*ordwid[x]*grad/(1.0+grad**2)
                # Calculate the equation of a line to the positive of this pixel
                grad = 0.5*(tilts[x]+tilts[x+1])
                xstp = 0.5*(ordwid[x]+ordwid[x+1])/((1.0+grad**2)*down)
                ystp = 0.5*(ordwid[x]+ordwid[x+1])*grad/((1.0+grad**2)*down)
                polyarr[2,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + xstp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                polyarr[2,1] = 0.5*(ordycen[x]+ordycen[x+1]) + ystp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
                polyarr[3,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + xstp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                polyarr[3,1] = 0.5*(ordycen[x]+ordycen[x+1]) + ystp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
            elif x==sz_o-1:
                # Calculate the equation of a line to the negative of this pixel
                grad = 0.5*(tilts[x]+tilts[x-1])
                xstp = 0.5*(ordwid[x]+ordwid[x-1])/((1.0+grad**2)*down)
                ystp = 0.5*(ordwid[x]+ordwid[x-1])*grad/((1.0+grad**2)*down)
                polyarr[0,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + xstp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
                polyarr[0,1] = 0.5*(ordycen[x]+ordycen[x-1]) + ystp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
                polyarr[1,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + xstp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
                polyarr[1,1] = 0.5*(ordycen[x]+ordycen[x-1]) + ystp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
                # Calculate the equation of a line to the positive of this pixel
                grad = tilts[x]
                xstp = (ordwid[x])/((1.0+grad**2)*down)
                ystp = (ordwid[x])*grad/((1.0+grad**2)*down)
                polyarr[2,0] = ordxcen[x] + xstp*<double>(w+1) - 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[2,1] = ordycen[x] + ystp*<double>(w+1) - 0.5*ordwid[x]*grad/(1.0+grad**2)
                polyarr[3,0] = ordxcen[x] + xstp*<double>(w+0) - 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[3,1] = ordycen[x] + ystp*<double>(w+0) - 0.5*ordwid[x]*grad/(1.0+grad**2)
                rend = 1
            else:
                # Calculate the equation of a line to the negative of this pixel
                grad = 0.5*(tilts[x]+tilts[x-1])
                xstp = 0.5*(ordwid[x]+ordwid[x-1])/((1.0+grad**2)*down)
                ystp = 0.5*(ordwid[x]+ordwid[x-1])*grad/((1.0+grad**2)*down)
                polyarr[0,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + xstp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
                polyarr[0,1] = 0.5*(ordycen[x]+ordycen[x-1]) + ystp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
                polyarr[1,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + xstp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
                polyarr[1,1] = 0.5*(ordycen[x]+ordycen[x-1]) + ystp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
                # Calculate the equation of a line to the positive of this pixel
                if ordxcen[x+1] == -999999.9:
                    grad = tilts[x]
                    xstp = (ordwid[x])/((1.0+grad**2)*down)
                    ystp = (ordwid[x])*grad/((1.0+grad**2)*down)
                    polyarr[2,0] = ordxcen[x] + xstp*<double>(w+1) - 0.5*ordwid[x]/(1.0+grad**2)
                    polyarr[2,1] = ordycen[x] + ystp*<double>(w+1) - 0.5*ordwid[x]*grad/(1.0+grad**2)
                    polyarr[3,0] = ordxcen[x] + xstp*<double>(w+0) - 0.5*ordwid[x]/(1.0+grad**2)
                    polyarr[3,1] = ordycen[x] + ystp*<double>(w+0) - 0.5*ordwid[x]*grad/(1.0+grad**2)
                    rend = 1
                else:
                    grad = 0.5*(tilts[x]+tilts[x+1])
                    xstp = 0.5*(ordwid[x]+ordwid[x+1])/((1.0+grad**2)*down)
                    ystp = 0.5*(ordwid[x]+ordwid[x+1])*grad/((1.0+grad**2)*down)
                    polyarr[2,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + xstp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                    polyarr[2,1] = 0.5*(ordycen[x]+ordycen[x+1]) + ystp*<double>(w+1) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
                    polyarr[3,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + xstp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                    polyarr[3,1] = 0.5*(ordycen[x]+ordycen[x+1]) + ystp*<double>(w+0) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
            # Define the pixel array, integrate, and extract
            totarea = 0.0
            totflux = 0.0
            toterror = 0.0
            j=pixs
            donesome=0
            if dispaxis == 0:
                # yloc = slf._pixlocn[:,0,0]
                # xloc = slf._pixlocn[0,:,1]
                # ysiz = slf._pixlocn[:,0,2]
                # xsiz = slf._pixlocn[0,:,3]
                while True:
                    if j >= sz_x:
                        j = x+1
                        break
                    # Sum up over pixels area
                    uarea = 0.0
                    i=0
                    pixu = <int>(0.5+ordwid[x]/2.0)
                    while True:
                        if (pixcen[x]+i >= sz_y) or (pixcen[x]+i < 0):
                            if i>pixu+1:
#								print x, j, i, "BREAK 1"
                                break
                            i += 1
                            continue
                        if (bpix[j,pixcen[x]+i] == 1):
                            if i>pixu+1: break
                            i += 1
                            continue
                        pixlarr[0,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[0,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                        pixlarr[1,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[1,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                        pixlarr[2,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[2,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                        pixlarr[3,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[3,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[0,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[0,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[1,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[1,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[2,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[2,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[3,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[3,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                        # Calculate the overlapping polygon area
                        area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#						polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#						dme += 1
#						if dme >= 2000: return pixlarr, polyarr
                        #area = poly_area(polyout)
                        if area == np.nan:
                            print "AREA=NAN (up)"
                            print outlist
                            print inlist
                            print pixlarr
                        uarea += area
                        totflux += area*frame[j,pixcen[x]+i]
                        toterror += cpow(area*error[j,pixcen[x]+i],2.0)
                        if area == 0.0:
                            if i>pixu+1:
#								print x, j, i, "BREAK 2"
                                break
                        else:
                            if i>pixu: pixu=i
                        i += 1
                    # Sum up over pixels area
                    darea = 0.0
                    i=-1
                    pixd = <int>(0.5+ordwid[x]/2.0)
                    while True:
                        if (pixcen[x]+i >= sz_y) or (pixcen[x]+i < 0):
                            if -i>pixd+1:
#								print x, j, i, "BREAK 3"
                                break
                            i -= 1
                            continue
                        if (bpix[j,pixcen[x]+i] == 1):
                            if -i>pixd+1: break
                            i -= 1
                            continue
                        pixlarr[0,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[0,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                        pixlarr[1,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[1,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                        pixlarr[2,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[2,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                        pixlarr[3,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                        pixlarr[3,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[0,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[0,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[1,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[1,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[2,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[2,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#						pixlarr[3,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#						pixlarr[3,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                        # Calculate the overlapping polygon area
                        area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#						polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#						dme += 1
#						if dme >= 2000: return pixlarr, polyarr
                        #area = poly_area(polyout)
                        if area == np.nan:
                            print "AREA=NAN (down)"
                            print outlist
                            print inlist
                            print pixlarr
                        darea += area
                        totflux += area*frame[j,pixcen[x]+i]
                        toterror += cpow(area*error[j,pixcen[x]+i],2.0)
                        if area == 0.0:
                            if -i>pixd+1:
#								print x, j, i, "BREAK 4"
                                break
                        else:
                            if -i>pixd: pixd=-i
                        i -= 1
                    # Break if we've already found some overlap
                    if uarea+darea == 0.0:
                        if donesome == 1:
#							print x, j, i, "BREAK 5"
                            break
                    else:
                        donesome = 1
                    # If we haven't found overlap yet, then we must be off the chip
                    if j-pixll>=5 and donesome==0:
                        j = x+1
                        break
                    j += 1
                    totarea += (uarea+darea)
                pixll = pixl
                pixl  = j-1
            else:
                # yloc = slf._pixlocn[0,:,0]
                # xloc = slf._pixlocn[:,0,1]
                # ysiz = slf._pixlocn[0,:,2]
                # xsiz = slf._pixlocn[:,0,3]
                while True:
                    if j >= sz_y:
                        j = x+1
                        break
                    # Sum up over pixels area
                    uarea = 0.0
                    i=0
                    pixu = <int>(0.5+ordwid[x]/2.0)
                    while True:
                        if (pixcen[x]+i >= sz_x) or (pixcen[x]+i < 0):
                            if i>pixu+1: break
                            i += 1
                            continue
                        if (bpix[pixcen[x]+i,j] == 1):
                            if i>pixu+1: break
                            i += 1
                            continue
                        pixlarr[0,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[0,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                        pixlarr[1,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[1,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                        pixlarr[2,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[2,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                        pixlarr[3,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[3,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[0,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[0,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[1,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[1,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[2,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[2,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[3,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[3,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                        # Calculate the overlapping polygon area
                        area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#						polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#						dme += 1
#						if dme >= 2000: return pixlarr, polyarr
                        #area = poly_area(polyout)
                        if area == np.nan:
                            print "AREA=NAN"
                            print outlist
                            print inlist
                            print pixlarr
                        uarea += area
                        totflux += area*frame[pixcen[x]+i,j]
                        toterror += cpow(area*error[pixcen[x]+i,j],2.0)
                        if area == 0.0:
                            if i>pixu+1: break
                        else:
                            if i>pixu: pixu=i
                        i += 1
                    # Sum up over pixels area
                    darea = 0.0
                    i=-1
                    pixd = <int>(0.5+ordwid[x]/2.0)
                    while True:
                        if (pixcen[x]+i >= sz_x) or (pixcen[x]+i < 0):
                            if -i>pixd+1: break
                            i -= 1
                            continue
                        if (bpix[pixcen[x]+i,j] == 1):
                            if -i>pixd+1: break
                            i -= 1
                            continue
                        pixlarr[0,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[0,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                        pixlarr[1,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[1,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                        pixlarr[2,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[2,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                        pixlarr[3,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                        pixlarr[3,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[0,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[0,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[1,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[1,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[2,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[2,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#						pixlarr[3,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#						pixlarr[3,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                        # Calculate the overlapping polygon area
                        area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#						polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#						dme += 1
#						if dme >= 2000: return pixlarr, polyarr
                        #area = poly_area(polyout)
                        if area == np.nan:
                            print "AREA=NAN"
                            print outlist
                            print inlist
                            print pixlarr
                        darea += area
                        totflux += area*frame[pixcen[x]+i,j]
                        toterror += cpow(area*error[pixcen[x]+i,j],2.0)
                        if area == 0.0:
                            if -i>pixd+1: break
                        else:
                            if -i>pixd: pixd=-i
                        i -= 1
                    # Break if we've already found some overlap
                    if uarea+darea == 0.0:
                        if donesome == 1: break
                    else:
                        donesome = 1
                    # If we haven't found overlap yet, then we must be off the chip
                    if j-pixll>=5 and donesome==0:
                        j = x+1
                        break
                    j += 1
                    totarea += (uarea+darea)
                pixll = pixl
                pixl  = j-1
            if totarea == np.nan:
                print x, totarea, totflux, uarea, darea, donesome, darea==np.nan
                print x, polyarr
            if totarea != 0.0:
                extspec[x,w] = totflux/totarea
                errspec[x,w] = csqrt(toterror)/totarea
            else:
                extspec[x,w] = -999999.9
                errspec[x,w] = -999999.9
    return extspec, errspec


@cython.boundscheck(False)
def extract_mean(np.ndarray[DTYPE_t, ndim=2] frame not None,
                np.ndarray[DTYPE_t, ndim=2] bpix not None,
                np.ndarray[ITYPE_t, ndim=1] pixcen not None,
                np.ndarray[ITYPE_t, ndim=1] piycen not None,
                np.ndarray[DTYPE_t, ndim=3] pixloc not None,
                np.ndarray[DTYPE_t, ndim=1] tilts not None,
                np.ndarray[DTYPE_t, ndim=1] ordxcen not None,
                np.ndarray[DTYPE_t, ndim=1] ordycen not None,
                np.ndarray[DTYPE_t, ndim=1] ordwid not None,
                int ordlen, int ordnum, int dispaxis):
    """
    Sometimes there is a zero (or a few zeros) padded at the end of
    extspec. This might be because of a difference between ordnum and
    ordxcen[x+1] == -999999.9: rend=1
    ---> I might have solved this issue by defining ext1d to have a
    series of -999999.9 values rather than 0's in arextract.py, for
    when the array doesn't get completely filled in.
    """

    cdef int x, sz_x, sz_y, sz_o
    cdef double area, uarea, darea, totarea, totflux
    cdef double grad, intc
    cdef int pixll, pixl, pixu, pixd, i, j, donesome
#	cdef int dme = 0

    sz_o = ordxcen.shape[0]
    sz_x = frame.shape[dispaxis]
    sz_y = frame.shape[1-dispaxis]

    cdef np.ndarray[DTYPE_t, ndim=1] extspec = np.zeros(ordnum, dtype=DTYPE) # The extracted spectrum
    cdef np.ndarray[DTYPE_t, ndim=2] pixlarr = np.zeros((4,2), dtype=DTYPE) # The pixel coordinates
#	cdef np.ndarray[DTYPE_t, ndim=3] pixlarr = np.zeros((4,2,2000), dtype=DTYPE) # The pixel coordinates
    cdef np.ndarray[DTYPE_t, ndim=2] polyarr = np.zeros((4,2), dtype=DTYPE) # The polygon coordinates
#	cdef np.ndarray[DTYPE_t, ndim=2] polyout
    cdef int nsiz = 11
    cdef np.ndarray[DTYPE_t, ndim=2] outlist = np.zeros((nsiz,2), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] inlist = np.zeros((nsiz,2), dtype=DTYPE)

    cdef int rend = 0
    pixll = 0
    pixl = 0
    for x in range(sz_o):
        # Set the polygon extraction array
        if rend == 1: break # Reached the end of the order, the remaining values are masked
        if x == 0:
            # Calculate the equation of a line to the negative of this pixel
            grad = tilts[x]
            polyarr[0,0] = ordxcen[x] - 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[0,1] = ordycen[x] - 0.5*ordwid[x]*grad/(1.0+grad**2)
            polyarr[1,0] = ordxcen[x] + 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[1,1] = ordycen[x] + 0.5*ordwid[x]*grad/(1.0+grad**2)
            # Calculate the equation of a line to the positive of this pixel
            grad = 0.5*(tilts[x]+tilts[x+1])
            polyarr[2,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
            polyarr[2,1] = 0.5*(ordycen[x]+ordycen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
            polyarr[3,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
            polyarr[3,1] = 0.5*(ordycen[x]+ordycen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
        elif x==sz_o-1:
            # Calculate the equation of a line to the negative of this pixel
            grad = 0.5*(tilts[x]+tilts[x-1])
            polyarr[0,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[0,1] = 0.5*(ordycen[x]+ordycen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            polyarr[1,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[1,1] = 0.5*(ordycen[x]+ordycen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            # Calculate the equation of a line to the positive of this pixel
            grad = tilts[x]
            polyarr[2,0] = ordxcen[x] + 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[2,1] = ordycen[x] + 0.5*ordwid[x]*grad/(1.0+grad**2)
            polyarr[3,0] = ordxcen[x] - 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[3,1] = ordycen[x] - 0.5*ordwid[x]*grad/(1.0+grad**2)
            rend = 1
        else:
            # Calculate the equation of a line to the negative of this pixel
            grad = 0.5*(tilts[x]+tilts[x-1])
            polyarr[0,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[0,1] = 0.5*(ordycen[x]+ordycen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            polyarr[1,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[1,1] = 0.5*(ordycen[x]+ordycen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            # Calculate the equation of a line to the positive of this pixel
            if ordxcen[x+1] == -999999.9:
                grad = tilts[x]
                polyarr[2,0] = ordxcen[x] + 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[2,1] = ordycen[x] + 0.5*ordwid[x]*grad/(1.0+grad**2)
                polyarr[3,0] = ordxcen[x] - 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[3,1] = ordycen[x] - 0.5*ordwid[x]*grad/(1.0+grad**2)
                rend = 1
            else:
                grad = 0.5*(tilts[x]+tilts[x+1])
                polyarr[2,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                polyarr[2,1] = 0.5*(ordycen[x]+ordycen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
                polyarr[3,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                polyarr[3,1] = 0.5*(ordycen[x]+ordycen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
        # Define the pixel array, integrate, and extract
        totarea = 0.0
        totflux = 0.0
        j=pixll
        donesome=0
        if dispaxis == 0:
            # yloc = slf._pixlocn[:,0,0]
            # xloc = slf._pixlocn[0,:,1]
            # ysiz = slf._pixlocn[:,0,2]
            # xsiz = slf._pixlocn[0,:,3]
            while True:
                if j >= sz_x:
                    j = x+1
                    break
                # Sum up over pixels area
                uarea = 0.0
                i=0
                pixu = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_y) or (pixcen[x]+i < 0):
                        if i>pixu+1:
#							print x, j, i, "BREAK 1"
                            break
                        i += 1
                        continue
                    if (bpix[j,pixcen[x]+i] == 1):
                        if i>pixu+1: break
                        i += 1
                        continue
                    pixlarr[0,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[0,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[1,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[1,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[2,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[2,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[3,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[3,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[0,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[0,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[1,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[1,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[2,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[2,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[3,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[3,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#					polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#					dme += 1
#					if dme >= 2000: return pixlarr, polyarr
                    #area = poly_area(polyout)
                    if area == np.nan:
                        print "AREA=NAN (up)"
                        print outlist
                        print inlist
                        print pixlarr
                    uarea += area
                    totflux += area*frame[j,pixcen[x]+i]
                    if area == 0.0:
                        if i>pixu+1:
#							print x, j, i, "BREAK 2"
                            break
                    else:
                        if i>pixu: pixu=i
                    i += 1
                # Sum up over pixels area
                darea = 0.0
                i=-1
                pixd = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_y) or (pixcen[x]+i < 0):
                        if -i>pixd+1:
#							print x, j, i, "BREAK 3"
                            break
                        i -= 1
                        continue
                    if (bpix[j,pixcen[x]+i] == 1):
                        if -i>pixd+1: break
                        i -= 1
                        continue
                    pixlarr[0,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[0,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[1,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[1,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[2,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[2,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[3,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[3,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[0,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[0,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[1,0,dme] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[1,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[2,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[2,1,dme] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
#					pixlarr[3,0,dme] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
#					pixlarr[3,1,dme] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#					polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#					dme += 1
#					if dme >= 2000: return pixlarr, polyarr
                    #area = poly_area(polyout)
                    if area == np.nan:
                        print "AREA=NAN (down)"
                        print outlist
                        print inlist
                        print pixlarr
                    darea += area
                    totflux += area*frame[j,pixcen[x]+i]
                    if area == 0.0:
                        if -i>pixd+1:
#							print x, j, i, "BREAK 4"
                            break
                    else:
                        if -i>pixd: pixd=-i
                    i -= 1
                # Break if we've already found some overlap
                if uarea+darea == 0.0:
                    if donesome == 1:
#						print x, j, i, "BREAK 5"
                        break
                else:
                    donesome = 1
                # If we haven't found overlap yet, then we must be off the chip
                if j-pixll>=5 and donesome==0:
                    j = x+1
                    break
                j += 1
                totarea += (uarea+darea)
            pixll = pixl
            pixl  = j-1
        else:
            # yloc = slf._pixlocn[0,:,0]
            # xloc = slf._pixlocn[:,0,1]
            # ysiz = slf._pixlocn[0,:,2]
            # xsiz = slf._pixlocn[:,0,3]
            while True:
                if j >= sz_y:
                    j = x+1
                    break
                # Sum up over pixels area
                uarea = 0.0
                i=0
                pixu = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_x) or (pixcen[x]+i < 0):
                        if i>pixu+1: break
                        i += 1
                        continue
                    if (bpix[pixcen[x]+i,j] == 1):
                        if i>pixu+1: break
                        i += 1
                        continue
                    pixlarr[0,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[0,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[1,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[1,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[2,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[2,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[3,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[3,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[0,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[0,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[1,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[1,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[2,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[2,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[3,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[3,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#					polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#					dme += 1
#					if dme >= 2000: return pixlarr, polyarr
                    #area = poly_area(polyout)
                    if area == np.nan:
                        print "AREA=NAN"
                        print outlist
                        print inlist
                        print pixlarr
                    uarea += area
                    totflux += area*frame[pixcen[x]+i,j]
                    if area == 0.0:
                        if i>pixu+1: break
                    else:
                        if i>pixu: pixu=i
                    i += 1
                # Sum up over pixels area
                darea = 0.0
                i=-1
                pixd = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_x) or (pixcen[x]+i < 0):
                        if -i>pixd+1: break
                        i -= 1
                        continue
                    if (bpix[pixcen[x]+i,j] == 1):
                        if -i>pixd+1: break
                        i -= 1
                        continue
                    pixlarr[0,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[0,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[1,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[1,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[2,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[2,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[3,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[3,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[0,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[0,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[1,0,dme] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[1,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[2,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[2,1,dme] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
#					pixlarr[3,0,dme] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
#					pixlarr[3,1,dme] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
#					polyout = SH_poly_clip(pixlarr[:,:,dme],polyarr,outlist,inlist)
#					dme += 1
#					if dme >= 2000: return pixlarr, polyarr
                    #area = poly_area(polyout)
                    if area == np.nan:
                        print "AREA=NAN"
                        print outlist
                        print inlist
                        print pixlarr
                    darea += area
                    totflux += area*frame[pixcen[x]+i,j]
                    if area == 0.0:
                        if -i>pixd+1: break
                    else:
                        if -i>pixd: pixd=-i
                    i -= 1
                # Break if we've already found some overlap
                if uarea+darea == 0.0:
                    if donesome == 1: break
                else:
                    donesome = 1
                # If we haven't found overlap yet, then we must be off the chip
                if j-pixll>=5 and donesome==0:
                    j = x+1
                    break
                j += 1
                totarea += (uarea+darea)
            pixll = pixl
            pixl  = j-1
        if totarea == np.nan:
            print x, totarea, totflux, uarea, darea, donesome, darea==np.nan
            print x, polyarr
        if totarea != 0.0:
            extspec[x] = totflux/totarea
        else:
            extspec[x] = -999999.9
    return extspec


def extract_weighted(np.ndarray[DTYPE_t, ndim=2] frame not None,
                    np.ndarray[DTYPE_t, ndim=2] error not None,
                    np.ndarray[DTYPE_t, ndim=1] profile not None,
                    np.ndarray[DTYPE_t, ndim=2] bpix not None,
                    np.ndarray[ITYPE_t, ndim=1] pixcen not None,
                    np.ndarray[ITYPE_t, ndim=1] piycen not None,
                    np.ndarray[DTYPE_t, ndim=3] pixloc not None,
                    np.ndarray[DTYPE_t, ndim=1] tilts not None,
                    np.ndarray[DTYPE_t, ndim=1] ordxcen not None,
                    np.ndarray[DTYPE_t, ndim=1] ordycen not None,
                    np.ndarray[DTYPE_t, ndim=1] ordwid not None,
                    int ordlen, int ordnum, int intrpnum, int dispaxis):
    """
    Sometimes there is a zero (or a few zeros) padded at the end of
    extspec. This might be because of a difference between ordnum and
    ordxcen[x+1] == -999999.9: rend=1
    ---> I might have solved this issue by defining ext1d to have a
    series of -999999.9 values rather than 0's in arextract.py, for
    when the array doesn't get completely filled in.
    """

    cdef int x, sz_x, sz_y, sz_o, sz_p, xx, yy
    cdef double subarea, xpos, ypos, pweight
    cdef double sza, szb, szc, szd
    cdef double area, uarea, darea, totarea, totflux
    cdef double grad, intc
    cdef int pixll, pixl, pixu, pixd, i, j, donesome

    sz_o = ordxcen.shape[0]
    sz_x = frame.shape[dispaxis]
    sz_y = frame.shape[1-dispaxis]
    sz_p = profile.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] extspec = np.zeros(ordnum, dtype=DTYPE) # The extracted spectrum
    cdef np.ndarray[DTYPE_t, ndim=2] pixlarr = np.zeros((4,2), dtype=DTYPE) # The pixel coordinates
    cdef np.ndarray[DTYPE_t, ndim=2] polyarr = np.zeros((4,2), dtype=DTYPE) # The polygon coordinates
    cdef int nsiz = 11
    cdef np.ndarray[DTYPE_t, ndim=2] outlist = np.zeros((nsiz,2), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] inlist = np.zeros((nsiz,2), dtype=DTYPE)

    cdef int rend = 0
    pixll = 0
    pixl = 0
    for x in range(sz_o):
        # Set the polygon extraction array
        if rend == 1: break # Reached the end of the order, the remaining values are masked
        if x == 0:
            # Calculate the equation of a line to the negative of this pixel
            grad = tilts[x]
            polyarr[0,0] = ordxcen[x] - 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[0,1] = ordycen[x] - 0.5*ordwid[x]*grad/(1.0+grad**2)
            polyarr[1,0] = ordxcen[x] + 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[1,1] = ordycen[x] + 0.5*ordwid[x]*grad/(1.0+grad**2)
            # Calculate the equation of a line to the positive of this pixel
            grad = 0.5*(tilts[x]+tilts[x+1])
            polyarr[2,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
            polyarr[2,1] = 0.5*(ordycen[x]+ordycen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
            polyarr[3,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
            polyarr[3,1] = 0.5*(ordycen[x]+ordycen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
        elif x==sz_o-1:
            # Calculate the equation of a line to the negative of this pixel
            grad = 0.5*(tilts[x]+tilts[x-1])
            polyarr[0,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[0,1] = 0.5*(ordycen[x]+ordycen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            polyarr[1,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[1,1] = 0.5*(ordycen[x]+ordycen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            # Calculate the equation of a line to the positive of this pixel
            grad = tilts[x]
            polyarr[2,0] = ordxcen[x] + 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[2,1] = ordycen[x] + 0.5*ordwid[x]*grad/(1.0+grad**2)
            polyarr[3,0] = ordxcen[x] - 0.5*ordwid[x]/(1.0+grad**2)
            polyarr[3,1] = ordycen[x] - 0.5*ordwid[x]*grad/(1.0+grad**2)
            rend = 1
        else:
            # Calculate the equation of a line to the negative of this pixel
            grad = 0.5*(tilts[x]+tilts[x-1])
            polyarr[0,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[0,1] = 0.5*(ordycen[x]+ordycen[x-1]) - 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            polyarr[1,0] = 0.5*(ordxcen[x]+ordxcen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])/(1.0+grad**2)
            polyarr[1,1] = 0.5*(ordycen[x]+ordycen[x-1]) + 0.25*(ordwid[x]+ordwid[x-1])*grad/(1.0+grad**2)
            # Calculate the equation of a line to the positive of this pixel
            if ordxcen[x+1] == -999999.9:
                grad = tilts[x]
                polyarr[2,0] = ordxcen[x] + 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[2,1] = ordycen[x] + 0.5*ordwid[x]*grad/(1.0+grad**2)
                polyarr[3,0] = ordxcen[x] - 0.5*ordwid[x]/(1.0+grad**2)
                polyarr[3,1] = ordycen[x] - 0.5*ordwid[x]*grad/(1.0+grad**2)
                rend = 1
            else:
                grad = 0.5*(tilts[x]+tilts[x+1])
                polyarr[2,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                polyarr[2,1] = 0.5*(ordycen[x]+ordycen[x+1]) + 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
                polyarr[3,0] = 0.5*(ordxcen[x]+ordxcen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])/(1.0+grad**2)
                polyarr[3,1] = 0.5*(ordycen[x]+ordycen[x+1]) - 0.25*(ordwid[x]+ordwid[x+1])*grad/(1.0+grad**2)
        # Define the pixel array, integrate, and extract
        totarea = 0.0
        totflux = 0.0
        j=pixll
        donesome=0
        if dispaxis == 0:
            # yloc = slf._pixlocn[:,0,0]
            # xloc = slf._pixlocn[0,:,1]
            # ysiz = slf._pixlocn[:,0,2]
            # xsiz = slf._pixlocn[0,:,3]
            while True:
                if j >= sz_x:
                    j = x+1
                    break
                # Sum up over pixels area
                uarea = 0.0
                i=0
                pixu = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_y) or (pixcen[x]+i < 0):
                        if i>pixu+1:
#							print x, j, i, "BREAK 1"
                            break
                        i += 1
                        continue
                    if (bpix[j,pixcen[x]+i] == 1):
                        if i>pixu+1: break
                        i += 1
                        continue
                    pixlarr[0,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[0,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[1,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[1,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[2,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[2,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[3,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[3,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
                    # Calculate the profile weights
                    if area != 0.0:
                        subarea = (pixloc[j,pixcen[x]+i,3]*pixloc[j,pixcen[x]+i,2])/<double>(intrpnum*intrpnum)
                        for xx in range(0,intrpnum):
                            xpos = pixloc[j,pixcen[x]+i,1] + (-0.5 + (<double>(xx)+0.5)/<double>(intrpnum))*pixloc[j,pixcen[x]+i,3]
                            for yy in range(0,intrpnum):
                                ypos = pixloc[j,pixcen[x]+i,0] + (-0.5 + (<double>(yy)+0.5)/<double>(intrpnum))*pixloc[j,pixcen[x]+i,2]
                                pip = point_in_poly(xpos,ypos,inlist)
                                if pip == 1:
                                    sza = csqrt( (polyarr[0,0]-polyarr[3,0])*(polyarr[0,0]-polyarr[3,0]) + (polyarr[0,1]-polyarr[3,1])*(polyarr[0,1]-polyarr[3,1]) )
                                    szb = csqrt( (polyarr[0,0]-xpos)*(polyarr[0,0]-xpos) + (polyarr[0,1]-ypos)*(polyarr[0,1]-ypos) )
                                    szc = csqrt( (polyarr[3,0]-xpos)*(polyarr[3,0]-xpos) + (polyarr[3,1]-ypos)*(polyarr[3,1]-ypos) )
                                    print "Is there a negative sign needed here -- check the cosine rule?"
                                    szd = csqrt( (4.0*sza*sza*szb*szb - cpow((szc*szc - sza*sza - szb*szb),2.0))/(4.0*sza*sza) )
                                    ppix = <int>( 0.5 + <double>(sz_p)*szd/ordwid[x]) # 0.5 is to round, rather than floor.
                                    pweight = profile[ppix]
                                    print "check... is the subarea and/or pweight needed both times?"
                                    uarea += subarea*pweight
                                    totflux += subarea*pweight*frame[j,pixcen[x]+i]
                    if area == 0.0:
                        if i>pixu+1:
#							print x, j, i, "BREAK 2"
                            break
                    else:
                        if i>pixu: pixu=i
                    i += 1
                # Sum up over pixels area
                darea = 0.0
                i=-1
                pixd = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_y) or (pixcen[x]+i < 0):
                        if -i>pixd+1:
#							print x, j, i, "BREAK 3"
                            break
                        i -= 1
                        continue
                    if (bpix[j,pixcen[x]+i] == 1):
                        if -i>pixd+1: break
                        i -= 1
                        continue
                    pixlarr[0,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[0,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[1,0] = pixloc[j,pixcen[x]+i,1]-0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[1,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[2,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[2,1] = pixloc[j,pixcen[x]+i,0]+0.5*pixloc[j,pixcen[x]+i,2]
                    pixlarr[3,0] = pixloc[j,pixcen[x]+i,1]+0.5*pixloc[j,pixcen[x]+i,3]
                    pixlarr[3,1] = pixloc[j,pixcen[x]+i,0]-0.5*pixloc[j,pixcen[x]+i,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
                    darea += area
                    totflux += area*frame[j,pixcen[x]+i]
                    if area == 0.0:
                        if -i>pixd+1:
                            break
                    else:
                        if -i>pixd: pixd=-i
                    i -= 1
                # Break if we've already found some overlap
                if uarea+darea == 0.0:
                    if donesome == 1:
                        break
                else:
                    donesome = 1
                # If we haven't found overlap yet, then we must be off the chip
                if j-pixll>=5 and donesome==0:
                    j = x+1
                    break
                j += 1
                totarea += (uarea+darea)
            pixll = pixl
            pixl  = j-1
        else:
            # yloc = slf._pixlocn[0,:,0]
            # xloc = slf._pixlocn[:,0,1]
            # ysiz = slf._pixlocn[0,:,2]
            # xsiz = slf._pixlocn[:,0,3]
            while True:
                if j >= sz_y:
                    j = x+1
                    break
                # Sum up over pixels area
                uarea = 0.0
                i=0
                pixu = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_x) or (pixcen[x]+i < 0):
                        if i>pixu+1: break
                        i += 1
                        continue
                    if (bpix[pixcen[x]+i,j] == 1):
                        if i>pixu+1: break
                        i += 1
                        continue
                    pixlarr[0,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[0,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[1,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[1,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[2,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[2,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[3,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[3,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
                    if area == np.nan:
                        print "AREA=NAN"
                        print outlist
                        print inlist
                        print pixlarr
                    uarea += area
                    totflux += area*frame[pixcen[x]+i,j]
                    if area == 0.0:
                        if i>pixu+1: break
                    else:
                        if i>pixu: pixu=i
                    i += 1
                # Sum up over pixels area
                darea = 0.0
                i=-1
                pixd = <int>(0.5+ordwid[x]/2.0)
                while True:
                    if (pixcen[x]+i >= sz_x) or (pixcen[x]+i < 0):
                        if -i>pixd+1: break
                        i -= 1
                        continue
                    if (bpix[pixcen[x]+i,j] == 1):
                        if -i>pixd+1: break
                        i -= 1
                        continue
                    pixlarr[0,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[0,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[1,0] = pixloc[pixcen[x]+i,j,1]-0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[1,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[2,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[2,1] = pixloc[pixcen[x]+i,j,0]+0.5*pixloc[pixcen[x]+i,j,2]
                    pixlarr[3,0] = pixloc[pixcen[x]+i,j,1]+0.5*pixloc[pixcen[x]+i,j,3]
                    pixlarr[3,1] = pixloc[pixcen[x]+i,j,0]-0.5*pixloc[pixcen[x]+i,j,2]
                    # Calculate the overlapping polygon area
                    area = SH_poly_clip_area(polyarr,pixlarr,outlist,inlist)
                    if area == np.nan:
                        print "AREA=NAN"
                        print outlist
                        print inlist
                        print pixlarr
                    darea += area
                    totflux += area*frame[pixcen[x]+i,j]
                    if area == 0.0:
                        if -i>pixd+1: break
                    else:
                        if -i>pixd: pixd=-i
                    i -= 1
                # Break if we've already found some overlap
                if uarea+darea == 0.0:
                    if donesome == 1: break
                else:
                    donesome = 1
                # If we haven't found overlap yet, then we must be off the chip
                if j-pixll>=5 and donesome==0:
                    j = x+1
                    break
                j += 1
                totarea += (uarea+darea)
            pixll = pixl
            pixl  = j-1
        if totarea == np.nan:
            print x, totarea, totflux, uarea, darea, donesome, darea==np.nan
            print x, polyarr
        if totarea != 0.0:
            extspec[x] = totflux/totarea
        else:
            extspec[x] = -999999.9
    return extspec



@cython.boundscheck(False)
def extract_ritter_gauss(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                        np.ndarray[DTYPE_t, ndim=2] errfr not None,
                        np.ndarray[DTYPE_t, ndim=1] centint not None,
                        np.ndarray[DTYPE_t, ndim=1] widthint not None,
                        int sz_x):
    # Implements the extraction algorithm discussed by Ritter et al. 1311.4755, assuming a Gaussian profile
    cdef int x, l, sz_l
    cdef double sumW, sumWPD, sumWP, sumWD, sumWPP, Wlx, delL
    cdef double scierrval, skyerrval
    cdef double tstS, tstE

    sz_l = scifr.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] sciext = np.zeros(sz_l, dtype=DTYPE) # The extracted science spectrum
    cdef np.ndarray[DTYPE_t, ndim=1] scierr = np.zeros(sz_l, dtype=DTYPE) # The extracted science error spectrum
    cdef np.ndarray[DTYPE_t, ndim=1] skyext = np.zeros(sz_l, dtype=DTYPE) # The extracted sky spectrum
    cdef np.ndarray[DTYPE_t, ndim=1] skyerr = np.zeros(sz_l, dtype=DTYPE) # The extracted sky error spectrum
    cdef np.ndarray[DTYPE_t, ndim=1] profil = np.zeros(sz_x, dtype=DTYPE) # The profile shape

    for l in range(sz_l):
        sumW   = 0.0
        sumWPD = 0.0
        sumWP  = 0.0
        sumWD  = 0.0
        sumWPP = 0.0
        # Calculate the Profile shape
        gauss_prof(profil,centint[l],widthint[l])
        for x in range(1,sz_x-1):
            tstS = scifr[l,x]+999999.9
            tstE = errfr[l,x]+999999.9
            if tstS < 0.0: tstS *= -1.0
            if tstE < 0.0: tstE *= -1.0
            if tstS < 0.1 or tstE < 0.1 or errfr[l,x] == 0.0: continue
            Wlx = 1.0/(errfr[l,x]*errfr[l,x])
            sumW   += Wlx
            sumWPD += Wlx * profil[x] * scifr[l,x]
            sumWP  += Wlx * profil[x]
            sumWD  += Wlx * scifr[l,x]
            sumWPP += Wlx * profil[x] * profil[x]
        # Calculate the normalization
        delL = sumW*sumWPP - sumWP*sumWP
        if delL == 0.0:
            sciext[l] = -999999.9
            scierr[l] = -999999.9
            skyext[l] = -999999.9
            skyerr[l] = -999999.9
        else:
            # Calculate the errors
            scierrval = sumW/delL
            skyerrval = sumWPP/delL
            if scierrval < 0.0: scierrval = csqrt(-scierrval)
            else: scierrval = csqrt(scierrval)
            if skyerrval < 0.0: skyerrval = csqrt(-skyerrval)
            else: skyerrval = csqrt(skyerrval)
            # Assign these values to the arrays
            sciext[l] = (sumW*sumWPD - sumWP*sumWD)/delL
            scierr[l] = scierrval
            skyext[l] = (sumWPP*sumWD - sumWP*sumWPD)/delL
            skyerr[l] = skyerrval
    # Return the extracted arrays
    return sciext, scierr, skyext, skyerr


@cython.boundscheck(False)
def gauss_prof(np.ndarray[DTYPE_t, ndim=1] profile not None,
                double centroid, double width):
    cdef int x, sz_x
    cdef double xmin, xmax, abl, abr, erfabl, erfabr

    sz_x = profile.shape[0]

    cdef double pxsz = 1.0/(<double>(sz_x-1))
    xmin = -pxsz/2.0
    xmax = pxsz/2.0
    abl = (xmin-centroid)/(csqrt(2.0)*width)
    erfabl = cerf(abl)
    for x in range(sz_x):
        abr = (xmax-centroid)/(csqrt(2.0)*width)
        erfabr = cerf(abr)
        profile[x] = width*csqrt(1.5707963267948966)*(erfabr-erfabl)/(xmax-xmin)
        xmin = xmax
        xmax += pxsz
        abl = abr
        erfabl = erfabr
    return


@cython.boundscheck(False)
def get_locations(np.ndarray[DTYPE_t, ndim=2] inxcen not None,
                    np.ndarray[DTYPE_t, ndim=2] inxwid not None,
                    np.ndarray[DTYPE_t, ndim=1] ypix not None,
                    np.ndarray[DTYPE_t, ndim=2] tilts not None,
                    double binby):
    cdef int o, sz_o, x, sz_x, n, sz_n
    cdef double vol, tsiz, nlen, segm, totl, olen, ang

    sz_x = inxcen.shape[0]
    sz_o = inxcen.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] ordlen = np.zeros(sz_o, dtype=DTYPE) # The pixel box
    cdef np.ndarray[ITYPE_t, ndim=1] ordnum = np.zeros(sz_o, dtype=ITYPE) # The pixel box

    # First derive the total "spectral" (wavelength) length of the order, and how many divisions to make
    # NOTE: The "spectral" length takes into account the tilts of the order and the trace of the order on the CCD.
    sz_n = 0
    for o in range(0,sz_o):
        vol = 0.0
        for x in range(0,sz_x-1):
            # Due to poor notation, this appears like x/y rather than y/x --- as written, it is y/x, which is correct
            ang = catan((inxcen[x+1,o]-inxcen[x,o])/(ypix[x+1]-ypix[x])) - catan(tilts[x,o])
            vol += ccos(ang) * csqrt( (inxcen[x+1,o]-inxcen[x,o])**2 + (ypix[x+1]-ypix[x])**2 )
        ordlen[o] = vol
        ordnum[o] = <int>(0.5 + vol/binby) # 0.5 causes int to round off rather than round down
        if sz_n < ordnum[o]: sz_n = ordnum[o]

    cdef np.ndarray[DTYPE_t, ndim=2] ordxcen = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] ordycen = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] ordtilt = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] ordwid  = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] cumlen  = np.zeros((sz_x), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] xint, yint, tint, wint

    # Now derive the x & y locations of these divisions
    for o in range(0,sz_o):
        vol = 0.0
        for x in range(0,sz_x-1):
            cumlen[x] = vol
            # Due to poor notation, this appears like x/y rather than y/x --- as written, it is y/x, which is correct
            ang = catan((inxcen[x+1,o]-inxcen[x,o])/(ypix[x+1]-ypix[x])) - catan(tilts[x,o])
            vol += ccos(ang) * csqrt( (inxcen[x+1,o]-inxcen[x,o])**2 + (ypix[x+1]-ypix[x])**2 )
        cumlen[sz_x-1] = vol
        spllen = np.linspace(0.0,vol,ordnum[o])
        xint = np.interp(spllen,cumlen,inxcen[:,o])
        yint = np.interp(spllen,cumlen,ypix)
        tint = np.interp(spllen,cumlen,tilts[:,o])
        wint = np.interp(spllen,cumlen,inxwid[:,o])
        n = 0
        while n < ordnum[o]:
            ordxcen[n,o] = xint[n]
            ordycen[n,o] = yint[n]
            ordtilt[n,o] = tint[n]
            ordwid[n,o]  = wint[n]
            n += 1
        # Fill in the blanks
        while n < sz_n:
            ordxcen[n,o] = -999999.9
            ordycen[n,o] = -999999.9
            ordtilt[n,o] = -999999.9
            ordwid[n,o]  = -999999.9
            n += 1
    return ordxcen, ordycen, ordtilt, ordlen, ordwid, ordnum


@cython.boundscheck(False)
def get_wavelocations(np.ndarray[DTYPE_t, ndim=2] ordxcen not None,
                    np.ndarray[DTYPE_t, ndim=2] ordycen not None,
                    np.ndarray[DTYPE_t, ndim=2] ordtilt not None,
                    np.ndarray[DTYPE_t, ndim=2] ordwid not None,
                    np.ndarray[DTYPE_t, ndim=1] wave not None,
                    np.ndarray[DTYPE_t, ndim=2] waveids not None,
                    np.ndarray[ITYPE_t, ndim=2] waveidx not None,
                    int sz_n):
    print "ERROR: THIS FUNCTION IS DEPRECATED"
    print "ERROR: THIS FUNCTION IS DEPRECATED"
    print "ERROR: THIS FUNCTION IS DEPRECATED"
    cdef int o, sz_o, w, sz_w, n, sz
    cdef double tval

    sz_w = wave.shape[0]
    sz_o = waveids.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] oordxcen = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] oordycen = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] oordtilt = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] oordwid  = np.zeros((sz_n,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] xint, yint, tint, wint

    # Derive the x & y locations of these divisions, as specified by the wavelength array
    for o in range(0,sz_o):
        sz = waveidx[1,o]-waveidx[0,o]
        xint = np.interp(wave[waveidx[0,o]:waveidx[1,o]],waveids[:,o],ordxcen[:,o])
        yint = np.interp(wave[waveidx[0,o]:waveidx[1,o]],waveids[:,o],ordycen[:,o])
        wint  = np.interp(wave[waveidx[0,o]:waveidx[1,o]],waveids[:,o],ordwid[:,o])
        tint = np.interp(wave[waveidx[0,o]:waveidx[1,o]],waveids[:,o],ordtilt[:,o])
        n = 0
        while n < sz:
            oordxcen[n,o] = xint[n]
            oordycen[n,o] = yint[n]
            oordtilt[n,o] = tint[n]
            oordwid[n,o]  = wint[n]
            n += 1
        # Fill in the blanks
        while n < sz_n:
            oordxcen[n,o] = -999999.9
            oordycen[n,o] = -999999.9
            oordtilt[n,o] = -999999.9
            oordwid[n,o]  = -999999.9
            n += 1
    return oordxcen, oordycen, oordtilt, oordwid


@cython.boundscheck(False)
def maskedmean_aperture(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                        np.ndarray[DTYPE_t, ndim=2] varfr not None,
                        np.ndarray[ITYPE_t, ndim=1] bckidx not None,
                        np.ndarray[DTYPE_t, ndim=2] maskcr not None,
                        np.ndarray[DTYPE_t, ndim=1] masksky not None,
                        double maskval, int axis):

    cdef int r, sz_r, c, sz_c
    cdef double avprof, wtprof, tst

    sz_r = scifr.shape[0]
    sz_c = scifr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] meanaxA = np.zeros((sz_c), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] meanaxB = np.zeros((sz_r), dtype=DTYPE)

    if axis == 0:
        # Calculate the profile for each row
        for c in range(0,sz_c):
            avprof = 0.0
            wtprof = 0.0
            # Calculate weighted mean profile along this row
            for r in range(0,sz_r):
                tst = scifr[r,c] - maskval
                if tst < 0.0: tst *= -1.0
                if maskcr[r,c] == 1.0 and masksky[r] == 1.0 and tst > 0.1 and varfr[r,c] != 0.0 and bckidx[r] == 1:
                    avprof += scifr[r,c]/varfr[r,c]
                    wtprof += 1.0/varfr[r,c]
            if wtprof == 0.0:
                meanaxA[c] = maskval
            else:
                meanaxA[c] = avprof/wtprof
        return meanaxA
    else:
        # Calculate the profile for each column
        for r in range(0,sz_r):
            avprof = 0.0
            wtprof = 0.0
            # Calculate weighted mean profile along this column
            for c in range(0,sz_c):
                tst = scifr[r,c] - maskval
                if tst < 0.0: tst *= -1.0
                if maskcr[r,c] == 1.0 and masksky[r] == 1.0 and tst > 0.1 and varfr[r,c] != 0.0 and bckidx[c] == 1:
                    avprof += scifr[r,c]/varfr[r,c]
                    wtprof += 1.0/varfr[r,c]
            if wtprof == 0.0:
                meanaxB[r] = maskval
            else:
                meanaxB[r] = avprof/wtprof
        return meanaxB


@cython.boundscheck(False)
def masked_mean(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                np.ndarray[DTYPE_t, ndim=2] errfr not None,
                np.ndarray[DTYPE_t, ndim=2] maskcr not None,
                np.ndarray[DTYPE_t, ndim=1] masksky not None,
                double maskval):

    cdef int r, sz_r, c, sz_c
    cdef double avprof, wtprof, tst

    sz_r = scifr.shape[0]
    sz_c = scifr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] profile = np.zeros((sz_c), dtype=DTYPE)

    # Calculate the profile for each column
    for c in range(0,sz_c):
        avprof = 0.0
        wtprof = 0.0
        # Calculate weighted mean profile along this column
        for r in range(0,sz_r):
            tst = scifr[r,c] - maskval
            if tst < 0.0: tst *= -1.0
            if maskcr[r,c] == 1.0 and masksky[r] == 1.0 and tst > 0.1 and errfr[r,c] != 0.0:
                avprof += scifr[r,c]#/(errfr[r,c]*errfr[r,c])
                wtprof += 1.0#/(errfr[r,c]*errfr[r,c])
        if wtprof == 0.0:
            profile[c] = maskval
        else:
            profile[c] = avprof/wtprof
    return profile

@cython.boundscheck(False)
def maskedaverage_order(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                        np.ndarray[DTYPE_t, ndim=2] varfr not None,
                        double maskval):

    cdef int x, y, sz_x, sz_y
    cdef double sum, sumw

    sz_x = scifr.shape[0]
    sz_y = scifr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] avarr = np.zeros((sz_y), dtype=DTYPE)

    for y in range(sz_y):
        sum  = 0.0
        sumw = 0.0
        for x in range(sz_x):
            if (scifr[x,y] != maskval) and (varfr[x,y] > 0.0):
                sum  += scifr[x,y]/varfr[x,y]
                sumw += 1.0/varfr[x,y]
        if sumw == 0.0: avarr[y] = maskval # Weights cannot be = 0.0 if any value is unmasked.
        else: avarr[y] = sum/sumw
    return avarr


@cython.boundscheck(False)
def maskedmedian_order(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                        double maskval):

    cdef int x, y, j, sz_x, sz_y, cnt
    cdef double temp, numb

    sz_x = scifr.shape[0]
    sz_y = scifr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] medarr = np.zeros((sz_y), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] arrt = np.zeros((sz_x), dtype=DTYPE)

    for y in range(sz_y):
        cnt = 0
        for x in range(sz_x):
            # Fill in the array
            if scifr[x,y] != maskval:
                arrt[cnt] = scifr[x,y]
                cnt += 1
        if cnt == 0:
            medarr[y] = maskval
            continue
        # Sort the array
        for x in range(cnt-1):
            for j in range(x+1,cnt):
                if arrt[j] < arrt[x]:
                    temp = arrt[x]
                    arrt[x] = arrt[j]
                    arrt[j] = temp
        # Consider only the innermost 50 percent of pixels, and find the mean value
#		for x in range(cnt/4,(3*cnt)/4): medarr[y] += arrt[x]
#		temp = <double> ( (3*cnt)/4 - cnt/4 )
#		medarr[y] /= temp
        # Calculate the median
        if cnt%2==0: medarr[y] = 0.5*(arrt[cnt/2] + arrt[cnt/2 - 1])
        else: medarr[y] = arrt[(cnt-1)/2]
    return medarr


@cython.boundscheck(False)
def optimal_extract(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                        np.ndarray[DTYPE_t, ndim=2] errfr not None,
                        np.ndarray[DTYPE_t, ndim=2] maskcr not None,
                        np.ndarray[DTYPE_t, ndim=2] profile not None,
                        double maskval):

    cdef int r, sz_r, c, sz_c
    cdef double sumW, sumWPD, sumWP, sumWD, sumWPP, wgt, delL
    cdef double obj, bck, tst

    sz_r = scifr.shape[0]
    sz_c = scifr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] objarr = np.zeros((sz_r), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] bckarr = np.zeros((sz_r), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] objerr = np.zeros((sz_r), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] bckerr = np.zeros((sz_r), dtype=DTYPE)

    for r in range(0,sz_r):
        sumW   = 0.0
        sumWPD = 0.0
        sumWP  = 0.0
        sumWD  = 0.0
        sumWPP = 0.0
        for c in range(0,sz_c):
            tst = scifr[r,c] - maskval
            if tst < 0.0: tst *= -1.0
            if tst < 0.1 or errfr[r,c] == 0.0: continue
            wgt = maskcr[r,c]/(errfr[r,c]*errfr[r,c])
            sumW   += wgt
            sumWPD += wgt * profile[r,c] * scifr[r,c]
            sumWP  += wgt * profile[r,c]
            sumWD  += wgt * scifr[r,c]
            sumWPP += wgt * profile[r,c] * profile[r,c]
        delL = sumW*sumWPP - sumWP*sumWP
        obj  = sumW*sumWPD - sumWP*sumWD
        bck  = sumWPP*sumWD - sumWP*sumWPD
        if delL == 0.0:
            objarr[r] = maskval
            bckarr[r] = maskval
            objerr[r] = maskval
            bckerr[r] = maskval
        else:
            objarr[r] = obj/delL
            bckarr[r] = bck/delL
            objarr[r] = sumW/delL
            bckarr[r] = sumWPP/delL
    return objarr, objerr, bckarr, bckerr


@cython.boundscheck(False)
def optimal_getprofile(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                        np.ndarray[DTYPE_t, ndim=2] errfr not None,
                        np.ndarray[DTYPE_t, ndim=1] backgr not None,
                        np.ndarray[DTYPE_t, ndim=2] maskcr not None,
                        np.ndarray[DTYPE_t, ndim=1] masksky not None,
                        double maskval, int flag):

    cdef int r, sz_r, c, sz_c
    cdef double avprof, wtprof
    cdef double minA, minB, temp, tst

    sz_r = scifr.shape[0]
    sz_c = scifr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] profile = np.zeros((sz_r,sz_c), dtype=DTYPE)

    if flag == 1: # Make a first guess at the background by averaging the two lowest pixels
        for r in range(0,sz_r):
            c = 0
            minA = maskval
            minB = maskval
            while c < sz_c:
                if scifr[r,c] == maskval:
                    c += 1
                    continue
                if minA == maskval:
                    minA = scifr[r,c]
                    c += 1
                    continue
                if minB == maskval:
                    if scifr[r,c] < minA:
                        minB = minA
                        minA = scifr[r,c]
                    else:
                        minB = scifr[r,c]
                    c += 1
                    continue
                if scifr[r,c] < minA:
                    minB = minA
                    minA = scifr[r,c]
                elif scifr[r,c] < minB:
                    minB = scifr[r,c]
                c += 1
            if minB == maskval:
                if minA != maskval:
                    backgr[r] = minA
            else:
                backgr[r] = 0.5*(minA+minB)

    # Calculate the profile for each column
    for c in range(0,sz_c):
        avprof = 0.0
        wtprof = 0.0
        # Calculate weighted mean profile along this column
        for r in range(0,sz_r):
            tst = scifr[r,c] - maskval
            if tst < 0.0: tst *= -1.0
            if maskcr[r,c] == 1.0 and masksky[r] == 1.0 and tst > 0.1 and errfr[r,c] != 0.0:
                avprof += (scifr[r,c]-backgr[r])#/(errfr[r,c]*errfr[r,c])
                wtprof += 1.0#/(errfr[r,c]*errfr[r,c])
        # Assign the same profile to all rows
        for r in range(0,sz_r):
            profile[r,c] = avprof/wtprof
    # Normalize the profile
#	for r in range(0,sz_r):
#		avprof = profile[r,c]
#		for c in range(1,sz_c): avprof += profile[r,c]
#		for c in range(0,sz_c): profile[r,c] /= avprof
    # Return the normalized profile
    return profile

@cython.boundscheck(False)
def optimal_normalize(np.ndarray[DTYPE_t, ndim=2] profile not None):

    cdef int r, sz_r, c, sz_c
    cdef double avprof

    sz_r = profile.shape[0]
    sz_c = profile.shape[1]

    # Normalize the profile
    for r in range(0,sz_r):
        avprof = profile[r,c]
        for c in range(1,sz_c): avprof += profile[r,c]
        for c in range(0,sz_c): profile[r,c] /= avprof
    return profile

@cython.boundscheck(False)
def optimal_scaleerr(np.ndarray[DTYPE_t, ndim=2] scifr not None,
                        np.ndarray[DTYPE_t, ndim=2] errfr not None,
                        np.ndarray[DTYPE_t, ndim=2] maskcr not None,
                        np.ndarray[DTYPE_t, ndim=1] object not None,
                        np.ndarray[DTYPE_t, ndim=1] backgr not None,
                        np.ndarray[DTYPE_t, ndim=2] profile not None,
                        double maskval):

    cdef int r, sz_r, c, sz_c
    cdef int i, j, cnt
    cdef double scl, temp, meds, mads, mede

    sz_r = scifr.shape[0]
    sz_c = scifr.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] scaleerr = np.zeros((sz_r,sz_c), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] devarr = np.zeros((sz_r,sz_c), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] marr = np.zeros((sz_c), dtype=DTYPE)

    # Calculate the deviations from the model
    for r in range(0,sz_r):
        for c in range(0,sz_c):
            devarr[r,c] = scifr[r,c] - backgr[r] - profile[r,c]*object[r]

    # Calculate the scaling factor to be applied to the error spectrum
    for r in range(0,sz_r):
        ### Calculate the median deviation
        cnt = 0
        for c in range(0,sz_c):
            if scifr[r,c] != maskval and maskcr[r,c] != 0.0 and errfr[r,c] != 0.0:
                marr[cnt] = devarr[r,c]
                cnt += 1
        if cnt == 0:
            for c in range(0,sz_c):
                scaleerr[r,c] = errfr[r,c]
            continue
        # Sort the array
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if marr[j] < marr[i]:
                    temp = marr[i]
                    marr[i] = marr[j]
                    marr[j] = temp
        # Find the median value
        if cnt%2==0:
            meds = 0.5*(marr[cnt/2] + marr[cnt/2 - 1])
        else:
            meds = marr[(cnt-1)/2]
        ### Calculate the median absolute deviation --> standard deviation
        for i in range(0,cnt):
            temp = marr[i]-meds
            if temp < 0.0:
                marr[i] = -temp
            else:
                marr[i] = temp
        # Now get the median of madarr
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if marr[j] < marr[i]:
                    temp = marr[i]
                    marr[i] = marr[j]
                    marr[j] = temp
        # Find the median value (and multiply by 1.4826 to get an estimate of the standard deviation)
        if cnt%2==0:
            mads = 1.4826*0.5*(marr[cnt/2] + marr[cnt/2 - 1])
        else:
            mads = 1.4826*marr[(cnt-1)/2]
        ### Calculate the median error value
        cnt = 0
        for c in range(0,sz_c):
            if scifr[r,c] != maskval and maskcr[r,c] != 0.0 and errfr[r,c] != 0.0:
                marr[cnt] = errfr[r,c]
                cnt += 1
        # Sort the array
        for i in range(cnt-1):
            for j in range(i+1,cnt):
                if marr[j] < marr[i]:
                    temp = marr[i]
                    marr[i] = marr[j]
                    marr[j] = temp
        # Find the median value
        if cnt%2==0:
            mede = 0.5*(marr[cnt/2] + marr[cnt/2 - 1])
        else:
            mede = marr[(cnt-1)/2]

        if mede == 0.0: scl = 1.0
        else: scl = mads/mede
        # Scale the error spectrum for each column of this row
        for c in range(0,sz_c):
            scaleerr[r,c] = scl * errfr[r,c]
    # Return the scaled error spectrum
    return scaleerr


# @cython.boundscheck(False)
# def optimal(np.ndarray[DTYPE_t, ndim=2] prof not None,
# 			np.ndarray[DTYPE_t, ndim=2] xloc not None,
# 			np.ndarray[DTYPE_t, ndim=2] xsiz not None,
# 			int pixwid, int argf_profsamp):
# 
# 	cdef int i, j
# 
# 	shft = pixwid/2
# 	icen = None
# 	ynew = -999999.9*np.ones((frame.shape[dispaxis],2*shft + 1))
# 	for i in range(frame.shape[dispaxis]):
# 		minv = pixcen[i,o]-shft
# 		maxv = pixcen[i,o]+shft+1
# 		if minv < 0: minv = 0
# 		elif minv > frame.shape[1-dispaxis]: continue
# 		if maxv > frame.shape[1-dispaxis]: maxv = frame.shape[1-dispaxis]
# 		elif maxv < 0: continue
# 		xfit = xloc[minv:maxv]
# 		xprof = np.linspace(xloc[minv]-xsiz[minv]/2.0,xloc[maxv-1]+xloc[maxv-1]/2.0,argf_profsamp*np.max(pixwid))
# 		for j in range(0,maxv-minv):
# 			xmin = xloc[minv+j]-xsiz[minv+j]/2.0
# 			xmax = xloc[minv+j]+xsiz[minv+j]/2.0
# 			
# 			sum all pixels between xloc[minv]-xsiz[minv]/2.0 and xloc[minv]+xsiz[minv]/2.0
# 			the number of pixels between these two values indicates the number of profile points to be summed.




# @cython.boundscheck(False)
# def optimal(np.ndarray[DTYPE_t, ndim=2] p not None):
# 	"""
# 	This function performs an iterative optimal extraction, fitting
# 	a constant background, object flux, and object profile. The object's
# 	profile is subpixellated (i.e. oversampled).
# 	"""
# 	
# 	return blah


# @cython.boundscheck(False)
# def optimal_profileweights():
# 
# 	return weights


@cython.boundscheck(False)
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


@cython.boundscheck(False)
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


@cython.boundscheck(False)
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


@cython.boundscheck(False)
def rectify_withcrr(np.ndarray[DTYPE_t, ndim=2] frame not None,
                    np.ndarray[ITYPE_t, ndim=2] bpixmask not None,
                    np.ndarray[ITYPE_t, ndim=1] pixcen not None,
                    np.ndarray[ITYPE_t, ndim=1] pixledg not None,
                    np.ndarray[ITYPE_t, ndim=1] pixredg not None,
                    int pixwid, double sigmadet, double sigmarep,
                    double maskval, int dispaxis):

    cdef int x, y, sz_x, sz_y, sz_xx, sz_yy
    cdef int yidx, ymin, ymax, sz_p
    cdef int shift = pixwid/2

    sz_x  = frame.shape[dispaxis]
    sz_y  = frame.shape[1-dispaxis]
    sz_xx = frame.shape[0]
    sz_yy = frame.shape[1]
    sz_p = 2*shift+1

    cdef np.ndarray[DTYPE_t, ndim=2] rectify = np.zeros((sz_x,sz_p), dtype=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] maskcr = np.ones((sz_x,sz_p), dtype=ITYPE)

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
                if dispaxis == 0:
                    rectify[x,y] = frame[x,yidx]
                else:
                    rectify[x,y] = frame[yidx,x]

    # Identify the bad pixels
    cr_maskmedian(rectify, maskval, sigmadet, sigmarep, maskcr)

    # un-rectify the bad pixel mask
    if dispaxis == 0:
        for x in range(0,sz_xx):
            if pixcen[x] <= 0:
                ymin = (pixredg[x]-shift)-shift
                ymax = (pixredg[x]-shift)+shift+1
            elif pixcen[x] >= sz_yy-1:
                ymin = (pixledg[x]+shift)-shift
                ymax = (pixledg[x]+shift)+shift+1
            else:
                ymin = pixcen[x]-shift
                ymax = pixcen[x]+shift+1
            if ymin < 0: ymin = 0
            elif pixledg[x] >= sz_yy-1: continue
            if pixredg[x] <= 0: continue
            elif ymax > sz_yy: ymax = sz_yy
            for y in range(ymin,ymax):
                bpixmask[x,y] = maskcr[x,y-ymin]
    else:
        for y in range(0,sz_yy):
            if pixcen[y] <= 0:
                ymin = (pixredg[y]-shift)-shift
                ymax = (pixredg[y]-shift)+shift+1
            elif pixcen[y] >= sz_xx-1:
                ymin = (pixledg[y]+shift)-shift
                ymax = (pixledg[y]+shift)+shift+1
            else:
                ymin = pixcen[y]-shift
                ymax = pixcen[y]+shift+1
            if ymin < 0: ymin = 0
            elif pixledg[y] >= sz_xx-1: continue
            if pixredg[y] <= 0: continue
            elif ymax > sz_xx: ymax = sz_xx
            for x in range(ymin,ymax):
                bpixmask[x,y] = maskcr[x,y-ymin]

    return bpixmask

@cython.boundscheck(False)
def rectify(np.ndarray[DTYPE_t, ndim=2] frame not None,
            np.ndarray[ITYPE_t, ndim=2] mask not None,
            np.ndarray[ITYPE_t, ndim=1] pixcen not None,
            np.ndarray[ITYPE_t, ndim=1] pixledg not None,
            np.ndarray[ITYPE_t, ndim=1] pixredg not None,
            int pixwid, double maskval, int dispaxis):

    cdef int x, y, sz_x, sz_y, sz_xx, sz_yy
    cdef int yidx, ymin, ymax, sz_p
    cdef int shift = pixwid/2

    sz_x  = frame.shape[dispaxis]
    sz_y  = frame.shape[1-dispaxis]
    sz_xx = frame.shape[0]
    sz_yy = frame.shape[1]
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
                if dispaxis == 0:
                    if mask[x,yidx] == 0: rectify[x,y] = maskval
                    else: rectify[x,y] = frame[x,yidx]
                else:
                    if mask[yidx,x] == 0: rectify[x,y] = maskval
                    else: rectify[x,y] = frame[yidx,x]
    return rectify


@cython.boundscheck(False)
def rectify_undo(np.ndarray[DTYPE_t, ndim=2] recframe not None,
                np.ndarray[ITYPE_t, ndim=1] pixcen not None,
                np.ndarray[ITYPE_t, ndim=1] pixledg not None,
                np.ndarray[ITYPE_t, ndim=1] pixredg not None,
                int pixwid, double maskval, int sz_xx, int sz_yy,
                int dispaxis):
    """
    To be used in conjunction with the rectify function directly above.
    sz_xx and sz_yy correspond to the frame.shape[0] and frame.shape[1]
    respectively.
    """
    cdef int x, y, ymin, ymax
    cdef int shift = pixwid/2

    cdef np.ndarray[DTYPE_t, ndim=2] unrec = np.zeros((sz_xx,sz_yy), dtype=DTYPE)

    if dispaxis == 0:
        for x in range(0,sz_xx):
            if pixcen[x] <= 0:
                ymin = (pixredg[x]-shift)-shift
                ymax = (pixredg[x]-shift)+shift+1
            elif pixcen[x] >= sz_yy-1:
                ymin = (pixledg[x]+shift)-shift
                ymax = (pixledg[x]+shift)+shift+1
            else:
                ymin = pixcen[x]-shift
                ymax = pixcen[x]+shift+1
            if ymin < 0: ymin = 0
            elif pixledg[x] >= sz_yy-1: continue
            if pixredg[x] <= 0: continue
            elif ymax > sz_yy: ymax = sz_yy
            for y in range(ymin,ymax):
                unrec[x,y] = recframe[x,y-ymin]
    else:
        for y in range(0,sz_yy):
            if pixcen[y] <= 0:
                ymin = (pixredg[y]-shift)-shift
                ymax = (pixredg[y]-shift)+shift+1
            elif pixcen[y] >= sz_xx-1:
                ymin = (pixledg[y]+shift)-shift
                ymax = (pixledg[y]+shift)+shift+1
            else:
                ymin = pixcen[y]-shift
                ymax = pixcen[y]+shift+1
            if ymin < 0: ymin = 0
            elif pixledg[y] >= sz_xx-1: continue
            if pixredg[y] <= 0: continue
            elif ymax > sz_xx: ymax = sz_xx
            for x in range(ymin,ymax):
                unrec[y,x] = recframe[x,y-ymin]
    return unrec


##################################
#  Sutherland-Hodgman Algorithm  #
##################################

@cython.boundscheck(False)
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

########
# OLD SLOW VERSION
########

@cython.boundscheck(False)
def OLD_SH_inside(double px, double py, double cp1x, double cp1y, double cp2x, double cp2y):
    if (cp2x-cp1x)*(py-cp1y) > (cp2y-cp1y)*(px-cp1x): return 1
    else: return 0


@cython.boundscheck(False)
def OLD_SH_intersection(double sx, double sy,
                double ex, double ey,
                double cp1x, double cp1y,
                double cp2x, double cp2y):
    cdef double dcx, dcy, dpx, dpy, n1, n2, n3, rx, ry
    dcx = cp1x - cp2x
    dcy = cp1y - cp2y
    dpx = sx - ex
    dpy = sy - ey
    n1 = cp1x * cp2y - cp1y * cp2x
    n2 = sx * ey - sy * ex
    n3 = 1.0 / (dcx * dpy - dcy * dpx)
    rx = (n1*dpx - n2*dcx) * n3
    ry = (n1*dpy - n2*dcy) * n3
    return rx, ry

@cython.boundscheck(False)
def OLD_SH_poly_clip(np.ndarray[DTYPE_t, ndim=2] poly not None,
              np.ndarray[DTYPE_t, ndim=2] pixl not None,
              np.ndarray[DTYPE_t, ndim=2] outlist not None,
              np.ndarray[DTYPE_t, ndim=2] inlist not None):

    cdef int i, j, o
    cdef int sz_p, sz_e
    cdef double cp1x, cp1y, cp2x, cp2y
    cdef double sx, sy, ex, ey, rx, ry

    sz_p = pixl.shape[0]
    o = poly.shape[0]

    cp1x = pixl[sz_p-1,0]
    cp1y = pixl[sz_p-1,1]

    OLD_SH_empty(outlist)
    OLD_SH_empty(inlist)
    OLD_SH_set(outlist,poly)
    for i in range(0, sz_p):
        cp2x = pixl[i,0]
        cp2y = pixl[i,1]
        OLD_SH_update(inlist,outlist)
        OLD_SH_empty(outlist)
        sz_e = o
        sx = inlist[o-1,0]
        sy = inlist[o-1,1]
        o = 0
        for j in range(0,sz_e):
            ex = inlist[j,0]
            ey = inlist[j,1]
            if OLD_SH_inside(ex, ey, cp1x, cp1y, cp2x, cp2y) == 1:
                if OLD_SH_inside(sx, sy, cp1x, cp1y, cp2x, cp2y) == 0:
                    rx, ry = OLD_SH_intersection(sx, sy, ex, ey, cp1x, cp1y, cp2x, cp2y)
                    outlist[o,0] = rx
                    outlist[o,1] = ry
                    o += 1
                outlist[o,0] = ex
                outlist[o,1] = ey
                o += 1
            elif OLD_SH_inside(sx, sy, cp1x, cp1y, cp2x, cp2y) == 1:
                rx, ry = OLD_SH_intersection(sx, sy, ex, ey, cp1x, cp1y, cp2x, cp2y)
                outlist[o,0] = rx
                outlist[o,1] = ry
                o += 1
            sx = ex
            sy = ey
        cp1x = cp2x
        cp1y = cp2y
    return clip_arr(outlist,-999999.9)


@cython.boundscheck(False)
def OLD_SH_empty(np.ndarray[DTYPE_t, ndim=2] outlist not None):
    cdef int x, sz_x
    sz_x = outlist.shape[0]
    for x in range(0,sz_x):
        if outlist[x,0] == -999999.9: break
        outlist[x,0] = -999999.9
        outlist[x,1] = -999999.9
    return

@cython.boundscheck(False)
def OLD_SH_update(np.ndarray[DTYPE_t, ndim=2] inlist not None,
            np.ndarray[DTYPE_t, ndim=2] outlist not None):
    cdef int x, sz_x
    sz_x = outlist.shape[0]
    for x in range(0,sz_x):
        if inlist[x,0] == -999999.9 and outlist[x,0] == -999999.9: break
        inlist[x,0] = outlist[x,0]
        inlist[x,1] = outlist[x,1]
    return

@cython.boundscheck(False)
def OLD_SH_set(np.ndarray[DTYPE_t, ndim=2] outlist not None,
            np.ndarray[DTYPE_t, ndim=2] poly not None):
    cdef int x, sz_x
    sz_x = poly.shape[0]
    for x in range(0,sz_x):
        outlist[x,0] = poly[x,0]
        outlist[x,1] = poly[x,1]
    return
