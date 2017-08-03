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
    void gsl_multifit_linear_free(gsl_multifit_linear_workspace *w)
    void gsl_matrix_free(gsl_matrix *A)
    void gsl_vector_free(gsl_vector *v)


#@cython.boundscheck(False)
def clean_edges(np.ndarray[DTYPE_t, ndim=2] diff not None,
                np.ndarray[DTYPE_t, ndim=2] tedges not None):
    cdef int sz_x, sz_y
    cdef int x, y
    cdef int ldct, rdct, ymax
    cdef double lmax, rmax

    sz_x = tedges.shape[0]
    sz_y = tedges.shape[1]

    # Set up the array which marks the edge detections
    cdef np.ndarray[ITYPE_t, ndim=2] edgdet = np.zeros((sz_x, sz_y), dtype=ITYPE)

    for x in range(sz_x):
        ldct = 0
        rdct = 0
        for y in range(sz_y):
            if tedges[x, y] == -1.0:
                if ldct == 1:
                    if diff[x, y] > lmax:  # Recall, diff is very positive for significant left edge detections
                        lmax = diff[x, y]
                        ymax = y
                else:
                    lmax = diff[x, y]
                    ymax = y
                    ldct = 1
            elif tedges[x, y] == 1.0:
                if rdct == 1:
                    if diff[x, y] < rmax:  # Recall, diff is very negative for significant right edge detections
                        rmax = diff[x, y]
                        ymax = y
                else:
                    rmax = diff[x, y]
                    ymax = y
                    rdct = 1
            else:
                if ldct == 1:
                    edgdet[x, ymax] = -1
                    ldct = 0
                elif rdct == 1:
                    edgdet[x, ymax] = +1
                    rdct = 0
    return edgdet


#@cython.boundscheck(False)
def close_edges(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                np.ndarray[ITYPE_t, ndim=1] dets not None,
                int npix):

    cdef int sz_x, sz_y, sz_d
    cdef int x, y, d, s, mgap
#    cdef int tmp, tix

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]
    sz_d = dets.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=1] hasedge = np.zeros(sz_d, dtype=ITYPE)

    for d in range(0, sz_d):
#        tmp = sz_y
        for x in range(0, sz_x):
            for y in range(0, sz_y):
                if edgdet[x, y] != dets[d]:
                    continue
                else:
                    # Check if there's an edge nearby
                    mgap = y+npix+1
                    # Check limits
                    if mgap > sz_y:
                        mgap = sz_y
                    for s in range(y+1, mgap):
                        if edgdet[x, s] == dets[d]:
                            hasedge[d] = 1
                            break
                if hasedge[d] != 0:
                    break
            if hasedge[d] != 0:
                break
#        if tmp != sz_y:
#            hasedge[d] = tix
    return hasedge


#@cython.boundscheck(False)
def close_slits(np.ndarray[DTYPE_t, ndim=2] trframe not None,
                np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                np.ndarray[ITYPE_t, ndim=1] dets not None,
                int npix, int ednum):

    cdef int sz_x, sz_y, sz_d
    cdef int x, y, d, s, mgap, lgap, rgap, enum
    cdef double lminv, lmaxv, rminv, rmaxv
    cdef int tmp, tix, flg, diff

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]
    sz_d = dets.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=2] edgearr = np.zeros((sz_x, sz_y), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] hasedge = np.zeros(sz_d, dtype=ITYPE)

    for d in range(0, sz_d):
        tmp = sz_y
        for x in range(0, sz_x):
            for y in range(0, sz_y):
                if edgdet[x, y] != dets[d]:
                    continue
                else:
                    # Check if there's an edge nearby
                    mgap = y+npix+1
                    # Check limits
                    if mgap > sz_y:
                        mgap = sz_y
                    for s in range(y+1, mgap):
                        if edgdet[x, s] != 0:
                            if s-y < tmp:
                                tmp = s-y
                                tix = edgdet[x, s]
                            hasedge[d] = edgdet[x, s]
                            break
#                if hasedge[d] != 0:
#                    break
#            if hasedge[d] != 0:
#                break
        if tmp != sz_y:
            hasedge[d] = tix

    # Now, if there's an edge in hasedge, mark the corresponding index in hasedge with -1
    for d in range(0, sz_d):
        if hasedge[d] == dets[d]:
            # Close slits have caused a left/right edge to be labelled as one edge
            # Find only instances where there is a left and right edge detection
            # Then, take their average and set hadedge to be zero
            tmp = 0
            diff = 0
            for x in range(0, sz_x):
                for y in range(0, sz_y):
                    if edgdet[x, y] != dets[d]:
                        continue
                    else:
                        # Check if there's an edge nearby
                        mgap = y+npix+1
                        # Check limits
                        if mgap > sz_y:
                            mgap = sz_y
                        flg = 0
                        for s in range(y+1, mgap):
                            if edgdet[x, s] == edgdet[x, y]:
                                edgdet[x, s] = 0
                                edgdet[x, y] = 0
                                tix = y + <int>(0.5*<double>(s-y) + 0.5)  # +0.5 for rounding
                                edgdet[x, tix] = dets[d]
                                flg = 1
                                tmp += 1
                                diff += (s-y)
                                break
                        if flg == 0:
                            # If there isn't a common left/right edge for this pixel, ignore this single edge detection
                            edgdet[x, y] = 0
            hasedge[d] = diff/tmp
            continue
        if hasedge[d] > 0:
            for s in range(0, sz_d):
                if hasedge[d] == dets[s]:
                    hasedge[s] = -1
                    break

    # Introduce an edge in cases where no edge exists,
    # and redefine an edge where one does exist.
    enum = ednum
    for d in range(0, sz_d):
        tmp = 0
        for x in range(0, sz_x):
            for y in range(0, sz_y):
                if edgdet[x, y] != dets[d]:
                    continue
                if hasedge[d] >= ednum:
                    edgearr[x, y] = enum
                    # Relabel the appropriate hasedge
                    if tmp == 0:
                        for s in range(0, sz_d):
                            if hasedge[d] == dets[s]:
                                # Label hasedge as negative, to avoid confusion
                                # with the positive hasedge numbers
                                hasedge[s] = -enum
                                tmp = 1
                                break
                elif hasedge[d] < -1:
                    edgearr[x, y] = hasedge[d]
                elif hasedge[d] >= 0:
                    # Create a new edge
                    edgearr[x, y-(1+hasedge[d])] = enum
                    edgearr[x, y+(1+hasedge[d])] = -enum
                else:
                    print "      --->> BUG in arcytrace.close_edges(), check slit traces!!"
        if hasedge[d] >= 0:
            enum += 1
    # Finally return the new slit edges array
    return edgearr


#@cython.boundscheck(False)
def dual_edge(np.ndarray[ITYPE_t, ndim=2] edgearr not None,
              np.ndarray[ITYPE_t, ndim=2] edgearrcp not None,
              np.ndarray[ITYPE_t, ndim=1] wx not None,
              np.ndarray[ITYPE_t, ndim=1] wy not None,
              np.ndarray[ITYPE_t, ndim=1] wl not None,
              np.ndarray[ITYPE_t, ndim=1] wr not None,
              int shft, int npix, int newval):

    cdef int x, y, ee
    cdef int sz_x, sz_y, sz_a, sz_b, sz_e
    cdef int maxy, flg
    sz_x = edgearr.shape[0]
    sz_y = edgearr.shape[1]
    sz_a = wl.shape[0]
    sz_b = wr.shape[0]
    sz_e = wx.shape[0]

    # First go through the leftmost edge (suffix a)
    for x in range(sz_a):
        for ee in range(0, sz_e):
            if edgearr[wx[ee], wy[ee]] == wl[x]:
                # Update the value given to this edge
                edgearrcp[wx[ee], wy[ee]] = newval
                # Determine if an edge can be placed in this row
                maxy = npix
                if wy[ee] + maxy >= sz_y:
                    maxy = sz_y - wy[ee] - 1
                flg = 0
                for y in range(1, maxy):
                    if edgearrcp[wx[ee], wy[ee]+y] != 0:
                        flg = 1
                if flg == 0:
                    edgearrcp[wx[ee], wy[ee]+shft] = newval+1

    # Now go through the rightmost edge (suffix b)
    for x in range(sz_b):
        for ee in range(0, sz_e):
            if edgearr[wx[ee], wy[ee]] == wr[x]:
                # Update the value given to this edge
                edgearrcp[wx[ee], wy[ee]] = newval + 1
                # Determine if an edge can be placed in this row
                maxy = npix
                if wy[ee] - maxy < 0:
                    maxy = wy[ee] + 1
                flg = 0
                for y in range(1, maxy):
                    if edgearrcp[wx[ee], wy[ee]-y] != 0:
                        flg = 1
                if flg == 0:
                    edgearrcp[wx[ee], wy[ee]-shft] = newval

    return


#@cython.boundscheck(False)
def expand_slits(np.ndarray[DTYPE_t, ndim=2] msedge not None,
                 np.ndarray[ITYPE_t, ndim=2] ordcen not None,
                 np.ndarray[ITYPE_t, ndim=1] extord not None):

    cdef int x, sz_x, o, sz_o, y, sz_y, ymax
    cdef int mwid, pwid
    cdef double cenv, edgv

    sz_x = msedge.shape[0]
    sz_y = msedge.shape[1]
    sz_o = ordcen.shape[1]

    cdef np.ndarray[ITYPE_t, ndim=2] pordwid = np.zeros((ordcen.shape[0], sz_o), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] mordwid = np.zeros((ordcen.shape[0], sz_o), dtype=ITYPE)

    # Find the separation between orders
    mwid = -1
    pwid = -1
    for o in range(sz_o):
        for x in range(0, sz_x):
            if extord[o] == 1:
                mwid = -1
                pwid = -1
            else:
                # Trace from centre to left
                if o == 0:
                    # Don't worry about the left edge
                    mwid = -1
                else:
                    ymax = (ordcen[x, o-1] + ordcen[x, o])/2
                    if ymax < 0:
                        mwid = -1
                    else:
                        if msedge[x, ordcen[x, o]] < msedge[x, ymax]:
                            # Slits overlap
                            mwid = -1
                        else:
                            # There's a gap between echelle orders
                            # Find the slit edge
                            edgv = 0.5*(msedge[x, ordcen[x, o]] + msedge[x, ymax])
                            for y in range(ymax, ordcen[x, o]):
                                if msedge[x, y] > edgv:
                                    mwid = ordcen[x, o]-y
                                    break
                # Trace from centre to right
                if o == sz_o-1:
                    # Don't worry about the right edge
                    pwid = -1
                else:
                    ymax = (ordcen[x, o] + ordcen[x, o+1])/2
                    if ymax >= sz_y:
                        pwid = -1
                    else:
                        if msedge[x, ordcen[x, o]] < msedge[x, ymax]:
                            # Slits overlap
                            pwid = -1
                        else:
                            # There's a gap between echelle orders
                            # Find the slit edge
                            edgv = 0.5*(msedge[x, ordcen[x, o]] + msedge[x, ymax])
                            for y in range(0, ymax-ordcen[x, o]):
                                if msedge[x, ymax-y] > edgv:
                                    pwid = (ymax-ordcen[x, o])-y
                                    break
            mordwid[x, o] = mwid
            pordwid[x, o] = pwid

    return mordwid, pordwid


#@cython.boundscheck(False)
def find_peak_limits(np.ndarray[ITYPE_t, ndim=1] hist not None,
                    np.ndarray[ITYPE_t, ndim=1] pks not None):
    """
    Find all values between the zeros of hist
    """

    cdef int ii, sz_i, sz_h, lim

    sz_i = pks.shape[0]
    sz_h = hist.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=2] edges = np.zeros((sz_i,2), dtype=ITYPE)

    for ii in range(0, sz_i):
        # Search below the peak
        lim = pks[ii]
        while True:
            if lim < 0:
                break
            if hist[lim] == 0:
                break
            lim -= 1
        # Store the limit
        edges[ii, 0] = lim
        # Search above the peak
        lim = pks[ii]
        while True:
            if lim > sz_h-1:
                break
            if hist[lim] == 0:
                break
            lim += 1
        # Store the limit
        edges[ii, 1] = lim
    return edges


#@cython.boundscheck(False)
def find_between(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                np.ndarray[ITYPE_t, ndim=1] ledgem not None,
                np.ndarray[ITYPE_t, ndim=1] ledgep not None,
                int dirc):

    cdef int sz_x, sz_y
    cdef int x, y, ymn, ymx, ystrt

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]

    # Setup the coefficient arrays
    cdef np.ndarray[ITYPE_t, ndim=1] edgbtwn = np.zeros((3), dtype=ITYPE)
    edgbtwn[0] = -1
    edgbtwn[1] = -1
    edgbtwn[2] = -1

    for x in range(0,sz_x):
        if ledgem[x]<ledgep[x]:
            ymn = ledgem[x]
            ymx = ledgep[x]
        else:
            ymn = ledgep[x]
            ymx = ledgem[x]
        for y in range(ymn,ymx):
            if edgdet[x,y] > 0:
                if edgbtwn[0]==-1:
                    edgbtwn[0] = edgdet[x,y]
                elif edgdet[x,y] == edgbtwn[0]:
                    continue
                elif edgbtwn[1]==-1:
                    edgbtwn[1] = edgdet[x,y]

    # If no right order edges were found between these two left order edges, find the next right order edge
    if edgbtwn[0] == -1 and edgbtwn[1] == -1:
        for x in range(0,sz_x):
            if ledgem[x]<ledgep[x]:
                ystrt = ledgep[x]
            else:
                ystrt = ledgem[x]
            if dirc == 1:
                while ystrt < sz_y:
                    if edgdet[x,y] > 0:
                        if edgbtwn[2] == -1:
                            edgbtwn[2] = edgdet[x,y]
                        elif edgdet[x,y] < edgbtwn[2]:
                            edgbtwn[2] = edgdet[x,y]
                    ystrt+=1
            else:
                while ystrt >= 0:
                    if edgdet[x,y] > 0:
                        if edgbtwn[2] == -1:
                            edgbtwn[2] = edgdet[x,y]
                        elif edgdet[x,y] > edgbtwn[2]:
                            edgbtwn[2] = edgdet[x,y]
                    ystrt-=1
    # Now return the array
    return edgbtwn


#@cython.boundscheck(False)
def find_objects(np.ndarray[DTYPE_t, ndim=1] profile not None,
                int bgreg, double stddev):
    """
    Find significantly detected objects in the profile array
    For all objects found, the background regions will be defined.
    """
    cdef int o, x, sz_x
    cdef int cntr, imax
    cdef double maxv

    sz_x = profile.shape[0]

    # Define the object centroids array
    cdef np.ndarray[ITYPE_t, ndim=1] objl = -1*np.ones((sz_x), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] objr = -1*np.ones((sz_x), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] msk = np.zeros((sz_x), dtype=ITYPE)

    cntr = 0
    while True:
        # Find maximum flux point
        maxv = -1.0
        imax = -1
        for x in range(sz_x):
            if msk[x] == 1: continue
            if imax == -1:
                imax = x
                maxv = profile[x]
                continue
            if profile[x]>maxv:
                maxv = profile[x]
                imax = x
        if maxv < 5.0*stddev:
            # No more objects left to be found
            break
        msk[imax] = 1
        objl[cntr] = imax
        objr[cntr] = imax
        # Find the left edge of this object
        for x in range(1,imax):
            if profile[imax-x] > 3.0*stddev:
                objl[cntr] -= 1
                msk[imax-x] = 1
            else:
                objl[cntr] -= 1
                msk[imax-x] = 1
                break
        # Find the right edge of this object
        for x in range(imax+1,sz_x):
            if profile[x] > 3.0*stddev:
                objr[cntr] += 1
                msk[x] = 1
            else:
                objr[cntr] += 1
                msk[x] = 1
                break
        cntr += 1

    # Determine the background pixels for each object
    cdef np.ndarray[ITYPE_t, ndim=2] bckl = np.zeros((sz_x,cntr), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] bckr = np.zeros((sz_x,cntr), dtype=ITYPE)
    for o in range(cntr):
        for x in range(1,bgreg+1):
            if objl[o]-x >= 0:
                if msk[objl[o]-x]==0:
                    bckl[objl[o]-x,o]=1
            if objr[o]+x <= sz_x-1:
                if msk[objr[o]+x]==0:
                    bckr[objr[o]+x,o]=1
    return objl[:cntr], objr[:cntr], bckl, bckr


#@cython.boundscheck(False)
def find_shift(np.ndarray[DTYPE_t, ndim=2] mstrace not None,
                np.ndarray[DTYPE_t, ndim=1] minarr not None,
                np.ndarray[ITYPE_t, ndim=1] lopos not None,
                np.ndarray[ITYPE_t, ndim=1] diffarr not None,
                int numsrch):

    cdef int sz_x, sz_y
    cdef int x, y, s
    cdef int shift = 0
    cdef double cnts, npts, maxcnts

    sz_x = mstrace.shape[0]
    sz_y = mstrace.shape[1]

    maxcnts = -999999.9
    shift = 0
    for s in range(0,numsrch):
        cnts = 0.0
        npts = 0.0
        for x in range(0,sz_x):
            ymin = lopos[x]+s
            ymax = ymin + diffarr[x]
            if ymin < 0: ymin = 0
            elif ymax > sz_y: ymax = sz_y
            for y in range(ymin,ymax):
                cnts += (mstrace[x,y]-minarr[x])
                npts += 1.0
        if npts == 0.0:
            continue
        else:
            cnts = cnts/npts
        if cnts > maxcnts:
            maxcnts = cnts
            shift = s
    return shift


#@cython.boundscheck(False)
def ignore_orders(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                int fracpix, int lmin, int lmax, int rmin, int rmax):
    cdef int sz_x, sz_y
    cdef int x, y
    cdef int lnc, lxc, rnc, rxc

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]

    cdef np.ndarray[ITYPE_t, ndim=2] larr = np.zeros((2,lmax-lmin+1), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] rarr = np.zeros((2,rmax-rmin+1), dtype=ITYPE)

    for x in range(lmax-lmin+1):
        larr[0,x] = sz_x
    for x in range(rmax-rmin+1):
        rarr[0,x] = sz_x

    for x in range(sz_x):
        for y in range(sz_y):
            if edgdet[x,y] < 0:
                if x < larr[0,-edgdet[x,y]-lmin]:
                    larr[0,-edgdet[x,y]-lmin] = x
                elif x > larr[1,-edgdet[x,y]-lmin]:
                    larr[1,-edgdet[x,y]-lmin] = x
            elif edgdet[x,y] > 0:
                if x < rarr[0,edgdet[x,y]-rmin]:
                    rarr[0,edgdet[x,y]-rmin] = x
                elif x > rarr[1,edgdet[x,y]-rmin]:
                    rarr[1,edgdet[x,y]-rmin] = x

    # Go through the array once more to remove pixels that do not cover fracpix
    for x in range(sz_x):
        for y in range(sz_y):
            if edgdet[x,y] < 0:
                if larr[1,-edgdet[x,y]-lmin]-larr[0,-edgdet[x,y]-lmin] < fracpix:
                    edgdet[x,y]=0
            elif edgdet[x,y] > 0:
                if rarr[1,edgdet[x,y]-rmin]-rarr[0,edgdet[x,y]-rmin] < fracpix:
                    edgdet[x,y]=0

    # Check if lmin, lmax, rmin, and rmax need to be changed
    lnc = 0
    while True:
        if larr[1,lnc]-larr[0,lnc] < fracpix:
            lnc += 1
        else:
            break
    lxc = 0
    while True:
        if larr[1,lmax-lmin-lxc]-larr[0,lmax-lmin-lxc] < fracpix:
            lxc += 1
        else:
            break
    rnc = 0
    while True:
        if rarr[1,rnc]-rarr[0,rnc] < fracpix:
            rnc += 1
        else:
            break
    rxc = 0
    while True:
        if rarr[1,rmax-rmin-rxc]-rarr[0,rmax-rmin-rxc] < fracpix:
            rxc += 1
        else:
            break
    return lnc, lxc, rnc, rxc, larr, rarr


#@cython.boundscheck(False)
def limit_yval(int yc, int maxv):
    cdef int yn, yx
    if yc == 0:
        yn = 0
        yx = 4
    elif yc == 1:
        yn = -1
        yx = 4
    elif yc == 2:
        yn = -2
        yx = 4
    elif yc == maxv-3:
        yn = -3
        yx = 3
    elif yc == maxv-2:
        yn = -3
        yx = 2
    elif yc == maxv-1:
        yn = -3
        yx = 1
    else:
        yn = -3
        yx = 4
    return yn, yx


#@cython.boundscheck(False)
def locate_order(np.ndarray[DTYPE_t, ndim=1] lordloc not None,
                 np.ndarray[DTYPE_t, ndim=1] rordloc not None,
                 int sz_x, int sz_y, int pad):
    """ Generate a boolean image that identifies which pixels
    belong to the slit associated with the supplied left and
    right slit edges.

    Parameters
    ----------
    lordloc : ndarray
      Location of the left slit edges of 1 slit
    rordloc : ndarray
      Location of the right slit edges of 1 slit
    sz_x : int
      The size of an image in the spectral (0th) dimension
    sz_y : int
      The size of an image in the spatial (1st) dimension
    pad : int
      Additional pixels to pad the left and right slit edges

    Returns
    -------
    orderloc : ndarray
      An image the same size as the input frame, containing values from 0-1.
      0 = pixel is not in the specified slit
      1 = pixel is in the specified slit
    """

    cdef int x, y
    cdef int ymin, ymax
    cdef double ow, oc

    cdef np.ndarray[ITYPE_t, ndim=2] orderloc = np.zeros((sz_x,sz_y), dtype=ITYPE)

    for x in range(0, sz_x):
        ow = (rordloc[x]-lordloc[x])/2.0
        oc = (rordloc[x]+lordloc[x])/2.0
        ymin = <int>(oc-ow)-pad
        ymax = <int>(oc+ow)+1+pad
        # Check we are in bounds
        if ymin < 0:
            ymin = 0
        elif ymax < 0:
            continue
        elif ymax <= ymin:
            continue
        if ymax > sz_y-1:
            ymax = sz_y-1
        elif ymin > sz_y-1:
            continue
        elif ymin >= ymax:
            continue
        # Set these values to 1 for this slit
        for y in range(ymin, ymax):
            orderloc[x,y] = 1
    return orderloc


#@cython.boundscheck(False)
def match_edges(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                int ednum):
    cdef int sz_x, sz_y
    cdef int x, y, yn, yx, xr, xs, yt, s, t
    cdef int lcnt, rcnt, suc, anyt

    cdef int mr = 5   # This is the minimum number of acceptable pixels required to form the detection of an order edge
    cdef np.ndarray[ITYPE_t, ndim=1] mrxarr = np.zeros(mr, dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] mryarr = np.zeros(mr, dtype=ITYPE)

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]

    lcnt = 2*ednum
    rcnt = 2*ednum
    for y in range(sz_y):
        for x in range(sz_x):
            anyt = 0
            if edgdet[x, y] == -1:
                # Search upwards from x,y
                xs = x + 1
                yt = y
                while xs <= sz_x-1:
                    xr = 10
                    yn, yx = limit_yval(yt, sz_y)
                    if xs + xr >= sz_x:
                        xr = sz_x - xs - 1
                    suc = 0
                    for s in range(xs, xs+xr):
                        suc = 0
                        for t in range(yt + yn, yt + yx):
                            if edgdet[s, t] == -1:
                                edgdet[s, t] = -lcnt
                                suc = 1
                                if anyt < mr:
                                    mrxarr[anyt] = s
                                    mryarr[anyt] = t
                                anyt += 1
                                yt = t
                                break
                        if suc == 1:
                            xs = s + 1
                            break
                    if suc == 0: # The trace is lost!
                        break
                # Search downwards from x,y
                xs = x - 1
                yt = y
                while xs >= 0:
                    xr = 10
                    yn, yx = limit_yval(yt, sz_y)
                    if xs-xr < 0:
                        xr = xs
                    suc = 0
                    for s in range(0, xr):
                        suc = 0
                        for t in range(yt+yn, yt+yx):
                            if edgdet[xs-s, t] == -1:
                                edgdet[xs-s, t] = -lcnt
                                suc = 1
                                if anyt < mr:
                                    mrxarr[anyt] = xs-s
                                    mryarr[anyt] = t
                                anyt += 1
                                yt = t
                                break
                        if suc == 1:
                            xs = xs - s - 1
                            break
                    if suc == 0: # The trace is lost!
                        break
                if anyt > mr:
                    edgdet[x, y] = -lcnt
                    lcnt = lcnt + 1
                else:
                    edgdet[x, y] = 0
                    for s in range(anyt):
                        if mrxarr[s] != 0 and mryarr[s] != 0:
                            edgdet[mrxarr[s], mryarr[s]] = 0
            elif edgdet[x, y] == 1:
                # Search upwards from x,y
                xs = x+1
                yt = y
                while xs <= sz_x-1:
                    xr = 10
                    yn, yx = limit_yval(yt, sz_y)
                    if xs+xr >= sz_x:
                        xr = sz_x-xs-1
                    suc = 0
                    for s in range(xs, xs+xr):
                        suc = 0
                        for t in range(yt+yn, yt+yx):
                            if edgdet[s, t] == 1:
                                edgdet[s, t] = rcnt
                                suc = 1
                                if anyt < mr:
                                    mrxarr[anyt] = s
                                    mryarr[anyt] = t
                                anyt += 1
                                yt = t
                                break
                        if suc == 1:
                            xs = s + 1
                            break
                    if suc == 0: # The trace is lost!
                        break
                # Search downwards from x,y
                xs = x-1
                yt = y
                while xs >= 0:
                    xr = 10
                    yn, yx = limit_yval(yt, sz_y)
                    if xs-xr < 0:
                        xr = xs
                    suc = 0
                    for s in range(0, xr):
                        suc = 0
                        for t in range(yt+yn, yt+yx):
                            if edgdet[xs-s, t] == 1:
                                edgdet[xs-s, t] = rcnt
                                suc = 1
                                if anyt < mr:
                                    mrxarr[anyt] = xs-s
                                    mryarr[anyt] = t
                                anyt += 1
                                yt = t
                                break
                        if suc == 1:
                            xs = xs - s - 1
                            break
                    if suc == 0: # The trace is lost!
                        break
                if anyt > mr:
                    edgdet[x, y] = rcnt
                    rcnt = rcnt + 1
                else:
                    edgdet[x, y] = 0
                    for s in range(anyt):
                        if mrxarr[s] != 0 and mryarr[s] != 0:
                            edgdet[mrxarr[s], mryarr[s]] = 0
    return lcnt-2*ednum, rcnt-2*ednum


@cython.boundscheck(False)
def minbetween(np.ndarray[DTYPE_t, ndim=2] mstrace not None,
                np.ndarray[ITYPE_t, ndim=1] loord not None,
                np.ndarray[ITYPE_t, ndim=1] hiord not None):
    cdef int sz_x, sz_y
    cdef int x, y, ymin, ymax

    sz_x = mstrace.shape[0]
    sz_y = mstrace.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] minarr = np.zeros(sz_x, dtype=DTYPE)

    for x in range(sz_x):
        ymin = loord[x]
        ymax = hiord[x]
        if ymin < 0: ymin = 0
        elif ymax > sz_y: ymax = sz_y
        for y in range(ymin,ymax):
            if mstrace[x,y] < minarr[x] and mstrace[x,y] > 0.0:
                minarr[x] = mstrace[x,y]
            elif minarr[x] == 0.0:
                minarr[x] = mstrace[x,y]
    return minarr


#@cython.boundscheck(False)
def phys_to_pix(np.ndarray[DTYPE_t, ndim=2] array not None,
        np.ndarray[DTYPE_t, ndim=1] diff not None):

    cdef int sz_n, sz_a, sz_d
    cdef int a, d, n, mind
    cdef double mindv, test

    sz_a = array.shape[0]
    sz_n = array.shape[1]
    sz_d = diff.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=2] pixarr = np.zeros((sz_a,sz_n), dtype=ITYPE)

    for n in range(0,sz_n):
        for a in range(0,sz_a):
            mind = 0
            mindv = array[a,n]-diff[0]
            for d in range(1,sz_d):
                test = array[a,n]-diff[d]
                if test < 0.0: test *= -1.0
                if test < mindv:
                    mindv = test
                    mind = d
                if array[a,n]-diff[d] < 0.0: break
            pixarr[a,n] = mind
    return pixarr


#@cython.boundscheck(False)
def tilts_image(np.ndarray[DTYPE_t, ndim=2] tilts not None,
                np.ndarray[DTYPE_t, ndim=2] lordloc not None,
                np.ndarray[DTYPE_t, ndim=2] rordloc not None,
                int pad, int sz_y):
    """ Using the tilt (assumed to be fit with a first order polynomial)
    generate an image of the tilts for each slit.

    Parameters
    ----------
    tilts : ndarray
      An (m x n) 2D array specifying the tilt (i.e. gradient) at each
      pixel along the spectral direction (m) for each slit (n).
    lordloc : ndarray
      Location of the left slit edges
    rordloc : ndarray
      Location of the right slit edges
    pad : int
      Set the tilts within each slit, and extend a number of pixels
      outside the slit edges (this number is set by pad).
    sz_y : int
      Number of detector pixels in the spatial direction.

    Returns
    -------
    tiltsimg : ndarray
      An image the same size as the science frame, containing values from 0-1.
      0/1 corresponds to the bottom/top of the detector (in the spectral direction),
      and constant wavelength is represented by a single value from 0-1.
    """

    cdef int o, sz_o, x, sz_x, y
    cdef int ymin, ymax
    cdef double yv, ow, oc, dszx

    sz_x = tilts.shape[0]
    sz_o = tilts.shape[1]
    dszx = (sz_x-1.0)

    cdef np.ndarray[DTYPE_t, ndim=2] tiltsimg = np.zeros((sz_x,sz_y), dtype=DTYPE)

    for o in range(0, sz_o):
        for x in range(0, sz_x):
            ow = (rordloc[x,o]-lordloc[x,o])/2.0
            oc = (rordloc[x,o]+lordloc[x,o])/2.0
            ymin = <int>(oc-ow) - pad
            ymax = <int>(oc+ow) + 1 + pad
            # Check we are in bounds
            if ymin < 0:
                ymin = 0
            elif ymax < 0:
                continue
            if ymax > sz_y-1:
                ymax = sz_y-1
            elif ymin > sz_y-1:
                continue
            # Set the tilt value at each pixel in this row
            for y in range(ymin, ymax):
                yv = (<double>(y)-lordloc[x, o])/ow - 1.0
                tiltsimg[x,y] = (tilts[x,o]*yv + <double>(x))/dszx
    return tiltsimg
