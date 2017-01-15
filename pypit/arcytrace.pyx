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

#######
#  A  #
#######

@cython.boundscheck(False)
def assign_orders(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                int lcnt, int rcnt, int ednum):
    cdef int sz_x, sz_y, x, y, i, ni
    cdef int coml, comr, cntl, cntr
    cdef int nl, nr, tolp, tolm, tf, ymn, ymx
    cdef int idnum, nc
    cdef int rmin, rmax, lmin, lmax

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]

    cdef np.ndarray[ITYPE_t, ndim=1] larr = np.zeros((lcnt), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] rarr = np.zeros((rcnt), dtype=ITYPE)

    # Find which slit edge id is the most common
    for x in range(sz_x):
        for y in range(sz_y):
            if edgdet[x,y] == 0: continue
            elif edgdet[x,y] < 0:
                larr[-2*ednum-edgdet[x,y]] += 1
            elif edgdet[x,y] > 0:
                rarr[edgdet[x,y]-2*ednum] += 1
    coml = 0
    cntl = larr[0]
    for x in range(1,lcnt):
        if larr[x] > cntl:
            cntl = larr[x]
            coml = x
    comr = 0
    cntr = rarr[0]
    for x in range(1,rcnt):
        if rarr[x] > cntr:
            cntr = rarr[x]
            comr = x

    # Obtain the (x,y) values for the most common slit edge id
    cdef np.ndarray[ITYPE_t, ndim=2] lcarr = np.zeros((cntl,2), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] rcarr = np.zeros((cntr,2), dtype=ITYPE)
    nl = 0
    nr = 0
    for x in range(sz_x):
        for y in range(sz_y):
            if edgdet[x,y] == -2*ednum-coml:
                lcarr[nl,0] = x
                lcarr[nl,1] = y
                nl += 1
            elif edgdet[x,y] == 2*ednum+comr:
                rcarr[nr,0] = x
                rcarr[nr,1] = y
                nr += 1
            else:
                continue

    # Label this slit edge
    for x in range(nl):
        edgdet[lcarr[x,0],lcarr[x,1]] = -ednum
    for x in range(nr):
        edgdet[rcarr[x,0],rcarr[x,1]] = ednum

    # Find the closest set of points above and below this slit
    cdef np.ndarray[ITYPE_t, ndim=1] llab = np.zeros((lcnt), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] rlab = np.zeros((rcnt), dtype=ITYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] ldiffp = np.zeros((nl), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] ldiffm = np.zeros((nl), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] rdiffp = np.zeros((nr), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] rdiffm = np.zeros((nr), dtype=DTYPE)
    llab[coml] = ednum # Label the most common slit edge as ednum
    rlab[comr] = ednum # Label the most common slit edge as ednum
    for x in range(nl):
        for y in range(lcarr[x,1]+1,sz_y):
            if edgdet[lcarr[x,0],y] <= -2*ednum:
                ldiffp[x] = <double>(y - lcarr[x,1])
                break
        for y in range(0,lcarr[x,1]):
            if edgdet[lcarr[x,0],lcarr[x,1]-y-1] < 0:
                ldiffm[x] = <double>(y+1)
                break
    for x in range(nr):
        for y in range(rcarr[x,1]+1,sz_y):
            if edgdet[rcarr[x,0],y] >= 2*ednum:
                rdiffp[x] = <double>(y - rcarr[x,1])
                break
        for y in range(0,rcarr[x,1]):
            if edgdet[rcarr[x,0],rcarr[x,1]-y-1] > 0:
                rdiffm[x] = <double>(y+1)
                break
    # Calculate the median difference between the nearest identifications
    tolp = <int>(median(ldiffp))
    tolm = <int>(median(ldiffm))
    tf = 5
    ni = 2
    for x in range(nl):
        for i in range(1,ni+1):
            ymn = lcarr[x,1]+i*tolp-tolp/tf
            ymx = lcarr[x,1]+i*tolp+tolp/tf+1
            if ymn < 0: ymn = 0
            if ymx > sz_y: ymx = sz_y
            for y in range(ymn,ymx):
                if edgdet[lcarr[x,0],y] <= -2*ednum:
                    llab[-2*ednum-edgdet[lcarr[x,0],y]] = ednum+i
                    edgdet[lcarr[x,0],y] = -(ednum+i)
            ymn = lcarr[x,1]-i*tolm-tolm/tf
            ymx = lcarr[x,1]-i*tolm+tolm/tf+1
            if ymn < 0: ymn = 0
            if ymx > sz_y: ymx = sz_y
            for y in range(ymn,ymx):
                if edgdet[lcarr[x,0],y] <= -2*ednum:
                    llab[-2*ednum-edgdet[lcarr[x,0],y]] = ednum-i
                    edgdet[lcarr[x,0],y] = -(ednum-i)
    tolp = <int>(median(rdiffp))
    tolm = <int>(median(rdiffm))
    for x in range(nr):
        for i in range(1,ni+1):
            ymn = rcarr[x,1]+i*tolp-tolp/tf
            ymx = rcarr[x,1]+i*tolp+tolp/tf+1
            if ymn < 0: ymn = 0
            if ymx > sz_y: ymx = sz_y
            for y in range(ymn,ymx):
                if edgdet[rcarr[x,0],y] >= 2*ednum:
                    rlab[edgdet[rcarr[x,0],y]-2*ednum] = ednum+i
                    edgdet[rcarr[x,0],y] = ednum+i
            ymn = rcarr[x,1]-i*tolm-tolm/tf
            ymx = rcarr[x,1]-i*tolm+tolm/tf+1
            if ymn < 0: ymn = 0
            if ymx > sz_y: ymx = sz_y
            for y in range(ymn,ymx):
                if edgdet[rcarr[x,0],y] >= 2*ednum:
                    rlab[edgdet[rcarr[x,0],y]-2*ednum] = ednum-i
                    edgdet[rcarr[x,0],y] = ednum-i

    # Iterate through and label all of the left edges
    idnum = -ednum-1
    while True:
        nc, tolp, tolm = get_xy(edgdet, lcarr, idnum, -1) # -1 specifies the direction to move in
        if nc == 0: break
        if tolp == 0 and tolm == 0: break
        for x in range(nc):
            for i in range(1,ni+1):
                ymn = lcarr[x,1]+i*tolp-tolp/tf
                ymx = lcarr[x,1]+i*tolp+tolp/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[lcarr[x,0],y] <= -2*ednum:
                        llab[-2*ednum-edgdet[lcarr[x,0],y]] = -(idnum-i)
                        edgdet[lcarr[x,0],y] = idnum-i
                ymn = lcarr[x,1]-i*tolm-tolm/tf
                ymx = lcarr[x,1]-i*tolm+tolm/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[lcarr[x,0],y] <= -2*ednum:
                        llab[-2*ednum-edgdet[lcarr[x,0],y]] = -(idnum+i)
                        edgdet[lcarr[x,0],y] = idnum+i
        # Update the labels
        update_labels(edgdet, llab, idnum, ni, -1, ednum)
        idnum -= 1
    #
    # Now search and label the left edges, except now search backwards from edge ednum
    idnum = -ednum+1
    while True:
        nc, tolp, tolm = get_xy(edgdet, lcarr, idnum, +1) # +1 specifies the direction to move in
        if nc == 0: break
        if tolp == 0 and tolm == 0: break
        for x in range(nc):
            for i in range(1,ni+1):
                ymn = lcarr[x,1]+i*tolp-tolp/tf
                ymx = lcarr[x,1]+i*tolp+tolp/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[lcarr[x,0],y] <= -2*ednum:
                        llab[-2*ednum-edgdet[lcarr[x,0],y]] = -(idnum-i)
                        edgdet[lcarr[x,0],y] = idnum-i
                ymn = lcarr[x,1]-i*tolm-tolm/tf
                ymx = lcarr[x,1]-i*tolm+tolm/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[lcarr[x,0],y] <= -2*ednum:
                        llab[-2*ednum-edgdet[lcarr[x,0],y]] = -(idnum+i)
                        edgdet[lcarr[x,0],y] = idnum+i
        update_labels(edgdet, llab, idnum, ni, -1, ednum)
        idnum += 1
    #
    # Now iterate through and label all of the right edges
    idnum = ednum+1
    while True:
        nc, tolp, tolm = get_xy(edgdet, rcarr, idnum, +1) # +1 specifies the direction to move in
        if nc == 0: break
        if tolp == 0 and tolm == 0: break
        for x in range(nc):
            for i in range(1,ni+1):
                ymn = rcarr[x,1]+i*tolp-tolp/tf
                ymx = rcarr[x,1]+i*tolp+tolp/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[rcarr[x,0],y] >= 2*ednum:
                        rlab[edgdet[rcarr[x,0],y]-2*ednum] = idnum+i
                        edgdet[rcarr[x,0],y] = idnum+i
                ymn = rcarr[x,1]-i*tolm-tolm/tf
                ymx = rcarr[x,1]-i*tolm+tolm/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[rcarr[x,0],y] >= 2*ednum:
                        rlab[edgdet[rcarr[x,0],y]-2*ednum] = idnum-i
                        edgdet[rcarr[x,0],y] = idnum-i
        update_labels(edgdet, rlab, idnum, ni, +1, ednum)
        idnum += 1
    #
    # Now search and label the right edges, except now search backwards from edge ednum
    idnum = ednum-1
    while True:
        nc, tolp, tolm = get_xy(edgdet, rcarr, idnum, -1) # -1 specifies the direction to move in
        if nc == 0: break
        if tolp == 0 and tolm == 0: break
        for x in range(nc):
            for i in range(1,ni+1):
                ymn = rcarr[x,1]+i*tolp-tolp/tf
                ymx = rcarr[x,1]+i*tolp+tolp/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[rcarr[x,0],y] >= 2*ednum:
                        rlab[edgdet[rcarr[x,0],y]-2*ednum] = idnum+i
                        edgdet[rcarr[x,0],y] = idnum+i
                ymn = rcarr[x,1]-i*tolm-tolm/tf
                ymx = rcarr[x,1]-i*tolm+tolm/tf+1
                if ymn < 0: ymn = 0
                if ymx > sz_y: ymx = sz_y
                for y in range(ymn,ymx):
                    if edgdet[rcarr[x,0],y] >= 2*ednum:
                        rlab[edgdet[rcarr[x,0],y]-2*ednum] = idnum-i
                        edgdet[rcarr[x,0],y] = idnum-i
        update_labels(edgdet, rlab, idnum, ni, +1, ednum)
        idnum -= 1
    #
    # Search through edgdet and remove any unidentified (spurious) identifications
    rmin = ednum
    rmax = ednum
    lmin = ednum
    lmax = ednum
    for x in range(sz_x):
        for y in range(sz_y):
            if edgdet[x,y] <= -2*ednum: edgdet[x,y] = 0
            elif edgdet[x,y] >= 2*ednum: edgdet[x,y] = 0
            else:
                if edgdet[x,y] < 0:
                    if -edgdet[x,y] < lmin: lmin = -edgdet[x,y]
                    elif -edgdet[x,y] > lmax: lmax = -edgdet[x,y]
                elif edgdet[x,y] > 0:
                    if edgdet[x,y] < rmin: rmin = edgdet[x,y]
                    elif edgdet[x,y] > rmax: rmax = edgdet[x,y]
    return lmin, lmax, rmin, rmax

#######
#  B  #
#######


#######
#  C  #
#######

@cython.boundscheck(False)
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


@cython.boundscheck(False)
def clean_pix(np.ndarray[DTYPE_t, ndim=2] array not None,
            double threshold):

    cdef int sz_x
    cdef int x, y, k, l
    cdef double tt, tp

    sz_x = array.shape[0]

    cdef int clspix = 5 # This should be passed in as an input parameter --- it is the number of nearby edges the trace algorithm has to reject before the trace is considered lost

    cdef np.ndarray[DTYPE_t, ndim=1] meddist = np.zeros(sz_x, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] clsdist = np.zeros(clspix, dtype=DTYPE)

    # Find the clspix closest pixels
    for x in range(0,sz_x):
        # Reset the pixel distance array
        for y in range(0,clspix):
            clsdist[y] = -999
        # Now find the closest pixels to the x'th pixel
        for y in range(0,sz_x):
            if x==y: continue
            distt = csqrt( <double>((array[x,0]-array[y,0])**2) + <double>((array[x,1]-array[y,1])**2) )
            tp = distt
            for k in range(clspix):
                if clsdist[k,0] == -999 and clsdist[k,1] == -999:
                    clsdist[k] = distt
                else:
                    distp = csqrt( <double>((array[x,0]-array[k,0])**2) + <double>((array[x,1]-array[k,1])**2) )
                    l = k
                    while distt <= distp and l<clspix:
                        tt = clsdist[l]
                        clsdist[l] = tp
                        tp = tt
                        l+=1
        meddist[x] = median(clsdist)

    # Calculate the typical distance between the median of the nearest pixels
    medval, madval = medianmad(meddist)
    # Erase any pixels that are outliers
    nop = sz_x
    x = 0
    while x < nop:
        if meddist[nop] >= medval + threshold*madval:
            nop -= 1
            array[x,0] = array[nop,0]
            array[x,1] = array[nop,1]
            array[nop,0] = -999
            array[nop,1] = -999
        else:
            x += 1
    # Finally return the cleaned array
    return array


@cython.boundscheck(False)
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


@cython.boundscheck(False)
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


#######
#  D  #
#######

@cython.boundscheck(False)
def detect_edges(np.ndarray[DTYPE_t, ndim=2] array not None):
    # array     : trace frame

    cdef int sz_x, sz_y
    cdef int x, y
    cdef int cl, dl, cr, dr, tl, tr
    cdef double difvl, difvc, difvr
    cdef int lcnt, rcnt

    sz_x = array.shape[0]
    sz_y = array.shape[1]

    # Set up the array which marks detections
    cdef np.ndarray[ITYPE_t, ndim=2] edgdet = np.zeros((array.shape[0],array.shape[1]), dtype=ITYPE)

    cdef int medsum = 7  # This should be an odd number
    cdef int ms = (medsum-1)/2
    cdef np.ndarray[DTYPE_t, ndim=2] medarr = np.zeros((sz_y,medsum), dtype=DTYPE)

    for x in range(ms,sz_x-ms):
        median_ed(array,medarr,x,ms)
        cl = 0
        dl = 0
        cr = 0
        dr = 0
        for y in range(1,sz_y-2):
            if medarr[y,ms] == 0.0 or medarr[y-1,ms] == 0.0 or medarr[y+1,ms] == 0.0: continue
            difvl = medarr[y,ms]-medarr[y-1,ms]
            difvc = medarr[y+1,ms]-medarr[y,ms]
            difvr = medarr[y+2,ms]-medarr[y+1,ms]
            # Detect left edges
            if difvc-difvl >= 0.0:
                cl += 1
            else:
                cl = 0
            if (dl == 0) and (cl >= 4) and (difvr-difvc < 0.0): # Candidate left edge detection
                cl = 0
                dl += 1
                tl = y
            elif dl > 0:
                if (difvr-difvc < 0.0): dl += 1
                else: dl = 0
            # Detect right edges
            if difvl-difvc >= 0.0:
                cr += 1
            else:
                cr = 0
            if (dr == 0) and (cr >= 4) and (difvc-difvr < 0.0): # Candidate right edge detection
                cr = 0
                dr += 1
                tr = y
            elif dr > 0:
                if (difvc-difvr < 0.0):
                    dr += 1
                else:
                    dr = 0
            if dl == 4:
                # A left edge detection!
                edgdet[x,tl] = -1
                dl = 0
            if dr == 4:
                # A right edge detection!
                edgdet[x,tr] = 1
                dr = 0
    return edgdet


@cython.boundscheck(False)
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

#######
#  E  #
#######

#@cython.boundscheck(False)
def edge_sum(np.ndarray[ITYPE_t, ndim=1] edghist not None,
            np.ndarray[ITYPE_t, ndim=1] sumarr not None):
    cdef int s, sz_s

    sz_s = sumarr.shape[0]

    # Find which slit edge id is the most common
    for s in range(sz_s):
        edghist[sumarr[s]] += 1
    return edghist


@cython.boundscheck(False)
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


#######
#  F  #
#######

@cython.boundscheck(False)
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


@cython.boundscheck(False)
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


@cython.boundscheck(False)
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


@cython.boundscheck(False)
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



@cython.boundscheck(False)
def fit_edges(np.ndarray[DTYPE_t, ndim=2] array not None,
                int polyorder, int lcnt, int rcnt):
    # array     : edges frame
    # polyorder : The order of the polynomial to use in fitting
    # lcnt      : Number of left orders
    # rcnt      : Number of right orders

    cdef int sz_x, sz_y
    cdef int x, y, o
    cdef int nothing
    cdef double nothingalso

    sz_x = array.shape[0]
    sz_y = array.shape[1]

    # Setup the coefficient arrays
    cdef np.ndarray[DTYPE_t, ndim=2] lcoeff = np.zeros((polyorder+1, lcnt), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] rcoeff = np.zeros((polyorder+1, rcnt), dtype=DTYPE)

    # Fit the left orders
#	for o in range(lcnt):
#		for x in range(sz_x):
#			for y in range(sz_y):

    # Fit the right orders
#	for o in range(rcnt):
#		for x in range(sz_x):
#			for y in range(sz_y):

    return lcoeff, rcoeff

#######
#  G  #
#######

@cython.boundscheck(False)
def get_closest(np.ndarray[DTYPE_t, ndim=2] orderpx not None,
                  int swlin, int edx, int edy):

    cdef int sz_x
    cdef int y, k, l
    cdef double ttx, tty, tpx, tpy

    sz_x = orderpx.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=2] clsdist = np.zeros((swlin,2), dtype=DTYPE)

    # Reset the pixel distance array
    for y in range(0,swlin):
        clsdist[y] = -999
    # Now find the closest pixels to the x'th pixel
    for y in range(0,sz_x):
        distt = csqrt( <double>((edx-orderpx[y,0])**2) + <double>((edy-orderpx[y,1])**2) )
        tpx = orderpx[y,0]
        tpy = orderpx[y,1]
        for k in range(swlin):
            if clsdist[k,0] == -999 and clsdist[k,1] == -999:
                clsdist[k,0] = tpx
                clsdist[k,1] = tpy
            else:
                distp = csqrt( <double>((edx-orderpx[k,0])**2) + <double>((edy-orderpx[k,1])**2) )
                l = k
                while distt <= distp and l<swlin:
                    ttx = clsdist[l,0]
                    tty = clsdist[l,1]
                    clsdist[l,0] = tpx
                    clsdist[l,1] = tpy
                    tpx = ttx
                    tpy = tty
                    l+=1
    return clsdist

@cython.boundscheck(False)
def get_edges(np.ndarray[DTYPE_t, ndim=2] array not None,
                  double threshold, int srchspe, int srchspa):
    # array     : trace frame
    # threshold : sigma detection threshold
    # srchspe   : Number of pixels to search along the spectral direction
    # srchspa   : Number of pixels to search along the spatial direction (must be an odd number)

    cdef int sz_x, sz_y
    cdef int x, y, i, j, k, l, nx, ny
    cdef int trc, cmax, rmin, rmax
    cdef double difval, errval
    cdef int ival, cval, icnt
    cdef double dval, tval
    cdef int ttx, tty, tpx, tpy
    cdef int nop, isdev
    cdef double distt, distp
    cdef int pedge, medge
    sz_x = array.shape[0]
    sz_y = array.shape[1]

    cdef int trcign = 5 # This should be passed in as an input parameter --- it is the number of nearby edges the trace algorithm has to reject before the trace is considered lost
    cdef int swlin = 15 # This should (probably) be passed in as an input parameter --- it is the number of nearby edges the trace algorithm needs to accept before doing linear regression
    cdef double ordthres = 4.0 # This should be passed in as an input parameter --- it is the number of sigma a pixel needs to be away from the other pixels before it is rejected

    cdef np.ndarray[ITYPE_t, ndim=2] edgearr = np.zeros((sz_x,sz_y), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] edgenby = np.zeros((trcign,2), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] orderpx = np.zeros((1,2), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] closept = np.zeros((swlin,2), dtype=ITYPE)

    # start by masking the nearby edges array
    for i in range(trcign):
        edgenby[i,0] = -999
        edgenby[i,1] = -999

    for x in range(sz_x):
        ival = -1
        cval = -1
        icnt = 0
        dval = 0.0
        for y in range(1,sz_y):
            if array[x,y] == 0.0 or array[x,y-1] == 0.0: continue
            difval = array[x,y]-array[x,y-1]
            errval = array[x,y]+array[x,y-1]
            if errval > 0.0:
                errval = csqrt(errval)
            else:
                errval = csqrt(-errval)
            if difval > 0.0:
                if dval < 0.0:
                    if ival != -1 and icnt > 3: edgearr[x,ival] = -999
                    ival = -1
                    cval = -1
                    icnt = 0
                    dval = 0.0
                if difval/errval > threshold:
                    if cval == -1:
                        ival = y
                        icnt = 1
                        dval = difval
                    elif cval+1 == y:
                        if difval > dval:
                            dval = difval
                            ival = y
                        icnt += 1
                    else:
                        if icnt > 3: edgearr[x,ival] = 999
                        ival = y
                        icnt = 1
                        dval = difval
                    cval = y
                else:
                    if cval != -1 and ival != -1 and icnt > 3:
                        edgearr[x,ival] = 999
                    cval = -1
                    icnt = 0
                    dval = 0.0
            else:
                if dval > 0.0:
                    if ival != -1 and icnt > 3: edgearr[x,ival] = 999
                    ival = -1
                    cval = -1
                    icnt = 0
                    dval = 0.0
                if -difval/errval > threshold:
                    if cval == -1:
                        ival = y
                        icnt = 1
                        dval = difval
                    elif cval+1 == y:
                        if difval < dval:
                            dval = difval
                            ival = y
                        icnt += 1
                    else:
                        edgearr[x,ival] = -999
                        ival = y
                        icnt = 0
                        dval = difval
                    cval = y
                else:
                    if cval != -1 and ival != -1 and icnt > 3:
                        edgearr[x,ival] = -999
                    cval = -1
                    icnt = 0
                    dval = 0.0
    return edgearr
    # Match the edges to form groups of traces
    pedge = 1
    medge = -1
    for y in range(1,sz_y):
        for x in range(0,sz_x):
            if edgearr[x,y] == 0: continue
            elif edgearr[x,y] == 999:
                # Find nearby similar edges
                trc = 1 # Boolean, do I keep tracing?
                cmax = 1+srchspe # Go srchspe pixels beyond the current one
                nx = x
                while trc == 1:
                    if nx+cmax > sz_x: cmax = sz_x-nx
                    for i in range(1,cmax):
                        rmin = -(srchspa-1)/2
                        rmax = 1-rmin
                        if y+rmin < 1: rmin = 1 # 1 because there's nothing in the zeroth column
                        if y+rmax > sz_y: rmax = sz_y-y
                        bi = -1
                        bv = 0.0
                        for j in range(rmin,rmax):
                            if edgearr[nx+i,y+j] == 999:
                                if bi == -1:
                                    bi = y+j
                                    bv = array[nx+i,y+j]-array[nx,y]
                                    if bv < 0.0: bv = -bv
                                else:
                                    tval = array[nx+i,y+j]-array[nx,y]
                                    if tval < 0.0: tval = -tval
                                    if tval < bv:
                                        bi = y+j
                                        bv = tval
                        if bi != -1: # We're continuing the trace
                            edgearr[nx+i,bi] = pedge
                            nx += i
                            break
                    if bi == -1: # We've lost the trace!
                        #pedge += 1
                        trc = 0
            elif edgearr[x,y] == -999:
                #edgearr[x,y] = medge
                trc = 1 # Boolean, do I keep tracing?
                cmax = 1+srchspe # Go srchspe pixels beyond the current one
                nx = x
                while trc == 1:
                    if nx+cmax > sz_x: cmax = sz_x-nx
                    for i in range(1,cmax):
                        rmin = -(srchspa-1)/2
                        rmax = 1-rmin
                        if y+rmin < 1: rmin = 1 # 1 because there's nothing in the zeroth column
                        if y+rmax > sz_y: rmax = sz_y-y
                        bi = -1
                        bv = 0.0
                        for j in range(rmin,rmax):
                            if edgearr[nx+i,y+j] == -999:
                                if bi == -1:
                                    bi = y+j
                                    bv = array[nx+i,y+j]-array[nx,y]
                                    if bv < 0.0: bv = -bv
                                else:
                                    tval = array[nx+i,y+j]-array[nx,y]
                                    if tval < 0.0: tval = -tval
                                    if tval < bv:
                                        bi = y+j
                                        bv = tval
                        if bi != -1: # We're continuing the trace
                            edgearr[nx+i,bi] = pedge
                            nx += i
                            break
                    if bi == -1: # We've lost the trace!
                        #pedge += 1
                        trc = 0
            else: continue
        # Now trace in the other direction!
        for x in range(0,sz_x):
            if edgearr[sz_x-x-1,y] == 0: continue
            elif edgearr[sz_x-x-1,y] == 999:
                # Find nearby similar edges
                edgearr[sz_x-x-1,y] = pedge
                trc = 1 # Boolean, do I keep tracing?
                cmax = 1+srchspe # Go srchspe pixels beyond the current one
                nx = sz_x-x-1
                while trc == 1:
                    if nx-cmax < 0: cmax = 0
                    for i in range(cmax,sz_x-x-1):
                        rmin = -(srchspa-1)/2
                        rmax = 1-rmin
                        if y+rmin < 1: rmin = 1 # 1 because there's nothing in the zeroth column
                        if y+rmax > sz_y: rmax = sz_y-y
                        bi = -1
                        bv = 0.0
                        for j in range(rmin,rmax):
                            if edgearr[nx+i,y+j] == 999:
                                if bi == -1:
                                    bi = y+j
                                    bv = array[nx-i,y+j]-array[nx,y]
                                    if bv < 0.0: bv = -bv
                                else:
                                    tval = array[nx-i,y+j]-array[nx,y]
                                    if tval < 0.0: tval = -tval
                                    if tval < bv:
                                        bi = y+j
                                        bv = tval
                        if bi != -1: # We're continuing the trace
                            edgearr[nx-i,bi] = pedge
                            nx -= i
                            break
                    if bi == -1: # We've lost the trace!
                        pedge += 1
                        trc = 0
            elif edgearr[sz_x-x-1,y] == -999:
                edgearr[sz_x-x-1,y] = medge
                trc = 1 # Boolean, do I keep tracing?
                cmax = 1+srchspe # Go srchspe pixels beyond the current one
                nx = sz_x-x-1
                while trc == 1:
                    if nx-cmax < 0: cmax = 0
                    for i in range(cmax,sz_x-x-1):
                        rmin = -(srchspa-1)/2
                        rmax = 1-rmin
                        if y+rmin < 1: rmin = 1 # 1 because there's nothing in the zeroth column
                        if y+rmax > sz_y: rmax = sz_y-y
                        bi = -1
                        bv = 0.0
                        for j in range(rmin,rmax):
                            if edgearr[nx-i,y+j] == -999:
                                if bi == -1:
                                    bi = y+j
                                    bv = array[nx-i,y+j]-array[nx,y]
                                    if bv < 0.0: bv = -bv
                                else:
                                    tval = array[nx-i,y+j]-array[nx,y]
                                    if tval < 0.0: tval = -tval
                                    if tval < bv:
                                        bi = y+j
                                        bv = tval
                        if bi != -1: # We're continuing the trace
                            edgearr[nx-i,bi] = medge
                            nx -= i
                            break
                    if bi == -1: # We've lost the trace!
                        medge -= 1
                        trc = 0
            else: continue
    return edgearr


@cython.boundscheck(False)
def get_xy(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
            np.ndarray[ITYPE_t, ndim=2] coords not None,
            int idnum, int mvdir):
    # Obtain the (x,y) values for idnum
    cdef int x, y, sz_x, sz_y
    cdef int nc, fcm, fcp, tolp, tolm, flgp, flgm, yvm, yvp, yv

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=1] medarrp = np.zeros(coords.shape[0], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] medarrm = np.zeros(coords.shape[0], dtype=DTYPE)

    nc = 0
    fcm = 0
    fcp = 0
    for x in range(sz_x):
        flgm = 0
        flgp = 0
        for y in range(sz_y):
            if edgdet[x,y] == idnum+mvdir:
                yvp = y
                flgp += 1
            elif edgdet[x,y] == idnum-mvdir:
                yvm = y
                flgm += 1
            elif edgdet[x,y] == idnum and flgp==0 and flgm==0:
                coords[nc,0] = x
                coords[nc,1] = y
                yv = y
                flgp += 1
                flgm += 1
                nc += 1
        if flgp == 2:
            medarrp[fcp] = <double>(yvp-yv)
            if medarrp[fcp] < 0.0: medarrp[fcp] *= -1.0
            fcp += 1
        if flgm == 2:
            medarrm[fcm] = <double>(yvm-yv)
            if medarrm[fcm] < 0.0: medarrm[fcm] *= -1.0
            fcm += 1
    if fcp != 0: tolp = <int>(medianmaskmin(medarrp,0.0))
    else: tolp = 0
    if fcm != 0: tolm = <int>(medianmaskmin(medarrm,0.0))
    else: tolm = 0
    return nc, tolp, tolm


#######
#  H  #
#######


#######
#  I  #
#######

@cython.boundscheck(False)
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
def label_orders_two(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                int lcnt, int rcnt, int ednum):
    cdef int sz_x, sz_y, cnt
    cdef int x, y, c, j
    cdef double lvcnt, rvcnt, lncnt, rncnt
    cdef double tempa, tempb

    lcnt = lcnt - 2*ednum
    rcnt = rcnt - 2*ednum

    if lcnt > rcnt: cnt = lcnt
    else: cnt = rcnt

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] larr = np.zeros((2,lcnt), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] rarr = np.zeros((2,rcnt), dtype=DTYPE)

    for c in range(2*ednum,2*ednum+cnt):
        lvcnt = 0.0
        rvcnt = 0.0
        lncnt = 0.0
        rncnt = 0.0
        for x in range(sz_x):
            for y in range(sz_y):
                if edgdet[x,y] == c:
                    rvcnt += <double>(y)
                    rncnt += 1.0
                elif edgdet[x,y] == -c:
                    lvcnt += <double>(y)
                    lncnt += 1.0
        if lncnt != 0.0:
            larr[0,c-2*ednum] = <double>(c)
            larr[1,c-2*ednum] = lvcnt/lncnt
        if rncnt != 0.0:
            rarr[0,c-2*ednum] = <double>(c)
            rarr[1,c-2*ednum] = rvcnt/rncnt
    print "completed first loop"
    # Now sort according to average order distance from one edge
    for x in range(lcnt-1):
        for y in range(x+1,lcnt):
            if larr[1,y] < larr[1,x]:
                tempa = larr[0,x]
                tempb = larr[1,x]
                larr[0,x] = larr[0,y]
                larr[1,x] = larr[1,y]
                larr[0,y] = tempa
                larr[1,y] = tempb
    for x in range(rcnt-1):
        for y in range(x+1,rcnt):
            if rarr[1,y] < rarr[1,x]:
                tempa = rarr[0,x]
                tempb = rarr[1,x]
                rarr[0,x] = rarr[0,y]
                rarr[1,x] = rarr[1,y]
                rarr[0,y] = tempa
                rarr[1,y] = tempb
    print "applying orders"
    # Finally, apply these new numbers to the orders
    for c in range(2*ednum,2*ednum+cnt):
        for x in range(sz_x):
            for y in range(sz_y):
                if edgdet[x,y] == c:
                    for j in range(rcnt):
                        if rarr[0,j] == <double>(c):
                            edgdet[x,y] = j+1
                            break
                elif edgdet[x,y] == -c:
                    for j in range(lcnt):
                        if larr[0,j] == <double>(c):
                            edgdet[x,y] = -j-1
                            break
    return


@cython.boundscheck(False)
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




#(closept,edgenby[trc,0],edgenby[trc,1],ordthres)
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
def medianmaskmin(np.ndarray[DTYPE_t, ndim=1] array not None, double mask):
    # This routine assumes that the masked value is also the minimum in the array
    cdef int sz_x
    cdef int x, j, nm
    cdef double temp

    sz_x = array.shape[0]

    for x in range(sz_x-1):
        for j in range(x+1,sz_x):
            if array[j] < array[x]:
                temp = array[x]
                array[x] = array[j]
                array[j] = temp
    # Count the number of masked values in the array
    nm = 0
    for x in range(sz_x):
        if array[x] == mask:
            nm += 1
        else: break
    # Find and return the median value
    if (sz_x-nm)%2==0:
        return 0.5*(array[(sz_x+nm-2)/2] + array[(sz_x+nm)/2])
    else:
        return array[(sz_x+nm-1)/2]


@cython.boundscheck(False)
def median_ed(np.ndarray[DTYPE_t, ndim=2] array not None,
            np.ndarray[DTYPE_t, ndim=2] medarr not None,
            int x, int ms):
    cdef int sz_y
    cdef int y, i, j
    cdef double temp

    sz_y = array.shape[1]

    # Extract the relevant pieces of the array
    for y in range(sz_y):
        for i in range(2*ms+1):
            medarr[y,i] = array[x+i-ms,y]

    # Sort the array
    for y in range(sz_y-1):
        for i in range(2*ms):
            for j in range(i+1,2*ms+1):
                if medarr[y,j] < medarr[y,i]:
                    temp = medarr[y,i]
                    medarr[y,i] = medarr[y,j]
                    medarr[y,j] = temp
    return
#	# Find and return the median value
#	if sz_x%2==0:
#		return 0.5*(array[sz_x/2] + array[sz_x/2 - 1])
#	else:
#		return array[(sz_x-1)/2]


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


@cython.boundscheck(False)
def polyfit(np.ndarray[DTYPE_t, ndim=1] x not None,
        np.ndarray[DTYPE_t, ndim=1] y not None,
        int degree):

    cdef int sz_x
    cdef int i, j
    cdef double chisq

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat
    cdef gsl_vector *yvec
    cdef gsl_vector *c

    sz_x = x.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(degree, dtype=DTYPE)

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
    return coeffs


@cython.boundscheck(False)
def polyfit_mask(np.ndarray[DTYPE_t, ndim=1] xt not None,
                np.ndarray[DTYPE_t, ndim=1] yt not None,
                np.ndarray[DTYPE_t, ndim=1] coeffs not None,
                double mask):

    cdef int sz_x, sz_m, sz_c
    cdef int i, j
    cdef double chisq

    cdef gsl_multifit_linear_workspace *ws
    cdef gsl_matrix *cov
    cdef gsl_matrix *xmat
    cdef gsl_vector *yvec
    cdef gsl_vector *c

    # Prepare an array that masks out the appropriate values
    sz_m = xt.shape[0]
    sz_x = 0
    for i in range(0,sz_m):
        if xt[i] != mask:
            sz_x += 1
    # Create new arrays without the masked values
    cdef np.ndarray[DTYPE_t, ndim=1] x = np.zeros(sz_x, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] y = np.zeros(sz_x, dtype=DTYPE)
    # And fill them in
    sz_c = 0
    for i in range(0,sz_m):
        if xt[i] != mask:
            x[sz_c] = xt[i]
            y[sz_c] = yt[i]
            sz_c += 1

    # Define the coefficients array
    cdef int degree = coeffs.shape[0]

    # Check that the result is reasonable
    if sz_x <= degree:
        for i in range(degree):
            coeffs[i] = mask
        return

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
    return coeffs


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
        gsl_matrix_set(xmat, i, 0, 1.0)
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


@cython.boundscheck(False)
def prune_peaks(np.ndarray[ITYPE_t, ndim=1] hist not None,
                np.ndarray[ITYPE_t, ndim=1] pks not None,
                int pkidx):
    """
    Identify the most well defined peaks
    """

    cdef int ii, jj, sz_i, cnt, lgd

    sz_i = pks.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=1] msk = np.zeros(sz_i, dtype=ITYPE)

    lgd = 1  # Was the previously inspected peak a good one?
    for ii in range(0, sz_i-1):
        cnt = 0
        for jj in range(pks[ii], pks[ii+1]):
            if hist[jj] == 0:
                cnt += 1
        #if cnt < (pks[ii+1] - pks[ii])/2:
        if cnt < 2:
            # If the difference is unacceptable, both peaks are bad
            msk[ii] = 0
            msk[ii+1] = 0
            lgd = 0
        else:
            # If the difference is acceptable, the right peak is acceptable,
            # the left peak is acceptable if it was not previously labelled as unacceptable
            if lgd == 1:
                msk[ii] = 1
            msk[ii+1] = 1
            lgd = 1
    # Now only consider the peaks closest to the highest peak
    lgd = 1
    for ii in range(pkidx, sz_i):
        if msk[ii] == 0:
            lgd = 0
        elif lgd == 0:
            msk[ii] = 0
    lgd = 1
    for ii in range(0, pkidx):
        if msk[pkidx-ii] == 0:
            lgd = 0
        elif lgd == 0:
            msk[pkidx-ii] = 0
    return msk


#######
#  Q  #
#######

#######
#  R  #
#######

#######
#  S  #
#######

#######
#  T  #
#######

def trace_fweight(np.ndarray[DTYPE_t, ndim=2] fimage not None,
                  np.ndarray[DTYPE_t, ndim=2] invvar not None,
                  np.ndarray[DTYPE_t, ndim=1] xinit not None,
                  np.ndarray[DTYPE_t, ndim=1] ycen not None,
                  double radius):

    """
    Parameters
    ----------
    fimage: 2D ndarray
      Image for tracing
    invvar: ndarray, optional
      Inverse variance array for the image
    xinit: ndarray
      Initial guesses for x-trace
    ycen : ndarray
      Initial guesses for y-trace
    radius: float
      Radius for centroiding
    """

    cdef int nx, ny, fpix, spot
    cdef int y, x, sz_x, qbad, ix1, ix2
    cdef double xdiff, absxdiff, delx
    cdef double x1, x2
    cdef double sumw, sumwt, sumxw, var_term
    cdef double sumsx1, sumsx2

    # Init
    nx = fimage.shape[1]
    ny = fimage.shape[0]
    sz_x = xinit.shape[0]

    # Create xnew, xerr
    cdef np.ndarray[DTYPE_t, ndim=1] xnew = np.zeros(sz_x, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] xerr = -999999.9*np.ones(sz_x, dtype=DTYPE)

    for x in range(sz_x):
        x1 = xinit[x] - radius + 0.5
        x2 = xinit[x] + radius + 0.5
        ix1 = <int>(x1)
        ix2 = <int>(x2)
        # Apply a floor operation for negative values
        if x1 < 0.0:
            ix1 -= 1
        if x2 < 0.0:
            ix2 -= 1
        fpix = (ix2-ix1)-1
        if fpix < 0:
            fpix = 0
        # Zero the centering
        sumw = 0.0
        sumwt = 0.0
        sumxw = 0.0
        sumsx1 = 0.0
        sumsx2 = 0.0
        qbad = 0
        # Compute
        for y in range(0, fpix+3):
            spot = ix1 - 1 + y
            if spot < 0:
                ih = 0
            elif spot > nx-1:
                ih = nx-1
            else:
                ih = spot
            xdiff = <double>(spot) - xinit[x]
            if xdiff < 0.0:
                absxdiff = -1.0*xdiff
            else:
                absxdiff = xdiff
            wt = radius - absxdiff + 0.5
            if wt < 0.0:
                wt = 0.0
            elif wt > 1.0:
                wt = 1.0
            sumw += fimage[ycen[x], ih] * wt
            sumwt += wt
            sumxw += fimage[ycen[x], ih] * xdiff * wt
            if invvar[ycen[x], ih] == 0.0:
                var_term = wt**2
            else:
                var_term = wt**2 / invvar[ycen[x], ih]
            sumsx2 += var_term
            sumsx1 += xdiff**2 * var_term
            if invvar[ycen[x],ih] <= 0:
                qbad = 1

        # Fill up
        xnew[x] = xinit[x]
        if (sumw > 0.0) and (qbad == 0):
            delx = sumxw/sumw
            xnew[x] = delx + xinit
            xerr[x] = csqrt(sumsx1 + sumsx2*delx*delx)/sumw

        absxdiff = xnew[x]-xinit[x]
        if absxdiff < 0.0:
            absxdiff *= -1.0
        if (absxdiff > radius+0.5) or (xinit[x] < radius-0.5) or (xinit[x] > nx - 0.5 - radius):
            xnew[x] = xinit[x]
            xerr[x] = -999999.9

    # Return
    return xnew, xerr


@cython.boundscheck(False)
def tilts_image(np.ndarray[DTYPE_t, ndim=2] tilts not None,
                np.ndarray[DTYPE_t, ndim=2] lordloc not None,
                np.ndarray[DTYPE_t, ndim=2] rordloc not None,
                int sz_y):
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
            ymin = <int>(oc-ow)
            ymax = <int>(oc+ow)+1
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


@cython.boundscheck(False)
def trace_tilts(np.ndarray[DTYPE_t, ndim=2] msarc not None,
                np.ndarray[ITYPE_t, ndim=2] ordcen not None,
                np.ndarray[ITYPE_t, ndim=2] ordsiz not None,
                np.ndarray[ITYPE_t, ndim=2] arccen not None,
                int maxnum):
    """
    pseudo-maximum likelihood estimate for the trace of the arc lines.
    """
    cdef int sz_l, sz_o
    cdef int l, o, ll, x, y, b, bb, nfit
    cdef int xmin, xmax, ymin, ymax
    cdef double mval, wght, bval

    cdef double maskval = -999999.9

    sz_l = arccen.shape[0]
    sz_o = arccen.shape[1]

    cdef int sz_ax = msarc.shape[1]
    cdef int sz_ay = msarc.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=2] derv = np.zeros((sz_l,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] cent = np.zeros((sz_l,sz_o), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] xfit = np.zeros(maxnum, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] yfit = np.zeros(maxnum, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.zeros(2, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] back = np.zeros(2*maxnum, dtype=DTYPE)
    # Now loop through and update the array
    for o in range(sz_o):
        for l in range(sz_l):
            if arccen[l,o]==-1:
                for ll in range(l,sz_l):
                    derv[ll,o] = maskval
                    cent[ll,o] = maskval
                break
            # Find the local background level
            xmin = -ordsiz[arccen[l,o],o]
            xmax = 1+ordsiz[arccen[l,o],o]
            ymin = -4
            ymax = 1+4
            if xmin+ordcen[l,o] < 0: xmin = -ordcen[l,o]
            if xmax+ordcen[l,o] > sz_ax: xmax = sz_ax-ordcen[l,o]
            if ymin+arccen[l,o] < 0: ymin = -arccen[l,0]
            if ymax+arccen[l,o] > sz_ay: ymax = sz_ay-arccen[l,o]
            for x in range(xmin,xmax):
                for y in range(ymin,ymax):
                    if msarc[y,x] < back[2*maxnum-1] or back[2*maxnum-1] == -1.0:
                        for b in range(0,2*maxnum):
                            if back[b] == -1.0: back[b] = msarc[y,x]
                            elif back[b] > msarc[y,x]: # Rearrange
                                for bb in range(b+1,2*maxnum):
                                    back[2*maxnum+b-bb] = back[2*maxnum+b-bb-1]
                                back[b] = msarc[y,x]
                                break
            bval = 0.0
            for b in range(0,2*maxnum):
                bval += back[b]
            bval /= <double>(2*maxnum)
            # Find centre of this arc line
            mval = 0.0
            wght = 0.0
            nfit = 0
            ymin = arccen[l,o]-4
            ymax = arccen[l,o]+1+4
            if ymin < 0: ymin = 0
            if ymax > sz_ay: ymax = sz_ay
            for y in range(ymin,ymax):
                mval += <double>(y)*(msarc[y,ordcen[l,o]]-bval)
                wght += (msarc[y,ordcen[l,o]]-bval)
            xfit[0] = 0.0
            yfit[0] = mval/wght
            nfit += 1
            # Trace up
            xmin = 1
            xmax = 1+ordsiz[arccen[l,o],o]
            if ordcen[l,o]+xmin > sz_ax: xmin = sz_ax-ordcen[l,o]
            if ordcen[l,o]+xmax > sz_ax: xmax = sz_ax-ordcen[l,o]
            for x in range(xmin,xmax):
                mval = 0.0
                wght = 0.0
#				ymin = arccen[l,o]-ordsiz[arccen[l,o],o]
#				ymax = arccen[l,o]+1+ordsiz[arccen[l,o],o]
                ymin = <int>(yfit[nfit-1]-4.0+0.5)
                ymax = <int>(yfit[nfit-1]+1.0+4.0+0.5)
                if ymin < 0: ymin = 0
                if ymax > sz_ay: ymax = sz_ay
                for y in range(ymin,ymax):
                    mval += <double>(y)*(msarc[y,ordcen[l,o]+x]-bval)
                    wght += (msarc[y,ordcen[l,o]+x]-bval)
                xfit[nfit] = <double>(x)
                yfit[nfit] = mval/wght
                nfit += 1
            # Trace down
            xmin = 1
            xmax = 1+ordsiz[arccen[l,o],o]
            if ordcen[l,o]-xmin < 0: xmin = ordcen[l,o]
            for x in range(xmin,xmax):
                mval = 0.0
                wght = 0.0
                if x == xmin:
                    ymin = <int>(yfit[0]-4.0+0.5) # +/- 4 pixels, +0.5 so that int rounds
                    ymax = <int>(yfit[0]+1.0+4.0+0.5)
                else:
                    ymin = <int>(yfit[nfit-1]-4.0+0.5)
                    ymax = <int>(yfit[nfit-1]+1.0+4.0+0.5)
                if ymin < 0: ymin = 0
                if ymax > sz_ay: ymax = sz_ay
                for y in range(ymin,ymax):
                    mval += <double>(y)*(msarc[y,ordcen[l,o]-x]-bval)
                    wght += (msarc[y,ordcen[l,o]-x]-bval)
                xfit[nfit] = <double>(-x)
                yfit[nfit] = mval/wght
                nfit += 1
            # Mask out all the unused xfit/yfit
            for x in range(nfit,maxnum):
                xfit[x] = maskval
                yfit[x] = maskval
            # Fit the tilt for this arc line
            polyfit_mask(xfit,yfit,coeffs,maskval)
            # Assign the gradient to the output array
            derv[l,o] = coeffs[1]
            cent[l,o] = coeffs[0]
    return derv, cent


#######
#  U  #
#######

@cython.boundscheck(False)
def update_labels(np.ndarray[ITYPE_t, ndim=2] edgdet not None,
                np.ndarray[ITYPE_t, ndim=1] labels not None,
                int idnum, int ni, int sn, int ednum):
    cdef int sz_x, sz_y, sz_l
    cdef int x, y, l, nc

    sz_x = edgdet.shape[0]
    sz_y = edgdet.shape[1]
    sz_l = labels.shape[0]

    cdef np.ndarray[ITYPE_t, ndim=1] indx = np.zeros(sz_l, dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] indl = np.zeros(sz_l, dtype=ITYPE)

    # To be efficient, let's first single out what to look for
    nc = 0
    for l in range(sz_l):
        indx[l] = -1
        if labels[l] >= sn*idnum-ni and labels[l] <= sn*idnum+ni:
            indx[nc] = sn*(2*ednum+l)
            indl[nc] = sn*labels[l]
            nc += 1

    # Now loop through and update the array
    for x in range(sz_x):
        for y in range(sz_y):
            for l in range(nc):
                if edgdet[x,y] == indx[l]:
                    edgdet[x,y] = indl[l]
    return


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


