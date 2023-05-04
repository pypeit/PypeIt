""" Module for finding patterns in arc line spectra
"""
import numpy as np
import scipy.ndimage

#import numba as nb


def detect_2Dpeaks(image):
    """
    Takes a 2D image and returns 1 if a local maximum is found, and 0 otherwise

    Parameters
    ----------
    image : ndarray
        2D image

    Returns
    -------
    pimage : ndarray
        boolean mask of the peaks (1 when the pixel's value is the
        neighborhood maximum, 0 otherwise)

    """
    # Define an 8-connected neighborhood
    neighborhood = scipy.ndimage.generate_binary_structure(2, 2)

    # Apply the local maximum filter
    local_max = scipy.ndimage.maximum_filter(image, footprint=neighborhood) == image

    # Background mask
    background = (image == 0)

    # Remove artifacts from local maximum filter
    eroded_background = scipy.ndimage.binary_erosion(
        background, structure=neighborhood, border_value=1
    )

    # Remove the background from the local_max mask
    pimage = local_max ^ eroded_background

    return pimage


def match_quad_to_list(spec_lines, line_list, wv_guess, dwv_guess,
                  tol=2., dwv_uncertainty=0.2, min_ftol=0.005):
    """
    Parameters
    ----------
    spec_lines : ndarray
        pixel space
    line_list :
    tol :
    min_ftol : float, optional
        Minimum tolerance for matching 

    Returns
    -------
    possible_matches : list
        list of indices of matching quads

    """
    # Setup spec
    npix = spec_lines[-1]-spec_lines[0]
    spec_values = (spec_lines[1:-1]-spec_lines[0])/(
        spec_lines[-1]-spec_lines[0])
    ftol = max(tol/npix, min_ftol)
    #
    possible_start = np.where((line_list > wv_guess[0]) & (line_list < wv_guess[1]))[0]
    possible_matches = []
    for start in possible_start:
        #print("Trying {:g}".format(line_list[start]))
        # Find possible ends
        possible_ends = np.where( (line_list > line_list[start] + npix*dwv_guess*(1-dwv_uncertainty)) &
                                     (line_list < line_list[start] + npix*dwv_guess*(1+dwv_uncertainty)))[0]
        # Loop on ends
        for end in possible_ends:
            values = (line_list[start+1:end]-line_list[start]) / (
                line_list[end]-line_list[start])
            # Test
            diff0 = np.abs(values-spec_values[0])
            tst0 = diff0 < ftol
            diff1 = np.abs(values-spec_values[1])
            tst1 = diff1 < ftol
            #if np.abs(line_list[start]-6097.8) < 0.2:
            #    debugger.set_trace()
            if np.any(tst0) & np.any(tst1):
                i0 = np.argmin(diff0)
                i1 = np.argmin(diff1)
                #if np.sum(tst0) > 1:
                #    pdb.set_trace()
                #possible_matches.append([start, start+1+np.where(tst0)[0][0],
                #                         start+1+np.where(tst1)[0][0], end])
                possible_matches.append([start, start+1+i0, start+1+i1, end])
    # Return
    return possible_matches


def run_quad_match(tcent, twave, llist_wv, disp, swv_uncertainty=250., pix_tol=1.):
    """
    Parameters
    ----------
    tcent : ndarray
        Pixel positions of arc lines
    twave : ndarray
        Crude guess at wavelength solution, e.g. from wvcen, disp
    llist_wv : ndarray
        Lines to match against (from a line list)
    pix_tol : float
        Tolerance in units of pixels to match to

    Returns
    -------
    match_idx : dict
        Record of matches
    scores : ndarray
        str array of scores
    """

    # Init
    nlin = tcent.size
    match_idx = {}
    for ii in range(nlin):
        match_idx[ii] = {}
        match_idx[ii]['matches'] = []

    # Run -- quad
    for idx in range(nlin-4):
        for jj in range(4):
            sub_idx = idx + np.arange(5).astype(int)
            msk = np.array([True]*5)
            msk[jj+1] = False
            # Setup
            sidx = sub_idx[msk]
            spec_lines = np.array(tcent)[sidx]
            #
            widx = int(np.round(tcent[idx]))
            wvmnx = [twave[widx]-swv_uncertainty, twave[widx]+swv_uncertainty]
            # Run
            #import pdb; pdb.set_trace()
            matches = match_quad_to_list(spec_lines, llist_wv, wvmnx, disp, tol=pix_tol)
            # Save
            for match in matches:
                for ii in range(4):
                    match_idx[sidx[ii]]['matches'].append(match[ii])
            #pdb.set_trace()
    # Score
    scores = score_quad_matches(match_idx)
    scores = np.array(scores)

    # Return
    return match_idx, scores


def scan_for_matches(wvcen, disp, npix, cut_tcent, wvdata, best_dict=None,
                     swv_uncertainty=350., wvoff=1000., pix_tol=2., ampl=None):
    """

    .. warning::

        best_dict is updated in place

    Parameters
    ----------
    wvcen : float
        Guess at central wavelength
    disp : float
    npix
    cut_tcent
    wvdata
    best_dict
    swv_uncertainty
    wvoff
    pix_tol

    """

    # Setup
    #wvoff=10.
    #pdb.set_trace()
    dcen = swv_uncertainty*0.8
    wvcens = np.arange(wvcen-wvoff, wvcen+wvoff+dcen, dcen)
    # Best
    if best_dict is None:
        best_dict = dict(nmatch=0, ibest=-1, bwv=0.)

    # Scan on wv_cen
    #wvcens = [9400.] # DEBUGGING
    for ss,iwv_cen in enumerate(wvcens):
        # Wavelength array
        wave = iwv_cen + (np.arange(npix) - npix/2.)*disp
        match_idx, scores = run_quad_match(cut_tcent, wave, wvdata, disp,
                                           swv_uncertainty=swv_uncertainty,
                                           pix_tol=pix_tol)
        # Score
        mask = np.array([False]*len(cut_tcent))
        IDs = []
        for kk,score in enumerate(scores):
            if score in ['Perf', 'Good', 'Ok']:
                mask[kk] = True
                uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
                imx = np.argmax(counts)
                IDs.append(wvdata[uni[imx]])
            else:
                IDs.append(0.)
        ngd_match = np.sum(mask)
        #import pdb; pdb.set_trace()
        # Update in place
        if ngd_match > best_dict['nmatch']:
            best_dict['nmatch'] = ngd_match
            best_dict['midx'] = match_idx
            best_dict['mask'] = mask
            best_dict['scores'] = scores
            best_dict['ibest'] = ss
            best_dict['bwv'] = iwv_cen
            best_dict['IDs'] = IDs
            # Search parameters
            best_dict['swv_uncertainty'] = swv_uncertainty
            best_dict['wvoff'] = wvoff
            best_dict['pix_tol'] = pix_tol
            best_dict['ampl'] = ampl


def score_quad_matches(fidx):
    """
    Grades quad_match results

    Parameters
    ----------
    fidx

    Returns
    -------
    scores : list

    """
    # Loop on indices
    scores = []
    for key in fidx.keys():
        if len(fidx[key]['matches']) == 0:
            scores.append('None')
            continue
        matches = np.array(fidx[key]['matches'])
        nmatch = matches.size
        uni, counts = np.unique(matches, return_counts=True)
        nuni = len(uni)
        max_counts = max(counts)
        # Score
        if (nuni==1) & (nmatch >= 4):
            scores.append('Perf')
        elif nmatch == 0:
            scores.append('None')
        elif (max_counts == 4) & (nmatch == 5):
            scores.append('Good')
        elif (max_counts/nmatch >= 2./3) & (nmatch >= 6):
            scores.append('Good')
        elif (nuni == 1) & (nmatch == 3):
            scores.append('Good')
        elif (max_counts == 3) & (nmatch == 4):
            scores.append('OK')
        elif (nuni == 1) & (nmatch == 2):
            scores.append('Risk')
        else:
            scores.append('Amb')
    # Return
    return scores


#@nb.jit(nopython=True, cache=True)
def triangles(detlines, linelist, npixels, detsrch=5, lstsrch=10, pixtol=1.0):
    """

    Brute force pattern recognition using triangles. A triangle contains
    (for either detlines or linelist):

        1. a starting point (s),
        2. an end point (e), and
        3. something in between (b)

    Something like this::

        >                    |
        >                    |      |
        >    |               |      |
        >    |               |      |
        >    s               b      e

    Then, the value (b-s)/(e-s) is in the same coordinate system for
    both detlines and linelist.

    Parameters
    ----------
    detlines : ndarray
        list of detected lines in pixels (sorted, increasing)
    linelist : ndarray
        list of lines that should be detected (sorted, increasing)
    npixels : float
        Number of pixels along the dispersion direction
    detsrch : int
        Number of consecutive elements in detlines to use to create a pattern (-1 means all lines in detlines)
    lstsrch : int
        Number of consecutive elements in linelist to use to create a pattern (-1 means all lines in detlines)
    pixtol : float
        tolerance that is used to determine if a match is successful (in units of pixels)

    Returns
    -------
    dindex : ndarray
        Index array of all detlines used in each triangle
    lindex : ndarray
        Index array of the assigned line to each index in dindex
    wvcen : ndarray
        central wavelength of each triangle
    disps : ndarray
        Dispersion of each triangle (angstroms/pixel)
    """

    nptn = 3  # Number of lines used to create a pattern

    sz_d = detlines.size
    sz_l = linelist.size

    # Count the number of detlines patterns that will be created
    cntdet = 0
    dup = 0
    for d in range(detsrch-nptn+1):
        dup += d+1
        if d == detsrch-nptn:
            cntdet += dup*(sz_d-detsrch+1)
        else:
            cntdet += dup

    # Count the number of linelist patterns that will be created
    cntlst = 0
    lup = 0
    for l in range(lstsrch-nptn+1):
        lup += l+1
        if l == lstsrch-nptn:
            cntlst += lup*(sz_l-lstsrch+1)
        else:
            cntlst += lup

    lindex = np.zeros((cntdet*cntlst, nptn),dtype=np.uint64) #, dtype=nb.types.uint64)
    dindex = np.zeros((cntdet*cntlst, nptn),dtype=np.uint64) #, dtype=nb.types.uint64)
    wvcen = np.zeros((cntdet*cntlst))
    disps = np.zeros((cntdet*cntlst))

    # Test each detlines combination
    cntdet = 0
    for d in range(0, sz_d-nptn+1):      # d is the starting point of the pattern
        dup = d + detsrch
        if dup > sz_d:
            dup = sz_d
        if detsrch == -1:
            dup = sz_d
        for dd in range(d+nptn-1, dup):  # dd is the end point of the pattern
            for xd in range(d+1, dd):  # xd is the mid point of the pattern
                # Create the test pattern
                dval = (detlines[xd]-detlines[d])/(detlines[dd]-detlines[d])
                tol = pixtol/(detlines[dd]-detlines[d])
                # Search through all possible patterns in the linelist
                for l in range(0, sz_l-nptn+1):
                    lup = l + lstsrch
                    if lup > sz_l:
                        lup = sz_l
                    if lstsrch == -1:
                        lup = sz_l
                    for ll in range(l+nptn-1, lup):
                        for xl in range(l+1, ll):
                            lval = (linelist[xl]-linelist[l])/(linelist[ll]-linelist[l])
                            tst = lval-dval
                            if tst < 0.0:
                                tst *= -1.0
                            if tst <= tol:
                                lindex[cntdet, 0] = l
                                lindex[cntdet, 1] = xl
                                lindex[cntdet, 2] = ll
                                dindex[cntdet, 0] = d
                                dindex[cntdet, 1] = xd
                                dindex[cntdet, 2] = dd
                                tst = (linelist[ll]-linelist[l]) / (detlines[dd]-detlines[d])
                                wvcen[cntdet] = (npixels/2.0) * tst + (linelist[ll]-tst*detlines[dd])
                                disps[cntdet] = tst
                            cntdet += 1
    return dindex, lindex, wvcen, disps


#@nb.jit(nopython=True, cache=True)
def quadrangles(detlines, linelist, npixels, detsrch=5, lstsrch=10, pixtol=1.0):
    """

    Brute force pattern recognition using quadrangles. A quadrangle
    contains (for either detlines or linelist):

        1. a left line (l),
        2. a right line (r), and
        3. two lines in between (a, b)

    Something like this::

        >                   |
        >             |     |      |
        >   |         |     |      |
        >   |         |     |      |
        >   l         a     b      r

    Then, the values (a-ll)/(r-ll) and (b-ll)/(r-ll) are in the same
    coordinate system for both detlines and linelist.

    Parameters
    ----------
    detlines : ndarray
        list of detected lines in pixels (sorted, increasing)
    linelist : ndarray
        list of lines that should be detected (sorted, increasing)
    npixels : float
        Number of pixels along the dispersion direction
    detsrch : int
        Number of consecutive elements in detlines to use to create a pattern (-1 means all lines in detlines)
    lstsrch : int
        Number of consecutive elements in linelist to use to create a pattern (-1 means all lines in detlines)
    pixtol : float
        tolerance that is used to determine if a match is successful (in units of pixels)

    Returns
    -------
    dindex : ndarray
        Index array of all detlines used in each triangle
    lindex : ndarray
        Index array of the assigned line to each index in dindex
    wvcen : ndarray
        central wavelength of each triangle
    disps : ndarray
        Dispersion of each triangle (angstroms/pixel)
    """

    nptn = 4  # Number of lines used to create a pattern

    sz_d = detlines.size
    sz_l = linelist.size

    lindex = np.zeros((1, nptn),dtype=np.uint64) #, dtype=nb.types.uint64)
    dindex = np.zeros((1, nptn),dtype=np.uint64) #, dtype=nb.types.uint64)
    wvcen = np.zeros((1,), dtype=np.ulonglong) #, dtype=nb.types.ulong)
    disps = np.zeros((1,), dtype=np.ulonglong) #, dtype=nb.types.ulong)

    # Generate the patterns
    for dl in range(0, sz_d-nptn+1):  # dl is the starting point of the detlines pattern
        dup = dl + detsrch
        if dup > sz_d:
            dup = sz_d
        if detsrch == -1:
            dup = sz_d
        for dr in range(dl+nptn-1, dup):  # dr is the end point of the detlines pattern
            # Set the tolerance
            tol = pixtol / (detlines[dr] - detlines[dl])
            for da in range(dl+1, dr-1):  # da is the left mid point of the detlines pattern
                # Create the test pattern
                daval = (detlines[da]-detlines[dl])/(detlines[dr]-detlines[dl])
                for db in range(da+1, dr):  # db is the right mid point of the detlines pattern
                    # Create the test pattern
                    dbval = (detlines[db]-detlines[dl])/(detlines[dr]-detlines[dl])
                    # Search through all possible patterns in the linelist
                    for ll in range(0, sz_l-nptn+1):  # ll is the start point of the linelist pattern
                        lup = ll + lstsrch
                        if lup > sz_l:
                            lup = sz_l
                        if lstsrch == -1:
                            lup = sz_l
                        for lr in range(ll+nptn-1, lup):  # lr is the end point of the linelist pattern
                            for la in range(ll+1, lr-1):  # la is the end point of the linelist pattern
                                laval = (linelist[la] - linelist[ll]) / (linelist[lr] - linelist[ll])
                                tst = laval - daval
                                if tst < 0.0:
                                    tst *= -1.0
                                if tst <= tol:
                                    # The first pattern matches, check the second one.
                                    for lb in range(la+1, lr):  # la is the end point of the linelist pattern
                                        lbval = (linelist[lb] - linelist[ll]) / (linelist[lr] - linelist[ll])
                                        tst = lbval - dbval
                                        if tst < 0.0:
                                            tst *= -1.0
                                        if tst <= tol:
                                            # The second pattern matches, store the result!
                                            lindex = np.vstack((lindex, np.array([[ll, la, lb, lr]], dtype=np.uint64))) #, dtype=nb.types.uint64)))
                                            dindex = np.vstack((dindex, np.array([[dl, da, db, dr]], dtype=np.uint64))) #, dtype=nb.types.uint64)))
                                            tst = (linelist[lr] - linelist[ll]) / (detlines[dr] - detlines[dl])
                                            wvl = (npixels/2.)*tst + (linelist[lr]-tst*detlines[dr])
                                            wvcen = np.hstack((wvcen, np.array([wvl], dtype=np.ulonglong))) #, dtype=nb.types.ulong)))
                                            disps = np.hstack((disps, np.array([tst], dtype=np.ulonglong))) # , dtype=nb.types.ulong)))
    # Return, but first remove the spurious first entry due to array creation
    return dindex[1:, :], lindex[1:, :], wvcen[1:], disps[1:]


#@nb.jit(nopython=True, cache=True)
def curved_quadrangles(detlines, linelist, npixels, detsrch=5, lstsrch=10, pixtol=1.0):
    """

    Brute force pattern recognition using curved quadrangles.  A curved
    quadrangle contains (for either detlines or linelist):

        1. a left line (l),
        2. a right line (r),
        3. a mid line (m), and
        4. one line in between (c; c != m)

    Something like this:

        >                   |
        >             |     |      |
        >   |         |     |      |
        >   |         |     |      |
        >   l         c     m      r

    Then, the values (c-r)/(r-l) are in the same coordinate system for
    both detlines and linelist.

    Parameters
    ----------
    detlines : ndarray
        list of detected lines in pixels (sorted, increasing)
    linelist : ndarray
        list of lines that should be detected (sorted, increasing)
    npixels : float
        Number of pixels along the dispersion direction
    detsrch : int
        Number of consecutive elements in detlines to use to create a pattern (-1 means all lines in detlines)
    lstsrch : int
        Number of consecutive elements in linelist to use to create a pattern (-1 means all lines in detlines)
    pixtol : float
        tolerance that is used to determine if a match is successful (in units of pixels)

    Returns
    -------
    dindex : ndarray
        Index array of all detlines used in each triangle
    lindex : ndarray
        Index array of the assigned line to each index in dindex
    wvcen : ndarray
        central wavelength of each triangle
    disps : ndarray
        Dispersion of each triangle (angstroms/pixel)
    """

    nptn = 4  # Number of lines used to create a pattern

    sz_d = detlines.size
    sz_l = linelist.size

    lindex = np.zeros((1, nptn), dtype=np.uint64) #, dtype=nb.types.uint64)
    dindex = np.zeros((1, nptn), dtype=np.uint64) #, dtype=nb.types.uint64)
    wvcen = np.zeros((1,), dtype=np.ulonglong) #, dtype=nb.types.ulong)
    disps = np.zeros((1,), dtype=np.ulonglong) #, dtype=nb.types.ulong)

    # Generate the patterns
    for dl in range(0, sz_d-nptn+1):  # dl is the starting point of the detlines pattern
        dup = dl + detsrch
        if dup > sz_d:
            dup = sz_d
        if detsrch == -1:
            dup = sz_d
        for dr in range(dl+nptn-1, dup):  # dr is the end point of the detlines pattern
            # Set the tolerance
            tol = pixtol / (detlines[dr] - detlines[dl])
            for dm in range(dl+1, dr):  # dm is the left mid point of the detlines pattern
                # Prepare the pattern
                coeffs = np.linalg.lstsq(detlines[dm]-detlines[dl])/(detlines[dr]-detlines[dl])
                for da in range(dl+1, dr):  # da is the right mid point of the detlines pattern
                    if dm == da:
                        continue
                    # Create the test pattern
                    daval = (detlines[da]-detlines[dl])/(detlines[dr]-detlines[dl])
                    # Search through all possible patterns in the linelist
                    for ll in range(0, sz_l-nptn+1):  # ll is the start point of the linelist pattern
                        lup = ll + lstsrch
                        if lup > sz_l:
                            lup = sz_l
                        if lstsrch == -1:
                            lup = sz_l
                        for lr in range(ll+nptn-1, lup):  # lr is the end point of the linelist pattern
                            for la in range(ll+1, lr-1):  # la is the end point of the linelist pattern
                                laval = (linelist[la] - linelist[ll]) / (linelist[lr] - linelist[ll])
                                tst = laval - daval
                                if tst < 0.0:
                                    tst *= -1.0
                                if tst <= tol:
                                    # The first pattern matches, check the second one.
                                    for lb in range(la+1, lr):  # la is the end point of the linelist pattern
                                        lbval = (linelist[lb] - linelist[ll]) / (linelist[lr] - linelist[ll])
                                        tst = lbval - dbval
                                        if tst < 0.0:
                                            tst *= -1.0
                                        if tst <= tol:
                                            # The second pattern matches, store the result!
                                            lindex = np.vstack((lindex, np.array([[ll, la, lb, lr]], dtype=np.uint64)))#, dtype=nb.types.uint64)))
                                            dindex = np.vstack((dindex, np.array([[dl, dm, da, dr]], dtype=np.uint64)))#, dtype=nb.types.uint64)))
                                            tst = (linelist[lr] - linelist[ll]) / (detlines[dr] - detlines[dl])
                                            wvl = (npixels/2.)*tst + (linelist[lr]-tst*detlines[dr])
                                            wvcen = np.hstack((wvcen, np.array([wvl], dtype=np.ulonglong))) #, dtype=nb.types.ulong)))
                                            disps = np.hstack((disps, np.array([tst], dtype=np.ulonglong))) #, dtype=nb.types.ulong)))
    # Return, but first remove the spurious first entry due to array creation
    return dindex[1:, :], lindex[1:, :], wvcen[1:], disps[1:]


def empty_patt_dict(nlines):
    """ Return an empty patt_dict

    Parameters
    ----------
    nlines:
        Number of lines for creating the mask.

    Returns
    -------
    patt_dict: dict
        An empty pattern dictionary

    """
    patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., sign=1, mask=np.zeros(nlines, dtype=bool))
    return patt_dict


def solve_xcorr(detlines, linelist, dindex, lindex, line_cc, nreid_min = 4, cc_local_thresh = 0.8):
    """  Given a starting solution, find the best match for all detlines

    Parameters
    ----------
    detlines : ndarray
        list of detected lines in pixels (sorted, increasing)
    linelist : ndarray
        list of lines that should be detected (sorted, increasing)
    dindex : ndarray
        Index array of all detlines (pixels) used in each triangle
    lindex : ndarray
        Index array of the assigned line (wavelengths)to each index in dindex

    Returns
    -------
    patt_dict : dict
       Contains all relevant details of the IDs.  Keys are:

          - acceptable: bool: flag indicating success or failure
          - mask: ndarray, dtype =bool: mask indicating which lines are
            good
          - nmatch: int: Number of matching lines
          - scores: ndarray, str: Scores of the lines
          - IDs: ndarray, float: Wavelength IDs of the lines
          - cc_avg: ndarray, float: Average local zero-lag
            cross-correlation (over all the spectra for which a match
            was obtained) for the most often occuring wavlength ID

    """
    nlines = detlines.size

    # Find the best ID of each line
    detids = np.zeros(nlines)
    cc_avg = np.zeros(nlines)
    scores =np.zeros(nlines,dtype='<U50')
    scores[:] = 'None'
    mask = np.zeros(nlines, dtype=bool)
    ngd_match = 0
    for dd in range(nlines):
        # Grab all the instances of this detected line's pixel position index
        ww = (dindex == dd)
        if not np.any(ww):
            continue
        # Find the unique set of wavelength indices that this detected line has been matched to, and the number of times
        unq, cnts = np.unique(lindex[ww], return_counts=True)
        # Quantify the average xcorr for this line for each set of unique wavelength matches
        cc_per_match = np.zeros(unq.size,dtype=float)
        for iuniq, unq_val in enumerate(unq):
            cc_per_match[iuniq] = np.mean((line_cc[ww])[lindex[ww] == unq_val])
        unq = unq.astype(int)
        # Assign the ID of this line to be wavelength whose index appears the largest number of times
        detids[dd] = linelist[unq[np.argmax(cnts)]]
        # Assign the cross-correlation of this line to be the average of when it was matched to this most often occurring wavelength
        cc_avg[dd] = cc_per_match[np.argmax(cnts)]
        # Give this ID a score based on the number of occurences of the match to a particular wavelength,
        # and the average cross-correlation value for each of those wavelength matches
        scr = score_xcorr(cnts, cc_avg[dd], nreid_min=nreid_min, cc_local_thresh= cc_local_thresh)
        scores[dd] = scr
        if scr in ["Perfect", "Very Good", "Good", "OK"]:
            mask[dd] = True
            ngd_match += 1


    patt_dict = empty_patt_dict(nlines)
    # Iteratively fit this solution, and ID all lines.
    if ngd_match > patt_dict['nmatch']:
        patt_dict['acceptable'] = True
        patt_dict['mask'] = mask
        patt_dict['nmatch'] = ngd_match
        patt_dict['scores'] = scores
        patt_dict['IDs'] = detids
        patt_dict['cc_avg'] = cc_avg

    return patt_dict


def score_xcorr(counts, cc_avg, nreid_min = 4, cc_local_thresh = -1.0):
    """  Grades for the cross-correlation results

    Parameters
    ----------
    counts : ndarray
        Each element is a counter, representing a wavelength that is
        attributed to a given detected line.  The more times that a
        wavelength is attributed to a detected line, the higher the
        counts. The more different wavelengths that are attributed to
        the same detected line (i.e. not ideal) the longer the counts
        list will be.
    nmin_match: int, default = 4, optional
        Minimum number of slits/solutions that have to have been matched
        to receive a score of 'Perfect' or 'Very Good'

    Returns
    -------
    score : str
        A string indicating the relative quality of the ID
    """
    ncnt = counts.size
    max_counts = np.max(counts)
    sum_counts = np.sum(counts)
    # Score
    if (ncnt == 1) and (max_counts >= nreid_min) and (cc_avg > cc_local_thresh):
        score = 'Perfect'
    elif (sum_counts/max_counts >= 0.8) and (max_counts >= nreid_min) and (cc_avg > cc_local_thresh):
        score = 'Very Good'
    elif sum_counts/max_counts >= 0.65:
        score = 'Good'
    elif sum_counts/max_counts >= 0.5:
        score = 'OK'
    elif sum_counts/max_counts >= 0.3:
        score = 'Risky'
    else:
        score = 'Ambitious'
    # Return
    return score


def solve_triangles(detlines, linelist, dindex, lindex, patt_dict=None):
    """
    Given a starting solution, find the best match for all detlines

    Parameters
    ----------
    detlines : ndarray
        list of detected lines in pixels (sorted, increasing)
    linelist : ndarray
        list of lines that should be detected (sorted, increasing)
    dindex : ndarray
        Index array of all detlines (pixels) used in each triangle
    lindex : ndarray
        Index array of the assigned line (wavelengths)to each index in dindex
    patt_dict : dict
        Contains all relevant details of the fit

    """
    nlines = detlines.size
    if patt_dict is None:
        patt_dict = dict(acceptable=False, nmatch=0, ibest=-1, bwv=0., mask=np.zeros(nlines, dtype=bool))

    # Find the best ID of each line
    detids = np.zeros(nlines)
    scores = ['None' for xx in range(nlines)]
    mask = np.zeros(nlines, dtype=bool)
    ngd_match = 0
    for dd in range(nlines):
        # Grab all the instances of this detected line's pixel position index
        ww = (dindex == dd)
        if not np.any(ww):
            continue
        # Find the unique set of wavelength indices that this detected line has been matched to, and the number of times
        unq, cnts = np.unique(lindex[ww], return_counts=True)
        unq = unq.astype(int)
        # Assign the ID of this line to be wavelength whose index appears the largest number of times
        detids[dd] = linelist[unq[np.argmax(cnts)]]
        # Give this ID a score based on the number of occurences
        scr = score_triangles(cnts)
        scores[dd] = scr
        if scr in ["Perfect", "Very Good", "Good", "OK"]:
            mask[dd] = True
            ngd_match += 1

    # Iteratively fit this solution, and ID all lines.
    if ngd_match > patt_dict['nmatch']:
        patt_dict['acceptable'] = True
        patt_dict['mask'] = mask
        patt_dict['nmatch'] = ngd_match
        patt_dict['scores'] = scores
        patt_dict['IDs'] = detids


def score_triangles(counts):
    """
    Grades for the triangle results

    Parameters
    ----------
    counts : ndarray
        Each element is a counter, representing a wavelength that is
        attributed to a given detected line.  The more times that a
        wavelength is attributed to a detected line, the higher the
        counts. The more different wavelengths that are attributed to
        the same detected line (i.e. not ideal) the longer the counts
        list will be.

    Returns
    -------
    score : str
        A string indicating the relative quality of the ID

    """
    ncnt = counts.size
    max_counts = np.max(counts)
    sum_counts = np.sum(counts)
    # Score
    if (ncnt == 1) & (max_counts >= 4):
        score = 'Perfect'
    elif sum_counts/max_counts >= 0.8 and max_counts >= 4:
        score = 'Very Good'
    elif sum_counts/max_counts >= 0.65:
        score = 'Good'
    elif sum_counts/max_counts >= 0.5:
        score = 'OK'
    elif sum_counts/max_counts >= 0.3:
        score = 'Risky'
    else:
        score = 'Ambitious'
    # Return
    return score

