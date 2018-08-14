""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import numba as nb
import pdb


def detect_peaks(image):
    """
    Takes a 2D image and returns 1 if a local maximum is found, and 0 otherwise

    Parameters
    ----------
    image : ndarray
      2D image

    Returns
    -------
    pimage : ndarray
      boolean mask of the peaks (1 when the pixel's value is the neighborhood maximum, 0 otherwise)

    """
    # Define an 8-connected neighborhood
    neighborhood = generate_binary_structure(2, 2)

    # Apply the local maximum filter
    local_max = maximum_filter(image, footprint=neighborhood) == image

    # Background mask
    background = (image == 0)

    # Remove artifacts from local maximum filter
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)

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
    line_list
    tol
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

    Returns
    -------
    best_dict is updated in place
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
    """  Grades quad_match results
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


@nb.jit(nopython=True, cache=True)
def triangles(detlines, linelist, npixels, detsrch=5, lstsrch=10, pixtol=1.0):
    """
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

    lindex = np.zeros((cntdet*cntlst, nptn), dtype=nb.types.uint64)
    dindex = np.zeros((cntdet*cntlst, nptn), dtype=nb.types.uint64)
    wvcen = np.zeros((cntdet*cntlst))
    disps = np.zeros((cntdet*cntlst))

    # Test each detlines combination
    cntdet = 0
    for d in range(0, sz_d-nptn+1):
        dup = d + detsrch
        if dup > sz_d:
            dup = sz_d
        if detsrch == -1:
            dup = sz_d
        for dd in range(d+nptn-1, dup):
            for xd in range(d+1, dd):
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


def solve_triangles(detlines, linelist, dindex, lindex, patt_dict=None):
    """  Given a starting solution, find the best match for all detlines

    Parameters
    ----------
    detlines : ndarray
      list of detected lines in pixels (sorted, increasing)
    linelist : ndarray
      list of lines that should be detected (sorted, increasing)
    dindex : ndarray
      Index array of all detlines used in each triangle
    lindex : ndarray
      Index array of the assigned line to each index in dindex
    patt_dict : dict
      Contains all relevant details of the fit

    Returns
    -------

    """
    nlines = detlines.size
    if patt_dict is None:
        patt_dict = dict(nmatch=0, ibest=-1, bwv=0.)

    # Find the best ID of each line
    detids = np.zeros(nlines)
    scores = ['None' for xx in range(nlines)]
    mask = np.zeros(nlines, dtype=np.bool)
    ngd_match = 0
    for dd in range(nlines):
        ww = np.where(dindex == dd)
        if ww[0].size == 0:
            continue
        unq, cnts = np.unique(lindex[ww], return_counts=True)
        unq = unq.astype(np.int)
        detids[dd] = linelist[unq[np.argmax(cnts)]]
        scr = score_triangles(cnts)
        scores[dd] = scr
        if scr in ["Perfect", "Very Good", "Good", "OK"]:
            mask[dd] = True
            ngd_match += 1

    # Iteratively fit this solution, and ID all lines.
    if ngd_match > patt_dict['nmatch']:
        patt_dict['mask'] = mask
        patt_dict['nmatch'] = ngd_match
        patt_dict['scores'] = scores
        patt_dict['IDs'] = detids
    return


def score_triangles(counts):
    """  Grades for the triangle results
    Parameters
    ----------
    fidx

    Returns
    -------
    scores : list

    """
    ncnt = counts.size
    max_counts = np.max(counts)
    sum_counts = np.sum(counts)
    # Score
    if (ncnt == 1) & (max_counts >= 4):
        score = 'Perfect'
    elif sum_counts/max_counts >= 0.8:
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
