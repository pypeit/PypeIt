""" Module for core algorithms related to tracing slits/orders
These should primarily be called by the TraceSlits class
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect
import copy
from collections import Counter

import numpy as np

from scipy import ndimage

import matplotlib.pyplot as plt
from matplotlib import cm, font_manager

from pypit import msgs
from pypit import arqa
from pypit import arplot
from pypit import arutils
from pypit import arpca
from pypit import arpixels
from pypit import arproc
from pypit import arparse as settings
from pypit import arspecobj
from pypit import ardebug as debugger

try:
    from pypit import ginga
except ImportError:
    pass

from pypit import arcytrace

try:
    ustr = unicode
except NameError:
    ustr = str

# Testing
import time


def add_user_edges(edgearr, siglev, tc_dict, add_slits):
    """Add user-defined slit(s)

    Warning: There is no real error checking here.
    The user is assumed to know what they are doing!

    Parameters
    ----------
    edgearr
    siglev
    tc_dict
    add_slits

    Returns
    -------

    """

    # Indices
    lmin = np.min(edgearr)
    rmax = np.max(edgearr)
    new_l = lmin-1
    new_r = rmax+1

    # Grab the edge indexes and xval's
    left_idx, left_xval, right_idx, right_xval = tc_indices(tc_dict)

    # Loop me
    nrow = edgearr.shape[0]
    ycen = nrow//2
    for new_slit in add_slits:
        msgs.info("Adding a user-defined slit [x0, x1, yrow]:  {}".format(new_slit))
        # Parse
        xleft, xright, yrow = new_slit
        # Left or right
        for side in ['left','right']:
            # Trace crude and setup
            if side == 'left':
                xset, xerr = trace_crude_init(np.maximum(siglev, -0.1), np.array([xleft]), yrow)
                #
                new_i = new_l
                ref_x = left_xval
                ref_i = left_idx
            else:
                xset, xerr = trace_crude_init(np.maximum(-1*siglev, -0.1), np.array([xright]), yrow)
                #
                new_i = new_r
                ref_x = right_xval
                ref_i = right_idx
            # Was the trace good enough?
            ygd = np.where(xerr[:,0] != 999.)[0]
            new_xval = int(np.round(xset[ycen, 0]))  # Always defined at the 1/2 point
            if len(ygd) > nrow//2: # Use the trace if it was primarily successful
                xvals = np.round(xset[:, 0]).astype(int)
                edgearr[ygd, xvals[ygd]] = new_i
            else: # Otherwise, find the closest left edge and use that
                # Find the closest
                idx = np.argmin(np.abs(ref_x-new_xval))
                ref_slit = ref_i[idx]
                dx = ref_x[idx]-new_xval

                # Grab its pixels
                i_pix = np.where(edgearr == ref_slit)
                new_pix = (i_pix[0], i_pix[1]-dx)
                # And use them
                edgearr[new_pix] = new_i
            # Update
            tc_dict[side]['xval'][str(new_i)] = new_xval
            if side == 'left':
                new_l -= 1
            else:
                new_r += 1
    # Return
    return edgearr


def assign_slits(binarr, edgearr, ednum=100000, lor=-1, settings=None):
    """This routine will trace the locations of the slit edges
    Putative edges come in with |values| > 200000 (in the edgearr)
    and leave (ideally) with values near ednum

    Parameters
    ----------
    binarr : numpy ndarray
      Calibration frame that will be used to identify slit traces (in most cases, the slit edge)
      Typically previously processed by a uniform filter
    edgearr : numpy ndarray
      An array of negative/positive numbers (left/right edges respectively) and zeros (no edge)
    ednum : int
      A dummy number given to define slit edges
    lor : int (-1 or +1)
      A flag that indicates if the left edge (-1) or right edge (+1) should be assigned
    settings : dict (optional)
       Set polyorder, function
       Likely to be replaed by another approach

    Returns
    -------
    edgearr : numpy ndarray
      An array of negative/positive numbers (left/right edges respectively) and zeros (no edge)
      Modified in-place
    """
    if settings is None:
        settings = dict(trace={'slits': {'polyorder': 3, 'function': 'legendre'}})

    if lor == -1:
        lortxt = "left"
    else:
        lortxt = "right"
    outitnm = 1
    oldedgearr = edgearr.copy()
    prevedgearr = edgearr.copy()
    while True:
        msgs.prindent("Outer {0:s} edge loop, Iteration {1:d}".format(lortxt, outitnm))
        labnum = lor*ednum
        itnm = 0
        nslit = 0
        cmnold = None
        firstpass = True
        while True:
            # Array to hold edge information
            edgehist = np.zeros(binarr.shape[1]*2, dtype=np.int)
            itnm += 1
            # Locate edges relative to the most common edge
            if lor == -1:
                wl = np.where(edgearr <= -2*ednum)
            else:
                wl = np.where(edgearr >= 2*ednum)
            if wl[0].size == 0:
                break
            cl = Counter(edg for edg in edgearr[wl])
            comml = cl.most_common(1)
            # Pixels of the most common edge
            common_pix = np.where(edgearr == comml[0][0])
            if not firstpass:
                if (cmnold[0] == comml[0][0]) and (cmnold[1] == comml[0][1]):
                    # Nothing has changed since the previous iteration, so end the loop
                    break
                if comml[0][1] < binarr.shape[0]/100.0:
                    # Now considering an edge that spans less than 1 per cent of the detector ---> insignificant
                    break
            # Save
            cmnold = comml[0]
            # Extract just these elements
            tedgearr = edgearr[common_pix[0], :]
            # Set the offset
            offs = binarr.shape[1]
            # Add these into edgehist
            edgehist[offs] = np.sum(binarr[common_pix])#common_pix[0].size
            # And a fudge to give this edge detection some width (for peak finding, below)
            edgehist[offs-1] = 1 + common_pix[0].size/2
            edgehist[offs+1] = 1 + common_pix[0].size/2
            # Find the difference between unknown edges
            if lor == -1:  # Left edges
                www = np.where(tedgearr <= -2*ednum)
            else:
                www = np.where(tedgearr >= 2*ednum)
            if www[0].size == 0:
                break
            shft = www[1] - common_pix[1][www[0]]  # Calculate the shift between right edges
            shft += offs  # Apply the offset to the edgehist arr
            #np.add.at(edgehist, shft, 1)
            # Insert other edges into edgehist
            np.add.at(edgehist, shft, binarr[common_pix[0], :][www])
            # Smooth the histogram with a Gaussian of standard deviation 1 pixel to reduce noise
            smedgehist = ndimage.uniform_filter1d(edgehist, 3)
            # Identify peaks (which indicate the locations of the right slit edges)
            #   Might consider another peak-finding algorithm, e.g. the one used for arc lines
            arrlfr = smedgehist[0:-4]
            arrlft = smedgehist[1:-3]
            arrcen = smedgehist[2:-2]
            arrrgt = smedgehist[3:-1]
            arrrfr = smedgehist[4:]
            wpk = np.where((arrcen >= arrlft) & (arrcen > arrrgt) &  # Exactly one of these should be >=
                           ((arrlft > arrlfr) | (arrrgt > arrrfr)))[0]
            if wpk.size == 0:
                # No more peaks
                break
            if wpk.size != 1:
                try:
                    wpkmsk = prune_peaks(smedgehist, wpk, np.where(wpk+2 == offs)[0][0])#, debug=debug)
                except:
                    debugger.set_trace()
                wpk = wpk[np.where(wpkmsk == 1)]
            if wpk.size == 0:
                # After pruning, there are no more peaks
                break
            pks = wpk+2  # Shifted by 2 because of the peak finding algorithm above
#            print('calling find_peak_limits')
#            t = time.clock()
#            _pedges = arcytrace.find_peak_limits(smedgehist, pks)
#            print('Old find_peak_limits: {0} seconds'.format(time.clock() - t))
#            t = time.clock()
            pedges = new_find_peak_limits(smedgehist, pks)
#            print('New find_peak_limits: {0} seconds'.format(time.clock() - t))
#            assert np.sum(_pedges != pedges) == 0, \
#                    'Difference between old and new find_peak_limits'

            if np.all(pedges[:, 1]-pedges[:, 0] == 0):
                # Remaining peaks have no width
                break
            if msgs._debug['trace']:
                plt.clf()
                plt.plot(arrcen, 'k-', drawstyle='steps')
                plt.plot(wpk, np.zeros(wpk.size), 'ro')
                plt.show()
                debugger.set_trace()
            # Label all edge ids (in the original edgearr) that are located in each peak with the same number
            for ii in range(pks.size):
                shbad = np.zeros(edgearr.shape)
                wp = np.where((shft >= pedges[ii, 0]) & (shft <= pedges[ii, 1]))
                vals = np.unique(tedgearr[(www[0][wp], www[1][wp])])
                # Fit the edge detections in this edge and calculate the offsets
                strev = "np.where("
                for vv in vals:
                    strev += "(edgearr=={0:d})|".format(vv)
                strev = strev[:-1] + ")"
                widx = eval(strev)
                if widx[0].size < 2*settings['trace']['slits']['polyorder']:
                    continue
                badmsk, fitcof = arutils.robust_polyfit(widx[0], widx[1],
                                                        settings['trace']['slits']['polyorder'],
                                                        function=settings['trace']['slits']['function'],
                                                        minv=0, maxv=binarr.shape[0]-1)
                shbad[widx] = badmsk
                smallhist = np.zeros(101, dtype=np.int)
                meddiff = np.zeros(vals.size)
                for vv in range(vals.size):
                    widx = np.where((edgearr == vals[vv]) & (shbad == 0))
                    if widx[0].size == 0:
                        # These pixels were deemed to be bad
                        continue
                    diff = widx[1] - arutils.func_val(fitcof, widx[0],
                                                      settings['trace']['slits']['function'],
                                                      minv=0, maxv=binarr.shape[0]-1)
                    diff = 50 + np.round(diff).astype(np.int)
                    np.add.at(smallhist, diff, 1)
                    meddiff[vv] = np.median(diff)
                # Find the peaks of this distribution
                wspk = np.where((smallhist[1:-1] >= smallhist[2:]) & (smallhist[1:-1] > smallhist[:-2]))[0]
                wspk += 1  # Add one here to account for peak finding
                if msgs._debug['trace'] and False:
                    plt.clf()
                    plt.plot(smallhist, 'k-', drawstyle='steps')
                    plt.show()

                for pp in range(wspk.size):  # For all small peaks identified
                    for vv in range(vals.size):
                        if lor == -1 and vals[vv] > -2*ednum:
                            continue
                        elif lor == 1 and vals[vv] < 2*ednum:
                            continue
                        # Make sure this value is within 1 pixel of the peak
                        if meddiff[vv] < wspk[pp]-1:
                            continue
                        if meddiff[vv] > wspk[pp]+1:
                            continue
                        edgearr[np.where(edgearr == vals[vv])] = labnum
                        meddiff[vv] = -1  # Flag this val as done
                    labnum += lor*1
                # Find any vals that weren't applied
                for vv in range(vals.size):
                    if meddiff[vv] == -1:
                        continue
                    edgearr[np.where(edgearr == vals[vv])] = 0
            nslit += pks.size
            msgs.prindent("  Inner loop, Iteration {0:d}, {1:d} {2:s} edges assigned ({3:d} total)".format(itnm, pks.size, lortxt, nslit))
            firstpass = False
        outitnm += 1
        if lor == -1:
            edgearr[np.where(edgearr <= -2*ednum)] = 0
        else:
            edgearr[np.where(edgearr >= 2*ednum)] = 0
        if np.array_equal(edgearr, oldedgearr) or np.array_equal(edgearr, prevedgearr):
            break
        elif outitnm > 10:
            msgs.warn("Edge assignment may not have converged")
            msgs.info("Please check the slit edge traces")
            #debugger.set_trace()
            break
        else:
            oldedgearr = prevedgearr.copy()
            prevedgearr = edgearr.copy()
            if lor == -1:
                edgearr[np.where(edgearr <= -ednum)] -= ednum
            else:
                edgearr[np.where(edgearr >= ednum)] += ednum
    # Ignore any order detections that weren't identified in the loop
    if lor == -1:
        edgearr[np.where(edgearr <= -2*ednum)] = 0
    else:
        edgearr[np.where(edgearr >= 2*ednum)] = 0
    # Sort vals by increasing spatial position on the detector
    # First, determine the model for the most common slit edge
    if lor == -1:
        wcm = np.where(edgearr <= -ednum)
    else:
        wcm = np.where(edgearr >= ednum)
    if wcm[0].size != 0:
        cntr = Counter(edg for edg in edgearr[wcm])
        commn = cntr.most_common(1)
        wedx, wedy = np.where(edgearr == commn[0][0])
        msk, cf = arutils.robust_polyfit(wedx, wedy,
                                         settings['trace']['slits']['polyorder'],
                                         function=settings['trace']['slits']['function'],
                                         minv=0, maxv=binarr.shape[0]-1)
        cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                   settings['trace']['slits']['function'],
                                   minv=0, maxv=binarr.shape[0]-1)
        if lor == -1:
            vals = np.unique(edgearr[np.where(edgearr < 0)])
        else:
            vals = np.unique(edgearr[np.where(edgearr > 0)])
        diffarr = np.zeros(vals.size)
        diffstd = 0.0
        for jj in range(vals.size):
            wedx, wedy = np.where(edgearr == vals[jj])
            diffarr[jj] = np.mean(wedy-cenmodl[wedx])
            diffstd += np.std(wedy-cenmodl[wedx])
        diffstd /= vals.size
        dasrt = np.argsort(diffarr)
        # Relabel the edges from left to right
        edgearr[wcm] += lor*ednum
        labnum = lor*ednum
        diffarrsrt = diffarr[dasrt]
        diffs = diffarrsrt[1:] - diffarrsrt[:-1]
        for jj in range(vals.size):
            wrplc = np.where(edgearr == lor*ednum + vals[dasrt[jj]])
            edgearr[wrplc] = labnum
            if jj != vals.size-1:
                if diffs[jj] > 3.0*diffstd:
                    # The next edge must be a different edge
                    labnum += lor*1
    return







def add_edge(ref_slit, insert_offset, earr, t_dict, final_left, final_right, left=True):
    """  Add a new edge using a reference slit
    right_slit : Reference slit for the left one
    insert_offset : int
      Offset fromm the right slit for the new left slit
    """
    # New left edge index
    if left== 'left':
        new_e = np.min(earr)-1
    else:
        new_e = np.max(earr)+1
    # Use the reference edge for the shape
    i_pix = np.where(earr == ref_slit)
    new_pix = (i_pix[0], i_pix[1]+insert_offset)
    # Add it in
    earr[new_pix] = new_e
    if left:
        t_dict['left']['xval'][str(new_e)] = t_dict['right']['xval'][str(ref_slit)]+insert_offset
    else:
        t_dict['right']['xval'][str(new_e)] = t_dict['left']['xval'][str(ref_slit)]+insert_offset
    # Lists
    if left:
        final_left.append(new_e)
        final_right.append(ref_slit)
    else:
        final_right.append(new_e)
        final_left.append(ref_slit)


def edgearr_mslit_sync(edgearr, tc_dict, ednum, insert_buff=5, add_left_edge_slit=True, verbose=False):
    """ Method to synchronize the slit edges
    Adds in extra edges according to a few criteria

    Developed for ARMLSD

    Parameters
    ----------
    edgearr : ndarray
    tc_dict : dict
       For book-keeping
    ednum : int
    insert_buff : int, optional
       Offset from existing edge for any edge added in
    add_left_edge_slit : bool, optional
       Allow the method to add in a left slit at the edge of the detector

    Returns
    -------
    new_edgearr : ndarray
      Updated edgearr
    """
    # Init
    new_edgearr = np.zeros_like(edgearr, dtype=int)
    final_left = []
    final_right = []

    # Grab the edge indexes and xval's
    left_idx, left_xval, right_idx, right_xval = tc_indices(tc_dict)

    # Only one slit?
    if (len(left_xval) == 1) and (len(right_xval)==1):
        if left_xval[0] < right_xval[0]:  # Ok slit, otherwise continue
            return edgearr

    # First slit (often up against the detector)
    if (right_xval[0] < left_xval[0]) and add_left_edge_slit:
        right_pix = np.where(edgearr == right_idx[0])
        mn_rp = np.min(right_pix[1])
        if mn_rp <= insert_buff:
            msgs.warn("Partial or too small right edge at start of detector.  Skipping it.")
        else:
            ioff = -1*mn_rp + insert_buff
            msgs.warn("Adding in a left edge at start of detector which mirrors the first right edge")
            add_edge(right_idx[0], ioff, edgearr, tc_dict, final_left, final_right, left=True)

    # Loop on left edges
    for kk,left in enumerate(left_xval):

        # Grab location of the next left edge
        if kk < len(left_idx)-1:
            next_left = left_xval[kk+1]
        else:
            next_left = edgearr.shape[1]-1

        # Search for a proper right edge
        #  Should be right of the current left and left of the next left
        gd_right = np.where((right_xval < next_left) & (right_xval > left))[0]
        if len(gd_right) == 0:   # None found?
            # Last slit?
            if kk == len(left_idx)-1:
                msgs.warn("Last slit has no right edge.  Adding one in which will not touch the detector edge")
                left_pix = np.where(edgearr == left_idx[kk])
                mx_lp = np.max(left_pix[1])
                if mx_lp >= edgearr.shape[1]:
                    msgs.warn("Partial left edge at end of detector.  Skipping it.")
                else:
                    # Stay on the detector!
                    ioff = edgearr.shape[1] - mx_lp - insert_buff
                    # Add
                    add_edge(left_idx[kk], ioff, edgearr, tc_dict, final_left, final_right, left=False)
                continue
            else: # Not the last slit, add one in!
                msgs.warn("Missing a right edge for slit with left edge at {}".format(left))
                msgs.warn("Adding in a corresponding right edge!")
                # Offset from the next left edge
                ioff = next_left-left-insert_buff
                # Add
                add_edge(left_idx[kk], ioff, edgearr, tc_dict, final_left, final_right, left=False)
        else:
            # Add in the first right edge
            final_left.append(left_idx[kk])
            iright = np.min(gd_right[0])
            final_right.append(right_idx[iright])

            # Check for multiple right edges between the two lefts (i.e. missing Left)
            #     Will only add in one missing left
            if len(gd_right) > 1:
                msgs.warn("Missing a left edge for slit with right edge(s) at {}".format(
                    right_xval[gd_right[1:]]))
                msgs.warn("Adding one (and only one) in unless you turn off the setting [blah]")
                # Offset is difference between the two right slits + a buffer
                ioff = right_xval[gd_right[0]] - right_xval[gd_right[1]] + insert_buff
                # Add
                add_edge(right_idx[gd_right[1]], ioff, edgearr, tc_dict, final_left, final_right, left=True)

    # Finish by remaking the edgearr
    #  And update the book-keeping dict
    ldict, rdict = {}, {}
    # Left
    for ss, left_i in enumerate(final_left):
        newval = -1*ednum - ss
        # Edge
        pix = edgearr == left_i
        new_edgearr[pix] = newval
        # Dict
        ldict[str(newval)] = tc_dict['left']['xval'][str(left_i)]
    tc_dict['left']['xval'] = ldict
    # Right
    for ss, right_i in enumerate(final_right):
        newval = ednum + ss
        # Edge
        pix = edgearr == right_i
        new_edgearr[pix] = newval
        # Dict
        rdict[str(newval)] = tc_dict['right']['xval'][str(right_i)]
    tc_dict['right']['xval'] = rdict

    if verbose:
        print(tc_dict['left']['xval'])
        print(tc_dict['right']['xval'])

    # Return
    return new_edgearr


def edgearr_tcrude(edgearr, siglev, ednum, TOL=3., tfrac=0.33, verbose=False):
    """ Use trace_crude to refine slit edges
    It is also used to remove bad slit edges and merge slit edges

    Only recommended for ARMLSD

    Parameters
    ----------
    edgearr : ndarray
      Edge array
    siglev : ndarray
      Sigma level image
    ednum : int
    TOL : float (optional)
      Tolerance for matching 2 edges (and ignoring one)
    tfrac : float (optional)
      Fraction of the slit edge that must be traced to keep
      There are exceptions, however (e.g. single slits)

    Returns
    -------
    new_edgarr : ndarray
       A new version of the edgearr
    tc_dict : dict
       A dict that has book-keeping on the edges
         left/right
             xval -- Position of edge at ycen (1/2 point on the detector); most useful parameter recorded
             uni_idx -- Unique edge numbers
             flags -- Internal book-keeping on edge analysis
             xset -- trace set values
             xerr -- trace set errors
    """
    # Init
    ycen = edgearr.shape[0] // 2

    # Items to return
    new_edgarr = np.zeros_like(edgearr, dtype=int)
    tc_dict = {}

    # Loop on side
    for side in ['left', 'right']:
        tc_dict[side] = {}
        tc_dict[side]['xval'] = {}

        # Unique edge values
        if side == 'left':
            uni_e = np.unique(edgearr[edgearr < 0])
        else:
            uni_e = np.unique(edgearr[edgearr > 0])
        # Save sorted (in absolute value) -- This is necessary
        uni_e = uni_e[np.argsort(np.abs(uni_e))]
        tc_dict[side]['uni_idx'] = uni_e

        # Flag: 0=not traced; 1=traced; -1=duplicate
        tc_dict[side]['flags'] = np.zeros(len(uni_e), dtype=int)


        # Loop on edges to trace
        niter = 0
        while np.any(tc_dict[side]['flags'] == 0):

            # Most common yrow to speed up trace_crude
            if side == 'left':
                all_e = np.where(edgearr < 0)
            else:
                all_e = np.where(edgearr > 0)
            cnt = Counter(all_e[0])
            yrow = cnt.most_common(1)[0][0]

            # Grab the x values on that row
            xinit = all_e[1][all_e[0] == yrow]
            # Unique..
            _, uidx = np.unique(edgearr[yrow, xinit], return_index=True)
            xinit = xinit[uidx]  # Yes, this is needed..
            #
            msk = np.ones_like(xinit, dtype=bool)

            # If not first pass, look for duplicates
            if niter > 0:
                # Check if one of the existing slits is consistent
                #  with this new, offset edge
                for ii,xx in enumerate(xinit):
                    # Clip (this doesn't catch them all, but there is another check below)
                    if np.min(np.abs(xx-tc_dict[side]['xset'][yrow])) < TOL:
                        msk[ii] = False
                        eval = edgearr[yrow,xx]
                        msgs.warn("Duplicate {} edge at x={} and y={}.  Clipping..".format(side, xx,yrow))
                        tc_dict[side]['flags'][uni_e==eval] = -1  # Clip me
                        # Zero out edgearr
                        edgearr[edgearr==eval] = 0
                # Any of these new?  If not continue
                if not np.any(msk):
                    # Next iteration
                    niter += 1
                    continue
                else:
                    xinit = xinit[msk]
                    pass
            # Trace crude
            if side == 'left':
                xset, xerr = trace_crude_init(np.maximum(siglev, -0.1), np.array(xinit), yrow)
            else:
                xset, xerr = trace_crude_init(np.maximum(-1*siglev, -0.1), np.array(xinit), yrow)
            # Save
            if niter == 0:
                tc_dict[side]['xset'] = xset
                tc_dict[side]['xerr'] = xerr
            else: # Need to append
                tc_dict[side]['xset'] = np.append(tc_dict[side]['xset'], xset, axis=1)
                tc_dict[side]['xerr'] = np.append(tc_dict[side]['xerr'], xerr, axis=1)

            # Good values allowing for edge of detector
            goodx = np.any([(xerr != 999.), (xset==0.), (xset==edgearr.shape[1]-1.)], axis=0)
            # Fill in
            for kk, x in enumerate(xinit):
                yval = np.where(goodx[:,kk])[0]
                eval = edgearr[yrow,x]
                new_xval = int(np.round(xset[ycen, kk]))

                # Check whether the trace is well determined on tface of the detector
                #    If only one trace, this check is ignored
                if (len(yval) < int(tfrac*edgearr.shape[0])) and (len(uni_e) > 1):
                    msgs.warn("Edge at x={}, y={} traced less than {} of the detector.  Removing".format(
                        x,yrow,tfrac))
                    tc_dict[side]['flags'][uni_e==eval] = -1  # Clip me
                    # Zero out edgearr
                    edgearr[edgearr==eval] = 0
                    continue
                # Do not allow a right edge at x=0
                if (side == 'right') and (new_xval==0):
                    msgs.warn("Right edge at x=0 removed")
                    tc_dict[side]['flags'][uni_e==eval] = -1  # Clip me
                    edgearr[edgearr==eval] = 0
                    continue
                # or a left edge at x=end
                if (side == 'left') and (new_xval==edgearr.shape[1]-1):
                    msgs.warn("Left edge at detector right edge removed")
                    tc_dict[side]['flags'][uni_e==eval] = -1  # Clip me
                    edgearr[edgearr==eval] = 0
                    continue
                # Check it really is a new xval (within TOL)
                if niter > 0:
                    curr_xvals = np.array([tc_dict[side]['xval'][key] for key in tc_dict[side]['xval'].keys()])
                    if np.min(np.abs(new_xval-curr_xvals)) < TOL:
                        msgs.warn("Edge matched exiting xval within TOL, clipping")
                        tc_dict[side]['flags'][uni_e==eval] = -1  # Clip me
                        edgearr[edgearr==eval] = 0
                        continue

                # Edge is ok, keep it
                xvals = np.round(xset[:, kk]).astype(int)
                # Single edge requires a bit more care (it is so precious!)
                if len(uni_e) == 1:
                        if np.sum(edgearr==eval)>len(yval):
                            new_edgarr[edgearr==eval] = edgearr[edgearr==eval]
                        else:
                            new_edgarr[yval, xvals[yval]] = eval
                else:
                    # Traces can disappear and then the crude trace can wind up hitting a neighbor
                    # Therefore, take only the continuous good piece from the starting point
                    ybad_xerr = np.where(~goodx[:,kk])[0]
                    # Lower point
                    ylow = ybad_xerr < yrow
                    if np.any(ylow):
                        y0 = np.max(ybad_xerr[ylow])+1  # Avoid the bad one
                    else:
                        y0 = 0
                    # Upper
                    yhi = ybad_xerr > yrow
                    if np.any(yhi):
                        y1 = np.min(ybad_xerr[yhi])
                    else:
                        y1 = edgearr.shape[0]
                    new_yval = np.arange(y0,y1).astype(int)  # Yes, this is necessary;  slicing fails..
                    #
                    new_edgarr[new_yval, xvals[new_yval]] = eval
                    #new_edgarr[yval, xvals[yval]] = eval
                # Flag
                tc_dict[side]['flags'][uni_e == eval] = 1
                # Save new_xval
                tc_dict[side]['xval'][str(eval)] = new_xval
                # Zero out edgearr
                edgearr[edgearr==eval] = 0
            # Next pass
            niter += 1

    # Reset edgearr values to run sequentially and update the dict
    for side in ['left', 'right']:
        if np.any(tc_dict[side]['flags'] == -1):
            # Loop on good edges
            gde = np.where(tc_dict[side]['flags'] == 1)[0]
            for ss,igde in enumerate(gde):
                if side == 'left':
                    newval = -1*ednum - ss
                else:
                    newval = ednum + ss
                oldval = tc_dict[side]['uni_idx'][igde]
                pix = new_edgarr == oldval
                new_edgarr[pix] = newval
                # Reset the dict too..
                if newval == 0:
                    debugger.set_trace()
                tc_dict[side]['xval'][str(newval)] = tc_dict[side]['xval'].pop(str(oldval))
    # Remove uni_idx
    for side in ['left', 'right']:
        for key in ['uni_idx', 'xset', 'xerr']:
            tc_dict[side].pop(key)
    if verbose:
        print(tc_dict['left']['xval'])
        print(tc_dict['right']['xval'])
    # Return
    return new_edgarr, tc_dict.copy()


def edgearr_from_user(shape, ledge, redge, det):
    """ Add a user-defined slit?
     Syntax is a list of values, 2 per detector that define the slit
     according to column values.  The 2nd value (for the right edge)
     must be >0 to be applied.  Example for LRISr [-1, -1, 7, 295]
     which means the code skips user-definition for the first detector
     but adds one for the 2nd.

    Parameters
    ----------
    shape : tuple
    ledge : int
    redge : int
    det : int
      Only used for comment so could be removed..

    Returns
    -------
    edgearr : ndarray (int)
      -1 : left edge
      1 : right edge
    """
    edgearr = np.zeros(shape, dtype=np.int)
    if settings['trace']['slits']['single'][redge] > 0:
        msgs.warn("Using input slit edges on detector {:d}: [{:g},{:g}]".format(
            det, ledge, redge))
            #settings['trace']['slits']['single'][ledge],
            #settings['trace']['slits']['single'][redge]))
        msgs.warn("Better know what you are doing!")
        edgearr[:, ledge] = -1
        edgearr[:, redge] = +1
    return edgearr


def edgearr_from_binarr(binarr, binbpx, medrep=0, min_sqm=30., maskBadRows=False,
                        sobel_mode='nearest', sigdetect=30., number_slits=-1):
    """ Generate the edge array from an input, trace image (likely slightly filtered)

    Parameters
    ----------
    binarr : numpy ndarray
      Calibration frame that will be used to identify slit traces (in most cases, the slit edge)
      Lightly filtered
    binbpx : ndarray
    medrep : int, optional
        Number of times to perform median smoothing on the mstrace
        One uniform filter is always done
        medrep = 0 is recommended for ARMLSD
    sobel_mode : str, optional
        ndimage.sobel mode
    sigdetect : float, optional
        threshold for edge detection
    number_slits : int, optional
        if >0, restrict the number of slits to this value
        not well tested
    min_sqm : float, optional
        Minimum error used when detecting a slit edge
    maskBadRows : bool, optional
      Mostly useful for echelle data where the slit edges are bent relative to
      the pixel columns. Do not set this keyword to True if slit edges are
      almost aligned with the pixel columns.

    Returns
    -------
    siglev : ndarray
    edgearr : ndarray
    """
    # Specify how many times to repeat the median filter
    # Even better would be to fit the filt/sqrt(abs(binarr)) array with a Gaussian near the maximum in each column
    msgs.info("Detecting slit edges in the mstrace image")

    # Replace bad columns -- should put this method somewhere else
    bad_cols = np.sum(binbpx, axis=0) > (binbpx.shape[0]//2)
    if np.any(bad_cols):
        ms2 = arproc.replace_columns(binarr, bad_cols)
    else:
        ms2 = binarr.copy()

    # Generate sqrt image
    sqmstrace = np.sqrt(np.abs(ms2))

    # Median filter, as desired
    # TODO -- Try size=(7,3) to bring up the edges instead of (3,7)
    for ii in range(medrep):
        #sqmstrace = ndimage.median_filter(sqmstrace, size=(3, 7))
        sqmstrace = ndimage.median_filter(sqmstrace, size=(7, 3))

    # Make sure there are no spuriously low pixels
    sqmstrace[(sqmstrace < 1.0) & (sqmstrace >= 0.0)] = 1.0
    sqmstrace[(sqmstrace > -1.0) & (sqmstrace <= 0.0)] = -1.0

    # Filter with a Sobel
    filt = ndimage.sobel(sqmstrace, axis=1, mode=sobel_mode)
    filt *= (1.0 - binbpx)  # Apply to the bad pixel mask
    # siglev
    siglev = np.sign(filt)*(filt**2)/np.maximum(sqmstrace, min_sqm)
    # First edge move
    tedges = np.zeros(binarr.shape, dtype=np.float)
    wl = np.where(siglev > + sigdetect)  # A positive gradient is a left edge
    wr = np.where(siglev < - sigdetect)  # A negative gradient is a right edge
    tedges[wl] = -1.0
    tedges[wr] = +1.0
    #if False:
    #    import astropy.io.fits as pyfits
    #    hdu = pyfits.PrimaryHDU(filt)
    #    hdu.writeto("filt_{0:02d}.fits".format(det), overwrite=True)
    #    hdu = pyfits.PrimaryHDU(sqmstrace)
    #    hdu.writeto("sqmstrace_{0:02d}.fits".format(det), overwrite=True)
    #    hdu = pyfits.PrimaryHDU(binarr)
    #    hdu.writeto("binarr_{0:02d}.fits".format(det), overwrite=True)
    #    hdu = pyfits.PrimaryHDU(siglev)
    #    hdu.writeto("siglev_{0:02d}.fits".format(det), overwrite=True)
    # Clean the edges
    wcl = np.where((ndimage.maximum_filter1d(siglev, 10, axis=1) == siglev) & (tedges == -1))
    wcr = np.where((ndimage.minimum_filter1d(siglev, 10, axis=1) == siglev) & (tedges == +1))
    nedgear = np.zeros(siglev.shape, dtype=np.int)
    nedgear[wcl] = -1
    nedgear[wcr] = +1
    #nedgear = arcytrace.clean_edges(siglev, tedges)
    if maskBadRows:
        msgs.info("Searching for bad pixel rows")
        edgsum = np.sum(nedgear, axis=0)
        sigma = 1.4826*np.median(np.abs(edgsum-np.median(edgsum)))
        w = np.where(np.abs(edgsum) >= 1.5*sigma)[0]
        maskcols = np.unique(np.append(w, np.append(w+1, w-1)))
        msgs.info("Masking {0:d} bad pixel rows".format(maskcols.size))
        for i in range(maskcols.size):
            if maskcols[i] < 0 or maskcols[i] >= nedgear.shape[1]:
                continue
            nedgear[:, maskcols[i]] = 0
    ######
    msgs.info("Applying bad pixel mask")
    nedgear *= (1-binbpx.astype(np.int))  # Apply to the bad pixel mask

    # Cut down to user-desired number
    if number_slits > 0:
        sigedg = np.copy(siglev)
        sigedg[np.where(nedgear == 0)] = 0
        sigedg[np.where(np.isinf(sigedg) | np.isnan(sigedg))] = 0
        # Now identify the number of most significantly detected peaks (specified by the user)
        amnmx = np.argsort(sigedg, axis=1)
        edgearr = np.zeros(binarr.shape, dtype=np.int)
        xsm = np.arange(binarr.shape[0])
        for ii in range(0, number_slits):
            wset = np.where(sigedg[(xsm, amnmx[:, ii])] != 0)
            edgearr[(wset[0], amnmx[wset[0], ii])] = 1
            wset = np.where(sigedg[(xsm, amnmx[:, amnmx.shape[1] - 1 - ii])] != 0)
            edgearr[(wset[0], amnmx[wset[0], amnmx.shape[1] - 1 - ii])] = -1
    else:
        edgearr = np.copy(nedgear)

    # Return
    return siglev, edgearr


def edgearr_add_left_right(edgearr, binarr, binbpx, lcnt, rcnt, ednum):
    """ Add left/right edges in the event that none were found thus far
    This is especially useful for long slits that fill the full detector,
    e.g. Kast

    Parameters
    ----------
    edgearr : ndarray
    binarr : ndarray
    binbpx : ndarray
      Bad pixel mask
    lcnt : int
      Number of left edges
    rcnt : int
      Number of right edges
    ednum : int

    Returns
    -------
    edgearrcp : ndarray
      New edgearr
    lcnt : int
      Updated count
    rcnt : int
      Updated count
    """
    if lcnt == 1:
        letxt = "edge"
    else:
        letxt = "edges"
    if rcnt == 1:
        retxt = "edge"
    else:
        retxt = "edges"
    msgs.info("{0:d} left {1:s} and {2:d} right {3:s} were found in the trace".format(lcnt, letxt, rcnt, retxt))
    if (lcnt == 0) and (rcnt == 0):
        if np.median(binarr) > 500:
            msgs.warn("Found flux but no edges.  Assuming they go to the edge of the detector.")
            edgearr[:, -1] = 2*ednum
            rcnt = 1
            edgearr[:, 0] = -2*ednum
            lcnt = 1
        else:
            msgs.error("Unable to trace any edges"+msgs.newline()+"try a different method to trace the order edges")
    elif rcnt == 0:
        msgs.warn("Unable to find a right edge. Adding one in.")
        # Respecting the BPM (using first column where there is no mask)
        sum_bpm = np.sum(binbpx, axis=0)
        gdi1 = np.max(np.where(sum_bpm == 0)[0])
        # Apply
        edgearr[:, gdi1] = 2*ednum
        rcnt = 1
    elif lcnt == 0:
        msgs.warn("Unable to find a left edge. Adding one in.")
        # Respecting the BPM (using first column where there is no mask)
        sum_bpm = np.sum(binbpx, axis=0)
        gdi0 = np.min(np.where(sum_bpm == 0)[0])
        # Apply
        edgearr[:, gdi0] = -2*ednum
        lcnt = 1
    msgs.info("Assigning slit edge traces")
    # Find the most common set of edges
    edgearrcp = edgearr.copy()
    return edgearrcp, lcnt, rcnt


def edgearr_close_slits(binarr, edgearrcp, edgearr, ednum, settings):
    """ Fuss about with slits that are very close..

    Note that edgearrcp and edgearr are *intentionally*
    swapped from their names in driver_edge_slits.  Yes,
    this is painfully confusing..

    Not extenisvely tested
    NOT Recommended

    Return
    ------
    edgearrcp : ndarray
      Which was modified in place
    """
    vals = np.sort(np.unique(edgearrcp[np.where(edgearrcp != 0)]))
#        print('calling close_edges')
#        exit()
    hasedge = arcytrace.close_edges(edgearrcp, vals, int(settings['trace']['slits']['maxgap']))
    # Find all duplicate edges
    edgedup = vals[np.where(hasedge == 1)]
    if edgedup.size > 0:
        for jj in range(edgedup.size):
            # Raise all remaining edges by one
            if jj != edgedup.size-1:
                edgedup[jj+1:] += 1
            edgearrcp[np.where(edgearrcp > edgedup[jj])] += 1
            # Now investigate the duplicate
            wdup = np.where(edgearrcp == edgedup[jj])
            alldup = edgearr[wdup]
            alldupu = np.unique(alldup)
            cntr = Counter(edg for edg in alldup)
            commn = cntr.most_common(alldupu.size)
            shftsml = np.zeros(len(commn))
            shftarr = np.zeros(wdup[0].size, dtype=np.int)
            wghtarr = np.zeros(len(commn))
            duploc = [None for ii in range(len(commn))]
            for ii in range(len(commn)):
                wghtarr[ii] = commn[ii][1]
                duploc[ii] = np.where(edgearr[wdup] == commn[ii][0])
            changesmade = True
            while changesmade:
                # Keep shifting pixels until the best match is found
                changesmade = False
                # First calculate the old model
                cf = arutils.func_fit(wdup[0], wdup[1]+shftarr,
                                      settings['trace']['slits']['function'],
                                      settings['trace']['slits']['polyorder'],
                                      minv=0, maxv=binarr.shape[0]-1)
                cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                           settings['trace']['slits']['function'],
                                           minv=0, maxv=binarr.shape[0]-1)
                chisqold = np.abs(cenmodl[wdup[0]]-wdup[1]-shftarr).sum()
                for ii in range(1, len(commn)):
                    # Shift by +1
                    adj = np.zeros(wdup[0].size)
                    adj[duploc[ii]] += 1
                    cf = arutils.func_fit(wdup[0], wdup[1]+shftarr+adj,
                                          settings['trace']['slits']['function'],
                                          settings['trace']['slits']['polyorder'],
                                          minv=0, maxv=binarr.shape[0]-1)
                    cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                               settings['trace']['slits']['function'],
                                               minv=0, maxv=binarr.shape[0]-1)
                    chisqp = np.abs(cenmodl[wdup[0]]-wdup[1]-shftarr-adj).sum()
                    # Shift by -1
                    adj = np.zeros(wdup[0].size)
                    adj[duploc[ii]] -= 1
                    cf = arutils.func_fit(wdup[0], wdup[1]+shftarr+adj,
                                          settings['trace']['slits']['function'],
                                          settings['trace']['slits']['polyorder'],
                                          minv=0, maxv=binarr.shape[0]-1)
                    cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                               settings['trace']['slits']['function'],
                                               minv=0, maxv=binarr.shape[0]-1)
                    chisqm = np.abs(cenmodl[wdup[0]]-wdup[1]-shftarr-adj).sum()
                    # Test which solution is best:
                    if chisqold < chisqp and chisqold < chisqm:
                        # No changes are needed
                        continue
                    else:
                        changesmade = True
                        if chisqp < chisqm:
                            shftarr[duploc[ii]] += 1
                            shftsml[ii] += 1
                        else:
                            shftarr[duploc[ii]] -= 1
                            shftsml[ii] -= 1
            # Find the two most common edges
            cntr = Counter(sarr for sarr in shftarr)
            commn = cntr.most_common(2)
            if commn[0][0] > commn[1][0]:  # Make sure that suffix 'a' is assigned the leftmost edge
                wdda = np.where(shftarr == commn[0][0])
                wddb = np.where(shftarr == commn[1][0])
                shadj = commn[0][0] - commn[1][0]
            else:
                wdda = np.where(shftarr == commn[1][0])
                wddb = np.where(shftarr == commn[0][0])
                shadj = commn[1][0] - commn[0][0]
            wvla = np.unique(edgearr[wdup][wdda])
            wvlb = np.unique(edgearr[wdup][wddb])
            # Now generate the dual edge
#                print('calling dual_edge')
#                exit()
            arcytrace.dual_edge(edgearr, edgearrcp, wdup[0], wdup[1], wvla, wvlb, shadj,
                                int(settings['trace']['slits']['maxgap']), edgedup[jj])
    # Now introduce new edge locations
    vals = np.sort(np.unique(edgearrcp[np.where(edgearrcp != 0)]))
#        print('calling close_slits')
#        exit()
    edgearrcp = arcytrace.close_slits(binarr, edgearrcp, vals, int(settings['trace']['slits']['maxgap']),
                                      int(ednum))

    return edgearrcp



def edgearr_final_left_right(edgearr, ednum, siglev):
    """ Final fussing with left/right edges, as needed

    Adds in missing ones, truncates when there are too many of
    one type versus the other

    Parameters
    ----------
    edgearr : ndarray
    ednum : int
    siglev : ndarray

    Returns
    -------
    edgearr : ndarray
    lcnt : int
    rcnt: int
    """
    eaunq = np.unique(edgearr)
    lcnt = np.where(eaunq < 0)[0].size
    rcnt = np.where(eaunq > 0)[0].size
    if lcnt == 0:
        msgs.warn("Unable to find a left edge. Adding one in.")
        edgearr[:, 0] = -2 * ednum
        lcnt = 1
    if rcnt == 0:
        msgs.warn("Unable to find a right edge. Adding one in.")
        edgearr[:, -1] = 2 * ednum
        rcnt = 1
    if (lcnt == 1) & (rcnt > 1):
        msgs.warn("Only one left edge, and multiple right edges.")
        msgs.info("Restricting right edge detection to the most significantly detected edge.")
        wtst = np.where(eaunq > 0)[0]
        bval, bidx = -np.median(siglev[np.where(edgearr == eaunq[wtst[0]])]), 0
        for r in range(1, rcnt):
            wed = np.where(edgearr == eaunq[wtst[r]])
            tstv = -np.median(siglev[wed])
            if tstv > bval:
                bval = tstv
                bidx = r
        edgearr[np.where((edgearr > 0) & (edgearr != eaunq[wtst[bidx]]))] = 0
        edgearr[np.where(edgearr > 0)] = +ednum  # Reset the edge value
        rcnt = 1
    if (lcnt > 1) & (rcnt == 1):
        msgs.warn("Only one right edge, and multiple left edges.")
        msgs.info("Restricting left edge detection to the most significantly detected edge.")
        wtst = np.where(eaunq < 0)[0]
        bval, bidx = np.median(siglev[np.where(edgearr == eaunq[wtst[0]])]), 0
        for r in range(1, lcnt):
            wed = np.where(edgearr == eaunq[wtst[r]])
            tstv = np.median(siglev[wed])
            if tstv > bval:
                bval = tstv
                bidx = r
        edgearr[np.where((edgearr < 0) & (edgearr != eaunq[wtst[bidx]]))] = 0
        edgearr[np.where(edgearr < 0)] = -ednum  # Reset the edge value
        lcnt = 1
    if lcnt == 1:
        letxt = "edge"
    else:
        letxt = "edges"
    if rcnt == 1:
        retxt = "edge"
    else:
        retxt = "edges"
    msgs.info("{0:d} left {1:s} and {2:d} right {3:s} were found in the trace".format(lcnt, letxt, rcnt, retxt))
    return edgearr, lcnt, rcnt




def edgearr_ignore_orders(edgearr, fracignore, ednum):
    """ Ignore orders/slits that run off the edge of the detector
    Mainly used for echelle reductions

    Parameters
    ----------
    edgearr : ndarray
    fracignore : float

    Returns
    -------
    edgearr
    lmin
    lmax
    rmin
    rmax

    """

    iterate = True
    while iterate:
        # Calculate the minimum and maximum left/right edges
        ww = np.where(edgearr < 0)
        lmin, lmax = -np.max(edgearr[ww]), -np.min(edgearr[ww])  # min/max are switched because of the negative signs
        ww = np.where(edgearr > 0)
        rmin, rmax = np.min(edgearr[ww]), np.max(edgearr[ww])  # min/max are switched because of the negative signs
        # msgs.info("Ignoring any slit that spans < {0:3.2f}x{1:d} pixels on the detector".format(settings.argflag['trace']['slits']['fracignore'], int(edgearr.shape[0]*binby)))

        msgs.info("Ignoring any slit that spans < {0:3.2f}x{1:d} pixels on the detector".format(fracignore, int(edgearr.shape[0])))
        fracpix = int(fracignore*edgearr.shape[0])
        msgs.info("Ignoring any slit that spans < {0:3.2f}x{1:d} pixels on the detector".format(
            fracignore, int(edgearr.shape[1])))
        #fracpix = int(fracignore * edgearr.shape[1])
        #msgs.info("Ignoring any slit that spans < {0:3.2f}x{1:d} pixels on the detector".format(
        #    settings['trace']['slits']['fracignore'], int(edgearr.shape[1])))
        __edgearr = edgearr.copy()
        lnc, lxc, rnc, rxc, ldarr, rdarr = new_ignore_orders(__edgearr, fracpix, lmin, lmax, rmin, rmax)
        #        print('New ignore_orders: {0} seconds'.format(time.clock() - t))
        #        assert np.sum(_lnc != lnc) == 0, 'Difference between old and new ignore_orders, lnc'
        #        assert np.sum(_lxc != lxc) == 0, 'Difference between old and new ignore_orders, lxc'
        #        assert np.sum(_rnc != rnc) == 0, 'Difference between old and new ignore_orders, rnc'
        #        assert np.sum(_rxc != rxc) == 0, 'Difference between old and new ignore_orders, rxc'
        #        assert np.sum(_ldarr != ldarr) == 0, 'Difference between old and new ignore_orders, ldarr'
        #        assert np.sum(_rdarr != rdarr) == 0, 'Difference between old and new ignore_orders, rdarr'
        #        assert np.sum(__edgearr != _edgearr) == 0, 'Difference between old and new ignore_orders, edgearr'

        edgearr = __edgearr

        lmin += lnc
        rmin += rnc
        lmax -= lxc
        rmax -= rxc
        iterate = False
        if settings['trace']['slits']['number'] == 1:  # Another check on slits for singleSlit
            if lmax < lmin:
                msgs.warn("Unable to find a left edge2. Adding one in.")
                iterate = True
                edgearr[:, 0] = -2 * ednum
                lcnt = 1
            if rmax < rmin:
                msgs.warn("Unable to find a right edge2. Adding one in.")
                iterate = True
                edgearr[:, -1] = 2 * ednum
                rcnt = 1
    # Return
    return edgearr, lmin, lmax, rmin, rmax

def expand_slits(slf, mstrace, det, ordcen, extord):
    """
    This routine will traces the locations of the slit edges

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    mstrace: numpy ndarray
      Calibration frame that will be used to identify slit edges
    det : int
      Index of the detector
    ordcen : ndarray
      An array providing the physical pixel locations corresponding to the slit centres
    extord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)

    Returns
    -------
    lordloc : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    rordloc : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    """

    # Calculate the pixel locations of th eorder edges
    pixcen = phys_to_pix(ordcen, slf._pixlocn[det - 1], 1)
    msgs.info("Expanding slit traces to slit edges")
#    t = time.clock()
#    _mordwid, _pordwid = arcytrace.expand_slits(mstrace, pixcen, extord.astype(int))
#    print('Old expand_slits: {0} seconds'.format(time.clock() - t))
#    t = time.clock()
    mordwid, pordwid = new_expand_slits(mstrace, pixcen, extord.astype(int))
# TODO: old and new expand_slits do not produce the same result.  There
# was a bug in the old version, but we need to continue to check that
# this version gives good results.
#    print('New expand_slits: {0} seconds'.format(time.clock() - t))
#    assert np.sum(_mordwid != mordwid) == 0, 'Difference between old and new expand_slits, mordwid'
#    assert np.sum(_pordwid != pordwid) == 0, 'Difference between old and new expand_slits, pordwid'

    # Fit a function for the difference between left edge and the centre trace
    ldiff_coeff, ldiff_fit = arutils.polyfitter2d(mordwid, mask=-1,
                                                  order=settings.argflag['trace']['slits']['diffpolyorder'])
    # Fit a function for the difference between left edge and the centre trace
    rdiff_coeff, rdiff_fit = arutils.polyfitter2d(pordwid, mask=-1,
                                                  order=settings.argflag['trace']['slits']['diffpolyorder'])
    lordloc = ordcen - ldiff_fit.T
    rordloc = ordcen + rdiff_fit.T
    return lordloc, rordloc

def fit_edges(edgearr, lmin, lmax, plxbin, plybin, left=True,
              polyorder=3, function='ledgendre'):
    """ Fit the edges, either left or right ones

    Parameters
    ----------
    edgearr : ndarray
    lmin : int
      Minimum edge to fit
    lmax : int
      Maximum edge to fit
    plxbin : ndarray
    plybin : ndarray
    left : bool
      Fitting left (or right) edges?
    polyorder : int, optional
    function : str, optional

    Returns
    -------
    coeff : ndarray
      Fit coefficients (nedge x polyorder)
    nmbrarr : ndarray
      Indexing of the edges
    diffarr : ndarray
    wghtarr : ndarray
    """

    minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]

    # First, determine the model for the most common left slit edge
    if left:
        wcm = np.where(edgearr < 0)
    else:
        wcm = np.where(edgearr > 0)
    cntr = Counter(edg for edg in edgearr[wcm])
    commn = cntr.most_common(1)
    wedx, wedy = np.where(edgearr == commn[0][0])
    msk, cf = arutils.robust_polyfit(wedx, wedy,
                                     polyorder,
                                     function=function,
                                     minv=0, maxv=edgearr.shape[0] - 1)
    cenmodl = arutils.func_val(cf, np.arange(edgearr.shape[0]), function,
                               minv=0, maxv=edgearr.shape[0] - 1)

    if left:
        msgs.info("Fitting left slit traces")
    else:
        msgs.info("Fitting right slit traces")
    coeff = np.zeros((1 + polyorder, lmax - lmin + 1))
    diffarr = np.zeros(lmax - lmin + 1)
    wghtarr = np.zeros(lmax - lmin + 1)
    nmbrarr = np.zeros(lmax - lmin + 1)
    offs = cenmodl[int(edgearr.shape[0] / 2)]

    for i in range(lmin, lmax + 1):
        if left:
            w = np.where(edgearr == -i)
        else:
            w = np.where(edgearr == i)
        if np.size(w[0]) <= polyorder + 2:
            # lfail = np.append(lfail,i-lmin)
            continue
        tlfitx = plxbin[w]
        tlfity = plybin[w]
        diffarr[i - lmin] = np.mean(w[1] - cenmodl[w[0]]) + offs
        wghtarr[i - lmin] = np.size(w[0]) / float(edgearr.shape[0])
        if left:
            nmbrarr[i - lmin] = -i
        else:
            nmbrarr[i - lmin] = i
        msk, coeff[:, i - lmin] = arutils.robust_polyfit(tlfitx, tlfity, polyorder,
                                                          function=function,
                                                          minv=minvf, maxv=maxvf)
    # Return
    return coeff, nmbrarr, diffarr, wghtarr


def get_slitid(shape, lordloc, rordloc, islit, ypos=0.5):
    """ Convert slit position to a slitid
    Parameters
    ----------
    slf : SciExpObj or tuple
    det : int
    islit : int
    ypos : float, optional

    Returns
    -------
    slitid : int
      Slit center position on the detector normalized to range from 0-10000
    slitcen : float
      Slitcenter relative to the detector ranging from 0-1
    xslit : tuple
      left, right positions of the slit edges
    """
    #if isinstance(slf, tuple):
    #    shape, lordloc, rordloc = slf
    #else:
    #    shape = slf._mstrace[det-1].shape
    #    lordloc = slf._lordloc[det-1]
    #    rordloc = slf._rordloc[det-1]
    # Index at ypos
    yidx = int(np.round(ypos*lordloc.shape[0]))
    # Slit at yidx
    pixl_slit = lordloc[yidx, islit]
    pixr_slit = rordloc[yidx, islit]
    # Relative to full image
    xl_slit = pixl_slit/shape[1]
    xr_slit = pixr_slit/shape[1]
    # Center
    slitcen = np.mean([xl_slit, xr_slit])
    slitid = int(np.round(slitcen*1e4))
    # Return them all
    return slitid, slitcen, (xl_slit, xr_slit)


def new_expand_slits(msedge, ordcen, extord):

    t = time.clock()
    sz_x, sz_y = msedge.shape
    sz_o = ordcen.shape[1]

    # Get the pixels at the mid-point between orders
    mid_order = (ordcen[:,:-1] + ordcen[:,1:])//2

    # Instantiate the output
    pordwid = np.zeros(ordcen.shape, dtype=int)
    mordwid = np.zeros(ordcen.shape, dtype=int)

    # Ignore extracted orders
    mordwid[:,extord.astype(bool)] = -1
    pordwid[:,extord.astype(bool)] = -1

    # Set left edges to ignore
    lindx = (mid_order < 0) | (msedge[np.arange(sz_x)[:,None],ordcen[:,1:]] \
                                    < msedge[np.arange(sz_x)[:,None],mid_order])
    lindx = np.append(np.ones(sz_x, dtype=bool).reshape(-1,1), lindx, axis=1)
    mordwid[lindx] = -1

    # Set right edges to ignore
    rindx = (mid_order >= sz_y) | (msedge[np.arange(sz_x)[:,None],ordcen[:,:-1]] \
                                    < msedge[np.arange(sz_x)[:,None],mid_order])
    rindx = np.append(rindx, np.ones(sz_x, dtype=bool).reshape(-1,1), axis=1)
    pordwid[rindx] = -1

    # Find the separation between orders
    medgv = 0.5*(msedge[np.arange(sz_x)[:,None],ordcen[:,1:]] \
                    + msedge[np.arange(sz_x)[:,None],mid_order])
    pedgv = 0.5*(msedge[np.arange(sz_x)[:,None],ordcen[:,:-1]] \
                    + msedge[np.arange(sz_x)[:,None],mid_order])
    for o in range(sz_o):
        for x in range(sz_x):
            # Trace from centre to left
            if mordwid[x,o] != -1:
                mordwid[x,o] = -1
                for y in range(mid_order[x,o-1], ordcen[x, o]):
                    if msedge[x,y] > medgv[x,o-1]:
                        mordwid[x,o] = ordcen[x,o] - y
                        break

            # Trace from centre to right
            if pordwid[x,o] != -1:
                pordwid[x,o] = -1
                for y in range(mid_order[x,o], ordcen[x, o], -1):
                    if msedge[x,y] > pedgv[x,o]:
                        pordwid[x,o] = y-ordcen[x, o]
                        break

    return mordwid, pordwid


def new_find_between(edgdet, ledgem, ledgep, dirc):

    if len(edgdet.shape) != 2:
        msgs.error('Edge pixels array must be 2D.')
    if len(ledgem.shape) != 1 or len(ledgep.shape) !=1:
        msgs.error('Input must be 1D.')

    sz_x, sz_y = edgdet.shape

    # Setup the coefficient arrays
    edgbtwn = np.full(3, -1, dtype=int)

    for x in range(0,sz_x):
        rng = np.sort([ledgem[x],ledgep[x]])
        if not np.any(edgdet[x,slice(*rng)] > 0):
            continue
        e = edgdet[x,slice(*rng)][edgdet[x,slice(*rng)] > 0]
        if edgbtwn[0] == -1:
            edgbtwn[0] = e[0]
        indx = e != edgbtwn[0]
        if edgbtwn[1] == -1 and np.any(indx):
            edgbtwn[1] = e[indx][0]

    # If no right order edges were found between these two left order
    # edges, find the next right order edge
    if edgbtwn[0] == -1 and edgbtwn[1] == -1:
        for x in range(0,sz_x):
            ystrt = np.max([ledgem[x],ledgep[x]])
            if dirc == 1:
                emin = np.min(edgdet[x,ystrt:][edgdet[x,ystrt:] > 0])
                if edgbtwn[2] == -1 or emin < edgbtwn[2]:
                    edgbtwn[2] = emin
            else:
                emax = np.max(edgdet[x,:ystrt+1][edgdet[x,:ystrt+1] > 0])
                if edgbtwn[2] == -1 or emax > edgbtwn[2]:
                    edgbtwn[2] = emax

    # Now return the array
    return edgbtwn


def new_find_shift(mstrace, minarr, lopos, diffarr, numsrch):

    sz_y = mstrace.shape[1]
    maxcnts = -999999.9
    shift = 0
    d = mstrace - minarr[:,None]

    for s in range(0,numsrch):
        cnts = 0.0

        ymin = lopos + s
        ymin[ymin < 0] = 0
        ymax = ymin + diffarr
        ymax[ymax > sz_y] = sz_y

        indx = ymax > ymin

        if np.sum(indx) == 0:
            continue

        cnts = np.sum([ np.sum(t[l:h]) for t,l,h in zip(d[indx], ymin[indx], ymax[indx]) ]) \
                    / np.sum(ymax[indx]-ymin[indx])
        if cnts > maxcnts:
            maxcnts = cnts
            shift = s

    return shift


def new_ignore_orders(edgdet, fracpix, lmin, lmax, rmin, rmax):
    """
    .. warning::

        edgdet is alted by the function.
    """
    sz_x, sz_y = edgdet.shape

    lsize = lmax-lmin+1
    larr = np.zeros((2,lsize), dtype=int)
    larr[0,:] = sz_x

    rsize = rmax-rmin+1
    rarr = np.zeros((2,rsize), dtype=int)
    rarr[0,:] = sz_x

    # TODO: Can I remove the loop?  Or maybe just iterate through the
    # smallest dimension of edgdet?
    for x in range(sz_x):
        indx = edgdet[x,:] < 0
        if np.any(indx):
            larr[0,-edgdet[x,indx]-lmin] = np.clip(larr[0,-edgdet[x,indx]-lmin], None, x)
            larr[1,-edgdet[x,indx]-lmin] = np.clip(larr[1,-edgdet[x,indx]-lmin], x, None)
        indx = edgdet[x,:] > 0
        if np.any(indx):
            rarr[0,edgdet[x,indx]-rmin] = np.clip(rarr[0,edgdet[x,indx]-rmin], None, x)
            rarr[1,edgdet[x,indx]-rmin] = np.clip(rarr[1,edgdet[x,indx]-rmin], x, None)

    # Go through the array once more to remove pixels that do not cover fracpix
    edgdet = edgdet.ravel()
    lt_zero = np.arange(edgdet.size)[edgdet < 0]
    if len(lt_zero) > 0:
        edgdet[lt_zero[larr[1,-edgdet[lt_zero]-lmin]-larr[0,-edgdet[lt_zero]-lmin] < fracpix]] = 0
    gt_zero = np.arange(edgdet.size)[edgdet > 0]
    if len(gt_zero) > 0:
        edgdet[gt_zero[rarr[1,edgdet[gt_zero]-rmin]-rarr[0,edgdet[gt_zero]-rmin] < fracpix]] = 0
    edgdet = edgdet.reshape(sz_x,sz_y)

    # Check if lmin, lmax, rmin, and rmax need to be changed
    lindx = np.arange(lsize)[larr[1,:]-larr[0,:] > fracpix]
    lnc = lindx[0]
    lxc = lsize-1-lindx[-1]

    rindx = np.arange(rsize)[rarr[1,:]-rarr[0,:] > fracpix]
    rnc = rindx[0]
    rxc = rsize-1-rindx[-1]

    return lnc, lxc, rnc, rxc, larr, rarr


def new_limit_yval(yc, maxv):
    yn = 0 if yc == 0 else (-yc if yc < 3 else -3)
    yx = maxv-yc if yc > maxv-4 and yc < maxv else 4
    return yn, yx


def new_match_edges(edgdet, ednum, mr=50):
    """  Label groups of edge pixels and give them
    a unique identifier.

    Parameters
    ----------
    edgdet : ndarray
      Modified in place
    ednum : int
      a large dummy number used for slit edge assignment.
      ednum should be larger than the number of edges detected
    mr : int, optional
      minimum number of acceptable pixels required to form the detection of an order edge
      JXP increased the default value from 5 to 50
         50 is probably best for

    Returns
    -------
    lcnt-2*ednum
    rcnt-2*ednum
    """
    mrxarr = np.zeros(mr, dtype=int) -1  # -1 so as to be off the chip
    mryarr = np.zeros(mr, dtype=int) -1  # -1 so as to be off the chip

    sz_x, sz_y = edgdet.shape

    lcnt = 2*ednum
    rcnt = 2*ednum
    # TODO -- Consider starting at sz_x/2
    # Note:  x=rows and y=columns in the following
    for y in range(sz_y):
        for x in range(sz_x):
            if edgdet[x,y] != -1 and edgdet[x,y] != 1:
                # No edge at this pixel
                continue

            anyt = 0
            left = edgdet[x,y] == -1

            # Search upwards from x,y
            #xs = x + 1  (was x+1)
            xs = x
            yt = y
            while xs <= sz_x-1:
                xr = 10 if xs + 10 < sz_x else sz_x - xs - 1
                yn, yx = new_limit_yval(yt, sz_y)

                suc = 0
                for s in range(xs, xs+xr):
                    suc = 0
                    for t in range(yt + yn, yt + yx):
                        if edgdet[s, t] == -1 and left:
                            edgdet[s, t] = -lcnt
                        elif edgdet[s, t] == 1 and not left:
                            edgdet[s, t] = rcnt
                        else:
                            continue

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
                xr = xs if xs-10 < 0 else 10
                yn, yx = new_limit_yval(yt, sz_y)

                suc = 0
                for s in range(0, xr):
                    suc = 0
                    for t in range(yt+yn, yt+yx):
                        if edgdet[xs-s, t] == -1 and left:
                            edgdet[xs-s, t] = -lcnt
                        elif edgdet[xs-s, t] == 1 and not left:
                            edgdet[xs-s, t] = rcnt
                        else:
                            continue

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

            if anyt > mr and left:
                edgdet[x, y] = -lcnt
                lcnt = lcnt + 1
            elif anyt > mr and not left:
                edgdet[x, y] = rcnt
                rcnt = rcnt + 1
            else:
                edgdet[x, y] = 0
                for s in range(anyt):
                    if mrxarr[s] != -1 and mryarr[s] != -1:
                        edgdet[mrxarr[s], mryarr[s]] = 0

    return lcnt-2*ednum, rcnt-2*ednum




# TODO: Just a 1-1 mapping of the cython function to python.  Needs to
# be tested!!
def new_close_edges(edgdet, dets, npix):
    sz_x, sz_y = edgdet.shape
    sz_d = dets.size

    hasedge = np.zeros(sz_d, dtype=int)

    for d in range(sz_d):
        for x in range(sz_x):
            for y in range(sz_y):
                if edgdet[x,y] != dets[d]:
                    continue
                else:
                    # Check if there's an edge nearby
                    mgap = sz_y if y+npix+1 > sz_y else y+npix+1
                    for s in range(y+1, mgap):
                        if edgdet[x,s] == dets[d]:
                            hasedge[d] = 1
                            break
                if hasedge[d] != 0:
                    break
            if hasedge[d] != 0:
                break
    return hasedge


# TODO: Just a 1-1 mapping of the cython function to python.  Needs to
# be tested!!
def new_close_slits(trframe, edgdet, dets, npix, ednum):

    sz_x, sz_y = edgdet.shape
    sz_d = dets.size

    edgearr = np.zeros(edgdet.shape, dtype=int)
    hasedge = np.zeros(sz_d, dtype=int)

    for d in range(sz_d):
        tmp = sz_y
        for x in range(sz_x):
            for y in range(sz_y):
                if edgdet[x, y] != dets[d]:
                    continue
                else:
                    # Check if there's an edge nearby
                    mgap = sz_y if y+npix+1 > sz_y else y+npix+1
                    for s in range(y+1, mgap):
                        if edgdet[x,s] != 0:
                            if s-y < tmp:
                                tmp = s-y
                                tix = edgdet[x,s]
                            hasedge[d] = edgdet[x,s]
                            break
        if tmp != sz_y:
            hasedge[d] = tix

    # Now, if there's an edge in hasedge, mark the corresponding index
    # in hasedge with -1
    for d in range(sz_d):
        if hasedge[d] == dets[d]:
            # Close slits have caused a left/right edge to be labelled
            # as one edge. Find only instances where there is a left and
            # right edge detection. Then, take their average and set
            # hadedge to be zero
            tmp = 0
            diff = 0
            for x in range(sz_x):
                for y in range(sz_y):
                    if edgdet[x,y] != dets[d]:
                        continue
                    else:
                        # Check if there's an edge nearby
                        mgap = sz_y if y+npix+1 > sz_y else y+npix+1
                        flg = 0
                        for s in range(y+1, mgap):
                            if edgdet[x,s] == edgdet[x,y]:
                                edgdet[x,s] = 0
                                edgdet[x,y] = 0
                                # +0.5 for rounding
                                tix = y + int(0.5*(s-y) + 0.5)
                                edgdet[x,tix] = dets[d]
                                flg = 1
                                tmp += 1
                                diff += (s-y)
                                break
                        if flg == 0:
                            # If there isn't a common left/right edge
                            # for this pixel, ignore this single edge
                            # detection
                            edgdet[x,y] = 0
            hasedge[d] = diff/tmp
            continue
        if hasedge[d] > 0:
            for s in range(sz_d):
                if hasedge[d] == dets[s]:
                    hasedge[s] = -1
                    break

    # Introduce an edge in cases where no edge exists, and redefine an
    # edge where one does exist.
    enum = ednum
    for d in range(sz_d):
        tmp = 0
        for x in range(sz_x):
            for y in range(sz_y):
                if edgdet[x,y] != dets[d]:
                    continue
                if hasedge[d] >= ednum:
                    edgearr[x,y] = enum
                    # Relabel the appropriate hasedge
                    if tmp == 0:
                        for s in range(sz_d):
                            if hasedge[d] == dets[s]:
                                # Label hasedge as negative, to avoid
                                # confusion with the positive hasedge
                                # numbers
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
                    msgs.bug('Check slit traces in close_slits!')
        if hasedge[d] >= 0:
            enum += 1
    # Finally return the new slit edges array
    return edgearr


def new_dual_edge(edgearr, edgearrcp, wx, wy, wl, wr, shft, npix, newval):

    sz_x, sz_y = edgearr.shape
    sz_a = wl.shape[0]
    sz_b = wr.shape[0]
    sz_e = wx.shape[0]

    # First go through the leftmost edge (suffix a)
    for x in range(sz_a):
        for ee in range(sz_e):
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
        for ee in range(sz_e):
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


def new_minbetween(mstrace, loord, hiord):
    # TODO: Check shapes
    ymin = np.clip(loord, 0, mstrace.shape[1])
    ymax = np.clip(hiord, 0, mstrace.shape[1])
    minarr = np.zeros(mstrace.shape[0])
    indx = ymax > ymin
    minarr[indx] = np.array([ np.amin(t[l:h])
                                for t,l,h in zip(mstrace[indx], ymin[indx], ymax[indx]) ])
    return minarr


def new_find_peak_limits(hist, pks):
    """
    Find all values between the zeros of hist

    Parameters
    ----------
    hist : ndarray
      1D vector
    pks : ndarray
      1D vector
    """
    if len(hist.shape) != 1 or len(pks.shape) != 1:
        msgs.error('Arrays provided to find_peak_limits must be vectors.')
    # Pixel indices in hist for each peak
    hn = np.arange(hist.shape[0])
    indx = np.ma.MaskedArray(np.array([hn]*pks.shape[0]))
    # Instantiate output
    edges = np.zeros((pks.shape[0],2), dtype=int)
    # Find the left edges
    indx.mask = (hist != 0)[None,:] | (hn[None,:] > pks[:,None])
    edges[:,0] = np.ma.amax(indx, axis=1)
    # Find the right edges
    indx.mask = (hist != 0)[None,:] | (hn[None,:] < pks[:,None])
    edges[:,1] = np.ma.amin(indx, axis=1)
    return edges



def pca_order_slit_edges(binarr, edgearr, lcent, rcent, gord,
                   lcoeff, rcoeff, plxbin, slitcen, pixlocn, settings):
    """ Perform a PCA analyis on the order edges
    Primarily for extrapolation

    Parameters
    ----------
    binarr : ndarray
    edgearr : ndarray
    lcent : ndarray
      Left edges
    rcent : ndarray
      Right edges
    gord : ndarray
      Orders detected on both the left and right edge
    lcoeff : ndarray
      Fit coefficients for left edges
    rcoeff : ndarray
      Fit coefficients for right edges
    plxbin
    slitcen : ndarray
    pixlocn

    Returns
    -------

    """
    # Init
    wl = np.where(edgearr < 0)
    wr = np.where(edgearr > 0)

    ##############
    xv = plxbin[:, 0]
    minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]

    almin, almax = -np.max(edgearr[wl]), -np.min(edgearr[wl]) # min and max switched because left edges have negative values
    armin, armax = np.min(edgearr[wr]), np.max(edgearr[wr])

    # maskord = np.where((np.all(lcoeff[:,lg],axis=0)==False)|(np.all(rcoeff[:,rg],axis=0)==False))[0]
    maskord = np.where((np.all(lcoeff, axis=0) == False) | (np.all(rcoeff, axis=0) == False))[0]
    ordsnd = np.arange(min(almin, armin), max(almax, armax) + 1)
    totord = ordsnd[-1] + settings['trace']['slits']['pca']['extrapolate']['pos']
    # Identify the orders to be extrapolated during reconstruction
    extrapord = (1.0 - np.in1d(np.linspace(1.0, totord, totord), gord).astype(np.int)).astype(np.bool)
    msgs.info("Performing a PCA on the order edges")
    ofit = settings['trace']['slits']['pca']['params']
    lnpc = len(ofit) - 1
    msgs.work("May need to do a check here to make sure ofit is reasonable")
    coeffs = arutils.func_fit(xv, slitcen, settings['trace']['slits']['function'],
                              settings['trace']['slits']['polyorder'], minv=minvf, maxv=maxvf)
    for i in range(ordsnd.size):
        if i in maskord:
            if (i>=ordsnd[0]) and (i<ordsnd[-1]-1):  # JXP: Don't add orders that are already in there
                continue
            coeffs = np.insert(coeffs, i, 0.0, axis=1)
            slitcen = np.insert(slitcen, i, 0.0, axis=1)
            lcent = np.insert(lcent, i, 0.0, axis=0)
            rcent = np.insert(rcent, i, 0.0, axis=0)
    xcen = xv[:, np.newaxis].repeat(ordsnd.size, axis=1)
    fitted, outpar = arpca.basis(xcen, slitcen, coeffs, lnpc, ofit, x0in=ordsnd, mask=maskord,
                                 skipx0=False, function=settings['trace']['slits']['function'])
    if not msgs._debug['no_qa']:
        debugger.set_trace()  # NEED TO REMOVE slf
        # arpca.pca_plot(slf, outpar, ofit, "Slit_Trace", pcadesc=pcadesc)
    # Extrapolate the remaining orders requested
    orders = 1 + np.arange(totord)
    extrap_cent, outpar = arpca.extrapolate(outpar, orders, function=settings['trace']['slits']['function'])
    # Fit a function for the difference between left and right edges.
    diff_coeff, diff_fit = arutils.polyfitter2d(rcent - lcent, mask=maskord,
                                                order=settings['trace']['slits']['diffpolyorder'])
    # Now extrapolate the order difference
    ydet = np.linspace(0.0, 1.0, lcent.shape[0])
    ydetd = ydet[1] - ydet[0]
    lnum = ordsnd[0] - 1.0
    ydet = np.append(-ydetd * np.arange(1.0, 1.0 + lnum)[::-1], ydet)
    ydet = np.append(ydet,
                     1.0 + ydetd * np.arange(1.0, 1.0 + settings['trace']['slits']['pca']['extrapolate']['pos']))
    xde, yde = np.meshgrid(np.linspace(0.0, 1.0, lcent.shape[1]), ydet)
    extrap_diff = arutils.polyval2d(xde, yde, diff_coeff).T
    msgs.info("Refining the trace for reconstructed and predicted orders")
    # NOTE::  MIGHT NEED TO APPLY THE BAD PIXEL MASK HERE TO BINARR
    msgs.work("Should the bad pixel mask be applied to the frame here?")
    refine_cent, outpar = refine_traces(binarr, outpar, extrap_cent, extrap_diff,
                                        [gord[0] - orders[0], orders[-1] - gord[-1]], orders,
                                        ofit[0], pixlocn,
                                        function=settings['trace']['slits']['function'])
    # Generate the left and right edges
    lcen = refine_cent - 0.5 * extrap_diff
    rcen = refine_cent + 0.5 * extrap_diff

    # Return
    return lcen, rcen, extrapord


def pca_pixel_slit_edges(binarr, edgearr, lcoeff, rcoeff, ldiffarr, rdiffarr,
                         lnmbrarr, rnmbrarr, lwghtarr, rwghtarr, lcent, rcent, plxbin,
                         settings):
    """ PCA analysis for slit edges

    Parameters
    ----------
    binarr : ndarray
    edgearr : ndarray
    lcoeff : ndarray
      Fit coefficients for left edges
    rcoeff : ndarray
      Fit coefficients for right edges
    ldiffarr : ndarray
    rdiffarr : ndarray
    lnmbrarr : ndarray
    rnmbrarr : ndarray
    lwghtarr : ndarray
    rwghtarr : ndarray
    lcent : ndarray
    rcent : ndarray
    plxbin
    settings : dict-like object

    Returns
    -------
    lcen : ndarray
      Left edges
    rcen : ndarray
      Right edges
    extrapord

    """

    minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]
    maskord = np.where((np.all(lcoeff, axis=0) == False) | (np.all(rcoeff, axis=0) == False))[0]
    allord = np.arange(ldiffarr.shape[0])
    ww = np.where(np.in1d(allord, maskord) == False)[0]
    # Unmask where an order edge is located
    maskrows = np.ones(binarr.shape[1], dtype=np.int)
    # ldiffarr = np.round(ldiffarr[ww]).astype(np.int)
    # rdiffarr = np.round(rdiffarr[ww]).astype(np.int)
    ldiffarr = np.fmax(np.fmin(np.round(ldiffarr[ww]).astype(np.int), binarr.shape[1] - 1), 0)
    rdiffarr = np.fmax(np.fmin(np.round(rdiffarr[ww]).astype(np.int), binarr.shape[1] - 1), 0)
    maskrows[ldiffarr] = 0
    maskrows[rdiffarr] = 0
    # Extract the slit edge ID numbers associated with the acceptable traces
    lnmbrarr = lnmbrarr[ww]
    rnmbrarr = rnmbrarr[ww]
    # Fill in left/right coefficients
    tcoeff = np.ones((settings['trace']['slits']['polyorder'] + 1, binarr.shape[1]))
    tcoeff[:, ldiffarr] = lcoeff[:, ww]
    tcoeff[:, rdiffarr] = rcoeff[:, ww]
    # Weight the PCA fit by the number of detections in each slit edge
    pxwght = np.zeros(binarr.shape[1])
    pxwght[ldiffarr] = lwghtarr[ww]
    pxwght[rdiffarr] = rwghtarr[ww]
    maskrw = np.where(maskrows == 1)[0]
    maskrw.sort()
    extrap_row = maskrows.copy()
    xv = np.arange(binarr.shape[0])
    # trace values
    trcval = arutils.func_val(tcoeff, xv, settings['trace']['slits']['function'],
                              minv=minvf, maxv=maxvf).T
    msgs.work("May need to do a check here to make sure ofit is reasonable")
    ofit = settings['trace']['slits']['pca']['params']
    lnpc = len(ofit) - 1
    # Only do a PCA if there are enough good slits
    if np.sum(1.0 - extrap_row) > ofit[0] + 1:
        # Perform a PCA on the locations of the slits
        msgs.info("Performing a PCA on the slit traces")
        ordsnd = np.arange(binarr.shape[1])
        xcen = xv[:, np.newaxis].repeat(binarr.shape[1], axis=1)
        fitted, outpar = arpca.basis(xcen, trcval, tcoeff, lnpc, ofit, weights=pxwght,
                                     x0in=ordsnd, mask=maskrw, skipx0=False,
                                     function=settings['trace']['slits']['function'])
        if not msgs._debug['no_qa']:
            # JXP -- NEED TO REMOVE SLF FROM THE NEXT BIT
            msgs.warn("NEED TO REMOVE SLF FROM THE NEXT BIT")
            # arpca.pca_plot(slf, outpar, ofit, "Slit_Trace", pcadesc=pcadesc, addOne=False)
        # Now extrapolate to the whole detector
        pixpos = np.arange(binarr.shape[1])
        extrap_trc, outpar = arpca.extrapolate(outpar, pixpos,
                                               function=settings['trace']['slits']['function'])
        # Extract the resulting edge traces
        lcen = extrap_trc[:, ldiffarr]
        rcen = extrap_trc[:, rdiffarr]
        # Perform a final shift fit to ensure the traces closely follow the edge detections
        for ii in range(lnmbrarr.size):
            wedx, wedy = np.where(edgearr == lnmbrarr[ii])
            shft = np.mean(lcen[wedx, ii] - wedy)
            lcen[:, ii] -= shft
        for ii in range(rnmbrarr.size):
            wedx, wedy = np.where(edgearr == rnmbrarr[ii])
            shft = np.mean(rcen[wedx, ii] - wedy)
            rcen[:, ii] -= shft
    else:
        allord = np.arange(lcent.shape[0])
        maskord = np.where((np.all(lcent, axis=1) == False) | (np.all(rcent, axis=1) == False))[0]
        ww = np.where(np.in1d(allord, maskord) == False)[0]
        lcen = lcent[ww, :].T.copy()
        rcen = rcent[ww, :].T.copy()
    extrapord = np.zeros(lcen.shape[1], dtype=np.bool)

    # Return
    return lcen, rcen, extrapord




def prune_peaks(hist, pks, pkidx, debug=False):
    """ Identify the most well defined peaks
    And prune peaks too close to one another

    Parameters
    ----------
    hist : ndarray
      Histogram of detections
    pks : ndarray
      Indices of candidate peak locations
    pkidx : int
      Index of highest peak

    Returns
    -------
    msk : ndarray
      An mask of good peaks (1) and bad peaks (0)
    """

    sz_i = pks.shape[0]

    msk = np.zeros(sz_i, dtype=np.int)

    lgd = 1  # Was the previously inspected peak a good one?
    for ii in range(0, sz_i-1):
        cnt = 0
        for jj in range(pks[ii], pks[ii+1]):
            if hist[jj] == 0:
                cnt += 1
        if cnt < 2:  # Two peaks too close to each other.  Should we really eliminate both??
            msk[ii] = 1  # JXP modifies this from 0 to 1
            msk[ii+1] = 0
            lgd = 0
        else:
            # If the difference is acceptable, the right peak is acceptable,
            # the left peak is acceptable if it was not previously labelled as unacceptable
            if lgd == 1:
                msk[ii] = 1
            msk[ii+1] = 1
            lgd = 1
    if debug:
        debugger.set_trace()

    # Now only consider the peaks closest to the highest peak
    lgd = 1
    # If the highest peak was zeroed out, this will zero out everyone!
    for ii in range(pkidx, sz_i):
        if msk[ii] == 0:
            lgd = 0
        elif lgd == 0:
            msk[ii] = 0
    lgd = 1
    # If the highest peak was zeroed out, this will zero out everyone!
    for ii in range(0, pkidx):
        if msk[pkidx-ii] == 0:
            lgd = 0
        elif lgd == 0:
            msk[pkidx-ii] = 0

    if debug:
        debugger.set_trace()

    return msk


def refine_traces(binarr, outpar, extrap_cent, extrap_diff, extord, orders,
                  fitord, locations, function='polynomial'):
    """
    Parameters
    ----------
    binarr
    outpar
    extrap_cent
    extrap_diff
    extord
    orders
    fitord
    locations
    function : str, optional

    Returns
    -------

    """
    # Refine the orders in the positive direction
    i = extord[1]
    hiord = arpixels.phys_to_pix(extrap_cent[:, -i-2], locations, 1)
    nxord = arpixels.phys_to_pix(extrap_cent[:, -i-1], locations, 1)
    mask = np.ones(orders.size)
    mask[0:extord[0]] = 0.0
    mask[-extord[1]:] = 0.0
    extfit = extrap_cent.copy()
    outparcopy = copy.deepcopy(outpar)
    while i > 0:
        loord = hiord
        hiord = nxord
        nxord = arpixels.phys_to_pix(extrap_cent[:,-i], locations, 1)

        # Minimum counts between loord and hiord
#        print('calling minbetween')
#        t = time.clock()
#        _minarrL = arcytrace.minbetween(binarr, loord, hiord)
#        _minarrR = arcytrace.minbetween(binarr, hiord, nxord)
#        print('Old minbetween: {0} seconds'.format(time.clock() - t))
#        t = time.clock()
        minarrL = new_minbetween(binarr, loord, hiord)
        minarrR = new_minbetween(binarr, hiord, nxord)
#        print('New minbetween: {0} seconds'.format(time.clock() - t))
#        assert np.sum(_minarrL != minarrL) == 0, \
#                'Difference between old and new minbetween, minarrL'
#        assert np.sum(_minarrR != minarrR) == 0, \
#                'Difference between old and new minbetween, minarrR'

        minarr = 0.5*(minarrL+minarrR)
        srchz = np.abs(extfit[:,-i]-extfit[:,-i-1])/3.0
        lopos = arpixels.phys_to_pix(extfit[:,-i]-srchz, locations, 1)  # The pixel indices for the bottom of the search window
        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,-i])))
        diffarr = np.round(extrap_diff[:,-i]).astype(np.int)

#        print('calling find_shift')
#        t = time.clock()
#        _shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch)
#        print('Old find_shift: {0} seconds'.format(time.clock() - t))
#        t = time.clock()
        shift = new_find_shift(binarr, minarr, lopos, diffarr, numsrch)
#        print('New find_shift: {0} seconds'.format(time.clock() - t))
#        assert np.sum(_shift != shift) == 0, 'Difference between old and new find_shift, shift'

        relshift = np.mean(shift+extrap_diff[:,-i]/2-srchz)
        if shift == -1:
            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[-i])))
            relshift = 0.0
        else:
            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[-i]), relshift))
        # Renew guess for the next order
        mask[-i] = 1.0
        extfit, outpar, fail = arpca.refine_iter(outpar, orders, mask, -i, relshift, fitord, function=function)
        if fail:
            msgs.warn("Order refinement has large residuals -- check order traces")
            return extrap_cent, outparcopy
        i -= 1
    # Refine the orders in the negative direction
    i = extord[0]
    loord = arpixels.phys_to_pix(extrap_cent[:,i+1], locations, 1)
    extrap_cent = extfit.copy()
    outparcopy = copy.deepcopy(outpar)
    while i > 0:
        hiord = loord
        loord = arpixels.phys_to_pix(extfit[:,i], locations, 1)

#        print('calling minbetween')
#        t = time.clock()
#        _minarr = arcytrace.minbetween(binarr,loord, hiord)
#        print('Old minbetween: {0} seconds'.format(time.clock() - t))
#        t = time.clock()
        minarr = new_minbetween(binarr,loord, hiord)
#        print('New minbetween: {0} seconds'.format(time.clock() - t))
#        assert np.sum(_minarr != minarr) == 0, 'Difference between old and new minbetween, minarr'

        srchz = np.abs(extfit[:,i]-extfit[:,i-1])/3.0
        lopos = arpixels.phys_to_pix(extfit[:,i-1]-srchz, locations, 1)
        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,i-1])))
        diffarr = np.round(extrap_diff[:,i-1]).astype(np.int)

#        print('calling find_shift')
#        t = time.clock()
#        _shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch)
#        print('Old find_shift: {0} seconds'.format(time.clock() - t))
#        t = time.clock()
        shift = new_find_shift(binarr, minarr, lopos, diffarr, numsrch)
#        print('New find_shift: {0} seconds'.format(time.clock() - t))
#        assert np.sum(_shift != shift) == 0, 'Difference between old and new find_shift, shift'

        relshift = np.mean(shift+extrap_diff[:,i-1]/2-srchz)
        if shift == -1:
            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[i-1])))
            relshift = 0.0
        else:
            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[i-1]),relshift))
        # Renew guess for the next order
        mask[i-1] = 1.0
        extfit, outpar, fail = arpca.refine_iter(outpar, orders, mask, i-1, relshift, fitord, function=function)
        if fail:
            msgs.warn("Order refinement has large residuals -- check order traces")
            return extrap_cent, outparcopy
        i -= 1
    return extfit, outpar


def remove_slit(edgearr, lcen, rcen, tc_dict, rm_slits, TOL=3.):
    """ Remove slit

    Parameters
    ----------
    edgearr : ndarray
    lcen : ndarray
    rcen : ndarray
    tc_dict : dict
    rm_slits : list
      List of slits to remove
        [[left0, right0], [left1, right1]]
      Specified at ycen = nrows//2
    TOL : float
      Tolerance for specifying the left/right edge

    Returns
    -------
    edgearr : ndarray
    lcen
    rcen
    tc_dict

    """
    # Grab the edge indexes and xval's
    left_idx, left_xval, right_idx, right_xval = tc_indices(tc_dict)

    # Final edges
    ycen = lcen.shape[0] // 2
    lcen_yc = lcen[ycen,:]
    rcen_yc = rcen[ycen,:]
    msk = np.ones(lcen.shape[1], dtype=bool)

    # Loop on the slits to remove
    for rm_slit in rm_slits:
        left, right = rm_slit
        # Left check
        if (np.min(np.abs(left_xval-left)) < TOL) & (np.min(np.abs(lcen_yc-left)) < TOL):
            ileft = np.argmin(np.abs(left_xval-left))
            ilcen = np.argmin(np.abs(lcen_yc-left))
        else:
            msgs.warn("Could not find a left slit corresponding to {}".format(left))
            return
        # Right check
        if (np.min(np.abs(right_xval-right)) < TOL) & (np.min(np.abs(rcen_yc-right)) < TOL):
            iright = np.argmin(np.abs(right_xval-right))
            ircen = np.argmin(np.abs(rcen_yc-right))
            if ilcen != ircen:
                msgs.warn("lcen not in sync with rcen or you misdefined the slit to remove")
                return
        else:
            msgs.warn("Could not find a right slit corresponding to {}".format(right))
            return
        # Remove from final edges -- Am not sure these will be indexed identically to tc_dict..
        msk[ilcen] = False
        # Remove from edgearr
        edgearr[edgearr == left_idx[ileft]] = 0.
        edgearr[edgearr == right_idx[iright]] = 0.
        tc_dict['left']['xval'].pop(str(left_idx[ileft]))
        tc_dict['right']['xval'].pop(str(right_idx[iright]))
        msgs.info("Removed the slit at [left,right]: {}".format(rm_slit))

    # Do I need to reindex everything??
    lcen = lcen[:,msk]
    rcen = rcen[:,msk]

    # Return
    return edgearr, lcen, rcen, tc_dict

def synchronize_edges(binarr, edgearr, plxbin, lmin, lmax, lcoeff, rmin, rcoeff,
                      lnmbrarr, ldiffarr, lwghtarr, rnmbrarr, rdiffarr, rwghtarr, settings):
    """ Synchrnoizes the existing edges

    For ARMLSD, this step is largely unnecessary given multi_sync()

    @Ryan :: Could use help making this method less onerous..

    Parameters
    ----------
    binarr
    edgearr
    plxbin
    lmin
    lmax
    lcoeff
    rmin
    rcoeff
    lnmbrarr
    ldiffarr
    lwghtarr
    rnmbrarr
    rdiffarr
    rwghtarr
    settings

    Returns
    -------

    """
    minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]
    # Define the array of pixel values along the dispersion direction
    xv = plxbin[:, 0]
    num = (lmax - lmin) // 2
    lval = lmin + num  # Pick an order, somewhere in between lmin and lmax
    lv = (arutils.func_val(lcoeff[:, lval - lmin], xv, settings['trace']['slits']['function'], minv=minvf,
                           maxv=maxvf) + 0.5).astype(np.int)
    if np.any(lv < 0) or np.any(lv + 1 >= binarr.shape[1]):
        msgs.warn("At least one slit is poorly traced")
        msgs.info("Refer to the manual, and adjust the input trace parameters")
        msgs.error("Cannot continue without a successful trace")
    mnvalp = np.median(binarr[:, lv + 1])  # Go one row above and one row below an order edge,
    mnvalm = np.median(binarr[:, lv - 1])  # then see which mean value is greater.

    """
    lvp = (arutils.func_val(lcoeff[:,lval+1-lmin],xv,settings.argflag['trace']['slits']['function'],min=minvf,max=maxvf)+0.5).astype(np.int)
    edgbtwn = arcytrace.find_between(edgearr,lv,lvp,1)
    print lval, edgbtwn
    # edgbtwn is a 3 element array that determines what is between two adjacent left edges
    # edgbtwn[0] is the next right order along, from left order lval
    # edgbtwn[1] is only !=-1 when there's an order overlap.
    # edgebtwn[2] is only used when a left order is found before a right order
    if edgbtwn[0] == -1 and edgbtwn[1] == -1:
        rsub = edgbtwn[2]-(lval) # There's an order overlap
    elif edgbtwn[1] == -1: # No overlap
        rsub = edgbtwn[0]-(lval)
    else: # There's an order overlap
        rsub = edgbtwn[1]-(lval)
    """
    if msgs._debug['trace']:
        debugger.set_trace()
    if mnvalp > mnvalm:
        lvp = (arutils.func_val(lcoeff[:, lval + 1 - lmin], xv, settings['trace']['slits']['function'],
                                minv=minvf, maxv=maxvf) + 0.5).astype(np.int)
        #        t = time.clock()
        #        _edgbtwn = arcytrace.find_between(edgearr, lv, lvp, 1)
        #        print('Old find_between: {0} seconds'.format(time.clock() - t))
        #        t = time.clock()
        edgbtwn = new_find_between(edgearr, lv, lvp, 1)
        #        print('New find_between: {0} seconds'.format(time.clock() - t))
        #        assert np.sum(_edgbtwn != edgbtwn) == 0, 'Difference between old and new find_between'

        # edgbtwn is a 3 element array that determines what is between two adjacent left edges
        # edgbtwn[0] is the next right order along, from left order lval
        # edgbtwn[1] is only !=-1 when there's an order overlap.
        # edgebtwn[2] is only used when a left order is found before a right order
        if edgbtwn[0] == -1 and edgbtwn[1] == -1:
            rsub = edgbtwn[2] - lval  # There's an order overlap
        elif edgbtwn[1] == -1:  # No overlap
            rsub = edgbtwn[0] - lval
        else:  # There's an order overlap
            rsub = edgbtwn[1] - lval
    else:
        lvp = (arutils.func_val(lcoeff[:, lval - 1 - lmin], xv, settings['trace']['slits']['function'],
                                minv=minvf, maxv=maxvf) + 0.5).astype(np.int)
        #        t = time.clock()
        #        _edgbtwn = arcytrace.find_between(edgearr, lvp, lv, -1)
        #        print('Old find_between: {0} seconds'.format(time.clock() - t))
        #        t = time.clock()
        edgbtwn = new_find_between(edgearr, lvp, lv, -1)
        #        assert np.sum(_edgbtwn != edgbtwn) == 0, 'Difference between old and new find_between'
        #        print('New find_between: {0} seconds'.format(time.clock() - t))

        if edgbtwn[0] == -1 and edgbtwn[1] == -1:
            rsub = edgbtwn[2] - (lval - 1)  # There's an order overlap
        elif edgbtwn[1] == -1:  # No overlap
            rsub = edgbtwn[0] - (lval - 1)
        else:  # There's an order overlap
            rsub = edgbtwn[1] - (lval - 1)

    msgs.info("Relabelling slit edges")
    rsub = int(round(rsub))
    if lmin < rmin - rsub:
        esub = lmin - (settings['trace']['slits']['pca']['extrapolate']['neg'] + 1)
    else:
        esub = (rmin - rsub) - (settings['trace']['slits']['pca']['extrapolate']['neg'] + 1)

    wl = np.where(edgearr < 0)
    wr = np.where(edgearr > 0)
    edgearr[wl] += esub
    edgearr[wr] -= (esub + rsub)
    lnmbrarr += esub
    rnmbrarr -= (esub + rsub)

    # Insert new rows into coefficients arrays if rsub != 0 (if orders were not labelled correctly, there will be a mismatch for the lcoeff and rcoeff)
    almin, almax = -np.max(edgearr[wl]), -np.min(
        edgearr[wl])  # min and max switched because left edges have negative values
    armin, armax = np.min(edgearr[wr]), np.max(edgearr[wr])
    nmord = settings['trace']['slits']['polyorder'] + 1
    if armin != almin:
        if armin < almin:
            lcoeff = np.append(np.zeros((nmord, almin - armin)), lcoeff, axis=1)
            ldiffarr = np.append(np.zeros(almin - armin), ldiffarr)
            lnmbrarr = np.append(np.zeros(almin - armin), lnmbrarr)
            lwghtarr = np.append(np.zeros(almin - armin), lwghtarr)
        else:
            rcoeff = np.append(np.zeros((nmord, armin - almin)), rcoeff, axis=1)
            rdiffarr = np.append(np.zeros(armin - almin), rdiffarr)
            rnmbrarr = np.append(np.zeros(armin - almin), rnmbrarr)
            rwghtarr = np.append(np.zeros(armin - almin), rwghtarr)
    if armax != almax:
        if armax < almax:
            rcoeff = np.append(rcoeff, np.zeros((nmord, almax - armax)), axis=1)
            rdiffarr = np.append(rdiffarr, np.zeros(almax - armax))
            rnmbrarr = np.append(rnmbrarr, np.zeros(almax - armax))
            rwghtarr = np.append(rwghtarr, np.zeros(almax - armax))
        else:
            lcoeff = np.append(lcoeff, np.zeros((nmord, armax - almax)), axis=1)
            ldiffarr = np.append(ldiffarr, np.zeros(armax - almax))
            lnmbrarr = np.append(lnmbrarr, np.zeros(armax - almax))
            lwghtarr = np.append(lwghtarr, np.zeros(armax - almax))

    # import astropy.io.fits as pyfits        minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]
    # hdu = pyfits.PrimaryHDU(edgearr)
    # hdu.writeto("edgearr_{0:02d}.fits".format(det))

    # Now consider traces where both the left and right edges are detected
    ordunq = np.unique(edgearr)
    lunqt = ordunq[np.where(ordunq < 0)[0]]
    runqt = ordunq[np.where(ordunq > 0)[0]]
    lunq = np.arange(lunqt.min(), lunqt.max() + 1)
    runq = np.arange(runqt.min(), runqt.max() + 1)
    # Determine which orders are detected on both the left and right edge
    gord = np.intersect1d(-lunq, runq, assume_unique=True)
    # We need to ignore the orders labelled rfail and lfail.
    lg = np.where(np.in1d(-lunq, gord))[0]
    rg = np.where(np.in1d(runq, gord))[0]
    lgm = np.where(np.in1d(-lunq, gord, invert=True))[0]
    rgm = np.where(np.in1d(runq, gord, invert=True))[0]
    maxord = np.max(np.append(gord, np.append(-lunq[lgm], runq[rgm])))
    lcent = arutils.func_val(lcoeff[:, -lunq[lg][::-1] - 1 - settings['trace']['slits']['pca']['extrapolate']['neg']],
                             xv,
                             settings['trace']['slits']['function'], minv=minvf, maxv=maxvf)
    rcent = arutils.func_val(rcoeff[:, runq[rg] - 1 - settings['trace']['slits']['pca']['extrapolate']['neg']], xv,
                             settings['trace']['slits']['function'], minv=minvf, maxv=maxvf)
    # Return
    return lcent, rcent, gord, lcoeff, ldiffarr, lnmbrarr, lwghtarr, rcoeff, rdiffarr, rnmbrarr, rwghtarr


def trace_crude_init(image, xinit0, ypass, invvar=None, radius=2.,
    maxshift0=0.5, maxshift=0.15, maxerr=0.2):
    """Python port of trace_crude_idl.pro from IDLUTILS

    Modified for initial guess

    Parameters
    ----------
    image : 2D ndarray
      Image for tracing
    xinit : ndarray
      Initial guesses for trace peak at ypass
    ypass : int
      Row for initial guesses

    Returns
    -------
    xset : Trace for each fiber
    xerr : Estimated error in that trace
    """
    # Init
    xinit = xinit0.astype(float)
    #xinit = xinit[0:3]
    ntrace = xinit.size
    ny = image.shape[0]
    xset = np.zeros((ny,ntrace))
    xerr = np.zeros((ny,ntrace))
    if invvar is None:
        invvar = np.zeros_like(image) + 1.

    #
    #  Recenter INITIAL Row for all traces simultaneously
    #
    iy = ypass * np.ones(ntrace,dtype=int)
    xfit,xfiterr = trace_fweight(image, xinit, iy, invvar=invvar, radius=radius)
    # Shift
    xshift = np.clip(xfit-xinit, -1*maxshift0, maxshift0) * (xfiterr < maxerr)
    xset[ypass,:] = xinit + xshift
    xerr[ypass,:] = xfiterr * (xfiterr < maxerr)  + 999.0 * (xfiterr >= maxerr)

    #    /* LOOP FROM INITIAL (COL,ROW) NUMBER TO LARGER ROW NUMBERS */
    for iy in range(ypass+1, ny):
        xinit = xset[iy-1, :]
        ycen = iy * np.ones(ntrace,dtype=int)
        xfit,xfiterr = trace_fweight(image, xinit, ycen, invvar=invvar, radius=radius)
        # Shift
        xshift = np.clip(xfit-xinit, -1*maxshift, maxshift) * (xfiterr < maxerr)
        # Save
        xset[iy,:] = xinit + xshift
        xerr[iy,:] = xfiterr * (xfiterr < maxerr)  + 999.0 * (xfiterr >= maxerr)
    #      /* LOOP FROM INITIAL (COL,ROW) NUMBER TO SMALLER ROW NUMBERS */
    for iy in range(ypass-1, -1,-1):
        xinit = xset[iy+1, :]
        ycen = iy * np.ones(ntrace,dtype=int)
        xfit,xfiterr = trace_fweight(image, xinit, ycen, invvar=invvar, radius=radius)
        # Shift
        xshift = np.clip(xfit-xinit, -1*maxshift, maxshift) * (xfiterr < maxerr)
        # Save
        xset[iy,:] = xinit + xshift
        xerr[iy,:] = xfiterr * (xfiterr < maxerr) + 999.0 * (xfiterr >= maxerr)

    return xset, xerr


def trace_fweight(fimage, xinit, ycen=None, invvar=None, radius=2., debug=False):
    '''Python port of trace_fweight.pro from IDLUTILS

    Parameters
    ----------
    fimage: 2D ndarray
      Image for tracing
    xinit: ndarray
      Initial guesses for x-trace
    invvar: ndarray, optional
      Inverse variance array for the image
    radius: float, optional
      Radius for centroiding; default to 3.0
    '''
    # Definitions for Cython
    #cdef int nx,ny,ncen

    # Init
    nx = fimage.shape[1]
    ny = fimage.shape[0]
    ncen = len(xinit)
    # Create xnew, xerr
    xnew = xinit.astype(float)
    xerr = np.zeros(ncen) + 999.

    # ycen
    if ycen is None:
        if ncen != ny:
            raise ValueError('Bad input')
        ycen = np.arange(ny, dtype=int)
    else:
        if len(ycen) != ncen:
            raise ValueError('Bad ycen input.  Wrong length')
    x1 = xinit - radius + 0.5
    x2 = xinit + radius + 0.5
    ix1 = np.floor(x1).astype(int)
    ix2 = np.floor(x2).astype(int)

    fullpix = int(np.maximum(np.min(ix2-ix1)-1,0))
    sumw = np.zeros(ncen)
    sumxw = np.zeros(ncen)
    sumwt = np.zeros(ncen)
    sumsx1 = np.zeros(ncen)
    sumsx2 = np.zeros(ncen)
    qbad = np.array([False]*ncen)

    if invvar is None:
        invvar = np.zeros_like(fimage) + 1.

    # Compute
    for ii in range(0,fullpix+3):
        spot = ix1 - 1 + ii
        ih = np.clip(spot,0,nx-1)
        xdiff = spot - xinit
        #
        wt = np.clip(radius - np.abs(xdiff) + 0.5,0,1) * ((spot >= 0) & (spot < nx))
        sumw = sumw + fimage[ycen,ih] * wt
        sumwt = sumwt + wt
        sumxw = sumxw + fimage[ycen,ih] * xdiff * wt
        var_term = wt**2 / (invvar[ycen,ih] + (invvar[ycen,ih] == 0))
        sumsx2 = sumsx2 + var_term
        sumsx1 = sumsx1 + xdiff**2 * var_term
        #qbad = qbad or (invvar[ycen,ih] <= 0)
        qbad = np.any([qbad, invvar[ycen,ih] <= 0], axis=0)

    # Fill up
    good = (sumw > 0) &  (~qbad)
    if np.sum(good) > 0:
        delta_x = sumxw[good]/sumw[good]
        xnew[good] = delta_x + xinit[good]
        xerr[good] = np.sqrt(sumsx1[good] + sumsx2[good]*delta_x**2)/sumw[good]

    bad = np.any([np.abs(xnew-xinit) > radius + 0.5,xinit < radius - 0.5,xinit > nx - 0.5 - radius],axis=0)
    if np.sum(bad) > 0:
        xnew[bad] = xinit[bad]
        xerr[bad] = 999.0

    # Return
    return xnew, xerr

def tc_indices(tc_dict):
    """ Quick parser of tc_dict

    Parameters
    ----------
    tc_dict : dict

    Returns
    -------
    left_idx : list
    left_xval : ndarray
    right_idx : list
    right_xval : ndarray
    """
    # Grab the existing edges (code is duplicated in mslit_sync)
    left_idx = [int(key) for key in tc_dict['left']['xval']]  # These match to the edge values in edgearr
    left_idx.sort(reverse=True)
    left_xval = np.array([tc_dict['left']['xval'][str(idx)] for idx in left_idx])

    right_idx = [int(key) for key in tc_dict['right']['xval'].keys()]  # These match to the edge values in edgearr
    right_idx.sort()
    right_xval = np.array([tc_dict['right']['xval'][str(idx)] for idx in right_idx])

    # Return
    return left_idx, left_xval, right_idx, right_xval


def slit_trace_qa(frame, ltrace, rtrace, extslit, setup, desc="",
                  normalize=True, use_slitid=None):
    """ Generate a QA plot for the slit traces

    Parameters
    ----------
    frame : ndarray
      trace image
    ltrace : ndarray
      Left slit edge traces
    rtrace : ndarray
      Right slit edge traces
    extslit : ndarray
      Mask of extrapolated slits (True = extrapolated)
    setup : str
    desc : str, optional
      A description to be used as a title for the page
    root : str, optional
      Root name for generating output file, e.g. msflat_01blue_000.fits
    normalize: bool, optional
      Normalize the flat?  If not, use zscale for output
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'
    ticks_font = font_manager.FontProperties(family='times new roman', style='normal', size=16,
                                             weight='normal', stretch='normal')

    # Outfile
    method = inspect.stack()[0][3]
    outfile = arqa.set_qa_filename(setup, method)
    ntrc = ltrace.shape[1]
    ycen = np.arange(frame.shape[0])
    # Normalize flux in the traces
    if normalize:
        nrm_frame = np.zeros_like(frame)
        for ii in range(ntrc):
            xtrc = (ltrace[:, ii] + rtrace[:, ii])/2.
            ixtrc = np.round(xtrc).astype(int)
            # Simple 'extraction'
            dumi = np.zeros((frame.shape[0], 3))
            for jj in range(3):
                dumi[:, jj] = frame[ycen, ixtrc-1+jj]
            trc = np.median(dumi, axis=1)
            # Find portion of the image and normalize
            for yy in ycen:
                xi = max(0, int(ltrace[yy, ii])-3)
                xe = min(frame.shape[1], int(rtrace[yy, ii])+3)
                # Fill + normalize
                nrm_frame[yy, xi:xe] = frame[yy, xi:xe] / trc[yy]
        sclmin, sclmax = 0.4, 1.1
    else:
        nrm_frame = frame.copy()
        nrm_frame[frame > 0.0] = np.sqrt(nrm_frame[frame > 0.0])
        sclmin, sclmax = arplot.zscale(nrm_frame)

    # Plot
    plt.clf()

    ax = plt.gca()
#    set_fonts(ax)
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    cmm = cm.Greys_r
    mplt = plt.imshow(nrm_frame, origin='lower', cmap=cmm, interpolation=None,
                      extent=(0., frame.shape[1], 0., frame.shape[0]))
    mplt.set_clim(vmin=sclmin, vmax=sclmax)

    # Axes
    plt.xlim(0., frame.shape[1])
    plt.ylim(0., frame.shape[0])
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off',
                    labelbottom='off', labelleft='off')

    # Traces
    iy = int(frame.shape[0]/2.)
    for ii in range(ntrc):
        if extslit[ii] is True:
            ptyp = ':'
        else:
            ptyp = '--'
        # Left
        plt.plot(ltrace[:, ii]+0.5, ycen, 'r'+ptyp, linewidth=0.3, alpha=0.7)
        # Right
        plt.plot(rtrace[:, ii]+0.5, ycen, 'c'+ptyp, linewidth=0.3, alpha=0.7)
        # Label
        if use_slitid is not None:
            slitid, _, _ = get_slitid(frame.shape, ltrace, rtrace, ii, ypos=0.5)
            lbl = 'S{:04d}'.format(slitid)
        else:
            lbl = '{0:d}'.format(ii+1)
        plt.text(0.5*(ltrace[iy, ii]+rtrace[iy, ii]), ycen[iy], lbl, color='green', ha='center', size='small')
    # Title
    tstamp = arqa.gen_timestamp()
    if desc == "":
        plt.suptitle(tstamp)
    else:
        plt.suptitle(desc+'\n'+tstamp)

    # Write
    plt.savefig(outfile, dpi=800)
    plt.close()

    plt.rcdefaults()



