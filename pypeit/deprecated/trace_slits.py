""" Module for core algorithms related to tracing slits/orders
These should primarily be called by the TraceSlits class
"""
import inspect
import copy
from collections import Counter

import numpy as np

from scipy import ndimage
from scipy.special import erf
from scipy import signal
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib import cm, font_manager


from pypeit import msgs
from pypeit.core import qa
from pypeit.core import plot
from pypeit import utils
from pypeit.core import pca
from pypeit.core import pixels
from pypeit.core import procimg
from pypeit import debugger
from pypeit.core import extract
from pypeit.core import arc
from pypeit.core import pydl
from astropy.stats import sigma_clipped_stats
from IPython import embed
import copy

try:
    from pypeit import ginga
except ImportError:
    pass

# Testing
import time


def shift_slits(tslits_dict, tilts_dict, waveimg, shift):
    """
    Routine to shift slits in a tslits_dict. This routine is used to compensate for flexure.

    Args:
        tslits_dict (dict):
           Dictionary containing slit boundaries and other info

        shift (float):
           Shift to be added to the slit boundaries.

    Returns:
        tslits_shift (dict):
           A copy of the original tslits_dict but with the slit boundaries and the slitcen shifted.

    """

    tslits_shift = copy.deepcopy(tslits_dict)
    nspec, nslits = tslits_shift['nspec'], tslits_shift['nslits']

    for islit in range(nslits):
        tslits_shift['slit_left'][:,islit] += shift
        tslits_shift['slit_righ'][:,islit] += shift
        tslits_shift['slitcen'][:,islit] += shift

    # Not shifting the orignal slit_left_orig or slit_righ_orig
    return tslits_shift


def extrapolate_trace(traces_in, spec_min_max_in, fit_frac=0.2, npoly=1, method='poly'):
    """
    Extrapolates trace to fill in pixels that lie outside of the range spec_min, spec_max). This
    routine is useful for echelle spectrographs where the orders are shorter than the image by a signfiicant
    amount, since the polynomial trace fits often go wild.

    Args:
        traces (np.ndarray): shape = (nspec,) or (nspec, ntrace)
            Array containing object or slit boundary traces
        spec_min_max (np.ndarray):  shape = (2, ntrace)
            Array contaning the minimum  and maximum spectral region covered by each trace. If this is an array with
            ndim=1 array, the same numbers will be used for all traces in traces_in.  If a 2d array, then this must be
            an ndarray of shape (2, ntrace,) where the spec_min_max[0,:] are the minimua and spec_min_max[1,:] are the maxima.
        fit_frac (float):
            fraction of the good pixels to be used to fit when extrapolating traces. The upper fit_frac
            pixels are used to extrapolate to larger spectral position, and vice versa for lower spectral
            positions.
        npoly (int):
            Order of polynomial fit used for extrapolation
        method (str):
            Method used for extrapolation. Options are 'poly' or 'edge'. If 'poly' the code does a polynomial fit. If
            'edge' it just attaches the last good pixel everywhere. 'edge' is not currently used.
    Returns:
        trace_extrap (np.ndarray):
            Array with same size as trace containing the linearly extrapolated values for the bad spectral pixels.
    """

    #embed()
    # This little bit of code allows the input traces to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)
    spec_min_max_tmp = np.array(spec_min_max_in)
    if traces_in.ndim == 2:
        traces = traces_in
        nslits = traces.shape[1]
        if np.array(spec_min_max_in).ndim == 1:
            spec_min_max = np.outer(spec_min_max_tmp, np.ones(nslits))
        elif spec_min_max_tmp.ndim == 2:
            if (spec_min_max_tmp.shape[1] != nslits):
                msgs.error('If input as any arrays, spec_min_max needs to have dimensions (2,nslits)')
            spec_min_max = spec_min_max_tmp
        else:
            msgs.error('Invalid shapes for traces_min and traces_max')
    else:
        nslits = 1
        traces = traces_in.reshape(traces_in.size, 1)
        spec_min_max = spec_min_max_tmp

    nspec = traces.shape[0]
    spec_vec = np.arange(nspec,dtype=float)
    xnspecmin1 = spec_vec[-1]
    traces_extrap = traces.copy()

    # TODO should we be doing a more careful extrapolation here rather than just linearly using the nearest pixel
    # values??
    for islit in range(nslits):
        ibad_max = spec_vec > spec_min_max[1,islit]
        ibad_min = spec_vec < spec_min_max[0,islit]
        igood = (spec_vec >= spec_min_max[0,islit]) & (spec_vec <= spec_min_max[1,islit])
        nfit = int(np.round(fit_frac*np.sum(igood)))
        good_ind = np.where(igood)[0]
        igood_min = good_ind[0:nfit]
        igood_max = good_ind[-nfit:]
        if np.any(ibad_min):
            if 'poly' in method:
                coeff_min = utils.func_fit(spec_vec[igood_min], traces[igood_min, islit], 'legendre', npoly, minx=0.0, maxx=xnspecmin1)
                traces_extrap[ibad_min, islit] = utils.func_val(coeff_min, spec_vec[ibad_min], 'legendre', minx=0.0, maxx=xnspecmin1)
            elif 'edge' in method:
                traces_extrap[ibad_min, islit] = traces[good_ind[0], islit]
        if np.any(ibad_max):
            if 'poly' in method:
                coeff_max = utils.func_fit(spec_vec[igood_max], traces[igood_max, islit], 'legendre', npoly, minx=0.0,maxx=xnspecmin1)
                traces_extrap[ibad_max, islit] = utils.func_val(coeff_max, spec_vec[ibad_max], 'legendre', minx=0.0,maxx=xnspecmin1)
            elif 'edge' in method:
                traces_extrap[ibad_max, islit] = traces[good_ind[-1], islit]

        #ibad = np.invert(igood)
        #traces_extrap[ibad, islit] = interpolate.interp1d(spec_vec[igood], traces[igood, islit], kind='linear',
        #bounds_error=False, fill_value='extrapolate')(spec_vec[ibad])

    return traces_extrap


def add_user_edges(lcen, rcen, add_slits):
    """
    Add user-defined slit(s)

    Warning: There is no real error checking here.
    The user is assumed to know what they are doing!

    tc_dict is updated in place

    Args:
        lcen (np.ndarray): Left traces of slit/orders
        rcen (np.ndarray): Right traces of slit/orders
        add_slits (list):  List of slit info for adding
           y_spec, x_spat0, x_spat1 (int)

    Returns:
        np.ndarray, np.ndarray: new lcen, rcen arrays

    """
    nspec = lcen.shape[0]
    ycen = nspec//2

    # Loop me
    for new_slit in add_slits:
        msgs.info("Adding a user-defined slit [yrow, x0, x1]:  {}".format(new_slit))
        # Parse
        y_spec, x_spat0, x_spat1 = new_slit

        for xx, side in zip([x_spat0,x_spat1], ['left', 'right']):
            # Left
            ref_t = lcen if side == 'left' else rcen
            ref_x = ref_t[y_spec,:]
            # Find the closest
            idx = np.argmin(np.abs(xx-ref_x))
            dx = ref_x[idx]-xx
            # New trace
            new_trace = ref_t[:,idx] - dx
            if side == 'left':
                lcen = np.append(ref_t, new_trace.reshape(nspec,1), axis=1)
            else:
                rcen = np.append(ref_t, new_trace.reshape(nspec,1), axis=1)

    # Sort me
    for side in ['left', 'right']:
        ref_t = lcen if side == 'left' else rcen
        allx = ref_t[ycen,:]
        isrt = np.argsort(allx)
        # Do it
        if side == 'left':
            lcen = ref_t[:,isrt]
        else:
            rcen = ref_t[:,isrt]
    # Done
    return lcen, rcen

def rm_user_edges(lcen, rcen, rm_slits):
    """
    Remove one or more slits, as applicable

    Code compares exisiting slits (which must be sycnhronized)
    against the input request and removes any that match.

    Args:
        lcen (np.ndarray): Left traces of slit/orders
        rcen (np.ndarray): Right traces of slit/orders
        rm_slits: list
          y_spec, x_spat pairs

    Returns:
        np.ndarray, np.ndarray: new lcen, rcen arrays

    """
    # Mask me
    good_left = np.ones(lcen.shape[1], dtype=bool)
    good_right = np.ones(rcen.shape[1], dtype=bool)
    # Check
    if good_left.size != good_right.size:
        msgs.error("Slits need to be sync'd to use this method!")
    # Loop me
    for rm_slit in rm_slits:
        # Deconstruct
        y_spec, xcen = rm_slit
        # Edges
        lefts = lcen[y_spec,:]
        rights = rcen[y_spec,:]
        # Match?
        bad_slit = (lefts < xcen) & (rights > xcen)
        if np.any(bad_slit):
            # Double check
            if np.sum(bad_slit) != 1:
                msgs.error("Something went horribly wrong in edge tracing")
            #
            idx = np.where(bad_slit)[0]
            msgs.info("Removing user-supplied slit at {},{}".format(xcen,y_spec))
            good_right[idx] = False
            good_left[idx] = False

    # Finish
    lcen = lcen[:,good_left]
    rcen = rcen[:,good_right]

    return lcen, rcen


def orig_add_user_edges(edgearr, siglev, tc_dict, add_slits):
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
                xset, xerr = trace_crude_init(np.maximum(siglev, -0.1), np.array([xleft]), yrow, maxshift0=0.5, maxshift=0.15, maxerr=0.2)
                #
                new_i = new_l
                ref_x = left_xval
                ref_i = left_idx
            else:
                xset, xerr = trace_crude_init(np.maximum(-1*siglev, -0.1), np.array([xright]), yrow,maxshift0=0.5, maxshift=0.15, maxerr=0.2)
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

def assign_slits(binarr, edgearr, ednum=100000, lor=-1, function='legendre', polyorder=3):
    """
    This routine will trace the locations of the slit edges.  Putative
    edges come in with |values| > 200000 (in the edgearr) and leave
    (ideally) with values near ednum.

    Args:
        binarr (numpy.ndarray):
            Calibration frame that will be used to identify slit traces
            (in most cases, the slit edge).  Typically previously
            processed by a uniform filter

        edgearr (numpy.ndarray):
            An array of negative/positive numbers (left/right edges
            respectively) and zeros (no edge).

        ednum (:obj:`int`, optional):
            A dummy number given to define slit edges.

        lor (:obj:`int`, optional):
            A flag that indicates if the left edge (-1) or right edge
            (+1) should be assigned.

        function (:obj:`str`, optional):
            The type of function used to trace the slit edges.  Used by
            :func:`utils.robust_polyfit` and :func:`utils.func_val`.

        polyorder (:obj:`int`, optional):
            The order of the function used for the fit.  Used by
            :func:`utils.robust_polyfit`.

    Returns:
        numpy.ndarray: An array of negative/positive numbers (left/right
        edges respectively) and zeros (no edge).  The input
        :arg:`edgearr` is actually modified in place and returned.
    """
#    if settings is None:
#        settings = dict(trace={'slits': {'polyorder': 3, 'function': 'legendre'}})

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
                    wpkmsk = prune_peaks(smedgehist, wpk, np.where(wpk+2 == offs)[0][0])
                except:
                    debugger.set_trace()
                wpk = wpk[np.where(wpkmsk == 1)]
            if wpk.size == 0:
                # After pruning, there are no more peaks
                break
            pks = wpk+2  # Shifted by 2 because of the peak finding algorithm above
            pedges = find_peak_limits(smedgehist, pks)

            if np.all(pedges[:, 1]-pedges[:, 0] == 0):
                # Remaining peaks have no width
                break

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
                if widx[0].size < 2*polyorder:
                    continue
                badmsk, fitcof = utils.robust_polyfit(widx[0], widx[1], polyorder,
                                                        function=function, minx=0,
                                                        maxx=binarr.shape[0]-1)
                shbad[widx] = badmsk
                smallhist = np.zeros(101, dtype=np.int)
                meddiff = np.zeros(vals.size)
                for vv in range(vals.size):
                    widx = np.where((edgearr == vals[vv]) & (shbad == 0))
                    if widx[0].size == 0:
                        # These pixels were deemed to be bad
                        continue
                    diff = widx[1] - utils.func_val(fitcof, widx[0], function, minx=0,
                                                      maxx=binarr.shape[0]-1)
                    diff = 50 + np.round(diff).astype(np.int)
                    np.add.at(smallhist, diff, 1)
                    meddiff[vv] = np.median(diff)
                # Find the peaks of this distribution
                wspk = np.where((smallhist[1:-1] >= smallhist[2:]) & (smallhist[1:-1] > smallhist[:-2]))[0]
                wspk += 1  # Add one here to account for peak finding
#                if False:
#                    plt.clf()
#                    plt.plot(smallhist, 'k-', drawstyle='steps')
#                    plt.show()

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
        msk, cf = utils.robust_polyfit(wedx, wedy, polyorder, function=function,
                                         minx=0, maxx=binarr.shape[0]-1)
        cenmodl = utils.func_val(cf, np.arange(binarr.shape[0]), function,
                                   minx=0, maxx=binarr.shape[0]-1)
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


def new_add_edge(ref_slit, insert_offset, t_dict, left=True):
    """  Add a new edge using a reference slit

    Args:
        ref_slit: int
        insert_offset: int
          Offset from the right slit for the new left slit
          or vice-versa
        t_dict: dict
        left: bool, optional

    Returns:
        Fills tdict in-place

    """
    # Current indices (book-keeping)
    if left:
        use = 'right'
        fill = 'left'
    else:
        use = 'left'
        fill = 'right'
    traces = t_dict[use]['traces']

    # TODO - Use the PCA
    # Use the reference edge for the shape
    new_trace = traces[:,ref_slit] + insert_offset
    ypos = new_trace.shape[0]//2

    # Add it in
    t_dict[fill]['new_xval'].append(new_trace[ypos])
    t_dict[fill]['new_traces'].append(new_trace)

    # Return
    return


def sync_edges(tc_dict, nspat, insert_buff=5, verbose=False):
    """ Method to synchronize the slit edges
    Adds in extra edges according to a few criteria

    Developed for ARMLSD

    Parameters
    ----------
    tc_dict : dict
       For book-keeping
    ednum : int
    nspat : int
    insert_buff : int, optional
       Offset from existing edge for any edge added in

    Returns
    -------
    """
    # TODO - Should avoid adding a slit at the edge if those columns are masked in the BPM
    # Init
    for key in ['left', 'right']:
        tc_dict[key]['new_xval'] = []
        tc_dict[key]['new_traces'] = []

    # Grab the edge indexes and xval's
    #left_idx, left_xval, right_idx, right_xval = tc_indices(tc_dict)
    left_xval = tc_dict['left']['xval']
    right_xval = tc_dict['right']['xval']

    # Only one slit?
    if (len(left_xval) == 1) and (len(right_xval)==1):
        if left_xval[0] < right_xval[0]:  # Ok slit, otherwise continue
            return

    # Masks: True is a good edge, False is bad
    good_left = np.ones_like(left_xval, dtype=bool)
    good_right = np.ones_like(right_xval, dtype=bool)


    # Deal with missing left edges first (at left edge of detector)
    rights_missing_lefts = np.where(right_xval < left_xval[0])[0]

    for kk in rights_missing_lefts:
        # Grab the trace
        right_pix = tc_dict['right']['traces'][:,kk] #np.where(edgearr == right_idx[0])
        mn_rp = np.min(right_pix)
        if mn_rp <= insert_buff:
            good_right[kk] = False
            msgs.warn("Partial or too small right edge at start of detector.  Skipping it.")
        else:
            if kk == 0:
                ioff = -1*mn_rp + insert_buff
            else:
                ioff = right_xval[kk-1] - right_xval[kk] + insert_buff
            msgs.warn("Adding in a left edge near start of detector which mirrors the first right edge")
            new_add_edge(kk, ioff, tc_dict, left=True)

    # Loop on left edges
    for kk,left in enumerate(left_xval):

        # Grab location of the next left edge
        if kk < len(left_xval)-1:
            next_left = left_xval[kk+1]
        else:
            next_left = nspat-1

        # Search for a proper right edge
        #  Should be right of the current left and left of the next left
        gd_right = np.where((right_xval < next_left) & (right_xval > left))[0]
        if len(gd_right) == 0:   # None found?
            # Last slit?
            if kk == len(left_xval)-1:
                msgs.warn("Last slit has no right edge.  Adding one in which will not touch the detector edge")
                left_pix = tc_dict['left']['traces'][:, kk] #np.where(edgearr == left_idx[kk])
                mx_lp = np.max(left_pix[1])
                if mx_lp >= nspat-1:
                    msgs.warn("Partial left edge at end of detector.  Skipping it.")
                    good_left[kk] = False
                else:
                    # Stay on the detector!
                    ioff = nspat - mx_lp - insert_buff
                    # Add
                    new_add_edge(-1, ioff, tc_dict, left=False)
                continue
            else: # Not the last slit, add one in!
                msgs.warn("Missing a right edge for slit with left edge at {}".format(left))
                msgs.warn("Adding in a corresponding right edge!")
                # Offset from the next left edge
                ioff = next_left-left-insert_buff
                # Add
                new_add_edge(kk, ioff, tc_dict, left=False)
        else:
            # Check for multiple right edges between the two lefts (i.e. missing Left)
            #     Will only add in one missing left
            if len(gd_right) > 1:
                msgs.warn("Missing a left edge for slit with right edge(s) at {}".format(
                    right_xval[gd_right[1:]]))
                msgs.warn("Adding one (and only one)")
                # Offset is difference between the two right slits + a buffer
                ioff = right_xval[gd_right[0]] - right_xval[gd_right[1]] + insert_buff
                # Add
                new_add_edge(gd_right[1], ioff, tc_dict, left=True)
                #add_edge(right_idx[gd_right[1]], ioff, edgearr, tc_dict, final_left, final_right, left=True)

    # Deal with good
    tc_dict['right']['xval'] = tc_dict['right']['xval'][good_right]
    tc_dict['right']['traces'] = tc_dict['right']['traces'][:,good_right]
    tc_dict['left']['xval'] = tc_dict['left']['xval'][good_left]
    tc_dict['left']['traces'] = tc_dict['left']['traces'][:,good_left]


    # Add em in and then sort
    for side in ['left', 'right']:
        for kk, xval in enumerate(tc_dict[side]['new_xval']):
            tc_dict[side]['xval'] = np.append(tc_dict[side]['xval'], xval)
            tmp = tc_dict[side]['new_traces'][kk]
            tc_dict[side]['traces'] = np.append(tc_dict[side]['traces'], np.resize(tmp, (tmp.size,1)), axis=1)

        # Sort
        isrt = np.argsort(tc_dict[side]['xval'])
        tc_dict[side]['xval'] = tc_dict[side]['xval'][isrt]
        tc_dict[side]['traces'] = tc_dict[side]['traces'][:,isrt]

    # Return
    return

'''
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
'''


def edgearr_tcrude(edgearr, siglev, ednum, TOL=3., tfrac=0.33, verbose=False,
                   maxshift=0.15, bpm=None, skip_bad=True):
    """ Use trace_crude to refine slit edges
    It is also used to remove bad slit edges and merge slit edges

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
    maxshift : float
      Maximum shift in trace crude

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
    msgs.info("Crude tracing the edges")
    # Init
    nspec = edgearr.shape[0]
    ycen = nspec//2

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
        tc_dict[side]['xset'] = np.zeros((nspec,len(uni_e)))
        tc_dict[side]['xerr'] = np.zeros((nspec,len(uni_e))) + 999.

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
                xset, xerr = trace_crude_init(np.maximum(siglev, -0.1), np.array(xinit), yrow, maxshift=maxshift,maxshift0=0.5, maxerr=0.2)
            else:
                xset, xerr = trace_crude_init(np.maximum(-1*siglev, -0.1), np.array(xinit), yrow, maxshift=maxshift, maxshift0=0.5, maxerr=0.2)
            # Fill it up
            for kk,x in enumerate(xinit):
                # Annoying index
                idx = np.where(uni_e == edgearr[yrow,x])[0]
                #
                tc_dict[side]['xset'][:,idx[0]] = xset[:,kk]
                tc_dict[side]['xerr'][:,idx[0]] = xerr[:,kk]

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
                # All bad trace_crude?
                if not np.any(goodx[:,kk]):
                    msgs.warn("No good trace values. Rejecting")
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
                    if skip_bad:
                        ybad_xerr = np.array([])
                    else:
                        ybad_xerr = np.where(~goodx[:,kk])[0]
                    # Ignore bad pixels in BPM -- Somewhat kludgy
                    if bpm is not None:
                        keep_bad = []
                        for ibad in ybad_xerr:
                            xval = int(np.round(xset[ibad,kk]))
                            if bpm[ibad, xval] == 0:
                                keep_bad.append(ibad)
                        ybad_xerr = np.array(keep_bad)
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
            # Remove bad traces
            if len(gde) == 0:
                msgs.warn("Side {} had no good edges;  Keeping one!".format(side))
                # Keep 1 (this is mainly for Longslit)
                xval = tc_dict[side]['xset'][ycen,:]
                if side == 'left':
                    idx = np.argmin(xval)
                    tc_dict[side]['xval'][str(-1*ednum)] = xval[idx]
                else:
                    idx = np.argmax(xval)
                    tc_dict[side]['xval'][str(ednum)] = xval[idx]
                idx = np.array([idx]) # Needs to be an array for the shape
                #
                tc_dict[side]['xset'] = tc_dict[side]['xset'][:,idx]
                tc_dict[side]['xerr'] = tc_dict[side]['xerr'][:,idx]
            else:
                tc_dict[side]['xset'] = tc_dict[side]['xset'][:,gde]
                tc_dict[side]['xerr'] = tc_dict[side]['xerr'][:,gde]

    # Remove uni_idx
    for side in ['left', 'right']:
        for key in ['uni_idx']:
            tc_dict[side].pop(key)
    if verbose:
        print(tc_dict['left']['xval'])
        print(tc_dict['right']['xval'])
    # Return
    return new_edgarr, tc_dict.copy()


def edgearr_from_binarr(binarr, binbpx, medrep=0, min_sqm=30.,
                        sobel_mode='nearest', sigdetect=30.):
    """ Generate the edge array from an input, trace image (likely slightly filtered)
    Primary algorithm is to run a Sobolev filter on the image and then
    trigger on all significant features.

    The bad pixel mask is also used to fuss with bad columns, etc.

    Parameters
    ----------
    binarr : numpy ndarray
      Calibration frame that will be used to identify slit traces (in most cases, the slit edge)
      Lightly filtered
    binbpx : ndarray
      Bad pixel max image
    medrep : int, optional
        Number of times to perform median smoothing on the mstrace
        One uniform filter is always done
        medrep = 0 is recommended for ARMLSD
    sobel_mode : str, optional
        ndimage.sobel mode;  default is 'nearest'
    sigdetect : float, optional
        threshold for edge detection
    min_sqm : float, optional
        Minimum error used when detecting a slit edge

    Returns
    -------
    siglev : ndarray
    edgearr : ndarray
    """
    # Specify how many times to repeat the median filter
    # Even better would be to fit the filt/sqrt(abs(binarr)) array with a Gaussian near the maximum in each column
    msgs.info("Detecting slit edges in the mstrace image")

    # Replace bad columns
    #  TODO -- Should consider replacing bad 'rows' for rotated detectors (e.g. GMOS)
    bad_cols = np.sum(binbpx, axis=0) > (binbpx.shape[0]//2)
    if np.any(bad_cols):
        ms2 = procimg.replace_columns(binarr, bad_cols)
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

    # First edges assigned according to S/N
    tedges = np.zeros(binarr.shape, dtype=np.float)
    wl = np.where(siglev > + sigdetect)  # A positive gradient is a left edge
    wr = np.where(siglev < - sigdetect)  # A negative gradient is a right edge
    tedges[wl] = -1.0
    tedges[wr] = +1.0

    # Clean the edges
    wcl = np.where((ndimage.maximum_filter1d(siglev, 10, axis=1) == siglev) & (tedges == -1))
    wcr = np.where((ndimage.minimum_filter1d(siglev, 10, axis=1) == siglev) & (tedges == +1))
    nedgear = np.zeros(siglev.shape, dtype=np.int)
    nedgear[wcl] = -1
    nedgear[wcr] = +1

    ######
    msgs.info("Applying bad pixel mask")
    nedgear *= (1-binbpx.astype(np.int))  # Apply the bad pixel mask
    siglev *= (1-binbpx.astype(np.int))  # Apply the bad pixel mask

    # Return
    return siglev, nedgear


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
    If 0 is returned for both counts, this detector will be skipped
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
            msgs.warn("Unable to trace any edges"+msgs.newline()+"try a different method to trace the order edges")
            return None, 0, 0
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
    nspec, nspat = edgearr.shape
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
    if (lcnt == 1) & (rcnt > 1):  # This is mainly in here for LRISb which is a real pain..
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



def fit_edges(edgearr, lmin, lmax, plxbin, plybin, left=True, polyorder=3, function='ledgendre'):
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
    msk, cf = utils.robust_polyfit(wedx, wedy,
                                     polyorder,
                                     function=function,
                                     minx=0, maxx=edgearr.shape[0] - 1)
    cenmodl = utils.func_val(cf, np.arange(edgearr.shape[0]), function,
                               minx=0, maxx=edgearr.shape[0] - 1)

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
        msk, coeff[:, i - lmin] = utils.robust_polyfit(tlfitx, tlfity, polyorder,
                                                          function=function,
                                                          minx=minvf, maxx=maxvf)
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


# TODO: Need a better name for this function
def base_expand_slits(msedge, ordcen, extord):

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


def find_between(edgdet, ledgem, ledgep, dirc):
    """
    .. todo::
        Document this!
    """
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


def find_shift(mstrace, minarr, lopos, diffarr, numsrch):
    """
    .. todo::
        Document this!
    """
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


def ignore_orders(edgdet, fracpix, lmin, lmax, rmin, rmax):
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


def limit_yval(yc, maxv):
    yn = 0 if yc == 0 else (-yc if yc < 3 else -3)
    yx = maxv-yc if yc > maxv-4 and yc < maxv else 4
    return yn, yx


def match_edges(edgdet, ednum, mr=50):
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
                yn, yx = limit_yval(yt, sz_y)

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
                yn, yx = limit_yval(yt, sz_y)

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



# TODO: This pure python function was translated from the cython
# function above but was never fully tested; compare with function in
# pypeit/arcytrace.pyx!
def close_edges(edgdet, dets, npix):
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


# TODO: This pure python function was translated from the cython
# function above but was never fully tested; compare with function in
# pypeit/arcytrace.pyx!  The code is knarly and may need some attention!
def close_slits(trframe, edgdet, dets, npix, ednum):

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


# TODO: This pure python function was translated from the cython
# function above but was never fully tested; compare with function in
# pypeit/arcytrace.pyx!  This edits everything in place.  Do we want to
# do that?
def dual_edge(edgearr, edgearrcp, wx, wy, wl, wr, shft, npix, newval):

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


def minbetween(mstrace, loord, hiord):
    """
    .. todo::
        Document this!
    """
    # TODO: Check shapes
    ymin = np.clip(loord, 0, mstrace.shape[1])
    ymax = np.clip(hiord, 0, mstrace.shape[1])
    minarr = np.zeros(mstrace.shape[0])
    indx = ymax > ymin
    minarr[indx] = np.array([ np.amin(t[l:h])
                                for t,l,h in zip(mstrace[indx], ymin[indx], ymax[indx]) ])
    return minarr


def find_peak_limits(hist, pks):
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


# TODO -- REMOVE AS THIS IS NOW IN trace.py
def parse_user_slits(add_slits, this_det, rm=False):
    """
    Parse the parset syntax for adding slits

    Args:
        add_slits: str, list
          Taken from the parset
        this_det: int
          current detector
        rm: bool, optional
          Remove instead of add?

    Returns:
        user_slits: list or None
          if list,  [[x0,x1,yrow]] for add with one or more entries
          if list,  [[xcen,yrow]] for rm with one or more entries

    """
    # Might not be a list yet (only a str)
    if not isinstance(add_slits, list):
        add_slits = [add_slits]
    #
    user_slits = []
    for islit in add_slits:
        if not rm:
            det, x0, x1, yrow = [int(ii) for ii in islit.split(':')]
            if det == this_det:
                user_slits.append([x0,x1,yrow])
        else:
            det, xcen, yrow = [int(ii) for ii in islit.split(':')]
            if det == this_det:
                user_slits.append([xcen,yrow])
    # Finish
    if len(user_slits) == 0:
        return None
    else:
        return user_slits


#def pca_order_slit_edges(binarr, edgearr, lcent, rcent, gord, lcoeff, rcoeff, plxbin, slitcen,
#                         pixlocn, function='lengendre', polyorder=3, diffpolyorder=2,
#                         ofit=[3,2,1,0,0,0], extrapolate=[0,0], doqa=True):
#    """ Perform a PCA analyis on the order edges
#    Primarily for extrapolation
#
#    KBW: extrapolate neg was never used
#
#    Parameters
#    ----------
#    binarr : ndarray
#    edgearr : ndarray
#    lcent : ndarray
#      Left edges
#    rcent : ndarray
#      Right edges
#    gord : ndarray
#      Orders detected on both the left and right edge
#    lcoeff : ndarray
#      Fit coefficients for left edges
#    rcoeff : ndarray
#      Fit coefficients for right edges
#    plxbin
#    slitcen : ndarray
#    pixlocn
#
#    Returns
#    -------
#
#    """
#    # Init
#    wl = np.where(edgearr < 0)
#    wr = np.where(edgearr > 0)
#
#    ##############
#    xv = plxbin[:, 0]
#    minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]
#
#    # min and max switched because left edges have negative values
#    almin, almax = -np.max(edgearr[wl]), -np.min(edgearr[wl])
#    armin, armax = np.min(edgearr[wr]), np.max(edgearr[wr])
#
#    # maskord = np.where((np.all(lcoeff[:,lg],axis=0)==False) 
#    #                         | (np.all(rcoeff[:,rg],axis=0)==False))[0]
#    maskord = np.where((np.all(lcoeff, axis=0) == False) | (np.all(rcoeff, axis=0) == False))[0]
#    ordsnd = np.arange(min(almin, armin), max(almax, armax) + 1)
#    totord = ordsnd[-1] + extrapolate[1]
#
#    # Identify the orders to be extrapolated during reconstruction
#    extrapord = (1.0 - np.in1d(np.linspace(1.0, totord, totord), gord).astype(np.int)).astype(np.bool)
#    msgs.info("Performing a PCA on the order edges")
#    lnpc = len(ofit) - 1
#    msgs.work("May need to do a check here to make sure ofit is reasonable")
#    coeffs = utils.func_fit(xv, slitcen, function, polyorder, minx=minvf, maxx=maxvf)
#    for i in range(ordsnd.size):
#        if i in maskord:
#            if (i>=ordsnd[0]) and (i<ordsnd[-1]-1):  # JXP: Don't add orders that are already in there
#                continue
#            coeffs = np.insert(coeffs, i, 0.0, axis=1)
#            slitcen = np.insert(slitcen, i, 0.0, axis=1)
#            lcent = np.insert(lcent, i, 0.0, axis=0)
#            rcent = np.insert(rcent, i, 0.0, axis=0)
#    xcen = xv[:, np.newaxis].repeat(ordsnd.size, axis=1)
#    fitted, outpar = pca.basis(xcen, slitcen, coeffs, lnpc, ofit, x0in=ordsnd, mask=maskord,
#                                 skipx0=False, function=function)
##    if doqa:
##        debugger.set_trace()  # NEED TO REMOVE slf
##        # pca.pca_plot(slf, outpar, ofit, "Slit_Trace", pcadesc=pcadesc)
#
#    # Extrapolate the remaining orders requested
#    orders = 1 + np.arange(totord)
#    extrap_cent, outpar = pca.extrapolate(outpar, orders, function=function)
#
#    # Fit a function for the difference between left and right edges.
#    diff_coeff, diff_fit = utils.polyfitter2d(rcent - lcent, mask=maskord, order=diffpolyorder)
#
#    # Now extrapolate the order difference
#    ydet = np.linspace(0.0, 1.0, lcent.shape[0])
#    ydetd = ydet[1] - ydet[0]
#    lnum = ordsnd[0] - 1.0
#    ydet = np.append(-ydetd * np.arange(1.0, 1.0 + lnum)[::-1], ydet)
#    ydet = np.append(ydet, 1.0 + ydetd * np.arange(1.0, 1.0 + extrapolate[1]))
#    xde, yde = np.meshgrid(np.linspace(0.0, 1.0, lcent.shape[1]), ydet)
#    extrap_diff = utils.polyval2d(xde, yde, diff_coeff).T
#    msgs.info("Refining the trace for reconstructed and predicted orders")
#
#    # NOTE::  MIGHT NEED TO APPLY THE BAD PIXEL MASK HERE TO BINARR
#    msgs.work("Should the bad pixel mask be applied to the frame here?")
#    refine_cent, outpar = refine_traces(binarr, outpar, extrap_cent, extrap_diff,
#                                        [gord[0] - orders[0], orders[-1] - gord[-1]], orders,
#                                        ofit[0], pixlocn, function=function)
#
#    # Generate the left and right edges
#    lcen = refine_cent - 0.5 * extrap_diff
#    rcen = refine_cent + 0.5 * extrap_diff
#
#    # Return
#    return lcen, rcen, extrapord
#
#
#def pca_pixel_slit_edges(binarr, edgearr, lcoeff, rcoeff, ldiffarr, rdiffarr,
#                         lnmbrarr, rnmbrarr, lwghtarr, rwghtarr, lcent, rcent, plxbin,
#                         function='lengendre', polyorder=3, ofit=[3,2,1,0,0,0], doqa=True):
#    """ PCA analysis for slit edges
#
#    Parameters
#    ----------
#    binarr : ndarray
#    edgearr : ndarray
#    lcoeff : ndarray
#      Fit coefficients for left edges
#    rcoeff : ndarray
#      Fit coefficients for right edges
#    ldiffarr : ndarray
#    rdiffarr : ndarray
#    lnmbrarr : ndarray
#    rnmbrarr : ndarray
#    lwghtarr : ndarray
#    rwghtarr : ndarray
#    lcent : ndarray
#    rcent : ndarray
#    plxbin
#    settings : dict-like object
#
#    Returns
#    -------
#    lcen : ndarray
#      Left edges
#    rcen : ndarray
#      Right edges
#    extrapord
#
#    """
#
#    minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]
#    maskord = np.where((np.all(lcoeff, axis=0) == False) | (np.all(rcoeff, axis=0) == False))[0]
#    allord = np.arange(ldiffarr.shape[0])
#    ww = np.where(np.in1d(allord, maskord) == False)[0]
#
#    # Unmask where an order edge is located
#    maskrows = np.ones(binarr.shape[1], dtype=np.int)
#    # ldiffarr = np.round(ldiffarr[ww]).astype(np.int)
#    # rdiffarr = np.round(rdiffarr[ww]).astype(np.int)
#    ldiffarr = np.fmax(np.fmin(np.round(ldiffarr[ww]).astype(np.int), binarr.shape[1] - 1), 0)
#    rdiffarr = np.fmax(np.fmin(np.round(rdiffarr[ww]).astype(np.int), binarr.shape[1] - 1), 0)
#    maskrows[ldiffarr] = 0
#    maskrows[rdiffarr] = 0
#
#    # Extract the slit edge ID numbers associated with the acceptable traces
#    lnmbrarr = lnmbrarr[ww]
#    rnmbrarr = rnmbrarr[ww]
#
#    # Fill in left/right coefficients
#    tcoeff = np.ones((polyorder + 1, binarr.shape[1]))
#    tcoeff[:, ldiffarr] = lcoeff[:, ww]
#    tcoeff[:, rdiffarr] = rcoeff[:, ww]
#
#    # Weight the PCA fit by the number of detections in each slit edge
#    pxwght = np.zeros(binarr.shape[1])
#    pxwght[ldiffarr] = lwghtarr[ww]
#    pxwght[rdiffarr] = rwghtarr[ww]
#    maskrw = np.where(maskrows == 1)[0]
#    maskrw.sort()
#    extrap_row = maskrows.copy()
#    xv = np.arange(binarr.shape[0])
#
#    # trace values
#    trcval = utils.func_val(tcoeff, xv, function, minx=minvf, maxx=maxvf).T
#    msgs.work("May need to do a check here to make sure ofit is reasonable")
#    lnpc = len(ofit) - 1
#
#    # Only do a PCA if there are enough good slits
#    if np.sum(1.0 - extrap_row) > ofit[0] + 1:
#        # Perform a PCA on the locations of the slits
#        msgs.info("Performing a PCA on the slit traces")
#        ordsnd = np.arange(binarr.shape[1])
#        xcen = xv[:, np.newaxis].repeat(binarr.shape[1], axis=1)
#        fitted, outpar = pca.basis(xcen, trcval, tcoeff, lnpc, ofit, weights=pxwght, x0in=ordsnd, mask=maskrw, skipx0=False, function=function)
##        if doqa:
##            # JXP -- NEED TO REMOVE SLF FROM THE NEXT BIT
##            msgs.warn("NEED TO REMOVE SLF FROM THE NEXT BIT")
##            # pca.pca_plot(slf, outpar, ofit, "Slit_Trace", pcadesc=pcadesc, addOne=False)
#        # Now extrapolate to the whole detector
#        pixpos = np.arange(binarr.shape[1])
#        extrap_trc, outpar = pca.extrapolate(outpar, pixpos, function=function)
#        # Extract the resulting edge traces
#        lcen = extrap_trc[:, ldiffarr]
#        rcen = extrap_trc[:, rdiffarr]
#        # Perform a final shift fit to ensure the traces closely follow the edge detections
#        for ii in range(lnmbrarr.size):
#            wedx, wedy = np.where(edgearr == lnmbrarr[ii])
#            shft = np.mean(lcen[wedx, ii] - wedy)
#            lcen[:, ii] -= shft
#        for ii in range(rnmbrarr.size):
#            wedx, wedy = np.where(edgearr == rnmbrarr[ii])
#            shft = np.mean(rcen[wedx, ii] - wedy)
#            rcen[:, ii] -= shft
#    else:
#        allord = np.arange(lcent.shape[0])
#        maskord = np.where((np.all(lcent, axis=1) == False) | (np.all(rcent, axis=1) == False))[0]
#        ww = np.where(np.in1d(allord, maskord) == False)[0]
#        lcen = lcent[ww, :].T.copy()
#        rcen = rcent[ww, :].T.copy()
#    extrapord = np.zeros(lcen.shape[1], dtype=np.bool)
#
#    # Return
#    return lcen, rcen, extrapord




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


#def refine_traces(binarr, outpar, extrap_cent, extrap_diff, extord, orders,
#                  fitord, locations, function='polynomial'):
#    """
#    Parameters
#    ----------
#    binarr
#    outpar
#    extrap_cent
#    extrap_diff
#    extord
#    orders
#    fitord
#    locations
#    function : str, optional
#
#    Returns
#    -------
#
#    """
#    # Refine the orders in the positive direction
#    i = extord[1]
#    hiord = pixels.phys_to_pix(extrap_cent[:, -i-2], locations, 1)
#    nxord = pixels.phys_to_pix(extrap_cent[:, -i-1], locations, 1)
#    mask = np.ones(orders.size)
#    mask[0:extord[0]] = 0.0
#    mask[-extord[1]:] = 0.0
#    extfit = extrap_cent.copy()
#    outparcopy = copy.deepcopy(outpar)
#    while i > 0:
#        loord = hiord
#        hiord = nxord
#        nxord = pixels.phys_to_pix(extrap_cent[:,-i], locations, 1)
#
#        # Minimum counts between loord and hiord
#        minarrL = minbetween(binarr, loord, hiord)
#        minarrR = minbetween(binarr, hiord, nxord)
#
#        minarr = 0.5*(minarrL+minarrR)
#        srchz = np.abs(extfit[:,-i]-extfit[:,-i-1])/3.0
#        lopos = pixels.phys_to_pix(extfit[:,-i]-srchz, locations, 1)  # The pixel indices for the bottom of the search window
#        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,-i])))
#        diffarr = np.round(extrap_diff[:,-i]).astype(np.int)
#        shift = find_shift(binarr, minarr, lopos, diffarr, numsrch)
#
#        relshift = np.mean(shift+extrap_diff[:,-i]/2-srchz)
#        if shift == -1:
#            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[-i])))
#            relshift = 0.0
#        else:
#            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[-i]), relshift))
#        # Renew guess for the next order
#        mask[-i] = 1.0
#        extfit, outpar, fail = pca.refine_iter(outpar, orders, mask, -i, relshift, fitord, function=function)
#        if fail:
#            msgs.warn("Order refinement has large residuals -- check order traces")
#            return extrap_cent, outparcopy
#        i -= 1
#    # Refine the orders in the negative direction
#    i = extord[0]
#    loord = pixels.phys_to_pix(extrap_cent[:,i+1], locations, 1)
#    extrap_cent = extfit.copy()
#    outparcopy = copy.deepcopy(outpar)
#    while i > 0:
#        hiord = loord
#        loord = pixels.phys_to_pix(extfit[:,i], locations, 1)
#
#        minarr = minbetween(binarr,loord, hiord)
#        srchz = np.abs(extfit[:,i]-extfit[:,i-1])/3.0
#        lopos = pixels.phys_to_pix(extfit[:,i-1]-srchz, locations, 1)
#        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,i-1])))
#        diffarr = np.round(extrap_diff[:,i-1]).astype(np.int)
#        shift = find_shift(binarr, minarr, lopos, diffarr, numsrch)
#
#        relshift = np.mean(shift+extrap_diff[:,i-1]/2-srchz)
#        if shift == -1:
#            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[i-1])))
#            relshift = 0.0
#        else:
#            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[i-1]),relshift))
#        # Renew guess for the next order
#        mask[i-1] = 1.0
#        extfit, outpar, fail = pca.refine_iter(outpar, orders, mask, i-1, relshift, fitord, function=function)
#        if fail:
#            msgs.warn("Order refinement has large residuals -- check order traces")
#            return extrap_cent, outparcopy
#        i -= 1
#    return extfit, outpar


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


def synchronize_edges(binarr, edgearr, plxbin, lmin, lmax, lcoeff, rmin, rcoeff, lnmbrarr,
                      ldiffarr, lwghtarr, rnmbrarr, rdiffarr, rwghtarr, function='legendre',
                      polyorder=3, extrapolate=[0,0]):
    """ Synchrnoizes the existing edges

    For ARMLSD, this step is largely unnecessary given multi_sync()

    @Ryan :: Could use help making this method less onerous..

    
    KBW: extrapolate pos is never used


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
    lv = (utils.func_val(lcoeff[:, lval - lmin], xv, function, minx=minvf, maxx=maxvf) \
                + 0.5).astype(np.int)
    if np.any(lv < 0) or np.any(lv + 1 >= binarr.shape[1]):
        msgs.warn("At least one slit is poorly traced")
        msgs.info("Refer to the manual, and adjust the input trace parameters")
        debugger.set_trace()
        msgs.error("Cannot continue without a successful trace")
    mnvalp = np.median(binarr[:, lv + 1])  # Go one row above and one row below an order edge,
    mnvalm = np.median(binarr[:, lv - 1])  # then see which mean value is greater.

    """
    lvp = (utils.func_val(lcoeff[:,lval+1-lmin],xv,function,min=minvf,max=maxvf)+0.5).astype(np.int)
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
    if mnvalp > mnvalm:
        lvp = (utils.func_val(lcoeff[:, lval + 1 - lmin], xv, function, minx=minvf, maxx=maxvf) \
               + 0.5).astype(np.int)

        edgbtwn = find_between(edgearr, lv, lvp, 1)
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
        lvp = (utils.func_val(lcoeff[:, lval - 1 - lmin], xv, function,
                                minx=minvf, maxx=maxvf) + 0.5).astype(np.int)
        edgbtwn = find_between(edgearr, lvp, lv, -1)

        if edgbtwn[0] == -1 and edgbtwn[1] == -1:
            rsub = edgbtwn[2] - (lval - 1)  # There's an order overlap
        elif edgbtwn[1] == -1:  # No overlap
            rsub = edgbtwn[0] - (lval - 1)
        else:  # There's an order overlap
            rsub = edgbtwn[1] - (lval - 1)

    msgs.info("Relabelling slit edges")
    rsub = int(round(rsub))
    if lmin < rmin - rsub:
        esub = lmin - (extrapolate[0] + 1)
    else:
        esub = (rmin - rsub) - (extrapolate[0] + 1)

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
    nmord = polyorder + 1
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
    lcent = utils.func_val(lcoeff[:, -lunq[lg][::-1] - 1 - extrapolate[0]], xv, function,
                             minx=minvf, maxx=maxvf)
    rcent = utils.func_val(rcoeff[:, runq[rg] - 1 - extrapolate[0]], xv, function, minx=minvf,
                             maxx=maxvf)

    # Return
    return lcent, rcent, gord, lcoeff, ldiffarr, lnmbrarr, lwghtarr, \
                rcoeff, rdiffarr, rnmbrarr, rwghtarr

# TODO Make this a proper trace_crude, rename consistently with IDL
def trace_crude_init(image, xinit0, ypass, invvar=None, nave=5, radius=3.0,maxshift0=0.5, maxshift=0.1, maxerr=0.2):
    """Python port of trace_crude_idl.pro from IDLUTILS

    Modified for initial guess

    Parameters
    ----------
    image : 2D ndarray, shape (nspec, nspat)
      Image for tracing
    xinit : ndarray
      Initial guesses for trace peak at ypass
    ypass : int
      Row for initial guesses

    Optional Parameters
    -------------------
    radius: float, default = 3.0
        Radius for centroiding; default to 3.0
    nmed = int, default = None [NOT YET IMPLEMENTED!]
        Median filtering size down the nspec direction before performing trace
    nave = int, default = 5
        Boxcar averaging size down the nspec direction before performing trace. If set to None no averaging
        will be performed.
    maxerr: float, default = 0.2
        Maximum error in centroid allowed for valid recentering;
    maxshift: float, default = 0.1
        Maximum shift in centroid allowed for valid recentering.
    maxshift0: float, default 0.5
        Maximum shift in centroid allowed for initial row.

    Returns
    -------
    xset : Trace for each fiber
    xerr : Estimated error in that trace
    """
    # JFH TODO add error checking on input parameters

    # Init
    xinit = xinit0.astype(float)
    #xinit = xinit[0:3]
    ntrace = xinit.size
    ny = image.shape[0]
    xset = np.zeros((ny,ntrace))
    xerr = np.zeros((ny,ntrace))
    # Make copies of the image and the inverse variance image
    imgtemp = image.copy()
    if invvar is None:
        invtemp = np.zeros_like(image) + 1.
    else:
        invtemp = invvar.copy()

    # ToDo implement median filtering!

    # Boxcar-sum the entire image along columns by NAVE rows
    if nave is not None:
        nave = np.fmin(nave,ny)
        # Boxcar sum the entire image weighted by inverse variance over nave spectral pixels
        kernel = np.ones((nave, 1))/float(nave)
        imgconv = ndimage.convolve(imgtemp*invtemp, kernel, mode='nearest')
        # Add the weights
        invtemp = ndimage.convolve(invtemp, kernel, mode='nearest')
        # Look for pixels with infinite errors - replace with original values
        ibad = invtemp == 0.0
        invtemp[ibad] = 1.0
        imgconv[ibad] = imgtemp[ibad]
        # Renormalize the summed image by the weights
        imgtemp = imgconv/invtemp

    # JFH It seems odd to me that one is passing invtemp to trace_fweight, i.e. this is not correct
    # error propagation. While the image should be smoothed with inverse variance weights, the new noise
    # of the smoothed image has changed, and proper error propagation would then give:
    # var_convol = ndimage.convolve(1/invvar, kernel**2, mode='nearest')
    # invvar_convol = 1.0/var_convol
    # I have not implemented this for fear of breaking the behavior, and furthermore I think the desire was not
    # to have trace_fweight operate on formally correct errors.


    #  Recenter INITIAL Row for all traces simultaneously
    #
    iy = ypass * np.ones(ntrace,dtype=int)

    xfit,xfiterr = trace_fweight(imgtemp, xinit, ycen = iy, invvar=invtemp, radius=radius)
    # Shift
    xshift = np.clip(xfit-xinit, -1*maxshift0, maxshift0) * (xfiterr < maxerr)
    xset[ypass,:] = xinit + xshift
    xerr[ypass,:] = xfiterr * (xfiterr < maxerr)  + 999.0 * (xfiterr >= maxerr)

    #    /* LOOP FROM INITIAL (COL,ROW) NUMBER TO LARGER ROW NUMBERS */
    for iy in range(ypass+1, ny):
        xinit = xset[iy-1, :]
        ycen = iy * np.ones(ntrace,dtype=int)
        xfit,xfiterr = trace_fweight(imgtemp, xinit, ycen = ycen, invvar=invtemp, radius=radius)
        # Shift
        xshift = np.clip(xfit-xinit, -1*maxshift, maxshift) * (xfiterr < maxerr)
        # Save
        xset[iy,:] = xinit + xshift
        xerr[iy,:] = xfiterr * (xfiterr < maxerr)  + 999.0 * (xfiterr >= maxerr)
    #      /* LOOP FROM INITIAL (COL,ROW) NUMBER TO SMALLER ROW NUMBERS */
    for iy in range(ypass-1, -1,-1):
        xinit = xset[iy+1, :]
        ycen = iy * np.ones(ntrace,dtype=int)
        xfit,xfiterr = trace_fweight(imgtemp, xinit, ycen = ycen, invvar=invtemp, radius=radius)
        # Shift
        xshift = np.clip(xfit-xinit, -1*maxshift, maxshift) * (xfiterr < maxerr)
        # Save
        xset[iy,:] = xinit + xshift
        xerr[iy,:] = xfiterr * (xfiterr < maxerr) + 999.0 * (xfiterr >= maxerr)

    return xset, xerr


def trace_fweight(fimage, xinit_in, radius = 3.0, ycen=None, invvar=None):

    ''' Routine to recenter a trace using flux-weighted centroiding.

    Python port of trace_fweight.pro from IDLUTILS


    Parameters
    ----------
    fimage: 2D ndarray
      Image for tracing which shape (nspec, nspat)

    xinit: ndarray
      Initial guesses for spatial direction trace. This can either be an 2-d  array with shape
         (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.

    Optional Parameters:
    --------------------
    ycen: ndarray, default = None
      Optionally y-position of trace can be provided. It should be an integer array the same size as x-trace (nspec, nTrace). If
      not provided np.arange(nspec) will be assumed for each trace

    invvar: ndarray, default = None
         Inverse variance array for the image. Array with shape (nspec, nspat) matching fimage

    radius :  float or ndarray, default = 3.0
         Radius for centroiding in floating point pixels. This can be either be input as a scalar or as an array to perform
         centroiding with a varaible radius. If an array is input it must have the same size and shape as xinit_in, i.e.
         a 2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.

    Returns
    -------
    xnew:   ndarray
         Recentroided trace. The output will have the same shape as xinit i.e.  an 2-d  array with shape (nspec, nTrace)
         array if multiple traces were input, or a 1-d array with shape (nspec) for
         the case of a single trace.
    xerr: ndarray
         Formal propagated error on recentroided trace. These errors will only make sense if invvar is passed in, since
         otherwise invvar is set to 1.0.   The output will have the same shape as xinit i.e.  an 2-d  array with shape
         (nspec, nTrace) array if multiple traces were input, or a 1-d array with shape (nspec) for the case of a single
         trace. Locations will have this error set to 999 if any of the following conditions are true: 1) the flux weighted
         centroid deviates from the input guess by  > radius, or 2)  the centering window falls off the image, or 3) where any masked
         pixels (invvar == 0.0) contribute to the centroiding. The xnew values for the pixels which have xerr = 999 are reset to
         that of the input trace. These should thus be masked in any fit using the condition (xerr < 999)

         TODO we should probably either output a mask or set this 999 to something else, since I could imagine this causing problems.

     Revision History
     ----------------
     Python port of trace_fweight.pro from IDLUTILS
     24-Mar-1999  Written by David Schlegel, Princeton.
     27-Jun-2018  Ported to python by X. Prochaska and J. Hennawi
    """


    '''

    # Init
    nx = fimage.shape[1]

    # Checks on radius
    if isinstance(radius,(int, float)):
        radius_out = radius
    elif ((np.size(radius)==np.size(xinit_in)) & (np.shape(radius) == np.shape(xinit_in))):
        radius_out = radius
    else:
        raise ValueError('Boxcar radius must a be either an integer, a floating point number, or an ndarray '
                         'with the same shape and size as xinit_in')


    # Figure out dimensions of xinit
    dim = xinit_in.shape
    npix = dim[0]
    ndim = xinit_in.ndim
    if (ndim == 1):
        nTrace = 1
    else:
        nTrace = dim[1]

    ncen = xinit_in.size

    xinit = xinit_in.flatten()
    # Create xnew, xerr
    xnew = xinit.astype(float)
    xerr = np.full(ncen,999.0)

    if npix > fimage.shape[0]:
        raise ValueError('The number of pixels in xinit npix={:d} will run of the image nspec={:d}'.format(npix,fimage.shape[0]))

    if ycen is None:
        if ndim == 1:
            ycen = np.arange(npix, dtype='int')
        elif ndim == 2:
            ycen = np.outer(np.arange(npix, dtype='int'), np.ones(nTrace, dtype='int'))
        else:
            raise ValueError('xinit is not 1 or 2 dimensional')
    else: # check values of input ycen
        if (ycen.min() < 0) | (ycen.max() > (fimage.shape[0] - 1)):
            raise ValueError('Input ycen values will run off the fimage')

    ycen_out = ycen.astype(int)
    ycen_out = ycen_out.flatten()

    if np.size(xinit) != np.size(ycen_out):
        raise ValueError('Number of elements in xinit and ycen must be equal')

#    if npix != fimage.shape[0]:
#        raise ValueError('Number of elements in xinit npix = {:d} does not match spectral dimension of '
#                         'input image {:d}'.format(npix,fimage.shape[0]))

    if invvar is None:
        invvar = np.zeros_like(fimage) + 1.

    x1 = xinit - radius_out + 0.5
    x2 = xinit + radius_out + 0.5
    ix1 = np.floor(x1).astype(int)
    ix2 = np.floor(x2).astype(int)

    fullpix = int(np.maximum(np.min(ix2-ix1)-1,0))

    sumw = np.zeros_like(xinit)
    sumxw = np.zeros_like(xinit)
    sumwt = np.zeros_like(xinit)
    sumsx1 = np.zeros_like(xinit)
    sumsx2 = np.zeros_like(xinit)
    qbad = np.zeros_like(xinit,dtype=bool)

    # Compute
    for ii in range(0,fullpix+3):
        spot = ix1 - 1 + ii
        ih = np.clip(spot,0,nx-1)
        xdiff = spot - xinit
        #
        wt = np.clip(radius_out - np.abs(xdiff) + 0.5,0,1) * ((spot >= 0) & (spot < nx))
        sumw = sumw + fimage[ycen_out,ih] * wt
        sumwt = sumwt + wt
        sumxw = sumxw + fimage[ycen_out,ih] * xdiff * wt
        var_term = wt**2 / (invvar[ycen_out,ih] + (invvar[ycen_out,ih] == 0))
        sumsx2 = sumsx2 + var_term
        sumsx1 = sumsx1 + xdiff**2 * var_term
        #qbad = qbad or (invvar[ycen_out,ih] <= 0)
        #qbad = np.any([qbad, invvar[ycen_out,ih] <= 0], axis=0)
        qbad = qbad | (invvar[ycen_out,ih] <= 0)



    # Fill up
    good = (sumw > 0) &  (~qbad)
    if np.sum(good) > 0:
        delta_x = sumxw[good]/sumw[good]
        xnew[good] = delta_x + xinit[good]
        xerr[good] = np.sqrt(sumsx1[good] + sumsx2[good]*delta_x**2)/sumw[good]

    bad = np.any([np.abs(xnew-xinit) > radius_out + 0.5,xinit < radius_out - 0.5,xinit > nx - 0.5 - radius_out],axis=0)
    if np.sum(bad) > 0:
        xnew[bad] = xinit[bad]
        xerr[bad] = 999.0

    # Reshape to the right size for output if more than one trace was input
    if ndim > 1:
        xnew = xnew.reshape(npix,nTrace)
        xerr = xerr.reshape(npix,nTrace)

    # Return
    return xnew, xerr


def trace_gweight(fimage, xinit_in, sigma = 1.0, ycen = None, invvar=None, maskval=-999999.9):
    ''' Routine to recenter a trace using gaussian-weighted centroiding. Specifically the flux in the image is weighted
    by the integral of a Gaussian over a pixel. Port of idlutils trace_gweight.pro algorithm


    Parameters
    ----------
    fimage: 2D ndarray
      Image for tracing which shape (nspec, nspat)

    xinit: ndarray
      Initial guesses for spatial direction trace. This can either be an 2-d  array with shape
         (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.


    Optional Parameters:
    --------------------
    sigma :  float or ndarray, default = 1.0
         Sigma of Gaussian for centroiding in floating point pixels. This can be either be input as a scalar or as an array to perform
         centroiding with a varaible sigma. If an array is input it must have the same size and shape as xinit_in, i.e.
         a 2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.

    ycen: ndarray, default = None
      Optionally y-position of trace can be provided. It should be an integer array the same size as x-trace (nspec, nTrace). If
      not provided np.arange(nspec) will be assumed for each trace.

    invvar: ndarray, default = None
         Inverse variance array for the image. Array with shape (nspec, nspat) matching fimage

    Returns
    -------
    xnew:   ndarray
         Recentroided trace. The output will have the same shape as xinit i.e.  an 2-d  array with shape (nspec, nTrace)
         array if multiple traces were input, or a 1-d array with shape (nspec) for
         the case of a single trace.
    xerr: ndarray
         Formal propagated error on recentroided trace. These errors will only make sense if invvar is passed in, since
         otherwise invvar is set to 1.0.   The output will have the same shape as xinit i.e.  an 2-d  array with shape
         (nspec, nTrace) array if multiple traces were input, or a 1-d array with shape (nspec) for the case of a single
         trace. Locations where the gaussian weighted centroid falls off the image, will have this error set to 999 and
         will their xnew values set to that of the input trace. These should thus be masked in any fit via a condition like (xerr < 999)

     Revision History
     ----------------
     Python port of trace_fweight.pro from IDLUTILS
     17-Jan-2000  Written by Scott Burles, Chicago
     27-Jun-2018  Ported to python by X. Prochaska and J. Hennawi
    '''


    # Init
    nx = fimage.shape[1]

    # Checks on radius
    if isinstance(sigma,(int,float)):
        sigma_out = sigma
    elif ((np.size(sigma)==np.size(xinit_in)) & (np.shape(sigma) == np.shape(xinit_in))):
        sigma_out = sigma
    else:
        raise ValueError('Gaussian sigma must a be either an integer, a floating point number, or an ndarray '
                         'with the same shape and size as xinit_in')


    # Figure out dimensions of xinit
    dim = xinit_in.shape
    npix = dim[0]
    ndim = xinit_in.ndim
    if (ndim == 1):
        nTrace = 1
    else:
        nTrace = dim[1]

    ncen = xinit_in.size

    xinit = xinit_in.flatten()
    # Create xnew, xerr
    xnew = xinit.astype(float)
    xerr = np.full(ncen,999.)

    if npix > fimage.shape[0]:
        raise ValueError(
            'The number of pixels in xinit npix={:d} will run of the image nspec={:d}'.format(npix, fimage.shape[0]))

    if ycen is None:
        if ndim == 1:
            ycen = np.arange(npix, dtype=int)
        elif ndim == 2:
            ycen = np.outer(np.arange(npix, dtype='int'), np.ones(nTrace, dtype='int'))
        else:
            raise ValueError('xinit is not 1 or 2 dimensional')
    else: # check value of input ycen
        if (ycen.min() < 0) | (ycen.max() > (fimage.shape[0] - 1)):
            raise ValueError('Input ycen values will run off the fimage')

    ycen_out = ycen.astype(int)
    ycen_out = ycen_out.flatten()


    if np.size(xinit) != np.size(ycen_out):
        raise ValueError('Number of elements in xinit and ycen must be equal')


#    if npix != fimage.shape[0]:
#        raise ValueError('Number of elements in xinit npix = {:d} does not match spectral dimension of '
#                         'input image {:d}'.format(npix,fimage.shape[0]))

    if invvar is None:
        invvar = np.zeros_like(fimage) + 1.

    var = utils.inverse(invvar)
    # More setting up
    x_int = np.rint(xinit).astype(int)
    nstep = 2*int(3.0*np.max(sigma_out)) - 1

    weight = np.zeros_like(xinit)
    numer  = np.zeros_like(xinit)
    numer_var  = np.zeros_like(xinit)
    meanvar = np.zeros_like(xinit)
    qbad = np.zeros_like(xinit).astype(bool)

    nby2 = nstep//2
    for i in range(nstep):
        xh = x_int - nby2 + i
        xtemp = (xh - xinit - 0.5)/sigma_out/np.sqrt(2.0)
        g_int = (erf(xtemp+1./sigma_out/np.sqrt(2.0)) - erf(xtemp))/2.
        xs = np.fmin(np.fmax(xh,0),(nx-1))
        cur_weight = fimage[ycen_out, xs] * (invvar[ycen_out, xs] > 0) * g_int * ((xh >= 0) & (xh < nx))
        weight += cur_weight
        numer += cur_weight * xh
        numer_var += var[ycen_out,xs]*(invvar[ycen_out, xs] > 0) * (g_int**2) *((xh >= 0) & (xh < nx))
        # Below is Burles calculation of the error which I'm not following
        meanvar += cur_weight * cur_weight * (xinit-xh)**2/(invvar[ycen_out, xs] + (invvar[ycen_out, xs] == 0))
        qbad = qbad | (xh < 0) | (xh >= nx)
        # bad = np.any([bad, xh < 0, xh >= nx], axis=0)

    # Masking
    good = (~qbad) & (weight > 0)
    if np.sum(good) > 0:
        xnew[good] = numer[good]/weight[good]
        xerr[good] = np.sqrt(numer_var[good])/weight[good]
#        xerr[good] = np.sqrt(meanvar[good])/weight[good] # Burles error which I don't follow

    # For pixels with large deviations, simply reset to initial values and set large error as with trace_fweight
    bad = np.any([np.abs(xnew-xinit) > 2*sigma_out + 0.5,xinit < 2*sigma_out - 0.5,xinit > nx - 0.5 - 2*sigma_out],axis=0)
    if np.sum(bad) > 0:
        xnew[bad] = xinit[bad]
        xerr[bad] = 999.0

    # Reshape to the right size for output if more than one trace was input
    if ndim > 1:
        xnew = xnew.reshape(npix,nTrace)
        xerr = xerr.reshape(npix,nTrace)

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


# ToDo 1) Add an option where the user specifies the number of slits, and so it takes only the highest peaks
# from detect_lines
def trace_refine(filt_image, edges, edges_mask, ncoeff=5, npca = None, pca_explained_var = 99.8, coeff_npoly_pca = 3,
                 fwhm=3.0, sigthresh=100.0, trc_thresh=10.0, trc_median_frac=0.01, upper=2.0, lower=2.0, debug=False, fweight_boost=1.,
                 maxrej=1, smash_range=(0, 1)):
    """
    Refines input trace using a PCA analysis

    Args:
        filt_image: ndarray
          Filtered image (usually Sobolev)
        edges: ndarray
          Current set of edges
        edges_mask: ndarray
          Mask, of edges;  1 = Good
        ncoeff: int, optional
          Order of polynomial for fits
        npca: int, optional
          If provided, restrict the PCA analysis to this order
        pca_explained_var: float, optional
          If npca=None, the PCA will add coefficients until explaining this amount of the variance
        coeff_npoly_pca: int, optional
        fwhm: float, optional
          Size used for tracing (fweight and gweight)
        trc_thresh: float, optional
          Threshold for masking pixels when the tracing is done with iter_tracefit. Basically we extract the filt_image
          with a boxcar extraction, median filter it with a kernel that is trc_median_frac*nspec, and then mask all pixels
          which have an extracted value < trc_thresh in the fitting.
        fit: float, optional
          Threshold for an edge to be included

        upper: float, optional
        lower: float, optional
        debug:
        fweight_boost: float, optional
          Boost on fwhm for fweight.  This was 3.0 at one point (and that may be preferred for echelle instruments)
        maxrej: int, optional
          Rejection parameter for PCA.  1 makes the rejection go slowly (preferred)
        smash_range: tuple, optional
          Spectral range to smash (in fraction of nspec) when finding slit edges, e.g. (0.5, 1.0)
          If not provided, all rows are smashed

    Returns:
        trace_dict: dict
          dict containing the edge output

    """

    # edges_mask True = Good, Bad = False
    # filt image has left as positive, right as negative
    nedges = edges.shape[1]
    nspec = filt_image.shape[0]
    nspat = filt_image.shape[1]
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    # ToDo I would take this stuff out of this routine and put it in the calling routine. For the iterations, we've
    # already done fits, so need to re-fit
    edge_spec = np.outer(np.ones(nedges), spec_vec)
    tset = pydl.xy2traceset(edge_spec, edges.T, ncoeff=ncoeff, maxdev=5.0, maxiter=25, invvar=edges_mask.T.astype(float))
    edges_fit = tset.yfit.T

    # ToDO I think this is part is still okay
    spat_not_junk = np.sum(edges_mask, 1)
    iref = int(np.round(np.sum(spat_not_junk * spec_vec)/np.sum(spat_not_junk)))
    edges_ref = edges_fit[iref, :]

    msgs.info('PCA modeling {:d} slit edges'.format(nedges))
    pca_fit, poly_fit_dict, pca_mean, pca_vectors = extract.pca_trace(
        edges_fit, npca=npca, pca_explained_var = pca_explained_var,coeff_npoly=coeff_npoly_pca, order_vec=edges_ref,
        xinit_mean=edges_ref, upper = upper, lower = lower, minv = 0.0, maxv = float(nspat-1), debug= debug,
        maxrej=maxrej)

    # pca_poly_fit is list
    npca_out = len(poly_fit_dict)
    pca_coeff_spat = np.zeros((nspat, npca_out))
    for idim in range(npca_out):
        pca_coeff_spat[:, idim] = utils.func_val(poly_fit_dict[str(idim)]['coeffs'], spat_vec, 'polynomial',
                                                 minx=poly_fit_dict[str(idim)]['minv'],maxx=poly_fit_dict[str(idim)]['maxv'])

    trace_model = np.outer(pca_mean, np.ones(nspat)) + (np.dot(pca_coeff_spat, pca_vectors)).T + np.arange(nspat)
    # JFH What should this aperture size be? I think fwhm=3.0 since that is the width of the sobel filter
    trace_model_left = trace_model - fwhm/2.0
    trace_model_righ = trace_model + fwhm/2.0
#    trace_model_left = trace_model - 0.5 #fwhm/2.0
#    trace_model_righ = trace_model + 0.5 #fwhm/2.0

    msgs.info('Extracting filt_image along curved edge traces')
    filt_extract = extract.extract_asymbox2(filt_image, trace_model_left, trace_model_righ)
    if debug:
        ginga.show_image(filt_extract, chname ='rectified filt_image')

    # Smash the filtered image
    #   For instruments where light runs on only a portion of the detector,
    #   one is recommended to smash only that portion
    smash_spec = (int(smash_range[0]*nspec), int(smash_range[1]*nspec))
    filt_smash_mean, filt_smash_median, filt_smash_sig = sigma_clipped_stats(
        filt_extract[smash_spec[0]:smash_spec[1],:], axis=0, sigma=4.0)

    # Perform initial finding with a very liberal threshold
    # Put in Gaussian smoothing here?

    kernel_size = int(np.ceil(nspec*trc_median_frac)//2 * 2 + 1)  # This ensure kernel_size is odd
    trace_dict = {}
    for key,sign in zip(['left','right'], [1., -1.]):
        ypeak, _, edge_start, sigma_pk, _, igd, _, _ = arc.detect_lines(
            sign*filt_smash_mean, cont_subtract=False, fwhm=fwhm, input_thresh=sigthresh,
            max_frac_fwhm = 4.0, min_pkdist_frac_fwhm=5.0, debug=debug)
        # ToDO add error catching here if there are no peaks found!
        trace_dict[key] = {}
        trace_dict[key]['start'] = edge_start[igd]
        trace_dict[key]['nstart'] = len(edge_start[igd])
        msgs.info('Found {:d} {:s} slit edges'.format(len(edge_start[igd]),key))
        trace_crutch = trace_model[:, np.round(edge_start[igd]).astype(int)]
        msgs.info('Iteratively tracing {:s} edges'.format(key))
        # Extract a flux about the trace_crutch to mask out pixels that have no signal
        flux_fw = extract.extract_boxcar(np.fmax(sign*filt_image, -1.0*sign), trace_crutch,  fwhm)
        flux_fw_med = signal.medfilt(flux_fw, kernel_size=(1,kernel_size))
        trc_inmask_fw = (flux_fw_med.T > trc_thresh) & (trace_crutch > 0) & (trace_crutch < (nspat-1))
        trace_fweight, _, _, _ = extract.iter_tracefit(np.fmax(sign*filt_image, -1.0*sign), trace_crutch, ncoeff,
                                                       trc_inmask = trc_inmask_fw, fwhm=fweight_boost*fwhm, niter=9)
        # Extract a flux about the trace_fweight to mask out pixels that have no signal
        flux_gw = extract.extract_boxcar(np.fmax(sign*filt_image, -1.0*sign), trace_fweight,  fwhm)
        flux_gw_med = signal.medfilt(flux_gw, kernel_size=(1,kernel_size))
        trc_inmask_gw = (flux_gw_med.T > trc_thresh) & (trace_fweight > 0) & (trace_fweight < (nspat-1))
        trace_gweight, _, _, _ = extract.iter_tracefit(np.fmax(sign*filt_image, -1.0*sign), trace_fweight, ncoeff,
                                                       trc_inmask = trc_inmask_gw, fwhm=fwhm,gweight=True, niter=6)
        trace_dict[key]['trace'] = trace_gweight

    color = dict(left='green', right='red')
    if debug:
        viewer, ch = ginga.show_image(filt_image, chname='filt_image')
        for key in trace_dict.keys():
            for kk in range(trace_dict[key]['nstart']):
                ginga.show_trace(viewer, ch, trace_dict[key]['trace'][:,kk],trc_name = key + '_' + str(kk), color = color[key])


    return trace_dict


def slit_trace_qa(frame, ltrace, rtrace, slitmask, extslit, setup, desc="",
                  normalize=True, use_slitid=None, out_dir=None):
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
    nspec, nspat = frame.shape
    method = inspect.stack()[0][3]
    outfile = qa.set_qa_filename(setup, method, out_dir=out_dir)
    nslits = ltrace.shape[1]
    spec_vec = np.arange(nspec)
    slitcen = (ltrace + rtrace)/2.0
    # Normalize flux in the traces
    if normalize:
        sclmin, sclmax = 0.4, 1.1
        nrm_frame = np.zeros_like(frame)
        for islit in range(nslits):
            # Extract the flux down this trace
            flat_counts = extract.extract_boxcar(frame,slitcen[:,islit],1.5)/3.0
            trc_norm = np.outer(flat_counts,np.ones(nspat))
            slitind = slitmask == islit
            nrm_frame[slitind] = frame[slitind]/(trc_norm[slitind] + (trc_norm[slitind] <= 0.0))
    else:
        nrm_frame = frame.copy()
        nrm_frame[frame > 0.0] = np.sqrt(nrm_frame[frame > 0.0])
        sclmin, sclmax = plot.zscale(nrm_frame)

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
    for islit in range(nslits):
        if extslit[islit] is True:
            ptyp = ':'
        else:
            ptyp = '--'
        # Left
        plt.plot(ltrace[:, islit]+0.5, spec_vec, 'r'+ptyp, linewidth=0.3, alpha=0.7)
        # Right
        plt.plot(rtrace[:, islit]+0.5, spec_vec, 'c'+ptyp, linewidth=0.3, alpha=0.7)
        # Label
        if use_slitid is not None:
            slitid, _, _ = get_slitid(frame.shape, ltrace, rtrace, islit, ypos=0.5)
            lbl = 'S{:04d}'.format(slitid)
        else:
            lbl = '{0:d}'.format(islit+1)
        plt.text(0.5*(ltrace[iy, islit]+rtrace[iy, islit]), spec_vec[iy], lbl, color='green', ha='center', size='small')
    # Title
    tstamp = qa.gen_timestamp()
    if desc == "":
        plt.suptitle(tstamp)
    else:
        plt.suptitle(desc+'\n'+tstamp)

    # Write
    plt.savefig(outfile, dpi=800)
    plt.close()

    plt.rcdefaults()


def slit_spat_pos(tslits_dict):
    """
    Generate an array of the slit spat positions
    from the tslits_dict

    Parameters:
        tslits_dict (dict):  Trace slits dict

    Returns:
        np.ndarray

    """
    return (tslits_dict['slit_left'][tslits_dict['nspec']//2, :] +
            tslits_dict['slit_righ'][tslits_dict['nspec']//2,:]) /2/tslits_dict['nspat']
