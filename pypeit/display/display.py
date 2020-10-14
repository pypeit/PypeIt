"""
Module for ginga routines.  Mainly for debugging

.. include:: ../include/links.rst

"""
import os
import numpy as np
import time
from IPython import embed

import subprocess

# A note from ejeschke on how to use the canvas add command in ginga: https://github.com/ejeschke/ginga/issues/720
# c
# The add() command can add any of the shape types that are defined under ginga.canvas.types. The good part is that if you
# go to that directory in the ginga source tree (ginga/canvas/types) and browse the source, you will find a parameter
# description table at the beginning of each type definition, describing each parameter in the type and what it is for.
# Most of the standard geometric types are in basic.py and there are specialized ones in utils.py, astro.py and layer.py. Looking at
# the classes will also tell you which parameters are positional and which are keyword.

from astropy.io import fits

from ginga.util import grc

from pypeit import msgs

def connect_to_ginga(host='localhost', port=9000, raise_err=False, allow_new=False):
    """
    Connect to a RC Ginga.

    Args:
        host (:obj:`str`, optional):
            Host name.
        port (:obj:`int`, optional):
            Probably should remain at 9000
        raise_err (:obj:`bool`, optional):
            Raise an error if no connection is made, otherwise just
            raise a warning and continue
        allow_new (:obj:`bool`, optional):
            Allow a subprocess to be called to execute a new ginga
            viewer if one is not already running.

    Returns:
        RemoteClient: connection to ginga viewer.
    """
    # Start
    viewer = grc.RemoteClient(host, port)
    # Test
    sh = viewer.shell()
    try:
        tmp = sh.get_current_workspace()
    except:
        if allow_new:
            subprocess.Popen(['ginga', '--modules=RC'])

            # NOTE: time.sleep(3) is now insufficient. The loop below
            # continues to try to connect with the ginga viewer that
            # was just instantiated for a maximum number of iterations.
            # If the connection is remains unsuccessful, an error is
            # thrown stating that the connection timed out.
            maxiter = int(1e6)
            for i in range(maxiter):
                try:
                    viewer = grc.RemoteClient(host, port)
                    sh = viewer.shell()
                    tmp = sh.get_current_workspace()
                except:
                    continue
                else:
                    break
            if i == maxiter-1:
                msgs.error('Timeout waiting for ginga to start.  If window does not appear, type '
                           '`ginga --modules=RC` on the command line.  In either case, wait for '
                           'the ginga viewer to open and try the pypeit command again.')
            return viewer

        if raise_err:
            raise ValueError
        else:
            msgs.warn('Problem connecting to Ginga.  Launch an RC Ginga viewer and '
                      'then continue: \n    ginga --modules=RC')

    # Return
    return viewer


def show_image(inp, chname='Image', waveimg=None, bitmask=None, mask=None, exten=0, cuts=None,
               clear=False, wcs_match=False):
    """
    Display an image using Ginga.

    .. todo::
        - implement instrument specific reading
        - use the `mask` as a boolean mask if `bitmask` is not provided.

    Args:
        inp (:obj:`str`, numpy.ndarray):
            The image to view.  If a string is provided, it must be the
            name of a fits image that can be read by `astropy.io.fits`.
        chname (:obj:`str`, optional):
            The name of the ginga channel to use.
        waveimg (:obj:`str`, optional):
            The name of a FITS image with the relevant WCS coordinates
            in its header, mainly for wavelength array.  If None, no WCS
            is used.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            The object used to unpack the mask values.  If this is
            provided, mask must also be provided and the expectation is
            that a extraction image is being shown.
        mask (numpy.ndarray, optional):
            A boolean or bitmask array that designates a pixel as being
            masked.  Currently this is only used when displaying the
            spectral extraction result.
        exten (:obj:`int`, optional):
            The extension of the fits file with the image to show.  This
            is only used if the input is a file name.
        cuts (array-like, optional):
            Initial cut levels to apply when displaying the image.  This
            object must have a length of 2 with the lower and upper
            levels, respectively.
        clear (:obj:`bool`, optional):
            Clear any existing ginga viewer and its channels.
        wcs_match(:obj:`bool`, optional):
            Use this as a reference image for the WCS and match all
            image in other channels to it.

    Returns:
        ginga.util.grc.RemoteClient, ginga.util.grc._channel_proxy: The
        ginga remote client and the channel with the displayed image.

    Raises:
        ValueError:
            Raised if `cuts` is provided and does not have two elements
            or if bitmask is provided but mask is not.
    """
    # Input checks
    if cuts is not None and len(cuts) != 2:
        raise ValueError('Input cuts must only have two elements, the lower and upper cut.')
    if mask is not None and bitmask is None:
        raise ValueError('If providing a mask, must also provide the bitmask.')

    # Read or set the image data.  This will fail if the input is a
    # string and astropy.io.fits cannot read the image.
    img = fits.open(inp)[exten].data if isinstance(inp, str) else inp

    # Instantiate viewer
    viewer = connect_to_ginga()
    if clear:
        # Clear existing channels
        shell = viewer.shell()
        chnames = shell.get_channel_names()
        for ch in chnames:
            shell.delete_channel(ch)
    ch = viewer.channel(chname)
    # Header
    header = {}
    header['NAXIS1'] = img.shape[1]
    header['NAXIS2'] = img.shape[0]

    # Giddy up
#    waveimg = None
    if waveimg is not None:
        sh = viewer.shell()
        args = [chname, chname, grc.Blob(img.tobytes()), img.shape, img.dtype.name, header,
                grc.Blob(waveimg.tobytes()), waveimg.dtype.name, {}]
        sh.call_global_plugin_method('SlitWavelength', 'load_buffer', args, {})
    else:
        ch.load_np(chname, img, 'fits', header)
    canvas = viewer.canvas(ch._chname)

    # These commands set up the viewer. They can be found at
    # ginga/ginga/ImageView.py
    out = canvas.clear()
    if cuts is not None:
        out = ch.cut_levels(cuts[0], cuts[1])
    out = ch.set_color_map('ramp')
    out = ch.set_intensity_map('ramp')
    out = ch.set_color_algorithm('linear')
    out = ch.restore_contrast()
    out = ch.restore_cmap()

    # WCS Match this to other images with this as the reference image?
    if wcs_match:
        # After displaying all the images since up the images with WCS_MATCH
        shell = viewer.shell()
        out = shell.start_global_plugin('WCSMatch')
        out = shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [chname], {})

    # TODO: I would prefer to change the color map to indicate these
    # pixels rather than overplot points. Because for large numbers of
    # masked pixels, this is super slow. Need to ask ginga folks how to
    # do that.

    # If bitmask was passed in, assume this is an extraction qa image
    # and use the mask to identify why each pixel was masked
    if mask is not None:
        # Unpack the bitmask
        bpm, crmask, satmask, minmask, offslitmask, nanmask, ivar0mask, ivarnanmask, extractmask \
                = bitmask.unpack(mask)

        # These are the pixels that were masked by the bpm
        spec_bpm, spat_bpm = np.where(bpm & ~offslitmask)
        nbpm = len(spec_bpm)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_bpm = [dict(type='point', args=(float(spat_bpm[i]), float(spec_bpm[i]), 2),
                           kwargs=dict(style='plus', color='magenta')) for i in range(nbpm)]

        # These are the pixels that were masked by LACOSMICS
        spec_cr, spat_cr = np.where(crmask & ~offslitmask)
        ncr = len(spec_cr)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_cr = [dict(type='point', args=(float(spat_cr[i]), float(spec_cr[i]), 2),
                             kwargs=dict(style='plus', color='cyan')) for i in range(ncr)]

        # These are the pixels that were masked by the extraction
        spec_ext, spat_ext = np.where(extractmask & ~offslitmask)
        next = len(spec_ext)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_ext = [dict(type='point', args=(float(spat_ext[i]), float(spec_ext[i]), 2),
                            kwargs=dict(style='plus', color='red')) for i in range(next)]

        # These are the pixels that were masked for any other reason
        spec_oth, spat_oth = np.where(satmask | minmask | nanmask | ivar0mask | ivarnanmask
                                      & ~offslitmask)
        noth = len(spec_oth)
        # note: must cast numpy floats to regular python floats to pass
        # the remote interface
        points_oth = [dict(type='point', args=(float(spat_oth[i]), float(spec_oth[i]), 2),
                            kwargs=dict(style='plus', color='yellow')) for i in range(noth)]

        nspat = img.shape[1]
        nspec = img.shape[0]
        # Labels for the points
        text_bpm = [dict(type='text', args=(nspat / 2 -40, nspec / 2, 'BPM'),
                           kwargs=dict(color='magenta', fontsize=20))]

        text_cr = [dict(type='text', args=(nspat / 2 -40, nspec / 2 - 30, 'CR'),
                           kwargs=dict(color='cyan', fontsize=20))]

        text_ext = [dict(type='text', args=(nspat / 2 -40, nspec / 2 - 60, 'EXTRACT'),
                          kwargs=dict(color='red', fontsize=20))]

        text_oth = [dict(type='text', args=(nspat / 2 -40, nspec / 2 - 90, 'OTHER'),
                          kwargs=dict(color='yellow', fontsize=20))]

        canvas_list = points_bpm + points_cr + points_ext + points_oth + text_bpm + text_cr \
                        + text_ext + text_oth
        canvas.add('constructedcanvas', canvas_list)

    return viewer, ch

# TODO: Should we continue to allow rotate as an option?
def show_slits(viewer, ch, left, right, slit_ids=None, left_ids=None, right_ids=None, rotate=False,
               pstep=50, clear=False, synced=True):
    r"""
    Overplot slits on the image in Ginga in the given channel

    Args:
        viewer (ginga.util.grc.RemoteClient):
            Ginga RC viewer
        ch (ginga.util.grc._channel_proxy):
            Ginga channel
        left (`numpy.ndarray`_):
            Array with left slit edges. Shape must be :math:`(N_{\rm
            spec},)` or :math:`(N_{\rm spec}, N_{\rm l-edge})`, and
            can be different from ``right`` unless ``synced`` is
            True.
        right (`numpy.ndarray`_):
            Array with right slit edges. Shape must be :math:`(N_{\rm
            spec},)` or :math:`(N_{\rm spec}, N_{\rm r-edge})`, and
            can be different from ``left`` unless ``synced`` is True.
        slit_ids (:obj:`int`, array-like, optional):
            PypeIt ID numbers for the slits. If None, IDs run from -1 to
            :math:`-N_{\rm slits}`. If not None, shape must be
            :math:`(N_{\rm slits},)`. These are only used if
            ``synced`` is True.
        left_ids (:obj:`int`, array-like, optional):
            ID numbers for the left edges. If None, IDs run from -1
            to :math:`-N_{\rm l-edge}`. If not None, shape must be
            :math:`(N_{\rm l-edge},)`. These are only used if
            ``synced`` is False.
        right_ids (:obj:`int`, array-like, optional):
            ID numbers for the right edges. If None, IDs run from -1
            to :math:`-N_{\rm r-edge}`. If not None, shape must be
            :math:`(N_{\rm r-edge},)`. These are only used if
            ``synced`` is False.
        rotate (:obj:`bool`, optional):
            Rotate the image?
        pstep (:obj:`bool`, optional):
            Show every pstep point of the edges as opposed to *every*
            point, recommended for speed.
        clear (:obj:`bool`, optional):
            Clear the canvas?
        synced (:obj:`bool`, optional):
            Flag the left and right traces are synced into slits.
            Otherwise, the edges are treated separately. If True, the
            number of left and right edges must be the same and
            ``left_ids`` and ``right_ids`` are ignored.
    """
    # Setup the trace data and IDs
    _left = left.reshape(-1,1) if left.ndim == 1 else left
    nleft = _left.shape[1]

    _right = right.reshape(-1,1) if right.ndim == 1 else right
    nright = _right.shape[1]
    
    nspec = _left.shape[0]
    if _right.shape[0] != nspec:
        # TODO: Any reason to remove this restriction?
        msgs.error('Input left and right edges have different spectral lengths.')

    # Check input
    if synced:
        if left.shape != right.shape:
            msgs.error('Input left and right traces must have the same shape if they have been '
                       'synchronized into slits.')
        if left_ids is not None or right_ids is not None:
            msgs.warn('For showing synced edges, left and right ID numbers are ignored.')
        nslits = _left.shape[1]
        _left_ids = None
        _right_ids = None
        _slit_ids = np.arange(nslits) if slit_ids is None else np.atleast_1d(slit_ids)
        if len(_slit_ids) != nslits:
            msgs.error('Incorrect number of slit IDs provided.')
        _slit_id_loc = _left + 0.45*(_right - _left)
    else:
        _left_ids = -np.arange(nleft) if left_ids is None else np.atleast_1d(left_ids)
        if len(_left_ids) != nleft:
            msgs.error('Incorrect number of left IDs provided.')
        _left_id_loc = _left*1.05
        _right_ids = -np.arange(nright) if right_ids is None else np.atleast_1d(right_ids)
        if len(_right_ids) != nright:
            msgs.error('Incorrect number of right IDs provided.')
        _right_id_loc = _right*(1-0.05)

    # Canvas
    canvas = viewer.canvas(ch._chname)
    if clear:
        canvas.clear()

    # Spectral pixel location
    y = np.arange(nspec).astype(float)

    # Label positions
    top = int(2*nspec/3.)
    bot = int(nspec/2.)

    # Plot lefts. Points need to be int or float. Use of .tolist() on
    # each array insures this
    for i in range(nleft):
        points = list(zip(y[::pstep].tolist(), _left[::pstep,i].tolist())) if rotate \
                    else list(zip(_left[::pstep,i].tolist(), y[::pstep].tolist()))
        canvas.add(str('path'), points, color=str('green'))
        if not synced:
            # Add text
            xt, yt = float(_left_id_loc[top,i]), float(y[top])
            xb, yb = float(_left_id_loc[bot,i]), float(y[bot])
            if rotate:
                xt, yt = yt, xt
                xb, yb = yb, xb
            canvas.add(str('text'), xb, yb, str('S{0}'.format(_left_ids[i])), color=str('blue'),
                       fontsize=20.)
            #canvas.add(str('text'), xt, yt, str('{0}'.format(i)), color=str('green'), fontsize=20.)

    # Plot rights. Points need to be int or float. Use of .tolist() on
    # each array insures this
    for i in range(nright):
        points = list(zip(y[::pstep].tolist(), _right[::pstep,i].tolist())) if rotate \
                    else list(zip(_right[::pstep,i].tolist(), y[::pstep].tolist()))
        canvas.add(str('path'), points, color=str('red'))
        if not synced:
            # Add text
            xt, yt = float(_right_id_loc[top,i]), float(y[top])
            xb, yb = float(_right_id_loc[bot,i]), float(y[bot])
            if rotate:
                xt, yt = yt, xt
                xb, yb = yb, xb
            canvas.add(str('text'), xb, yb, str('S{0}'.format(_right_ids[i])), color=str('red'),
                       fontsize=20.)
            canvas.add(str('text'), xt, yt, str('{0}'.format(i)), color=str('red'),
                       fontsize=20.)

    # Plot slit labels, if synced
    if not synced:
        return
    for i in range(nslits):
        xt, yt = float(_slit_id_loc[top,i]), float(y[top])
        xb, yb = float(_slit_id_loc[bot,i]), float(y[bot])
        if rotate:
            xt, yt = yt, xt
            xb, yb = yb, xb
        # Slit IDs
        canvas.add(str('text'), xb, yb, str('S{0}'.format(_slit_ids[i])), color=str('blue'),
                   fontsize=20.)
        # TODO -- Fix indices if you really want to show them
        #canvas.add(str('text'), xt, yt, str('{0}'.format(i)), color=str('green'),
        #           fontsize=20.)


def show_trace(viewer, ch, trace, trc_name='Trace', color='blue', clear=False,
               rotate=False, pstep=50, yval=None):
    """

    Args:
        viewer (ginga.util.grc.RemoteClient):
            Ginga RC viewer
        ch (ginga.util.grc._channel_proxy):
            Ginga channel
        trace (np.ndarray):
            x-positions on the detector
        trc_name (str, optional):
            Trace name
        color (str, optional):
            Color for the trace
        clear (bool, optional):
            Clear the canvas?
        rotate (bool, optional):
            Rotate the image?
        pstep (int, optional):
            Show every pstep point of the edges as opposed to *every* point, recommended for speed
        yval (np.ndarray, optional):
            If not provided, it is assumed the input x values track y=0,1,2,3,etc.

    """
    # Canvas
    canvas = viewer.canvas(ch._chname)
    if clear:
        canvas.clear()
    # Show
    if yval is None:
        y = (np.arange(trace.size)[::pstep]).tolist()
    else:
        y = yval[::pstep].tolist()
    trace_list = trace[::pstep].tolist()
    xy = [trace_list, y]
    if rotate:
        xy[0], xy[1] = xy[1], xy[0]
    points = list(zip(xy[0], xy[1]))
    canvas.add(str('path'), points, color=str(color))
    # Text
    ohf = len(trace_list)//2
    xyt = [float(trace_list[ohf]), float(y[ohf])]
    if rotate:
        xyt[0], xyt[1] = xyt[1], xyt[0]
    # Do it
    canvas.add(str('text'), xyt[0], xyt[1], trc_name, rot_deg=90., color=str(color), fontsize=17.)


def clear_canvas(cname):
    """
    Clear the ginga canvas

    Args:
        cname (str):  Channel name

    """
    viewer = connect_to_ginga()
    ch = viewer.channel(cname)
    canvas = viewer.canvas(ch._chname)
    canvas.clear()


def clear_all():
    """
    Clear all of the ginga canvasses

    """
    viewer = connect_to_ginga()
    shell = viewer.shell()
    chnames = shell.get_channel_names()
    for ch in chnames:
        shell.delete_channel(ch)


def show_tilts(viewer, ch, trc_tilt_dict, sedges=None, yoff=0., xoff=0., pstep=1,
               points=True, clear_canvas=False):
    """
    Show the arc tilts on the input channel
      Not sure this is actually working correctly...

    Args:
        viewer (ginga.util.grc.RemoteClient):
            Ginga RC viewer
        ch (ginga.util.grc._channel_proxy):
            Ginga channel
        trc_tilt_dict (dict):
            Contains tilts info
        sedges (tuple, optional):
            Contains slit edges;  passed to show_slits()
        yoff (float, optional):
            Offset tilts by this amount
        xoff (float, optional):
            Offset tilts by this amount
        pstep (int, optional):
            Show every pstep point of the edges as opposed to *every* point, recommended for speed
        points (bool, optional):
            Plot the Gaussian-weighted tilt centers
        clear_canvas (bool, optional):
            Clear the canvas first?

    """
    canvas = viewer.canvas(ch._chname)
    if clear_canvas:
        canvas.clear()

    if sedges is not None:
        show_slits(viewer, ch, sedges[0], sedges[1])

    tilts = trc_tilt_dict['tilts']
    # Crutch is set plot the crutch instead of the tilt itself
    tilts_fit = trc_tilt_dict['tilts_fit']

    tilts_spat = trc_tilt_dict['tilts_spat']
    tilts_mask = trc_tilt_dict['tilts_mask']
    tilts_err = trc_tilt_dict['tilts_err']

    use_tilt = trc_tilt_dict['use_tilt']
    # Show a trace
    nspat = trc_tilt_dict['nspat']
    nspec = trc_tilt_dict['nspec']
    nlines = tilts.shape[1]
    for iline in range(nlines):
        x = tilts_spat[:,iline] + xoff # FOR IMAGING (Ginga offsets this value by 1 internally)
        this_mask = tilts_mask[:,iline]
        this_err = (tilts_err[:,iline] > 900)
        if np.sum(this_mask) > 0:
            if points: # Plot the gaussian weighted tilt centers
                y = tilts[:, iline] + yoff
                # Plot the actual flux weighted centroids of the arc lines that were traced
                goodpix = (this_mask == True) & (this_err == False)
                ngood = np.sum(goodpix)
                if ngood > 0:
                    xgood = x[goodpix]
                    ygood = y[goodpix]
                    # note: must cast numpy floats to regular python floats to pass the remote interface
                    points_good = [dict(type='squarebox',
                                        args=(float(xgood[i]), float(ygood[i]), 0.7),
                                        kwargs=dict(color='cyan',fill=True, fillalpha=0.5)) for i in range(ngood)]
                    canvas.add('constructedcanvas', points_good)
                badpix = (this_mask == True) & (this_err == True)
                nbad = np.sum(badpix)
                if nbad > 0:
                    xbad = x[badpix]
                    ybad = y[badpix]
                    # Now show stuff that had larger errors
                    # note: must cast numpy floats to regular python floats to pass the remote interface
                    points_bad = [dict(type='squarebox',
                                       args=(float(xbad[i]), float(ybad[i]), 0.7),
                                       kwargs=dict(color='red', fill=True,fillalpha=0.5)) for i in range(nbad)]
                    canvas.add('constructedcanvas', points_bad)
                # Now plot the polynomial fits to the the Gaussian weighted centroids
            y = tilts_fit[:, iline] + yoff
            points = list(zip(x[this_mask][::pstep].tolist(),y[this_mask][::pstep].tolist()))
            if use_tilt[iline]:
                clr = 'blue'  # Good line
            else:
                clr = 'yellow'  # Bad line
            canvas.add('path', points, color=clr, linewidth=3)


    canvas.add(str('text'), nspat//2 - 40, nspec//2,      'good tilt fit', color=str('blue'),fontsize=20.)
    canvas.add(str('text'), nspat//2 - 40, nspec//2 - 30, 'bad  tilt fit', color=str('yellow'),fontsize=20.)
    canvas.add(str('text'), nspat//2 - 40, nspec//2 - 60, 'trace good', color=str('cyan'),fontsize=20.)
    canvas.add(str('text'), nspat//2 - 40, nspec//2 - 90, 'trace masked', color=str('red'),fontsize=20.)

