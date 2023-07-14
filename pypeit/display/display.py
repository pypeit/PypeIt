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
from pypeit import io
from pypeit import utils

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
            subprocess.Popen(['ginga', '--modules=RC,SlitWavelength'])

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
                           '`ginga --modules=RC,SlitWavelength` on the command line.  In either case, wait for '
                           'the ginga viewer to open and try the pypeit command again.')
            return viewer

        if raise_err:
            raise ValueError
        else:
            msgs.warn('Problem connecting to Ginga.  Launch an RC Ginga viewer and '
                      'then continue: \n    ginga --modules=RC,SlitWavelength')

    # Return
    return viewer


def show_image(inp, chname='Image', waveimg=None, mask=None, exten=0, cuts=None, clear=False,
               wcs_match=False):
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
        waveimg (:obj:`numpy.ndarray`, optional):
            Wavelength image
        mask (:class:`~pypeit.images.ImageBitMaskArray`, optional):
            A bitmask array that designates a pixel as being masked.  Currently
            this is only used when displaying the spectral extraction result.
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
            Raised if `cuts` is provided and does not have two elements.
    """
    # Input checks
    if cuts is not None and len(cuts) != 2:
        raise ValueError('Input cuts must only have two elements, the lower and upper cut.')

    # Instantiate viewer
    viewer = connect_to_ginga()
    # Read or set the image data.  This will fail if the input is a
    # string and astropy.io.fits cannot read the image.
    img = io.fits_open(inp)[exten].data if isinstance(inp, str) else inp

    if clear:
        clear_all()

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

    # These commands set up the viewer. They can be found at
    # ginga/ginga/ImageView.py
    canvas = viewer.canvas(ch._chname)
    out = canvas.clear()
    out = ch.set_color_map('ramp')
    out = ch.set_intensity_map('ramp')
    out = ch.set_color_algorithm('linear')
    out = ch.restore_contrast()
    out = ch.restore_cmap()
    if cuts is not None:
        out = ch.cut_levels(float(cuts[0]), float(cuts[1]))

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
        # Select pixels on any slit
        onslit = mask.flagged('OFFSLITS', invert=True)

        # These are the pixels that were masked by the bpm
        spec_bpm, spat_bpm = np.where(mask.bpm & onslit)
        nbpm = len(spec_bpm)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_bpm = [dict(type='point', args=(float(spat_bpm[i]), float(spec_bpm[i]), 2),
                           kwargs=dict(style='plus', color='magenta')) for i in range(nbpm)]

        # These are the pixels that were masked by LACOSMICS
        spec_cr, spat_cr = np.where(mask.cr & onslit)
        ncr = len(spec_cr)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_cr = [dict(type='point', args=(float(spat_cr[i]), float(spec_cr[i]), 2),
                             kwargs=dict(style='plus', color='cyan')) for i in range(ncr)]

        # These are the pixels that were masked by the extraction
        spec_ext, spat_ext = np.where(mask.extract & onslit)
        next = len(spec_ext)
        # note: must cast numpy floats to regular python floats to pass the remote interface
        points_ext = [dict(type='point', args=(float(spat_ext[i]), float(spec_ext[i]), 2),
                            kwargs=dict(style='plus', color='red')) for i in range(next)]

        # Get the "rest" of the flags
        other_flags = list(mask.bit_keys())
        other_flags.remove('OFFSLITS')
        other_flags.remove('BPM')
        other_flags.remove('CR')
        other_flags.remove('EXTRACT')
        # Determine where any of them are flagged
        other_bpm = mask.flagged(flag=other_flags)

        # These are the pixels that were masked for any other reason (and on a slit)
        spec_oth, spat_oth = np.where(other_bpm & onslit)
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


def show_points(viewer, ch, spec, spat, color='cyan', legend=None, legend_spec=None, legend_spat=None):
    """
    Plot points in a ginga viewer

    Parameters
    ----------
    viewer (ginga.util.grc.RemoteClient):
        Ginga RC viewer
    ch (ginga.util.grc._channel_proxy):
        Ginga channel
    spec (list):
        List of spectral positions on image to plot
    spat (list):
        List of spatial positions on image to plot
    color (str):
        Color for points
    legend (str):
        Label for a legeng
    legend_spec (float):
        Spectral pixel loation for legend
    legend_spat (float):
        Pixel loation for legend

    """
    canvas = viewer.canvas(ch._chname)
    npoints = len(spec)
    canvas_list = [dict(type='point', args=(float(spat[i]), float(spec[i]), 2),
                         kwargs=dict(style='plus', color=color)) for i in range(npoints)]
    if legend is not None:
        spec_label = np.mean(np.array(spec)) if legend_spec is None else legend_spec
        spat_label = (np.mean(np.array(spat)) + 30) if legend_spat is None else legend_spat
        text = [dict(type='text', args=(spat_label, spec_label, legend), kwargs=dict(color=color, fontsize=20))]
        canvas_list += text

    canvas.add('constructedcanvas', canvas_list)


# TODO: Should we continue to allow rotate as an option?
def show_slits(viewer, ch, left, right, slit_ids=None, left_ids=None, right_ids=None, maskdef_ids=None, spec_vals=None,
               rotate=False, pstep=50, clear=False, synced=True):
    r"""
    Overplot slits on the image in Ginga in the given channel

    Args:
        viewer (ginga.util.grc.RemoteClient):
            Ginga RC viewer
        ch (ginga.util.grc._channel_proxy):
            Ginga channel
        left (`numpy.ndarray`_):
            Array with spatial position of left slit edges. Shape must be :math:`(N_{\rm
            spec},)` or :math:`(N_{\rm spec}, N_{\rm l-edge})`, and
            can be different from ``right`` unless ``synced`` is
            True.
        right (`numpy.ndarray`_):
            Array with spatial position of right slit edges. Shape must be :math:`(N_{\rm
            spec},)` or :math:`(N_{\rm spec}, N_{\rm r-edge})`, and
            can be different from ``left`` unless ``synced`` is True.
        spec_vals (`numpy.ndarray`_, optional):
            Array with spectral position of left and right slit edges. Shape must be :math:`(N_{\rm
            spec},)` or :math:`(N_{\rm spec}, N_{\rm r-edge})`. Currently it is only possible to input
            a single set of spec_vals for both ``left`` and ``right`` edges, but not possible to pass distinct
            spec_vals for left and right. If not passed in the default of np.arange(:math:`(N_{\rm spec},)`) will be used.
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
        maskdef_ids (:obj:`int`, array-like, optional):
            slitmask IDs assigned to each slits. If None, IDs will not
            be shown. If not None, shape must be
            :math:`(N_{\rm slits},)`. These are only used if
            ``synced`` is True.
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

    # Spectral pixel location
    if spec_vals is not None:
        y = spec_vals.reshape(-1,1) if spec_vals.ndim == 1 else spec_vals
    else:
        y = np.repeat(np.arange(nspec).astype(float)[:, np.newaxis], nright, axis=1)

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
        if maskdef_ids is not None and maskdef_ids.size == nslits:
            _maskdef_ids = np.atleast_1d(maskdef_ids)
        else:
            _maskdef_ids = None
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

    # Label positions
    top = int(2*nspec/3.)
    bot = int(nspec/2.)

    # Plot lefts. Points need to be int or float. Use of .tolist() on
    # each array insures this
    canvas_list = [dict(type=str('path'),
                        args=(list(zip(y[::pstep, i].tolist(), _left[::pstep,i].tolist())),) if rotate
                        else (list(zip(_left[::pstep,i].tolist(), y[::pstep, i].tolist())),),
                        kwargs=dict(color=str('green'))) for i in range(nleft)]
    if not synced:
        # Add text
        canvas_list += [dict(type='text',
                             args=(float(y[bot, i]), float(_left_id_loc[bot,i]), str('S{0}'.format(_left_ids[i]))) if rotate
                             else (float(_left_id_loc[bot,i]), float(y[bot, i]), str('S{0}'.format(_left_ids[i]))),
                             kwargs=dict(color=str('aquamarine'), fontsize=20., rot_deg=90.)) for i in range(nleft)]

    # Plot rights. Points need to be int or float. Use of .tolist() on
    # each array insures this
    canvas_list += [dict(type=str('path'),
                        args=(list(zip(y[::pstep, i].tolist(), _right[::pstep,i].tolist())),) if rotate
                        else (list(zip(_right[::pstep,i].tolist(), y[::pstep, i].tolist())),),
                        kwargs=dict(color=str('magenta'))) for i in range(nright)]
    if not synced:
        # Add text
        canvas_list += [dict(type='text',
                             args=(float(y[bot, i]), float(_right_id_loc[bot,i]), str('S{0}'.format(_right_ids[i]))) if rotate
                             else (float(_right_id_loc[bot,i]), float(y[bot, i]), str('S{0}'.format(_right_ids[i]))),
                             kwargs=dict(color=str('magenta'), fontsize=20., rot_deg=90.)) for i in range(nright)]

    canvas.add('constructedcanvas', canvas_list)

    # Plot slit labels, if synced
    if synced:
        # Slit IDs
        canvas_list += [dict(type='text',
                             args=(float(y[bot, i]), float(_slit_id_loc[bot,i])-400, str('S{0}'.format(_slit_ids[i]))) if rotate
                             else (float(_slit_id_loc[bot,i]), float(y[bot, i])-400, str('S{0}'.format(_slit_ids[i]))),
                             kwargs=dict(color=str('aquamarine'), fontsize=20., rot_deg=90.)) for i in range(nslits)]
        # maskdef_ids
        if _maskdef_ids is not None:
            canvas_list += [dict(type='text',
                                 args=(float(y[bot, i]), float(_slit_id_loc[bot,i])-250, str('{0}'.format(_maskdef_ids[i]))) if rotate
                                 else (float(_slit_id_loc[bot,i]), float(y[bot, i])-250, str('{0}'.format(_maskdef_ids[i]))),
                                 kwargs=dict(color=str('cyan'), fontsize=20., rot_deg=90.)) for i in range(nslits)]

    canvas.add('constructedcanvas', canvas_list)


def show_trace(viewer, ch, trace, trc_name=None, maskdef_extr=None, manual_extr=None, clear=False,
               rotate=False, pstep=50, yval=None):
    r"""

    Args:
        viewer (ginga.util.grc.RemoteClient):
            Ginga RC viewer.
        ch (ginga.util.grc._channel_proxy):
            Ginga channel.
        trace (`numpy.ndarray`_):
            Array with spatial position of the object traces on the detector.
            Shape must be :math:`(N_{\rm spec},)` or :math:`(N_{\rm spec}, N_{\rm trace})`.
        trc_name (`numpy.ndarray`_, optional):
            Array with Trace names. Shape must be :math:`(N_{\rm trace},)`.
        maskdef_extr (`numpy.ndarray`_, optional):
            Array with the maskdef extraction flags. Shape must be :math:`(N_{\rm trace},)`.
        manual_extr (`numpy.ndarray`_, optional):
            Array with the manual extraction flags. Shape must be :math:`(N_{\rm trace},)`.
        clear (:obj:`bool`, optional):
            Clear the canvas?
        rotate (:obj:`bool`, optional):
            Rotate the image?
        pstep (:obj:`bool`, optional):
            Show every pstep point of the edges as opposed to *every* point, recommended for speed.
        yval (`numpy.ndarray`_, optional):
            Array with spectral position of the object traces. Shape must be :math:`(N_{\rm spec},)`
            or :math:`(N_{\rm spec}, N_{\rm trace})`. If not passed in, the default of
            np.arange(:math:`(N_{\rm spec},)`) will be used.

    """
    # Canvas
    canvas = viewer.canvas(ch._chname)
    if clear:
        canvas.clear()

    if trace.ndim == 1:
        trace = trace.reshape(-1,1)

    # Show
    if yval is None:
        y = np.repeat(np.arange(trace.shape[0]).astype(float)[:, None], trace.shape[1], axis=1)
    else:
        y = yval.reshape(-1, 1) if yval.ndim == 1 else yval

    canvas_list = []
    for i in range(trace.shape[1]):
        if maskdef_extr[i]:
            color = '#f0e442'
        elif manual_extr[i]:
            color = '#33ccff'
        else:
            color = 'orange'
        canvas_list += [dict(type=str('path'),
                        args=(list(zip(y[::pstep,i].tolist(), trace[::pstep,i].tolist())),) if rotate
                        else (list(zip(trace[::pstep,i].tolist(), y[::pstep,i].tolist())),),
                        kwargs=dict(color=color))]
        # Text
        ohf = len(trace[:,i])//2
        # Do it
        canvas_list += [dict(type='text',args=(float(y[ohf,i]), float(trace[ohf,i]), str(trc_name[i])) if rotate
                             else (float(trace[ohf,i]), float(y[ohf,i]), str(trc_name[i])),
                             kwargs=dict(color=color, fontsize=17., rot_deg=90.))]

    canvas.add('constructedcanvas', canvas_list)


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


def clear_all(allow_new=False):
    """
    Clear all of the ginga canvasses.

    Args:
        allow_new (:obj:`bool`, optional):
            Allow a subprocess to be called to execute a new ginga viewer if one
            is not already running.  See :func:`connect_to_ginga`.
    """
    viewer = connect_to_ginga(allow_new=allow_new)
    shell = viewer.shell()
    chnames = shell.get_channel_names()
    for ch in chnames:
        shell.delete_channel(ch)


def show_tilts(viewer, ch, tilt_traces, yoff=0., xoff=0., points=True, nspec=None, pstep=3, clear_canvas=False):
    """
    Show the arc tilts on the input channel

    Args:
        viewer (ginga.util.grc.RemoteClient):
            Ginga RC viewer
        ch (ginga.util.grc._channel_proxy):
            Ginga channel
        tilt_traces (`astropy.table.Table`_):
            Table containing the traced and fitted tilts
        yoff (float, optional):
            Offset tilts by this amount
        xoff (float, optional):
            Offset tilts by this amount
        points (bool, optional):
            Plot the Gaussian-weighted tilt centers
        nspec (int, optional):
            Number of spectral pixels in the TiltImage
        pstep (int, optional):
            Show every pstep point of the tilts as opposed to *every*
            point, recommended for speed.
        clear_canvas (bool, optional):
            Clear the canvas first?

    """
    if tilt_traces is None:
        return msgs.error('No tilts have been traced or fitted')

    canvas = viewer.canvas(ch._chname)
    if clear_canvas:
        canvas.clear()

    canvas_list = []

    # Plot traced tilts
    # We just plot the points, so we do not need to loop over each slit/line
    # This makes the plotting much very slow, this is why we make it optional by using the points keyword
    if 'goodpix_tilt' in tilt_traces.keys() and tilt_traces['goodpix_tilt'][0].size > 0 and points:
        # note: must cast numpy floats to regular python floats to pass the remote interface
        goodpix_spat = tilt_traces['goodpix_spat'][0] + xoff
        goodpix_tilt = tilt_traces['goodpix_tilt'][0] + yoff
        canvas_list += [dict(type='squarebox', args=(float(goodpix_spat[i]), float(goodpix_tilt[i]), 1),
                             kwargs=dict(color='cyan', fill=False)) for i in range(goodpix_tilt.size)]

    # Plot the 2D fitted tilts
    # loop over each line, this allows to use type='path' and therefore a faster plotting
    if 'good2dfit_lid' in tilt_traces.keys():
        for iline in np.unique(tilt_traces['good2dfit_lid'][0]):
            # good fit
            this_line = tilt_traces['good2dfit_lid'][0] == iline
            if np.any(this_line):
                good2dfit_spat = tilt_traces['good2dfit_spat'][0][this_line] + xoff
                good2dfit_tilt = tilt_traces['good2dfit_tilt'][0][this_line] + yoff
                canvas_list += [dict(type=str('path'),
                                     args=(list(zip(good2dfit_spat[::pstep].tolist(), good2dfit_tilt[::pstep].tolist())),),
                                     kwargs=dict(color='blue', linewidth=1))]

    # Now plot the masked traces and the rejected 2D fits
    # We just plot the points, so we do not need to loop over each slit/line
    # masked traces
    if 'badpix_tilt' in tilt_traces.keys() and tilt_traces['badpix_tilt'][0].size > 0:
        # note: must cast numpy floats to regular python floats to pass the remote interface
        badpix_spat = tilt_traces['badpix_spat'][0] + xoff
        badpix_tilt = tilt_traces['badpix_tilt'][0] + yoff
        canvas_list += [dict(type='squarebox', args=(float(badpix_spat[i]), float(badpix_tilt[i]), 1),
                             kwargs=dict(color='red', fill=False)) for i in range(badpix_tilt.size)]
    # rejected fit
    if 'bad2dfit_tilt' in tilt_traces.keys() and tilt_traces['bad2dfit_tilt'][0].size > 0:
        # note: must cast numpy floats to regular python floats to pass the remote interface
        bad2dfit_spat = tilt_traces['bad2dfit_spat'][0] + xoff
        bad2dfit_tilt = tilt_traces['bad2dfit_tilt'][0] + yoff
        canvas_list += [dict(type='squarebox', args=(float(bad2dfit_spat[i]), float(bad2dfit_tilt[i]), 1),
                             kwargs=dict(color='yellow', fill=False)) for i in range(bad2dfit_tilt.size)]

    # Add text
    text_xpos = 20
    start_ypos = 20
    ypos_step = 0.03*nspec if nspec is not None else 50.
    text_ypos = [start_ypos, start_ypos + ypos_step, start_ypos + 2*ypos_step]
    text_str = ['Masked pixel', 'Rejected in fit', 'Good tilt fit']
    text_color = ['red', 'yellow', 'blue']
    if points:
        text_ypos += [start_ypos + 3*ypos_step]
        text_str += ['Good pixel']
        text_color += ['cyan']
    canvas_list += [dict(type='text', args=(float(text_xpos), float(text_ypos[i]), str(text_str[i])),
                    kwargs=dict(color=text_color[i], fontsize=20)) for i in range(len(text_str))]

    canvas.add('constructedcanvas', canvas_list)





