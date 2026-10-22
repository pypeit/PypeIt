"""
Module for ginga routines.  Mainly for debugging

All displays use a ginga viewer that runs in the same python process (see
:mod:`pypeit.display.viewer`).  The ``show_*`` functions add images and
overlays to the viewer; call :func:`show` to raise the viewer window and wait
for the user to close it before continuing.

.. include:: ../include/links.rst

"""
from pathlib import Path

import numpy as np
from IPython import embed

# A note from ejeschke on how to use the canvas add command in ginga: https://github.com/ejeschke/ginga/issues/720
# c
# The add() command can add any of the shape types that are defined under ginga.canvas.types. The good part is that if you
# go to that directory in the ginga source tree (ginga/canvas/types) and browse the source, you will find a parameter
# description table at the beginning of each type definition, describing each parameter in the type and what it is for.
# Most of the standard geometric types are in basic.py and there are specialized ones in utils.py, astro.py and layer.py. Looking at
# the classes will also tell you which parameters are positional and which are keyword.

from pypeit import log
from pypeit import PypeItError
from pypeit import io
from pypeit.display.viewer import get_viewer, show

#: Prefix for the tags of all canvas objects added by PypeIt
TAG_PREFIX = 'pypeit'


def connect_to_ginga(host=None, port=None, raise_err=False, allow_new=False):
    """
    Return the PypeIt ginga viewer.

    .. deprecated::
        PypeIt no longer connects to a separate ginga process.  This function
        is kept for backwards compatibility and simply returns the viewer from
        :func:`pypeit.display.viewer.get_viewer`; all arguments are ignored.

    Returns
    -------
    :class:`~pypeit.display.viewer.PypeItViewer`
        The in-process ginga viewer.
    """
    return get_viewer()


def _canvas(ch):
    """Return the canvas of a channel's main viewer."""
    return ch.fitsimage.get_canvas()


def _add(viewer, ch, obj, tagpfx=TAG_PREFIX):
    """Add a canvas object to a channel, tagged as a PypeIt object."""
    _canvas(ch).add(obj, tagpfx=tagpfx, redraw=False)
    # NOTE: Unlike a redraw of the canvas, this redraw is held while the
    # viewer's redraws are suppressed (see ``suppress_redraw``).
    ch.fitsimage.redraw(whence=3)
    viewer.pending = True


def _clear_channel(ch, tagpfx=TAG_PREFIX):
    """Delete the PypeIt objects with a given tag prefix from a channel canvas."""
    canvas = _canvas(ch)
    tags = canvas.get_tags_by_tag_pfx(tagpfx)
    if len(tags) > 0:
        canvas.delete_objects_by_tag(tags, redraw=True)


def show_image(inp, chname='Image', waveimg=None, mask=None, exten=0, cuts=None, clear=False,
               wcs_match=False, viewer=None):
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
        viewer (:class:`~pypeit.display.viewer.GingaShellViewer`, optional):
            Viewer to use.  If None, use PypeIt's in-process viewer (see
            :func:`~pypeit.display.viewer.get_viewer`).

    Returns:
        :obj:`tuple`: The viewer (:class:`~pypeit.display.viewer.PypeItViewer`)
        and the ginga channel (``ginga.rv.Channel.Channel``) with the displayed
        image.

    Raises:
        ValueError:
            Raised if `cuts` is provided and does not have two elements.
    """
    # Input checks
    if cuts is not None and len(cuts) != 2:
        raise ValueError('Input cuts must only have two elements, the lower and upper cut.')

    if viewer is None:
        viewer = get_viewer()
    # Read or set the image data.  This will fail if the input is a
    # string and astropy.io.fits cannot read the image.
    img = io.fits_open(inp)[exten].data if isinstance(inp, (str, Path)) else inp
    img = np.asarray(img)
    if img.dtype == bool:
        img = img.astype(np.uint8)

    if clear:
        viewer.clear_all()

    # NOTE: SlitImage requires the ginga widget toolkit to be chosen, which
    # is done when the viewer is created.
    from pypeit.display.slitwavelength import SlitImage
    ch = viewer.channel(chname)
    image = SlitImage(wav_np=waveimg, logger=viewer.shell.logger)
    image.load_data(img)
    image.set(name=chname, path=None)
    # These commands set up the viewer. They can be found at
    # ginga/ginga/ImageView.py
    fi = ch.fitsimage
    with fi.suppress_redraw:
        viewer.shell.add_image(chname, image, chname=chname)
        _clear_channel(ch)
        fi.set_color_map('ramp')
        fi.set_intensity_map('ramp')
        fi.set_color_algorithm('linear')
        fi.restore_contrast()
        fi.restore_cmap()
        if cuts is not None:
            fi.cut_levels(float(cuts[0]), float(cuts[1]))

    # WCS Match this to other images with this as the reference image?
    if wcs_match:
        viewer.set_wcs_match(chname)

    # If bitmask was passed in, assume this is an extraction qa image
    # and use the mask to identify why each pixel was masked
    if mask is not None:
        show_mask_flags(viewer, ch, mask)

    viewer.pending = True
    return viewer, ch


def show_mask_flags(viewer, ch, mask):
    """
    Overlay the pixels flagged by an image bitmask, colored by the reason.

    Pixels flagged as bad pixels are magenta, cosmic rays are cyan, pixels
    masked during extraction are red, and pixels masked for any other reason
    are yellow.  Pixels off all slits are not marked.

    Args:
        viewer (:class:`~pypeit.display.viewer.PypeItViewer`):
            PypeIt ginga viewer
        ch (``ginga.rv.Channel.Channel``):
            Ginga channel
        mask (:class:`~pypeit.images.ImageBitMaskArray`):
            Bitmask array with the reasons pixels are masked.
    """
    from ginga.RGBImage import RGBImage
    dc = viewer.canvas_types()

    # Select pixels on any slit
    onslit = mask.flagged('OFFSLITS', invert=True)
    other_flags = list(mask.bit_keys())
    for flag in ['OFFSLITS', 'BPM', 'CR', 'EXTRACT']:
        other_flags.remove(flag)

    # Later flags take precedence where pixels have more than one flag
    rgba = np.zeros(onslit.shape + (4,), dtype=np.uint8)
    layers = [(mask.flagged(flag=other_flags), (255, 255, 0), 'OTHER', 'yellow'),
              (mask.extract, (255, 0, 0), 'EXTRACT', 'red'),
              (mask.cr, (0, 255, 255), 'CR', 'cyan'),
              (mask.bpm, (255, 0, 255), 'BPM', 'magenta')]
    for flagged, rgb, _, _ in layers:
        rgba[flagged & onslit] = rgb + (255,)
    _add(viewer, ch, dc.Image(0, 0, RGBImage(data_np=rgba), flipy=False))

    # Legend
    labels = [dc.Text(10, 20 + 20*i, label, color=color, fontsize=14, coord='window',
                      bgcolor='black', bgalpha=0.6)
              for i, (_, _, label, color) in enumerate(layers[::-1])]
    _add(viewer, ch, dc.CompoundObject(*labels))


def show_points(viewer, ch, spec, spat, color='cyan', legend=None, legend_spec=None, legend_spat=None):
    """
    Plot points in a ginga viewer

    Parameters
    ----------
    viewer : :class:`~pypeit.display.viewer.PypeItViewer`
        PypeIt ginga viewer
    ch : ``ginga.rv.Channel.Channel``
        Ginga channel
    spec : list
        List of spectral positions on image to plot
    spat : list
        List of spatial positions on image to plot
    color : str
        Color for points
    legend : str
        Label for a legend
    legend_spec : float
        Spectral pixel location for legend
    legend_spat : float
        Pixel location for legend

    """
    dc = viewer.canvas_types()
    objs = [dc.Point(float(_spat), float(_spec), 2, style='plus', color=color)
            for _spat, _spec in zip(spat, spec)]
    if legend is not None:
        spec_label = np.mean(np.array(spec)) if legend_spec is None else legend_spec
        spat_label = (np.mean(np.array(spat)) + 30) if legend_spat is None else legend_spat
        objs += [dc.Text(float(spat_label), float(spec_label), legend, color=color, fontsize=20)]
    _add(viewer, ch, dc.CompoundObject(*objs))


def show_text(viewer, ch, spat, spec, text, color='cyan', fontsize=20):
    """
    Write text on an image in a ginga viewer.

    Parameters
    ----------
    viewer : :class:`~pypeit.display.viewer.PypeItViewer`
        PypeIt ginga viewer
    ch : ``ginga.rv.Channel.Channel``
        Ginga channel
    spat, spec : :obj:`float`
        Spatial (horizontal) and spectral (vertical) image coordinates of the
        start of the text.
    text : :obj:`str`
        Text to write.
    color : :obj:`str`, optional
        Text color.
    fontsize : :obj:`float`, optional
        Font size.
    """
    dc = viewer.canvas_types()
    _add(viewer, ch, dc.Text(float(spat), float(spec), text, color=color, fontsize=fontsize))


# TODO: Should we continue to allow rotate as an option?
def show_slits(viewer, ch, left, right, slit_ids=None, left_ids=None, right_ids=None, maskdef_ids=None, spec_vals=None,
               rotate=False, pstep=50, clear=False, synced=True, tagpfx=TAG_PREFIX):
    r"""
    Overplot slits on the image in Ginga in the given channel

    Args:
        viewer (:class:`~pypeit.display.viewer.PypeItViewer`):
            PypeIt ginga viewer
        ch (``ginga.rv.Channel.Channel``):
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
        pstep (:obj:`int`, optional):
            Minimum step between the edge points that are drawn.  At most
            100 points are drawn for the part of each edge in view.
        clear (:obj:`bool`, optional):
            Clear the canvas?
        synced (:obj:`bool`, optional):
            Flag the left and right traces are synced into slits.
            Otherwise, the edges are treated separately. If True, the
            number of left and right edges must be the same and
            ``left_ids`` and ``right_ids`` are ignored.
        tagpfx (:obj:`str`, optional):
            Prefix for the tags of the canvas objects; see
            :func:`clear_overlays`.  Must start with ``'pypeit'``.
    """
    from pypeit.display.overlays import TraceOverlay

    # Setup the trace data and IDs
    _left = left.reshape(-1,1) if left.ndim == 1 else left
    nleft = _left.shape[1]

    _right = right.reshape(-1,1) if right.ndim == 1 else right
    nright = _right.shape[1]

    nspec = _left.shape[0]
    if _right.shape[0] != nspec:
        # TODO: Any reason to remove this restriction?
        raise PypeItError('Input left and right edges have different spectral lengths.')

    # Spectral pixel location
    if spec_vals is not None:
        y = spec_vals.reshape(-1,1) if spec_vals.ndim == 1 else spec_vals
    else:
        y = np.arange(nspec, dtype=float)[:, np.newaxis]

    # Check input
    if synced:
        if left.shape != right.shape:
            raise PypeItError('Input left and right traces must have the same shape if they have been '
                       'synchronized into slits.')
        if left_ids is not None or right_ids is not None:
            log.warning('For showing synced edges, left and right ID numbers are ignored.')
        nslits = _left.shape[1]
        _slit_ids = np.arange(nslits) if slit_ids is None else np.atleast_1d(slit_ids)
        if len(_slit_ids) != nslits:
            raise PypeItError('Incorrect number of slit IDs provided.')
        _maskdef_ids = np.atleast_1d(maskdef_ids) \
                            if maskdef_ids is not None and maskdef_ids.size == nslits else None
    else:
        _left_ids = -np.arange(nleft) if left_ids is None else np.atleast_1d(left_ids)
        if len(_left_ids) != nleft:
            raise PypeItError('Incorrect number of left IDs provided.')
        _right_ids = -np.arange(nright) if right_ids is None else np.atleast_1d(right_ids)
        if len(_right_ids) != nright:
            raise PypeItError('Incorrect number of right IDs provided.')

    if clear:
        _clear_channel(ch)

    def _overlay(x, **kwargs):
        # Swap the image axes if requested
        _x, _y = np.broadcast_arrays(x, y)
        if rotate:
            _x, _y = _y, _x
        return TraceOverlay(_x, y=_y, min_stride=pstep, rot_deg=0. if rotate else 90.,
                            what='slit edges', owner=ch.fitsimage, **kwargs)

    # Edges
    _add(viewer, ch, tagpfx=tagpfx, obj=_overlay(_left, colors='green', fontsize=20.,
                              labels=None if synced else [f'S{i}' for i in _left_ids],
                              label_colors='aquamarine', label_offset=(12., 0.)))
    _add(viewer, ch, tagpfx=tagpfx, obj=_overlay(_right, colors='magenta', fontsize=20.,
                              labels=None if synced else [f'S{i}' for i in _right_ids],
                              label_offset=(-12., 0.)))

    # Plot slit labels, if synced
    if synced:
        _slit_id_loc = _left + 0.45*(_right - _left)
        # NOTE: The labels are rotated by 90 deg, so they extend upward from
        # their anchor.  The slit and slitmask IDs are placed below the center
        # of the visible part of the slit, where they do not overlap the
        # object names drawn by show_trace.
        _add(viewer, ch, tagpfx=tagpfx, obj=_overlay(_slit_id_loc, draw_lines=False, fontsize=20.,
                                  labels=[f'S{i}' for i in _slit_ids], label_colors='aquamarine',
                                  label_offset=(0., -100.)))
        if _maskdef_ids is not None:
            _add(viewer, ch, tagpfx=tagpfx, obj=_overlay(_slit_id_loc, draw_lines=False, fontsize=20.,
                                      labels=[f'{i}' for i in _maskdef_ids], label_colors='cyan',
                                      label_offset=(0., -220.)))


def show_trace(viewer, ch, trace, trc_name=None, maskdef_extr=None, manual_extr=None, clear=False,
               rotate=False, pstep=3, yval=None, color='blue', tagpfx=TAG_PREFIX):
    r"""
    Overplot object traces on the image in Ginga in the given channel.

    For more than 200 traces (e.g., fiber-fed IFUs), the traces are only drawn
    when the view is zoomed in enough to distinguish adjacent traces.

    Args:
        viewer (:class:`~pypeit.display.viewer.PypeItViewer`):
            PypeIt ginga viewer.
        ch (``ginga.rv.Channel.Channel``):
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
        pstep (:obj:`int`, optional):
            Minimum step between the trace points that are drawn.  At most
            100 points are drawn for the part of each trace in view.
        yval (`numpy.ndarray`_, optional):
            Array with spectral position of the object traces. Shape must be :math:`(N_{\rm spec},)`
            or :math:`(N_{\rm spec}, N_{\rm trace})`. If not passed in, the default of
            np.arange(:math:`(N_{\rm spec},)`) will be used.
        color (str, optional):
            If None, the traces are not drawn (only their names).  Otherwise,
            traces are orange, or yellow (maskdef extractions) or light blue
            (manual extractions).
        tagpfx (:obj:`str`, optional):
            Prefix for the tags of the canvas objects; see
            :func:`clear_overlays`.  Must start with ``'pypeit'``.

    """
    from pypeit.display.overlays import TraceOverlay

    if clear:
        _clear_channel(ch)

    if trace.ndim == 1:
        trace = trace.reshape(-1,1)

    ntrace = trace.shape[1]
    _maskdef_extr = ntrace*[False] if maskdef_extr is None else maskdef_extr
    _manual_extr = ntrace*[False] if manual_extr is None else manual_extr
    _trc_name = ntrace*[''] if trc_name is None else trc_name

    # Show
    if yval is None:
        y = np.arange(trace.shape[0], dtype=float)[:, None]
    else:
        y = yval.reshape(-1, 1) if yval.ndim == 1 else yval
    x, y = np.broadcast_arrays(trace, y)
    if rotate:
        x, y = y, x

    colors = []
    for i in range(ntrace):
        if _maskdef_extr[i]:
            _color = '#f0e442' if color is not None else color
        elif _manual_extr[i]:
            _color = '#33ccff' if color is not None else color
        else:
            _color = 'orange' if color is not None else color
        colors.append(_color)
    label_colors = [c if c is not None else 'orange' for c in colors]

    overlay = TraceOverlay(x, y=y, colors=colors, labels=_trc_name, label_colors=label_colors,
                           fontsize=17., rot_deg=0. if rotate else 90., min_stride=pstep,
                           what='object traces', note_index=1, owner=ch.fitsimage)
    _add(viewer, ch, overlay, tagpfx=tagpfx)
    if not rotate:
        overlay.enable_hover(ch.fitsimage)


def clear_overlays(ch, tagpfx):
    """
    Remove overlays added with a given tag prefix (e.g., by
    :func:`show_slits` or :func:`show_trace`) from a channel.

    Args:
        ch (``ginga.rv.Channel.Channel``):
            Ginga channel
        tagpfx (:obj:`str`):
            Tag prefix of the overlays to remove.
    """
    _clear_channel(ch, tagpfx=tagpfx)


def clear_canvas(cname):
    """
    Clear the ginga canvas

    Args:
        cname (str):  Channel name

    """
    viewer = get_viewer()
    if cname in viewer.channel_names():
        _clear_channel(viewer.channel(cname))


def clear_all(allow_new=False):
    """
    Clear all of the ginga canvasses.

    Args:
        allow_new (:obj:`bool`, optional):
            Ignored; kept for backwards compatibility.
    """
    get_viewer().clear_all()


def show_tilts(viewer, ch, tilt_traces, yoff=0., xoff=0., points=True, nspec=None, pstep=3, clear_canvas=False):
    """
    Show the arc tilts on the input channel

    Args:
        viewer (:class:`~pypeit.display.viewer.PypeItViewer`):
            PypeIt ginga viewer
        ch (``ginga.rv.Channel.Channel``):
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
    from pypeit.display.overlays import MarkerOverlay

    if tilt_traces is None:
        raise PypeItError('No tilts have been traced or fitted')

    if clear_canvas:
        _clear_channel(ch)

    dc = viewer.canvas_types()

    # Plot traced tilts
    if 'goodpix_tilt' in tilt_traces.keys() and tilt_traces['goodpix_tilt'][0].size > 0 and points:
        _add(viewer, ch, MarkerOverlay(tilt_traces['goodpix_spat'][0] + xoff,
                                       tilt_traces['goodpix_tilt'][0] + yoff, color='cyan'))

    # Plot the 2D fitted tilts
    # loop over each line, this allows to use type='path' and therefore a faster plotting
    if 'good2dfit_lid' in tilt_traces.keys():
        objs = []
        for iline in np.unique(tilt_traces['good2dfit_lid'][0]):
            # good fit
            this_line = tilt_traces['good2dfit_lid'][0] == iline
            if np.any(this_line):
                good2dfit_spat = tilt_traces['good2dfit_spat'][0][this_line] + xoff
                good2dfit_tilt = tilt_traces['good2dfit_tilt'][0][this_line] + yoff
                objs.append(dc.Path(np.column_stack([good2dfit_spat[::pstep],
                                                     good2dfit_tilt[::pstep]]),
                                    color='blue', linewidth=1))
        if len(objs) > 0:
            _add(viewer, ch, dc.CompoundObject(*objs))

    # Now plot the masked traces and the rejected 2D fits
    # masked traces
    if 'badpix_tilt' in tilt_traces.keys() and tilt_traces['badpix_tilt'][0].size > 0:
        _add(viewer, ch, MarkerOverlay(tilt_traces['badpix_spat'][0] + xoff,
                                       tilt_traces['badpix_tilt'][0] + yoff, color='red'))
    # rejected fit
    if 'bad2dfit_tilt' in tilt_traces.keys() and tilt_traces['bad2dfit_tilt'][0].size > 0:
        _add(viewer, ch, MarkerOverlay(tilt_traces['bad2dfit_spat'][0] + xoff,
                                       tilt_traces['bad2dfit_tilt'][0] + yoff, color='yellow'))

    # Add text
    text_str = ['Masked pixel', 'Rejected in fit', 'Good tilt fit']
    text_color = ['red', 'yellow', 'blue']
    if points:
        text_str += ['Good pixel']
        text_color += ['cyan']
    _add(viewer, ch, dc.CompoundObject(*[dc.Text(10, 20 + 20*i, s, color=c, fontsize=14,
                                                 coord='window', bgcolor='black', bgalpha=0.6)
                                         for i, (s, c) in enumerate(zip(text_str, text_color))]))


def show_scattered_light(image_list, slits=None, wcs_match=True):
    """
    Interface to ginga to show the quality of the Scattered Light subtraction

    Parameters
    ----------
    image_list : zip
        A zip of the images to show, their names, and the scales
    slits : :class:`~pypeit.slittrace.SlitTraceSet`, optional
        The current slit traces
    wcs_match : :obj:`bool`, optional
        Use a reference image for the WCS and match all image in other channels to it.
    """
    if slits is not None:
        left, right, mask = slits.select_edges()
        gpm = mask == 0
    # Loop me
    clear = True
    for img, name, cut in image_list:
        if img is None:
            continue
        viewer, ch = show_image(img, chname=name, cuts=cut, wcs_match=wcs_match, clear=clear)
        if slits is not None:
            show_slits(viewer, ch, left[:, gpm], right[:, gpm], slit_ids=slits.spat_id[gpm])
        # Turn off clear
        if clear:
            clear = False
    show()


def show_1dspec(filename, ext=0, masked=True, fluxed=False, extraction='OPT'):
    """
    Interface to ginga to show a 1dspec and manipulate with Spec1dView plugin.

    Parameters
    ----------
    filename : str
        spec1d FITS file to show in the viewer
    ext : int
        extension to show (which spectrum)
    """
    viewer = get_viewer()
    shell = viewer.shell

    chname, plname = "Spec1d", "Spec1dView"
    viewer.channel(chname)
    # start the plugin
    shell.start_local_plugin(chname, plname)
    viewer._update(0.2)
    # set up the options as passed
    shell.call_local_plugin_method(chname, plname, 'set_params', [],
                                   dict(ext=int(ext), extraction=extraction, masked=masked,
                                        fluxed=fluxed))
    # load the file
    shell.load_file(str(filename), chname=chname, wait=False)
    viewer.pending = True
