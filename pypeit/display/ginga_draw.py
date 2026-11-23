"""
In-process drawing helpers for the PypeIt ginga plugins.

The functions in this module run *inside* the ginga reference viewer process.
They build native ginga data and canvas objects, as opposed to the functions
in :mod:`pypeit.display.display`, which send data to a viewer over the remote
control (RC) interface.

Functions that build data objects (:func:`make_image`, :func:`slit_overlay`,
:func:`trace_overlay`) may be called from any thread.  Functions that change
the viewer (:func:`show_image`, :func:`add_overlay`, :func:`clear_channels`)
must be called on the ginga GUI thread, e.g., via ``fv.gui_call``.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np

from ginga.AstroImage import AstroImage
from ginga.canvas.CanvasObject import get_canvas_types

from pypeit import PypeItError


class SlitImage(AstroImage):
    """
    Custom image type for a 2D slit image.

    This is a child of ginga's `AstroImage`_ and primarily enables the display
    of the slit wavelengths (and, optionally, the RMS of the wavelength
    solution) in the ginga cursor readout.

    Parameters
    ----------
    wav_np : `numpy.ndarray`_, optional
        Wavelength map image.
    rms_np : `numpy.ndarray`_, optional
        Per-pixel image of the wavelength-fit RMS.  Each slit has a single
        RMS, so every pixel within a given slit holds the same value.  Must
        have the same shape as ``wav_np`` when provided.  Pixels outside all
        slits should be 0; those values are suppressed in the readout.
    **kwargs:
        Keyword arguments for `AstroImage`_
    """
    def __init__(self, wav_np=None, rms_np=None, **kwargs):
        super().__init__(**kwargs)
        self.wav_np = wav_np
        self.rms_np = rms_np

    def info_xy(self, data_x, data_y, settings):
        """
        Return the cursor readout information for an image coordinate.

        Parameters
        ----------
        data_x : float
            Image x coordinate
        data_y : float
            Image y coordinate
        settings : ``ginga.misc.Settings.SettingGroup``
            Ginga settings group

        Returns
        -------
        info : ``ginga.misc.Bunch.Bunch``
            Metadata for this coordinate
        """
        info = super().info_xy(data_x, data_y, settings)
        if self.wav_np is None:
            # Allow the wavelengths to be undefined and just return the
            # same result as would be provided by AstroImage.
            return info

        # We report the value across the pixel, even though the coords
        # change halfway across the pixel
        _d_x = int(np.floor(data_x + 0.5))
        _d_y = int(np.floor(data_y + 0.5))

        # Convert wavelength value to a string. Return an empty
        # string if the provided coordinate is off the image
        _ht, _wd = self.wav_np.shape
        in_bounds = 0 <= _d_y < _ht and 0 <= _d_x < _wd
        wav_s = '{:<14.6g}'.format(self.wav_np[_d_y, _d_x]) if in_bounds else ''

        # Build the RMS readout when an RMS image has been supplied and the
        # cursor is inside a slit (rms value > 0 and finite).
        rms_lbl = ''
        rms_s = ''
        if self.rms_np is not None and in_bounds:
            rms_val = self.rms_np[_d_y, _d_x]
            if np.isfinite(rms_val) and rms_val > 0:
                rms_lbl = 'RMS'
                rms_s = '{:<10.4g}'.format(rms_val)

        info.update(dict(ra_lbl="λ", ra_txt=wav_s, dec_lbl=rms_lbl, dec_txt=rms_s))
        return info


def make_image(img, name, waveimg=None, rms_img=None, logger=None):
    """
    Construct a ginga image object from a numpy array.

    This function does not touch the viewer and can be called from any thread.

    Parameters
    ----------
    img : `numpy.ndarray`_
        2D image to display.
    name : :obj:`str`
        Name for the image in ginga.  Must be unique within the channel the
        image is added to.
    waveimg : `numpy.ndarray`_, optional
        Wavelength image with the same shape as ``img``.  If provided, the
        returned object is a :class:`SlitImage` that shows the wavelength in
        the cursor readout.
    rms_img : `numpy.ndarray`_, optional
        Per-pixel image of the wavelength-fit RMS, with the same shape as
        ``img``.  Ignored if ``waveimg`` is None.
    logger : `logging.Logger`_, optional
        Logger passed to the ginga image object.

    Returns
    -------
    `AstroImage`_
        The image object; a :class:`SlitImage` if ``waveimg`` is provided.
    """
    if waveimg is None:
        image = AstroImage(logger=logger)
    else:
        if waveimg.shape != img.shape:
            raise PypeItError(
                f'Wavelength image shape {waveimg.shape} does not match image shape '
                f'{img.shape}.'
            )
        image = SlitImage(wav_np=waveimg, rms_np=rms_img, logger=logger)
    image.load_data(img)
    image.set(name=name, path=None, nothumb=True)
    return image


def show_image(fv, chname, image, cuts=None):
    """
    Add an image to a channel and set the PypeIt display defaults.

    The channel is created if it does not exist.  Must be called on the ginga
    GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    chname : :obj:`str`
        Name of the channel.
    image : `AstroImage`_
        Image to show, typically constructed with :func:`make_image`.
    cuts : array-like, optional
        Lower and upper cut levels.  If None, the viewer's automatic cut
        levels are used.

    Returns
    -------
    channel : ``ginga.rv.Channel.Channel``
        The channel showing the image.
    """
    if cuts is not None and len(cuts) != 2:
        raise PypeItError('Input cuts must only have two elements, the lower and upper cut.')
    channel = fv.get_channel_on_demand(chname)
    fv.add_image(image.get('name'), image, chname=channel.name)
    viewer = channel.fitsimage
    viewer.set_color_map('ramp')
    viewer.set_intensity_map('ramp')
    viewer.set_color_algorithm('linear')
    viewer.restore_contrast()
    viewer.restore_cmap()
    if cuts is not None:
        viewer.cut_levels(float(cuts[0]), float(cuts[1]))
    return channel


def _edge_points(spat, spec, pstep, rotate):
    """
    Return a decimated list of (x, y) display coordinates for one trace.

    Parameters
    ----------
    spat : `numpy.ndarray`_
        Spatial positions along the trace.
    spec : `numpy.ndarray`_
        Spectral positions along the trace.
    pstep : :obj:`int`
        Keep every ``pstep``-th point.
    rotate : :obj:`bool`
        If True, the spectral direction is along x.

    Returns
    -------
    :obj:`list`
        List of (x, y) tuples.
    """
    x, y = (spec, spat) if rotate else (spat, spec)
    return list(zip(x[::pstep].tolist(), y[::pstep].tolist()))


def _spec_coords(spec_vals, nspec, nedge):
    r"""
    Return the spectral coordinates of a set of edges or traces.

    Parameters
    ----------
    spec_vals : `numpy.ndarray`_
        Spectral coordinates, with shape :math:`(N_{\rm spec},)` or
        :math:`(N_{\rm spec}, N_{\rm edge})`.  If None, the spectral pixel
        index is used.
    nspec : :obj:`int`
        Number of spectral positions.
    nedge : :obj:`int`
        Number of edges or traces.

    Returns
    -------
    `numpy.ndarray`_
        Spectral coordinates with shape :math:`(N_{\rm spec}, N_{\rm edge})`.
    """
    y = np.arange(nspec, dtype=float) if spec_vals is None else np.asarray(spec_vals, dtype=float)
    if y.ndim == 1:
        return np.repeat(y[:,np.newaxis], nedge, axis=1)
    if y.shape != (nspec, nedge):
        raise PypeItError(
            f'Spectral coordinates have shape {y.shape}; expected {(nspec, nedge)}.'
        )
    return y


def _label_xy(spat, spec, rotate):
    """
    Return the display (x, y) of a label at a spatial and spectral position.

    Parameters
    ----------
    spat : :obj:`float`
        Spatial position.
    spec : :obj:`float`
        Spectral position.
    rotate : :obj:`bool`
        If True, the spectral direction is along x.

    Returns
    -------
    :obj:`tuple`
        The (x, y) display coordinates.
    """
    return (float(spec), float(spat)) if rotate else (float(spat), float(spec))


def slit_overlay(
    left, right, slit_ids=None, left_ids=None, right_ids=None, maskdef_ids=None, spec_vals=None,
    rotate=False, pstep=50, synced=True
):
    r"""
    Construct a ginga canvas object with slit edges and their labels.

    This is the in-process equivalent of :func:`pypeit.display.display.show_slits`.
    Left edges are green and right edges magenta; slit IDs are aquamarine and
    slit-mask design IDs cyan.  Can be called from any thread.

    Parameters
    ----------
    left : `numpy.ndarray`_
        Spatial position of the left slit edges.  Shape must be
        :math:`(N_{\rm spec},)` or :math:`(N_{\rm spec}, N_{\rm l-edge})`.
    right : `numpy.ndarray`_
        Spatial position of the right slit edges.  Shape must be
        :math:`(N_{\rm spec},)` or :math:`(N_{\rm spec}, N_{\rm r-edge})`, and
        must match ``left`` if ``synced`` is True.
    slit_ids : array-like, optional
        IDs for the slits, used if ``synced`` is True.  If None, IDs run from
        0 to :math:`N_{\rm slits}-1`.
    left_ids : array-like, optional
        IDs for the left edges, used if ``synced`` is False.  If None, IDs run
        from 0 to :math:`-(N_{\rm l-edge}-1)`.
    right_ids : array-like, optional
        IDs for the right edges, used if ``synced`` is False.  If None, IDs run
        from 0 to :math:`-(N_{\rm r-edge}-1)`.
    maskdef_ids : array-like, optional
        Slit-mask design IDs for the slits, used if ``synced`` is True.  Not
        shown if None or if the length does not match the number of slits.
    spec_vals : `numpy.ndarray`_, optional
        Spectral positions of the edges, with shape :math:`(N_{\rm spec},)`;
        or, if ``left`` and ``right`` have the same shape, the shape of the
        edge arrays.  If None, the spectral pixel index is used.
    rotate : :obj:`bool`, optional
        If True, the spectral direction is along x.
    pstep : :obj:`int`, optional
        Keep every ``pstep``-th point along the edges.
    synced : :obj:`bool`, optional
        If True, the left and right edges are synchronized into slits and
        labeled by slit; otherwise each edge is labeled separately.

    Returns
    -------
    ``ginga.canvas.types.layer.Canvas``
        Canvas holding all edges and labels.  The same canvas can be added
        to the canvases of several channels.
    """
    _left = left.reshape(-1,1) if left.ndim == 1 else left
    _right = right.reshape(-1,1) if right.ndim == 1 else right
    nspec, nleft = _left.shape
    nright = _right.shape[1]
    if _right.shape[0] != nspec:
        raise PypeItError('Input left and right edges have different spectral lengths.')

    if synced and _left.shape != _right.shape:
        raise PypeItError(
            'Input left and right traces must have the same shape if they have been '
            'synchronized into slits.'
        )
    y_left = _spec_coords(spec_vals, nspec, nleft)
    y = _spec_coords(spec_vals, nspec, nright)

    dc = get_canvas_types()
    objs = []
    for i in range(nleft):
        objs.append(
            dc.Path(_edge_points(_left[:,i], y_left[:,i], pstep, rotate), color='green')
        )
    for i in range(nright):
        objs.append(dc.Path(_edge_points(_right[:,i], y[:,i], pstep, rotate), color='magenta'))

    bot = nspec // 2
    if synced:
        _slit_ids = np.arange(nleft) if slit_ids is None else np.atleast_1d(slit_ids)
        if len(_slit_ids) != nleft:
            raise PypeItError('Incorrect number of slit IDs provided.')
        _maskdef_ids = None
        if maskdef_ids is not None and np.atleast_1d(maskdef_ids).size == nleft:
            _maskdef_ids = np.atleast_1d(maskdef_ids)
        id_loc = _left + 0.45*(_right - _left)
        for i in range(nleft):
            x, yy = _label_xy(id_loc[bot,i], y[bot,i] - 400, rotate)
            objs.append(
                dc.Text(x, yy, text=f'S{_slit_ids[i]}', color='aquamarine', fontsize=20.,
                        rot_deg=90.)
            )
            if _maskdef_ids is not None:
                x, yy = _label_xy(id_loc[bot,i], y[bot,i] - 250, rotate)
                objs.append(
                    dc.Text(x, yy, text=f'{_maskdef_ids[i]}', color='cyan', fontsize=20.,
                            rot_deg=90.)
                )
    else:
        _left_ids = -np.arange(nleft) if left_ids is None else np.atleast_1d(left_ids)
        if len(_left_ids) != nleft:
            raise PypeItError('Incorrect number of left IDs provided.')
        _right_ids = -np.arange(nright) if right_ids is None else np.atleast_1d(right_ids)
        if len(_right_ids) != nright:
            raise PypeItError('Incorrect number of right IDs provided.')
        for i in range(nleft):
            x, yy = _label_xy(_left[bot,i]*1.05, y_left[bot,i], rotate)
            objs.append(
                dc.Text(x, yy, text=f'S{_left_ids[i]}', color='aquamarine', fontsize=20.,
                        rot_deg=90.)
            )
        for i in range(nright):
            x, yy = _label_xy(_right[bot,i]*(1-0.05), y[bot,i], rotate)
            objs.append(
                dc.Text(x, yy, text=f'S{_right_ids[i]}', color='magenta', fontsize=20.,
                        rot_deg=90.)
            )

    return dc.Canvas(*objs)


def trace_overlay(
    trace, trc_name=None, maskdef_extr=None, manual_extr=None, rotate=False, pstep=3, yval=None
):
    r"""
    Construct a ginga canvas object with object traces and their names.

    This is the in-process equivalent of :func:`pypeit.display.display.show_trace`.
    Traces are orange, except those extracted at slit-mask design positions
    (yellow) or extracted manually (light blue).  Can be called from any
    thread.

    Parameters
    ----------
    trace : `numpy.ndarray`_
        Spatial position of the object traces.  Shape must be
        :math:`(N_{\rm spec},)` or :math:`(N_{\rm spec}, N_{\rm trace})`.
    trc_name : array-like, optional
        Names of the traces, shape :math:`(N_{\rm trace},)`.
    maskdef_extr : array-like, optional
        Boolean flags selecting traces extracted at slit-mask design
        positions, shape :math:`(N_{\rm trace},)`.
    manual_extr : array-like, optional
        Boolean flags selecting manually extracted traces, shape
        :math:`(N_{\rm trace},)`.
    rotate : :obj:`bool`, optional
        If True, the spectral direction is along x.
    pstep : :obj:`int`, optional
        Keep every ``pstep``-th point along the traces.
    yval : `numpy.ndarray`_, optional
        Spectral positions of the traces, with shape :math:`(N_{\rm spec},)`
        or the same shape as ``trace``.  If None, the spectral pixel index is
        used.

    Returns
    -------
    ``ginga.canvas.types.layer.Canvas``
        Canvas holding all traces and labels.  The same canvas can be added
        to the canvases of several channels.
    """
    _trace = trace.reshape(-1,1) if trace.ndim == 1 else trace
    nspec, ntrace = _trace.shape
    _maskdef_extr = np.zeros(ntrace, dtype=bool) if maskdef_extr is None \
                        else np.asarray(maskdef_extr, dtype=bool)
    _manual_extr = np.zeros(ntrace, dtype=bool) if manual_extr is None \
                        else np.asarray(manual_extr, dtype=bool)
    _trc_name = ntrace*[''] if trc_name is None else list(trc_name)
    y = _spec_coords(yval, nspec, ntrace)

    dc = get_canvas_types()
    objs = []
    half = nspec // 2
    for i in range(ntrace):
        if _maskdef_extr[i]:
            color = '#f0e442'
        elif _manual_extr[i]:
            color = '#33ccff'
        else:
            color = 'orange'
        objs.append(dc.Path(_edge_points(_trace[:,i], y[:,i], pstep, rotate), color=color))
        x, yy = _label_xy(_trace[half,i], y[half,i], rotate)
        objs.append(
            dc.Text(x, yy, text=str(_trc_name[i]), color=color, fontsize=17., rot_deg=90.)
        )
    return dc.Canvas(*objs)


def add_overlay(channel, overlay, tag):
    """
    Add an overlay to a channel's canvas, replacing any with the same tag.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    channel : ``ginga.rv.Channel.Channel``
        Channel to draw on.
    overlay : ``ginga.canvas.CanvasObject.CanvasObjectBase``
        Object to add, e.g. from :func:`slit_overlay` or :func:`trace_overlay`.
    tag : :obj:`str`
        Tag for the overlay, used to replace or remove it later.
    """
    canvas = channel.fitsimage.get_canvas()
    canvas.delete_object_by_tag(tag, redraw=False)
    canvas.add(overlay, tag=tag, redraw=True)


def clear_channels(fv, chnames):
    """
    Delete the named channels that exist in the viewer.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    chnames : :obj:`list`
        Names of the channels to delete.  Names of channels that do not exist
        are ignored.
    """
    for chname in chnames:
        if fv.has_channel(chname):
            fv.delete_channel(chname)
