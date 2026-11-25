"""
Command handlers run by the :class:`~pypeit.display.pypeit_remote.PypeItRemote`
ginga plugin.

Each handler runs one PypeIt display command inside the ginga process.  It is
called on a ginga worker thread as ``handler(fv, args, job)``, where ``fv`` is
the ginga reference viewer shell, ``args`` is the `argparse.Namespace`_ parsed
by the script that sent the command, and ``job`` is the
:class:`~pypeit.display.remote_jobs.RemoteJob` running it.

A handler must:

    - do file I/O and computation directly, on the worker thread;
    - make every change to the viewer through ``fv.gui_call``;
    - call ``job.check_cancel()`` between steps so that a newer command can
      supersede it;
    - raise :class:`~pypeit.PypeItError` for errors the user should see in the
      terminal.

Warnings logged with the PypeIt logger by the worker thread are reported to
the script that sent the command.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path

import numpy as np

from astropy.io import fits

from pypeit import edgetrace
from pypeit import flatfield
from pypeit import scattlight
from pypeit import log
from pypeit import PypeItError
from pypeit import slittrace
from pypeit import specobjs
from pypeit import wavetilts
from pypeit.images import buildimage
from pypeit.images import pypeitimage
from pypeit.display import ginga_draw
from pypeit.display import image_data
from pypeit.display import spec2d_data


def get_handlers():
    """
    Return the command handlers known to the ``PypeItRemote`` plugin.

    The command names match the names used by the scripts that send them
    (e.g., ``'show_2dspec'`` for ``pypeit_show_2dspec``).

    Returns
    -------
    :obj:`dict`
        Handler functions keyed by command name.
    """
    return {
        'chk_edges': chk_edges, 'chk_flats': chk_flats, 'show_1dspec': show_1dspec,
        'show_2dspec': show_2dspec, 'show_calib': show_calib, 'chk_tilts': chk_tilts,
        'view_fits': view_fits, 'show_pixflat': show_pixflat, 'chk_scattlight': chk_scattlight
    }


def clear_workspace(fv):
    """
    Delete all channels in the viewer.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    """
    ginga_draw.clear_channels(fv, list(fv.get_channel_names()))


def spec1d_extension(sobjs, exten=1, obj=None, extract='OPT'):
    """
    Select the spectrum to show from a spec1d file.

    Parameters
    ----------
    sobjs : :class:`~pypeit.specobjs.SpecObjs`
        Spectra read from the spec1d file.
    exten : :obj:`int`, optional
        1-based number of the spectrum to show.  Ignored if ``obj`` is
        provided.
    obj : :obj:`str`, optional
        Name of the spectrum to show, e.g. ``SPAT0424-SLIT0000-DET01``.
    extract : :obj:`str`, optional
        Extraction to show, ``'OPT'`` or ``'BOX'``.

    Returns
    -------
    :obj:`int`
        0-based index of the spectrum in ``sobjs``.

    Raises
    ------
    :class:`~pypeit.PypeItError`
        Raised if the spectrum does not exist or does not have the requested
        extraction.
    """
    nspec = len(sobjs)
    if obj is None:
        indx = exten - 1
        if indx < 0 or indx >= nspec:
            raise PypeItError(f'Extension {exten} is not available; the file has {nspec} spectra.')
    else:
        match = np.where(sobjs.NAME == obj)[0]
        if len(match) == 0:
            raise PypeItError(f'No spectrum named {obj}; use --list to see the names.')
        indx = int(match[0])
    if extract == 'OPT' and sobjs[indx]['OPT_WAVE'] is None:
        raise PypeItError('Spectrum not extracted with OPT.  Try --extract BOX')
    return indx


def _start_spec1dview(fv, chname, params):
    """
    Clear the viewer and start the Spec1dView plugin in a new channel.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    chname : :obj:`str`
        Name of the channel to create.
    params : :obj:`dict`
        Keyword arguments for ``Spec1dView.set_params``.
    """
    clear_workspace(fv)
    fv.add_channel(chname)
    fv.call_local_plugin_method(chname, 'Spec1dView', 'set_params', [], params)
    fv.start_local_plugin(chname, 'Spec1dView')
    fv.change_channel(chname)


def show_1dspec(fv, args, job):
    """
    Show a spectrum from a spec1d file with the Spec1dView plugin.

    Handler for ``pypeit_show_1dspec``; see
    :class:`~pypeit.scripts.show_1dspec.Show1DSpec` for the arguments.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_show_1dspec``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    sobjs = specobjs.SpecObjs.from_fitsfile(args.file, chk_version=False)
    indx = spec1d_extension(sobjs, exten=args.exten, obj=args.obj, extract=args.extract)
    job.check_cancel()
    chname = 'Spec1d'
    params = dict(ext=indx, extraction=args.extract, masked=args.masked, fluxed=args.flux)
    fv.gui_call(_start_spec1dview, fv, chname, params)
    job.check_cancel()
    fv.load_file(args.file, chname=chname, wait=True)


def find_slits_file(calib_key, calib_dir, directory):
    """
    Find the slits file associated with a calibration.

    The default slits file name is looked for first in the calibration
    directory recorded in the calibration file, and then in ``directory``,
    which allows for calibrations that were moved after they were made.

    Parameters
    ----------
    calib_key : :obj:`str`
        Calibration key, e.g. ``A_0_DET01``.  If None, no file is looked for.
    calib_dir : :obj:`str`
        Calibration directory recorded in the calibration file.  Can be None.
    directory : `Path`_
        Directory of the calibration file being shown.

    Returns
    -------
    `Path`_
        The slits file, or None if it was not found.  A warning is logged if
        it was not found.
    """
    if calib_key is None:
        log.warning('Calibration key unknown; cannot find the slits file.')
        return None
    directories = [directory] if calib_dir is None else [calib_dir, directory]
    candidates = [
        slittrace.SlitTraceSet.construct_file_name(calib_key, calib_dir=d) for d in directories
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    log.warning(
        'Could not find the slits file; looked for: '
        f'{", ".join(str(c) for c in candidates)}.  Slit edges will not be shown.'
    )
    return None


def read_slits(slits_file, chk_version=True):
    """
    Read a slits file to overlay on a display, if possible.

    Parameters
    ----------
    slits_file : `Path`_
        The slits file.  If None, nothing is read.
    chk_version : :obj:`bool`, optional
        Perform strict version checking when reading the file.

    Returns
    -------
    :class:`~pypeit.slittrace.SlitTraceSet`
        The slits, or None if ``slits_file`` is None or the file could not be
        read.  A warning is logged if the file could not be read.
    """
    if slits_file is None:
        return None
    try:
        return slittrace.SlitTraceSet.from_file(slits_file, chk_version=chk_version)
    except (PypeItError, OSError, TypeError, ValueError, KeyError) as e:
        # Files written by older versions of PypeIt can fail in several ways
        log.warning(f'Could not read {slits_file}; slit edges will not be shown.  ({e})')
        return None


class Overlay:
    """
    A canvas drawn on a set of channels, which a view plugin can show or hide.

    Parameters
    ----------
    name : :obj:`str`
        Name shown in the view plugin, e.g. ``'Slit edges'``.
    tag : :obj:`str`
        Tag of the canvas on each channel's canvas.
    canvas : ``ginga.canvas.types.layer.Canvas``
        The overlay; the same canvas is added to every channel that shows it.
    """
    def __init__(self, name, tag, canvas):
        self.name = name
        self.tag = tag
        self.canvas = canvas
        self.chnames = []
        self.visible = True


class ViewState:
    """
    What a command displayed, handed to the view plugin that controls it.

    Parameters
    ----------
    command : :obj:`str`
        Name of the command; used to re-run it.
    args : `argparse.Namespace`_
        Arguments of the command; used to re-run it.
    path : :obj:`str`
        File that was displayed.
    """
    def __init__(self, command, args, path):
        self.command = command
        self.args = args
        self.path = path
        self.chnames = []
        """Names of the channels shown, the first one being the reference."""
        self.overlays = []
        """List of :class:`Overlay` objects."""
        self.slits = []
        """List of dictionaries describing each slit; see :func:`slit_rows`."""
        self.options = {}
        """Additional information for the view plugin, e.g. the detectors."""


def make_overlays(slits=None, traces=None):
    """
    Build the overlays a view plugin can toggle.

    Parameters
    ----------
    slits : :obj:`dict`, optional
        Keyword arguments for :func:`~pypeit.display.ginga_draw.slit_overlay`.
    traces : :obj:`dict`, optional
        Keyword arguments for :func:`~pypeit.display.ginga_draw.trace_overlay`.

    Returns
    -------
    :obj:`dict`
        :class:`Overlay` objects keyed by ``'slit_edges'``, ``'slit_ids'``,
        ``'traces'`` and ``'trace_names'``; only those that can be built are
        included.
    """
    overlays = {}
    if slits is not None:
        overlays['slit_edges'] = Overlay(
            'Slit edges', 'pypeit_slit_edges', ginga_draw.slit_overlay(**slits, labels=False)
        )
        overlays['slit_ids'] = Overlay(
            'Slit IDs', 'pypeit_slit_ids', ginga_draw.slit_overlay(**slits, lines=False)
        )
    if traces is not None:
        overlays['traces'] = Overlay(
            'Object traces', 'pypeit_traces', ginga_draw.trace_overlay(**traces, labels=False)
        )
        overlays['trace_names'] = Overlay(
            'Object names', 'pypeit_trace_names',
            ginga_draw.trace_overlay(**traces, lines=False)
        )
    return overlays


def slit_rows(left, right, slit_ids=None, maskdef_ids=None):
    r"""
    Describe each slit for the slit navigator of a view plugin.

    Parameters
    ----------
    left : `numpy.ndarray`_
        Left edges, shape :math:`(N_{\rm spec}, N_{\rm slits})`.
    right : `numpy.ndarray`_
        Right edges, with the same shape as ``left``.
    slit_ids : array-like, optional
        Slit IDs.  If None, the slits are numbered from 0.
    maskdef_ids : array-like, optional
        Slit-mask design IDs.  Ignored if None or of the wrong length.

    Returns
    -------
    :obj:`list`
        One dictionary per slit, with the slit ``id``, ``maskdef_id`` (None if
        unknown), and the spatial center (``x``), spectral center (``y``) and
        width (``width``) at the middle of the spectral range, in pixels.
    """
    nspec, nslits = left.shape
    _ids = np.arange(nslits) if slit_ids is None else np.atleast_1d(slit_ids)
    _maskdef = None if maskdef_ids is None or np.atleast_1d(maskdef_ids).size != nslits \
                    else np.atleast_1d(maskdef_ids)
    mid = nspec // 2
    return [
        dict(id=int(_ids[i]), maskdef_id=None if _maskdef is None else int(_maskdef[i]),
             x=float(0.5*(left[mid,i] + right[mid,i])), y=float(mid),
             width=float(right[mid,i] - left[mid,i]))
        for i in range(nslits)
    ]


def _show_channel(fv, chname, image, cuts, overlays):
    """
    Show an image and its overlays in a channel.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    chname : :obj:`str`
        Name of the channel.
    image : `AstroImage`_
        Image to show.
    cuts : :obj:`tuple`
        Lower and upper cut levels, or None.
    overlays : :obj:`list`
        :class:`Overlay` objects to draw.
    """
    channel = ginga_draw.show_image(fv, chname, image, cuts=cuts)
    for overlay in overlays:
        ginga_draw.add_overlay(channel, overlay.canvas, overlay.tag)


def _link_views(fv, chnames):
    """
    Lock the pan, zoom, rotation and flips of a set of channels together.

    The images in the channels must have the same shape; the link is in
    pixel coordinates, so the images need no WCS.  The first channel is the
    reference.  Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    chnames : :obj:`list`
        Names of the channels to link.
    """
    reference = fv.get_channel(chnames[0]).fitsimage
    for chname in chnames[1:]:
        reference.copy_attributes(
            fv.get_channel(chname).fitsimage, ['transforms', 'rotation', 'zoom', 'pan'],
            share=True
        )


def display_channels(
    fv, job, channels, overlays=None, waveimg=None, rms_img=None, clear=True, link_views=True
):
    """
    Show a set of images, one per channel, with their overlays.

    Must be called on a ginga worker thread.  The images are built on the
    calling thread; only the display is done on the GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command; checked for cancellation before each
        channel.
    channels : :obj:`list`
        One dictionary per channel, with the ``image`` (a
        `numpy.ndarray`_), the channel name (``chname``), the cut levels
        (``cuts``, a tuple or None), and the keys of the ``overlays`` to draw
        on it (``overlays``, a list).
    overlays : :obj:`dict`, optional
        :class:`Overlay` objects, as returned by :func:`make_overlays`.  Each
        object's list of channel names is filled in.
    waveimg : `numpy.ndarray`_, optional
        Wavelength image shown in the cursor readout of every channel.
    rms_img : `numpy.ndarray`_, optional
        Wavelength-RMS image shown in the cursor readout of every channel.
    clear : :obj:`bool`, optional
        Delete all channels first.  If False, only the channels being shown
        are replaced.
    link_views : :obj:`bool`, optional
        Lock the pan, zoom and rotation of all channels together; see
        :func:`_link_views`.

    Returns
    -------
    :obj:`list`
        Names of the channels shown.
    """
    _overlays = {} if overlays is None else overlays
    chnames = [channel['chname'] for channel in channels]
    if clear:
        fv.gui_call(clear_workspace, fv)
    else:
        fv.gui_call(ginga_draw.clear_channels, fv, chnames)
    for channel in channels:
        job.check_cancel()
        image = ginga_draw.make_image(
            channel['image'], channel['chname'], waveimg=waveimg, rms_img=rms_img,
            logger=fv.logger
        )
        drawn = [_overlays[key] for key in channel['overlays'] if key in _overlays]
        for overlay in drawn:
            overlay.chnames.append(channel['chname'])
        fv.gui_call(_show_channel, fv, channel['chname'], image, channel['cuts'], drawn)
    if link_views and len(chnames) > 1:
        fv.gui_call(_link_views, fv, chnames)
    return chnames


def _start_view(fv, plugin_name, state):
    """
    Start a view plugin on the reference channel of a display.

    Must be called on the ginga GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    plugin_name : :obj:`str`
        Name of the local plugin, e.g. ``'Spec2dView'``.
    state : :class:`ViewState`
        What was displayed.
    """
    chname = state.chnames[0]
    fv.call_local_plugin_method(chname, plugin_name, 'set_state', [state], {})
    fv.start_local_plugin(chname, plugin_name)
    fv.change_channel(chname)


def chk_flats(fv, args, job):
    """
    Show the images in a flat-field calibration file.

    Handler for ``pypeit_chk_flats``; see
    :class:`~pypeit.scripts.chk_flats.ChkFlats` for the arguments.  The
    display is controlled by the ``FlatsView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_chk_flats``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    chk_version = not args.try_old
    flats = flatfield.FlatImages.from_file(args.file, chk_version=chk_version)
    slits_file = find_slits_file(flats.calib_key, flats.calib_dir, Path(args.file).parent)
    slits = read_slits(slits_file, chk_version=chk_version)
    job.check_cancel()
    products = flats.ginga_products(frametype=args.type, slits=slits, chk_version=chk_version)
    if len(products['images']) == 0:
        raise PypeItError(f'{args.file} has no {args.type} flat-field images to show.')

    state = ViewState('chk_flats', args, args.file)
    slit_kwargs = None
    if products['slits'] is not None:
        left, right, mask = products['slits'].select_edges()
        gpm = mask == 0
        slit_kwargs = dict(left=left[:,gpm], right=right[:,gpm],
                           slit_ids=products['slits'].spat_id[gpm])
        state.slits = slit_rows(slit_kwargs['left'], slit_kwargs['right'],
                                slit_ids=slit_kwargs['slit_ids'])
    overlays = make_overlays(slits=slit_kwargs)
    channels = [dict(image=img, chname=chname, cuts=cuts, overlays=['slit_edges', 'slit_ids'])
                for img, chname, cuts in products['images']]
    state.chnames = display_channels(
        fv, job, channels, overlays=overlays, waveimg=products['waveimg']
    )
    state.overlays = list(overlays.values())
    fv.gui_call(_start_view, fv, 'FlatsView', state)


def chk_edges(fv, args, job):
    """
    Show the trace image and the slit-edge traces in an edges calibration
    file.

    Handler for ``pypeit_chk_edges``; see
    :class:`~pypeit.scripts.chk_edges.ChkEdges` for the arguments.  The
    display is controlled by the ``EdgesView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_chk_edges``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    chk_version = not args.try_old
    edges = edgetrace.EdgeTraceSet.from_file(args.trace_file, chk_version=chk_version)

    # Use the slits file provided, or else look for the default one
    slits_file = None if args.slits_file is None else Path(args.slits_file)
    if slits_file is not None and not slits_file.exists():
        log.warning(f'{slits_file} does not exist!')
        slits_file = None
    if slits_file is None:
        slits_file = find_slits_file(
            edges.traceimg.calib_key, edges.traceimg.calib_dir, Path(args.trace_file).parent
        )
    slits = read_slits(slits_file, chk_version=chk_version)

    job.check_cancel()
    products = edges.ginga_products(thin=10, slits=slits)
    state = ViewState('chk_edges', args, args.trace_file)
    slit_kwargs = products['overlay']
    if slit_kwargs is not None and slit_kwargs['synced']:
        state.slits = slit_rows(
            slit_kwargs['left'], slit_kwargs['right'], slit_ids=slit_kwargs.get('slit_ids'),
            maskdef_ids=slit_kwargs.get('maskdef_ids')
        )
    overlays = make_overlays(slits=slit_kwargs)
    channels = [dict(image=img, chname=chname, cuts=None, overlays=['slit_edges', 'slit_ids'])
                for img, chname in products['images']]
    state.chnames = display_channels(fv, job, channels, overlays=overlays)
    state.overlays = list(overlays.values())
    fv.gui_call(_start_view, fv, 'EdgesView', state)


def show_2dspec(fv, args, job):
    """
    Show the images in a spec2d file, with the slit edges and object traces.

    Handler for ``pypeit_show_2dspec``; see
    :class:`~pypeit.scripts.show_2dspec.Show2DSpec` for the arguments.  The
    display is controlled by the ``Spec2dView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_show_2dspec``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    data = spec2d_data.load_spec2d(args)
    job.check_cancel()

    state = ViewState('show_2dspec', args, args.file)
    state.options = dict(detname=data['detname'], detectors=data['detectors'])
    if data['slits'] is not None:
        state.slits = slit_rows(
            data['slits']['left'], data['slits']['right'], slit_ids=data['slits']['slit_ids'],
            maskdef_ids=data['slits']['maskdef_ids']
        )
    overlays = make_overlays(slits=data['slits'], traces=data['traces'])
    channels = []
    for channel in data['channels']:
        keys = []
        if channel['slits']:
            keys += ['slit_edges', 'slit_ids']
        if channel['traces']:
            keys += ['traces', 'trace_names']
        channels.append(dict(image=channel['image'], chname=channel['chname'],
                             cuts=channel['cuts'], overlays=keys))
    state.chnames = display_channels(
        fv, job, channels, overlays=overlays, waveimg=data['waveimg'], rms_img=data['rms_img'],
        clear=args.clear
    )
    state.overlays = list(overlays.values())
    fv.gui_call(_start_view, fv, 'Spec2dView', state)


def read_calib_image(calib_file, chk_version=True):
    """
    Read a single-image calibration file, e.g. a Bias or Arc file.

    The calibration class is identified by the ``DMODCLS`` keyword of the
    first extension.  If the class is unknown, or the file cannot be parsed
    with the current datamodel and ``chk_version`` is False, the first image
    extension is read directly.

    Parameters
    ----------
    calib_file : :obj:`str`, `Path`_
        The calibration file.
    chk_version : :obj:`bool`, optional
        Require the datamodel version of the file to match the current one.

    Returns
    -------
    :obj:`dict`
        Dictionary with the ``image``, the inverse variance (``ivar``) and
        mask bits (``mask``), which are None if not available, and the
        calibration key and directory (``calib_key``, ``calib_dir``), which
        are None if unknown.
    """
    with fits.open(calib_file) as hdu:
        clsname = hdu[1].header.get('DMODCLS') if len(hdu) > 1 else None
    cls = getattr(buildimage, clsname, None) if isinstance(clsname, str) else None
    if cls is not None and issubclass(cls, pypeitimage.PypeItCalibrationImage):
        try:
            calib = cls.from_file(calib_file, chk_version=chk_version)
        except (PypeItError, TypeError, ValueError, KeyError):
            # Files written by older versions of PypeIt can fail in several ways
            if chk_version:
                raise
            log.warning(f'Could not parse {calib_file} as a {clsname}; showing its first image.')
        else:
            return dict(
                image=calib.image, ivar=calib.ivar,
                mask=None if calib.fullmask is None else calib.fullmask.mask,
                calib_key=calib.calib_key, calib_dir=calib.calib_dir
            )
    with fits.open(calib_file) as hdu:
        images = [h.data for h in hdu if h.data is not None and h.data.ndim == 2]
    if len(images) == 0:
        raise PypeItError(f'{calib_file} has no image to show.')
    return dict(image=images[0], ivar=None, mask=None, calib_key=None, calib_dir=None)


def show_calib(fv, args, job):
    """
    Show a single-image calibration file, e.g. a Bias, Arc or Tiltimg file.

    Handler for ``pypeit_show_calib``; see
    :class:`~pypeit.scripts.show_calib.ShowCalib` for the arguments.  The
    display is controlled by the ``CalibView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_show_calib``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    chk_version = not args.try_old
    calib = read_calib_image(args.file, chk_version=chk_version)
    job.check_cancel()

    state = ViewState('show_calib', args, args.file)
    slit_kwargs = None
    if args.slits:
        slits_file = find_slits_file(calib['calib_key'], calib['calib_dir'],
                                     Path(args.file).parent)
        slits = read_slits(slits_file, chk_version=chk_version)
        if slits is not None:
            left, right, mask = slits.select_edges()
            gpm = mask == 0
            slit_kwargs = dict(left=left[:,gpm], right=right[:,gpm],
                               slit_ids=slits.slitord_id[gpm],
                               maskdef_ids=None if slits.maskdef_id is None
                                           else slits.maskdef_id[gpm])
            state.slits = slit_rows(slit_kwargs['left'], slit_kwargs['right'],
                                    slit_ids=slit_kwargs['slit_ids'],
                                    maskdef_ids=slit_kwargs['maskdef_ids'])
    overlays = make_overlays(slits=slit_kwargs)

    stem = Path(args.file).name.split('.')[0]
    keys = ['slit_edges', 'slit_ids']
    channels = [dict(image=calib['image'], chname=stem, cuts=None, overlays=keys)]
    if args.ivar:
        if calib['ivar'] is None:
            log.warning(f'{args.file} has no inverse variance image.')
        else:
            channels.append(dict(image=calib['ivar'], chname=f'{stem}-ivar', cuts=None,
                                 overlays=keys))
    if args.mask:
        if calib['mask'] is None:
            log.warning(f'{args.file} has no mask.')
        else:
            channels.append(dict(image=calib['mask'], chname=f'{stem}-mask', cuts=None,
                                 overlays=keys))
    state.chnames = display_channels(fv, job, channels, overlays=overlays)
    state.overlays = list(overlays.values())
    fv.gui_call(_start_view, fv, 'CalibView', state)


def _slit_overlay_args(slits):
    """
    Return the slit-overlay arguments and slit-table rows for a set of slits.

    Parameters
    ----------
    slits : :class:`~pypeit.slittrace.SlitTraceSet`
        The slits.  Only the unmasked slits are used.

    Returns
    -------
    :obj:`tuple`
        Keyword arguments for :func:`~pypeit.display.ginga_draw.slit_overlay`
        and the rows for the slit navigator (see :func:`slit_rows`).
    """
    left, right, mask = slits.select_edges()
    gpm = mask == 0
    kwargs = dict(left=left[:,gpm], right=right[:,gpm], slit_ids=slits.spat_id[gpm],
                  maskdef_ids=None if slits.maskdef_id is None else slits.maskdef_id[gpm])
    return kwargs, slit_rows(kwargs['left'], kwargs['right'], slit_ids=kwargs['slit_ids'],
                             maskdef_ids=kwargs['maskdef_ids'])


def chk_tilts(fv, args, job):
    """
    Show the tilt image with the traced and fitted arc-line tilts.

    Handler for ``pypeit_chk_tilts``; see
    :class:`~pypeit.scripts.chk_tilts.ChkTilts` for the arguments.  The
    display is controlled by the ``TiltsView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_chk_tilts``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    chk_version = not args.try_old
    tilts = wavetilts.WaveTilts.from_file(args.file, chk_version=chk_version)
    products = tilts.ginga_products(calib_dir=Path(args.file).parent, chk_version=chk_version)
    job.check_cancel()

    state = ViewState('chk_tilts', args, args.file)
    slit_kwargs = None
    if products['left'] is not None:
        slit_kwargs = dict(left=products['left'], right=products['right'],
                           slit_ids=products['slit_ids'])
        state.slits = slit_rows(products['left'], products['right'],
                                slit_ids=products['slit_ids'])
    overlays = make_overlays(slits=slit_kwargs)
    for key, (name, canvas) in ginga_draw.tilt_overlays(
            tilts.tilt_traces, products['image'].shape[0], show_traces=args.show_traces).items():
        overlays[key] = Overlay(name, f'pypeit_{key}', canvas)
    channels = [dict(image=products['image'], chname='Tilts', cuts=products['cuts'],
                     overlays=list(overlays.keys()))]
    state.chnames = display_channels(
        fv, job, channels, overlays=overlays, waveimg=products['waveimg']
    )
    state.overlays = list(overlays.values())
    fv.gui_call(_start_view, fv, 'TiltsView', state)


def view_fits(fv, args, job):
    """
    Show a raw or processed image.

    Handler for ``pypeit_view_fits``; see
    :class:`~pypeit.scripts.view_fits.ViewFits` for the arguments.  The
    display is controlled by the ``FileView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_view_fits``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    images = image_data.view_fits_images(args)
    job.check_cancel()
    state = ViewState('view_fits', args, args.file)
    channels = [dict(image=img, chname=chname, cuts=None, overlays=[])
                for img, chname in images]
    state.chnames = display_channels(fv, job, channels, clear=args.clear)
    fv.gui_call(_start_view, fv, 'FileView', state)


def show_pixflat(fv, args, job):
    """
    Show the detector images in an archived pixel-flat file.

    Handler for ``pypeit_show_pixflat``; see
    :class:`~pypeit.scripts.show_pixflat.ShowPixFlat` for the arguments.
    The display is controlled by the ``FileView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_show_pixflat``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    images = image_data.pixflat_images(args.file, det=args.det)
    job.check_cancel()
    state = ViewState('show_pixflat', args, args.file)
    channels = [dict(image=img, chname=chname, cuts=(0.9, 1.1), overlays=[])
                for img, chname in images]
    # The detectors can have different shapes, so the channels are not linked
    state.chnames = display_channels(fv, job, channels, link_views=False)
    fv.gui_call(_start_view, fv, 'FileView', state)


def chk_scattlight(fv, args, job):
    """
    Show the scattered-light data, model and residuals.

    Handler for ``pypeit_chk_scattlight``; see
    :class:`~pypeit.scripts.chk_scattlight.ChkScattLight` for the arguments.
    The display is controlled by the ``FileView`` plugin.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_chk_scattlight``.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command.
    """
    if args.spec2d is not None:
        raise PypeItError('Displaying the spec2d scattered light is not currently supported.')
    chk_version = not args.try_old
    model = scattlight.ScatteredLight.from_file(args.file, chk_version=chk_version)
    slits = slittrace.SlitTraceSet.from_file(args.slits, chk_version=chk_version)
    images = model.ginga_products(slits=slits, mask=args.mask)
    job.check_cancel()

    state = ViewState('chk_scattlight', args, args.file)
    slit_kwargs, state.slits = _slit_overlay_args(slits)
    overlays = make_overlays(slits=slit_kwargs)
    channels = [dict(image=img, chname=chname, cuts=cuts, overlays=['slit_edges', 'slit_ids'])
                for img, chname, cuts in images]
    state.chnames = display_channels(fv, job, channels, overlays=overlays)
    state.overlays = list(overlays.values())
    fv.gui_call(_start_view, fv, 'FileView', state)
