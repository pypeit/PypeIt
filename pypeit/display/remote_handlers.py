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

from pypeit import edgetrace
from pypeit import flatfield
from pypeit import log
from pypeit import PypeItError
from pypeit import slittrace
from pypeit import specobjs
from pypeit.display import ginga_draw
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
        'show_2dspec': show_2dspec
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


def _show_channel(fv, chname, image, cuts, overlay):
    """
    Show an image, and optionally an overlay, in a channel.

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
    overlay : ``ginga.canvas.CanvasObject.CanvasObjectBase``
        Overlay to draw, or None.
    """
    channel = ginga_draw.show_image(fv, chname, image, cuts=cuts)
    if overlay is not None:
        ginga_draw.add_overlay(channel, overlay, 'pypeit_slits')


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


def show_images(fv, job, images, overlay=None, waveimg=None, link_views=False):
    """
    Clear the viewer and show a set of images, one per channel.

    Must be called on a ginga worker thread.  The images are built on the
    calling thread; only the display is done on the GUI thread.

    Parameters
    ----------
    fv : ``ginga.rv.Control.GingaShell``
        The ginga reference viewer shell.
    job : :class:`~pypeit.display.remote_jobs.RemoteJob`
        The job running the command; checked for cancellation before each
        channel.
    images : :obj:`list`
        List of ``(image, chname, cuts)`` tuples, where ``image`` is a
        `numpy.ndarray`_, ``chname`` is the channel name, and ``cuts`` is a
        tuple with the lower and upper cut levels or None.
    overlay : :obj:`dict`, optional
        Keyword arguments for :func:`~pypeit.display.ginga_draw.slit_overlay`.
        The overlay is built once and shared by every channel.  If None, no
        overlay is drawn.
    waveimg : `numpy.ndarray`_, optional
        Wavelength image shown in the cursor readout of every channel.
    link_views : :obj:`bool`, optional
        Lock the pan, zoom and rotation of all channels together; see
        :func:`_link_views`.
    """
    # One overlay canvas is built and shared by all channels
    slits = None if overlay is None else ginga_draw.slit_overlay(**overlay)
    fv.gui_call(clear_workspace, fv)
    for img, chname, cuts in images:
        job.check_cancel()
        image = ginga_draw.make_image(img, chname, waveimg=waveimg, logger=fv.logger)
        fv.gui_call(_show_channel, fv, chname, image, cuts, slits)
    if link_views and len(images) > 1:
        fv.gui_call(_link_views, fv, [chname for _, chname, _ in images])


def chk_flats(fv, args, job):
    """
    Show the images in a flat-field calibration file.

    Handler for ``pypeit_chk_flats``; see
    :class:`~pypeit.scripts.chk_flats.ChkFlats` for the arguments.

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
    overlay = None
    if products['slits'] is not None:
        left, right, mask = products['slits'].select_edges()
        gpm = mask == 0
        overlay = dict(left=left[:,gpm], right=right[:,gpm],
                       slit_ids=products['slits'].spat_id[gpm])
    show_images(
        fv, job, products['images'], overlay=overlay, waveimg=products['waveimg'],
        link_views=True
    )


def chk_edges(fv, args, job):
    """
    Show the trace image and the slit-edge traces in an edges calibration
    file.

    Handler for ``pypeit_chk_edges``; see
    :class:`~pypeit.scripts.chk_edges.ChkEdges` for the arguments.

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
    images = [(img, chname, None) for img, chname in products['images']]
    show_images(fv, job, images, overlay=products['overlay'], link_views=True)


def _show_spec2d_channel(fv, chname, image, cuts, overlays):
    """
    Show one spec2d image and its overlays.

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
    overlays : :obj:`dict`
        Overlay canvases to draw, keyed by their tag.
    """
    channel = ginga_draw.show_image(fv, chname, image, cuts=cuts)
    for tag, overlay in overlays.items():
        ginga_draw.add_overlay(channel, overlay, tag)


def show_2dspec(fv, args, job):
    """
    Show the images in a spec2d file, with the slit edges and object traces.

    Handler for ``pypeit_show_2dspec``; see
    :class:`~pypeit.scripts.show_2dspec.Show2DSpec` for the arguments.

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

    # Each overlay is built once and shared by the channels that show it
    slits = None if data['slits'] is None else ginga_draw.slit_overlay(**data['slits'])
    traces = None if data['traces'] is None else ginga_draw.trace_overlay(**data['traces'])

    chnames = [channel['chname'] for channel in data['channels']]
    if args.clear:
        fv.gui_call(clear_workspace, fv)
    else:
        # Replace only the channels being shown again
        fv.gui_call(ginga_draw.clear_channels, fv, chnames)

    for channel in data['channels']:
        job.check_cancel()
        image = ginga_draw.make_image(
            channel['image'], channel['chname'], waveimg=data['waveimg'],
            rms_img=data['rms_img'], logger=fv.logger
        )
        overlays = {}
        if channel['slits'] and slits is not None:
            overlays['pypeit_slits'] = slits
        if channel['traces'] and traces is not None:
            overlays['pypeit_traces'] = traces
        fv.gui_call(
            _show_spec2d_channel, fv, channel['chname'], image, channel['cuts'], overlays
        )
    if len(chnames) > 1:
        fv.gui_call(_link_views, fv, chnames)
