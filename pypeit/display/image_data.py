"""
Prepare raw or processed images for display.

The functions here read the images shown by ``pypeit_view_fits`` and
``pypeit_show_pixflat``, without interacting with ginga.  They are used by the
:class:`~pypeit.display.pypeit_remote.PypeItRemote` ginga plugin; see
:mod:`pypeit.display.remote_handlers`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import io
from pypeit import log
from pypeit import PypeItError
from pypeit.images import buildimage
from pypeit.images import pypeitimage
from pypeit.spectrographs import util


def _raw_detector(spectrograph, det):
    """
    Parse the detectors requested for ``pypeit_view_fits``.

    Parameters
    ----------
    spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        The spectrograph.
    det : :obj:`list`
        Detector numbers as strings, or ``['mosaic']``.

    Returns
    -------
    :obj:`tuple`
        The detector (an :obj:`int`, or a :obj:`tuple` for a mosaic) and
        whether it is a mosaic.
    """
    if 'mosaic' in det:
        _det = spectrograph.default_mosaic
        if _det is None:
            raise PypeItError(f'{spectrograph.name} does not have a known mosaic')
        return _det, True
    try:
        _det = tuple(int(d) for d in det)
    except ValueError as e:
        raise PypeItError('Could not convert detector input to integer.') from e
    return (_det, True) if len(_det) > 1 else (_det[0], False)


def view_fits_images(args):
    """
    Read the images shown by ``pypeit_view_fits``.

    Parameters
    ----------
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_view_fits``; see
        :class:`~pypeit.scripts.view_fits.ViewFits`.

    Returns
    -------
    :obj:`list`
        List of ``(image, chname)`` tuples: the image, and the bad-pixel mask
        if requested with ``--showmask`` and ``--proc``.
    """
    if args.proc and args.exten is not None:
        raise PypeItError('You cannot specify --proc and --exten, since --exten shows the raw '
                          'image')
    if args.exten is not None and 'mosaic' in args.det:
        raise PypeItError('You cannot specify --exten and --det mosaic, since --mosaic '
                          'displays multiple extensions by definition')

    # Intermediate file?
    if args.inter:
        return [(pypeitimage.PypeItImage.from_file(args.file).image, args.chname)]
    if args.exten is not None:
        with io.fits_open(args.file) as hdu:
            return [(hdu[args.exten].data, args.chname)]

    spectrograph = util.load_spectrograph(args.spectrograph)
    bad_read_message = 'Unable to construct image due to a read or image processing ' \
                       'error.  Use case interpreted from command-line inputs requires ' \
                       'a raw image, not an output image product from pypeit.  To show ' \
                       'a pypeit output image, specify the extension using --exten.  ' \
                       'Use --list to show the extension names.'
    det, mosaic = _raw_detector(spectrograph, args.det)
    if not args.proc:
        if args.showmask:
            log.info('You need to use --proc with --showmask to show the mask.  Ignoring your '
                     'argument')
        try:
            return [(spectrograph.get_rawimage(args.file, det)[1], args.chname)]
        except Exception as e:
            raise PypeItError(bad_read_message + f'  Original exception -- {type(e).__name__}: '
                              f'{str(e)}') from e

    # Use the biasframe processing parameters because processing these frames
    # is independent of any other frames (ie., does not perform bias
    # subtraction or flat-fielding)
    par = spectrograph.default_pypeit_par()['calibrations']['biasframe']
    images = []
    for f in [args.file] if args.bkg_file is None else [args.file, args.bkg_file]:
        try:
            images.append(buildimage.buildimage_fromlist(spectrograph, det, par, [f],
                                                         mosaic=mosaic))
        except Exception as e:
            raise PypeItError(bad_read_message + f'  Original exception -- {type(e).__name__}: '
                              f'{str(e)}') from e
    img = images[0] if len(images) == 1 else images[0].sub(images[1])
    shown = [(img.image, args.chname)]
    if args.showmask:
        shown.append((img.fullmask.bpm.astype(np.uint8), 'BPM'))
    return shown


def pixflat_images(pixflat_file, det=None):
    """
    Read the detector images in an archived pixel-flat file.

    Parameters
    ----------
    pixflat_file : :obj:`str`, `Path`_
        The pixel-flat file.
    det : :obj:`list`, optional
        Detector numbers to read.  If None, all detectors are read.

    Returns
    -------
    :obj:`list`
        List of ``(image, chname)`` tuples, one per detector.
    """
    with io.fits_open(pixflat_file) as hdu:
        # get all the available detectors in the file
        file_dets = [int(h.name.split('-')[0].split('DET')[1]) for h in hdu[1:]]
        if det is not None:
            in_file = np.isin(det, file_dets)
            if not np.any(in_file):
                raise PypeItError('Provided detector(s) not found in the file. Available '
                                  f'detectors are {file_dets}')
            if np.any(np.logical_not(in_file)):
                log.warning(f'Detector(s) {np.array(det)[np.logical_not(in_file)]} not found '
                            f'in the file. Available detectors are {file_dets}')
        return [(h.data, h.name) for h, d in zip(hdu[1:], file_dets)
                if det is None or d in det]
