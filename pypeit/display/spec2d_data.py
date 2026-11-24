"""
Prepare the images and overlays of a spec2d file for display.

The functions here read a spec2d file and build everything needed to display
it, without interacting with ginga.  They are used by the
``show_2dspec`` command of the
:class:`~pypeit.display.pypeit_remote.PypeItRemote` ginga plugin; see
:func:`pypeit.display.remote_handlers.show_2dspec`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path

import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import __version__
from pypeit import io
from pypeit import log
from pypeit import PypeItBitMaskError
from pypeit import PypeItDataModelError
from pypeit import PypeItError
from pypeit import spec2dobj
from pypeit import specobjs
from pypeit.images.detector_container import DetectorContainer


def wave_rms_image(slits, wavesol, flexure=None):
    """
    Construct an image with the RMS of the wavelength solution of each slit.

    Each wavelength solution has a single RMS, measured in pixels, so every
    pixel in a given slit (or echelle order) is assigned the same value.

    Parameters
    ----------
    slits : :class:`~pypeit.slittrace.SlitTraceSet`
        Slit traces used to associate each pixel with a slit.
    wavesol : `astropy.table.Table`_
        Wavelength-calibration diagnostics table, as constructed by
        :func:`~pypeit.wavecalib.WaveCalib.wave_diagnostics`.  Can be None.
    flexure : :obj:`float`, optional
        Spatial flexure shift, in pixels, applied to the slit edges.

    Returns
    -------
    `numpy.ndarray`_
        Image with the RMS, in pixels, of the wavelength solution for the slit
        that each pixel falls in.  Pixels that are not in an unmasked slit, or
        that are in a slit without a wavelength solution, are set to 0.  None is
        returned if the diagnostics table is unavailable.
    """
    if wavesol is None or any(c not in wavesol.colnames for c in ['SpatOrderID', 'RMS']):
        return None

    # RMS for each slit, matched by spatial/order ID; 0 means no solution
    wave_id = np.asarray(wavesol['SpatOrderID'])
    wave_rms = np.asarray(wavesol['RMS'], dtype=float)
    slit_rms = np.zeros(slits.nslits, dtype=float)
    for i, slitord_id in enumerate(slits.slitord_id):
        indx = np.where(wave_id == slitord_id)[0]
        if indx.size > 0:
            slit_rms[i] = wave_rms[indx[0]]

    # Pixels outside of all (unmasked) slits have an index of -1
    slitidx_img = slits.slit_img(flexure=flexure, use_spatial=False)
    return np.where(slitidx_img >= 0, slit_rms[slitidx_img], 0.)


def trace_overlay_args(sobjs, detname):
    """
    Collect the object traces on one detector for display.

    Parameters
    ----------
    sobjs : :class:`~pypeit.specobjs.SpecObjs`
        Objects extracted from the spec2d file.  Can be None.
    detname : :obj:`str`
        Name of the detector or mosaic, e.g. ``DET01``.

    Returns
    -------
    :obj:`dict`
        Keyword arguments for :func:`~pypeit.display.ginga_draw.trace_overlay`,
        or None if there are no objects on the detector.
    """
    if sobjs is None:
        return None
    in_det = np.where(sobjs.DET == detname)[0]
    if len(in_det) == 0:
        log.warning('spec1d file found, but no objects were extracted for this detector.')
        return None
    traces = []
    names = []
    maskdef_extr = []
    manual_extr = []
    for kk in in_det:
        maskdef_objname = sobjs[kk].MASKDEF_OBJNAME
        traces.append(sobjs[kk]['TRACE_SPAT'])
        names.append(sobjs[kk].NAME if maskdef_objname is None
                     else f'{sobjs[kk].NAME}     OBJNAME:{maskdef_objname}')
        maskdef_extr.append(sobjs[kk].MASKDEF_EXTRACT is True)
        manual_extr.append(sobjs[kk].hand_extract_flag is True)
    return dict(
        trace=np.swapaxes(np.array(traces), 1, 0), trc_name=np.array(names),
        maskdef_extr=np.array(maskdef_extr), manual_extr=np.array(manual_extr)
    )


def spec2d_detectors(spec2d_file):
    """
    Return the names of the detectors or mosaics in a spec2d file.

    Parameters
    ----------
    spec2d_file : :obj:`str`, `Path`_
        The spec2d file.

    Returns
    -------
    :obj:`list`
        Detector or mosaic names, e.g. ``['DET01', 'DET02']``.
    """
    return fits.getval(spec2d_file, 'HIERARCH ALLSPEC2D_DETS').split(',')


def spec2d_detname(spec2d_file, det=None):
    """
    Return the name of the detector or mosaic to show.

    Parameters
    ----------
    spec2d_file : :obj:`str`, `Path`_
        The spec2d file.
    det : :obj:`str`, optional
        Detector name or number.  If a number, the name is constructed
        assuming the reduction is for a single detector.  If None, the first
        detector in the file is used.

    Returns
    -------
    :obj:`str`
        Detector or mosaic name, e.g. ``DET01`` or ``MSC01``.
    """
    if det is None:
        return spec2d_detectors(spec2d_file)[0]
    try:
        return DetectorContainer.get_name(int(det))
    except ValueError:
        return det


def _read_spec2d_hdus(spec2d_file, detname, show_wave_rms=False):
    """
    Read the components of a spec2d file directly from its HDUs.

    This is the fallback used when the file cannot be parsed by the current
    :class:`~pypeit.spec2dobj.Spec2DObj` datamodel.

    Parameters
    ----------
    spec2d_file : :obj:`str`, `Path`_
        The spec2d file.
    detname : :obj:`str`
        Name of the detector or mosaic.
    show_wave_rms : :obj:`bool`, optional
        Whether the wavelength RMS was requested; only used to issue a
        warning, because it cannot be shown in this case.

    Returns
    -------
    :obj:`dict`
        The images, masks and slit edges; see :func:`_read_spec2d_object`.
    """
    with io.fits_open(spec2d_file) as hdu:
        names = [h.name for h in hdu]
        if not any([detname in n for n in names]):
            raise PypeItError(f'Provided file has no extensions including {detname}.')
        for ext in ['SCIIMG', 'SKYMODEL', 'OBJMODEL', 'IVARMODEL']:
            if f'{detname}-{ext}' not in names:
                raise PypeItError(f'{spec2d_file} missing extension {detname}-{ext}.')

        data = dict(
            sciimg=hdu[f'{detname}-SCIIMG'].data, skymodel=hdu[f'{detname}-SKYMODEL'].data,
            objmodel=hdu[f'{detname}-OBJMODEL'].data, ivarmodel=hdu[f'{detname}-IVARMODEL'].data,
            rms_img=None, sobjs=None, slit_left=None, slit_right=None, slit_mask=None,
            slit_spat_id=None, slit_mask_id=None, slit_ids=None
        )
        _ext = f'{detname}-BPMMASK'
        data['bpmmask'] = hdu[_ext].data if _ext in names else None
        _ext = f'{detname}-WAVEIMG'
        data['waveimg'] = hdu[_ext].data if _ext in names else None
        if show_wave_rms:
            log.warning('Cannot show the wavelength RMS when the spec2d file cannot be '
                        'parsed directly.')

        data['img_gpm'] = np.ones(data['sciimg'].shape, dtype=bool) \
                            if data['bpmmask'] is None else data['bpmmask'] == 0
        data['model_gpm'] = data['img_gpm'].copy()

        _ext = f'{detname}-SLITS'
        if _ext not in names:
            log.warning(f'{spec2d_file} missing extension {_ext}; cannot show slit edges.')
            return data

        cols = hdu[_ext].columns.names
        slits = hdu[_ext].data
        data['slit_spat_id'] = slits['spat_id'] if 'spat_id' in cols else None
        data['slit_left'] = slits['left_tweak'].T if 'left_tweak' in cols \
                                else (slits['left_init'].T if 'left_init' in cols else None)
        data['slit_right'] = slits['right_tweak'].T if 'right_tweak' in cols \
                                else (slits['right_init'].T if 'right_init' in cols else None)
        data['slit_mask'] = slits['mask'] if 'mask' in cols else None
        data['slit_mask_id'] = slits['maskdef_id'] if 'maskdef_id' in cols else None
        hdr = hdu[f'{detname}-SCIIMG'].header
        if 'SCI_SPAT_FLEXURE' in hdr and data['slit_left'] is not None \
                and data['slit_right'] is not None:
            flexure = float(hdr['SCI_SPAT_FLEXURE'])
            data['slit_left'] = data['slit_left'] + flexure
            data['slit_right'] = data['slit_right'] + flexure
            log.info(f'Offseting slits by {flexure} pixels.')
        pypeline = hdr['PYPELINE'] if 'PYPELINE' in hdr else None
        if pypeline in ['MultiSlit', 'SlicerIFU', 'Fiber']:
            data['slit_ids'] = data['slit_spat_id']
        elif pypeline == 'Echelle':
            data['slit_ids'] = slits['ech_order'] if 'ech_order' in cols else None
    return data


def _read_spec2d_object(spec2DObj, spec2d_file, detname, show_wave_rms=False,
                        ignore_extract_mask=False):
    """
    Collect the components of a parsed spec2d file.

    Parameters
    ----------
    spec2DObj : :class:`~pypeit.spec2dobj.Spec2DObj`
        The parsed spec2d data for one detector.
    spec2d_file : :obj:`str`, `Path`_
        The spec2d file; used to find the matching spec1d file.
    detname : :obj:`str`
        Name of the detector or mosaic.
    show_wave_rms : :obj:`bool`, optional
        Construct the wavelength RMS image.
    ignore_extract_mask : :obj:`bool`, optional
        Ignore the extraction mask in the model images.

    Returns
    -------
    :obj:`dict`
        Dictionary with the images (``sciimg``, ``skymodel``, ``objmodel``,
        ``ivarmodel``, ``waveimg``, ``rms_img``), masks (``bpmmask``,
        ``img_gpm``, ``model_gpm``), slit edges and IDs (``slit_left``,
        ``slit_right``, ``slit_mask``, ``slit_spat_id``, ``slit_mask_id``,
        ``slit_ids``), and the extracted objects (``sobjs``).
    """
    data = dict(
        sciimg=spec2DObj.sciimg, skymodel=spec2DObj.skymodel, objmodel=spec2DObj.objmodel,
        ivarmodel=spec2DObj.ivarmodel, bpmmask=spec2DObj.bpmmask, waveimg=spec2DObj.waveimg,
        rms_img=None
    )
    if show_wave_rms and data['waveimg'] is not None:
        data['rms_img'] = wave_rms_image(spec2DObj.slits, spec2DObj.wavesol,
                                         flexure=spec2DObj.sci_spat_flexure)
        if data['rms_img'] is None:
            log.warning(f'{spec2d_file} does not include the wavelength-calibration '
                        'diagnostics; cannot show the wavelength RMS.')

    data['img_gpm'] = spec2DObj.select_flag(invert=True)
    if not np.any(data['img_gpm']):
        log.warning('The full science image is masked!')
    data['model_gpm'] = data['img_gpm'].copy()
    if ignore_extract_mask:
        data['model_gpm'] |= spec2DObj.select_flag(flag='EXTRACT')

    if spec2DObj.sci_spat_flexure is not None:
        log.info(f'Offseting slits by {spec2DObj.sci_spat_flexure}')
    data['slit_left'], data['slit_right'], data['slit_mask'] \
            = spec2DObj.slits.select_edges(flexure=spec2DObj.sci_spat_flexure)
    data['slit_spat_id'] = spec2DObj.slits.spat_id
    data['slit_mask_id'] = spec2DObj.slits.maskdef_id
    data['slit_ids'] = spec2DObj.slits.slitord_id

    # Object traces from the spec1d file
    spec1d_file = Path(str(spec2d_file).replace('spec2d', 'spec1d'))
    if spec1d_file.suffix == '.gz':
        spec1d_file = spec1d_file.with_suffix('')
    if spec1d_file.is_file():
        data['sobjs'] = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=False)
    else:
        data['sobjs'] = None
        log.warning(f'Could not find spec1d file: {spec1d_file}\nNo objects were extracted.')
    return data


def _read_spec2d(spec2d_file, detname, chk_version=True, show_wave_rms=False,
                 ignore_extract_mask=False):
    """
    Read the components of a spec2d file for one detector.

    The file is parsed with :class:`~pypeit.spec2dobj.Spec2DObj` if possible.
    If its datamodel is incompatible and ``chk_version`` is False, the
    components are read directly from the HDUs.

    Parameters
    ----------
    spec2d_file : :obj:`str`, `Path`_
        The spec2d file.
    detname : :obj:`str`
        Name of the detector or mosaic.
    chk_version : :obj:`bool`, optional
        Require the datamodel version of the file to match the current one.
    show_wave_rms : :obj:`bool`, optional
        Construct the wavelength RMS image.
    ignore_extract_mask : :obj:`bool`, optional
        Ignore the extraction mask in the model images.

    Returns
    -------
    :obj:`dict`
        The images, masks, slit edges and objects; see
        :func:`_read_spec2d_object`.
    """
    try:
        spec2DObj = spec2dobj.Spec2DObj.from_file(spec2d_file, detname, chk_version=chk_version)
    except (PypeItDataModelError, PypeItBitMaskError):
        try:
            file_pypeit_version = fits.getval(spec2d_file, 'VERSPYP', 0)
        except KeyError:
            file_pypeit_version = '*unknown*'
        if chk_version:
            addendum = 'To allow the script to attempt to read the data anyway, use the ' \
                       '--try_old command-line option.  This will first try to simply ' \
                       'ignore the version number.  If the datamodels are incompatible ' \
                       '(e.g., the new datamodel contains components not in a previous ' \
                       'version), this may not be enough and the script will continue by ' \
                       'trying to parse only the components necessary for use by this ' \
                       'script. In either case, BEWARE that the displayed data may be in ' \
                       'error!'
        else:
            addendum = 'The datamodels are sufficiently different that the script will now ' \
                       'try to parse only the components necessary for use by this ' \
                       'script.  BEWARE that the displayed data may be in error!'
        message = (
            f'Your installed version of PypeIt ({__version__}) cannot be used to parse '
            f'{spec2d_file}, which was reduced using version {file_pypeit_version}.  You '
            'are strongly encouraged to re-reduce your data using this (or, better yet, '
            'the most recent) version of PypeIt.  ' + addendum
        )
        if chk_version:
            raise PypeItError(message)
        log.warning(message)
        return _read_spec2d_hdus(spec2d_file, detname, show_wave_rms=show_wave_rms)
    return _read_spec2d_object(
        spec2DObj, spec2d_file, detname, show_wave_rms=show_wave_rms,
        ignore_extract_mask=ignore_extract_mask
    )


def _clipped_cuts(image, gpm):
    """
    Return cut levels of mean - 1 sigma and mean + 4 sigma of the good pixels.

    Parameters
    ----------
    image : `numpy.ndarray`_
        Image.
    gpm : `numpy.ndarray`_
        Good-pixel mask.

    Returns
    -------
    :obj:`tuple`
        Lower and upper cut levels.
    """
    mean, _, sigma = sigma_clipped_stats(image[gpm], sigma_lower=5.0, sigma_upper=5.0)
    return (mean - 1.0*sigma, mean + 4.0*sigma)


def load_spec2d(args):
    """
    Read a spec2d file and build the images and overlays to display.

    Parameters
    ----------
    args : `argparse.Namespace`_
        Arguments parsed by ``pypeit_show_2dspec``; see
        :class:`~pypeit.scripts.show_2dspec.Show2DSpec`.

    Returns
    -------
    :obj:`dict`
        Dictionary with:

            - ``channels``: list of dictionaries, one per channel to show, in
              order, each with the ``image`` to show, the channel name
              (``chname``), the cut levels (``cuts``), and whether to draw
              the slit edges (``slits``) and object traces (``traces``);
            - ``waveimg`` and ``rms_img``: the wavelength and wavelength-RMS
              images shown in the cursor readout (either can be None);
            - ``slits``: keyword arguments for
              :func:`~pypeit.display.ginga_draw.slit_overlay`, or None;
            - ``traces``: keyword arguments for
              :func:`~pypeit.display.ginga_draw.trace_overlay`, or None;
            - ``detname``: the detector or mosaic shown;
            - ``detectors``: all detectors or mosaics in the file.
    """
    detname = spec2d_detname(args.file, det=args.det)
    show_channels = [0, 1, 2, 3] if args.channels is None \
                        else [int(item) for item in args.channels.split(',')]

    # The MMT Binospec IFU produces >700 fiber traces per exposure;
    # overplotting them all makes the display unusably slow (and
    # uninformative), so they are skipped.
    try:
        pyp_spec = fits.getval(args.file, 'PYP_SPEC', 0)
    except (KeyError, OSError):
        pyp_spec = None
    show_overlays = pyp_spec != 'mmt_binospec_ifu'
    if not show_overlays:
        log.info('Detected MMT Binospec IFU spec2d; skipping object-trace '
                 'and slit-edge overlays for performance (>700 fibers).')

    data = _read_spec2d(
        args.file, detname, chk_version=not args.try_old, show_wave_rms=args.show_wave_rms,
        ignore_extract_mask=args.ignore_extract_mask
    )

    # Select the slits to show
    slits = None
    if data['slit_left'] is not None and data['slit_right'] is not None:
        # TODO: This may be too restrictive, i.e. ignore BADFLTCALIB??
        slit_gpm = np.ones(data['slit_left'].shape[1], dtype=bool) \
                        if data['slit_mask'] is None else data['slit_mask'] == 0
        if data['slit_spat_id'] is not None and args.spat_id is not None:
            slit_gpm &= data['slit_spat_id'] == args.spat_id
        elif data['slit_mask_id'] is not None and args.maskID is not None:
            slit_gpm &= data['slit_mask_id'] == args.maskID
        slits = dict(
            left=data['slit_left'][:,slit_gpm], right=data['slit_right'][:,slit_gpm],
            slit_ids=None if data['slit_ids'] is None else data['slit_ids'][slit_gpm],
            maskdef_ids=None if data['slit_mask_id'] is None
                        else data['slit_mask_id'][slit_gpm]
        )
    traces = trace_overlay_args(data['sobjs'], detname)
    if not show_overlays:
        slits = None
        traces = None

    sciimg, skymodel, objmodel = data['sciimg'], data['skymodel'], data['objmodel']
    ivarmodel, img_gpm, model_gpm = data['ivarmodel'], data['img_gpm'], data['model_gpm']
    show_resid_traces = not args.removetrace
    channels = []
    if args.showmask is not None and data['bpmmask'] is not None:
        mask = data['bpmmask'].mask if len(args.showmask) == 0 \
                    else data['bpmmask'].flagged(flag=args.showmask).astype(np.uint8)
        channels.append(dict(image=mask, chname='MASK', cuts=None, slits=False, traces=False))
    if 0 in show_channels:
        channels.append(dict(
            image=sciimg, chname=f'{args.prefix}sciimg-{detname}',
            cuts=_clipped_cuts(sciimg, img_gpm), slits=True, traces=True
        ))
    if 1 in show_channels and skymodel is not None:
        image = (sciimg - skymodel) * model_gpm.astype(float)
        channels.append(dict(
            image=image, chname=f'{args.prefix}skysub-{detname}',
            cuts=_clipped_cuts(image, model_gpm), slits=True, traces=show_resid_traces
        ))
    if 2 in show_channels and skymodel is not None:
        image = (sciimg - skymodel) * np.sqrt(ivarmodel) * model_gpm.astype(float)
        channels.append(dict(
            image=image, chname=f'{args.prefix}sky_resid-{detname}', cuts=(-5.0, 5.0),
            slits=True, traces=show_resid_traces
        ))
    if 3 in show_channels and objmodel is not None:
        image = (sciimg - skymodel - objmodel) * np.sqrt(ivarmodel) * img_gpm.astype(float)
        channels.append(dict(
            image=image, chname=f'{args.prefix}resid-{detname}', cuts=(-5.0, 5.0),
            slits=True, traces=show_resid_traces
        ))
    if len(channels) == 0:
        raise PypeItError('No images to show; check the --channels option.')

    return dict(
        channels=channels, waveimg=data['waveimg'], rms_img=data['rms_img'], slits=slits,
        traces=traces, detname=detname, detectors=spec2d_detectors(args.file)
    )
