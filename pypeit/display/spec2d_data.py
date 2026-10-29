"""
Prepare the contents of a spec2d file for display.

This module reads a PypeIt spec2d file (and its spec1d file, if present) and
constructs the images, cut levels, slit edges, and object traces shown by
``pypeit_show_2dspec`` and the ``Spec2dView`` ginga plugin.  It does not
depend on ginga.

.. include:: ../include/links.rst
"""
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import log
from pypeit import PypeItError, PypeItDataModelError, PypeItBitMaskError
from pypeit import __version__
from pypeit import io
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit.images.detector_container import DetectorContainer

#: Names of the channels shown for each detector, in order; the index is the
#: channel number used by the ``--channels`` option of ``pypeit_show_2dspec``.
CHANNEL_NAMES = ['sciimg', 'skysub', 'sky_resid', 'resid']


@dataclass
class ChannelImage:
    """
    An image to show in one ginga channel.

    Attributes
    ----------
    name : :obj:`str`
        Channel name.
    image : `numpy.ndarray`_
        Image to show.
    cuts : :obj:`tuple`, optional
        Lower and upper cut levels.
    show_traces : :obj:`bool`
        Overplot the object traces in this channel.
    show_slits : :obj:`bool`
        Overplot the slit edges in this channel.
    key : :obj:`int`, :obj:`str`
        Index of the channel in :data:`CHANNEL_NAMES`, or ``'mask'`` for the
        mask channel.
    """
    name: str
    image: np.ndarray
    cuts: tuple = None
    show_traces: bool = True
    show_slits: bool = True
    key: object = None


@dataclass
class Spec2DDisplay:
    """
    Everything needed to display one detector of a spec2d file.

    Attributes
    ----------
    file : `Path`_
        The spec2d file.
    detname : :obj:`str`
        Name of the detector or mosaic shown.
    detnames : :obj:`list`
        Names of all detectors/mosaics in the file.
    pyp_spec : :obj:`str`
        PypeIt name of the spectrograph, if available.
    waveimg : `numpy.ndarray`_
        Wavelength image, if available.
    channels : :obj:`list`
        :class:`ChannelImage` objects for each channel to show.  Other
        channels can be constructed using :meth:`make_channel`.
    left, right : `numpy.ndarray`_
        Left and right edges of the selected slits, shape
        :math:`(N_{\\rm spec}, N_{\\rm slit})`.  None if not available.
    slit_ids : `numpy.ndarray`_
        Slit or order IDs of the selected slits.
    maskdef_ids : `numpy.ndarray`_
        Slit-mask IDs of the selected slits.
    traces : `numpy.ndarray`_
        Spatial positions of the object traces for this detector, shape
        :math:`(N_{\\rm spec}, N_{\\rm obj})`.  None if there are none.
    trace_names : `numpy.ndarray`_
        Names of the traced objects.
    maskdef_extract : `numpy.ndarray`_
        Flags for objects extracted using the slit-mask design.
    manual_extract : `numpy.ndarray`_
        Flags for manually extracted objects.
    prefix : :obj:`str`
        Prefix for the channel names.
    removetrace : :obj:`bool`
        Do not overplot the object traces in the skysub, sky_resid, and resid
        channels.
    showmask : :obj:`list`
        Mask flags used to construct the mask channel; see
        :func:`load_spec2d_display`.
    sciimg, skymodel, objmodel, ivarmodel : `numpy.ndarray`_
        The science image, sky and object models, and the inverse variance of
        the model, used to construct the channel images.
    img_gpm, model_gpm : `numpy.ndarray`_
        Good-pixel masks for the science and model images.
    bpmmask : :class:`~pypeit.images.imagebitmask.ImageBitMaskArray`, `numpy.ndarray`_
        Bit mask of the science image.
    """
    file: Path
    detname: str
    detnames: list
    pyp_spec: str = None
    waveimg: np.ndarray = None
    channels: list = field(default_factory=list)
    left: np.ndarray = None
    right: np.ndarray = None
    slit_ids: np.ndarray = None
    maskdef_ids: np.ndarray = None
    traces: np.ndarray = None
    trace_names: np.ndarray = None
    maskdef_extract: np.ndarray = None
    manual_extract: np.ndarray = None
    prefix: str = ''
    removetrace: bool = False
    showmask: list = None
    sciimg: np.ndarray = None
    skymodel: np.ndarray = None
    objmodel: np.ndarray = None
    ivarmodel: np.ndarray = None
    img_gpm: np.ndarray = None
    model_gpm: np.ndarray = None
    bpmmask: object = None

    def channel_name(self, key):
        """
        Return the name of a channel.

        Parameters
        ----------
        key : :obj:`int`, :obj:`str`
            Index of the channel in :data:`CHANNEL_NAMES`, or ``'mask'``.

        Returns
        -------
        :obj:`str`
            Channel name.
        """
        if key == 'mask':
            return f'{self.prefix}MASK'
        return f'{self.prefix}{CHANNEL_NAMES[key]}-{self.detname}'

    def make_channel(self, key):
        """
        Construct the image for a channel.

        Parameters
        ----------
        key : :obj:`int`, :obj:`str`
            Index of the channel in :data:`CHANNEL_NAMES`, or ``'mask'`` for
            the mask channel.

        Returns
        -------
        :class:`ChannelImage`
            The channel image, or None if the data needed are not available.
        """
        name = self.channel_name(key)
        if key == 'mask':
            if self.bpmmask is None:
                return None
            if isinstance(self.bpmmask, np.ndarray):
                _mask = self.bpmmask
            else:
                _mask = self.bpmmask.mask if self.showmask is None or len(self.showmask) == 0 \
                            else self.bpmmask.flagged(flag=self.showmask)
            return ChannelImage(name=name, image=_mask, show_traces=False, show_slits=False,
                                key=key)
        if key == 0:
            # SCIIMG
            return ChannelImage(name=name, image=self.sciimg,
                                cuts=_cuts(self.sciimg, self.img_gpm), key=key)
        if key == 1 and self.skymodel is not None:
            # SKYSUB
            image = (self.sciimg - self.skymodel) * self.model_gpm.astype(float)
            return ChannelImage(name=name, image=image, cuts=_cuts(image, self.model_gpm),
                                show_traces=not self.removetrace, key=key)
        if key == 2 and self.skymodel is not None:
            # SKYRESIDS
            image = (self.sciimg - self.skymodel) * np.sqrt(self.ivarmodel) \
                        * self.model_gpm.astype(float)
            return ChannelImage(name=name, image=image, cuts=(-5.0, 5.0),
                                show_traces=not self.removetrace, key=key)
        if key == 3 and self.objmodel is not None:
            # RESIDS: full model residual map
            image = (self.sciimg - self.skymodel - self.objmodel) * np.sqrt(self.ivarmodel) \
                        * self.img_gpm.astype(float)
            return ChannelImage(name=name, image=image, cuts=(-5.0, 5.0),
                                show_traces=not self.removetrace, key=key)
        return None


def get_detnames(file):
    """
    Return the names of the detectors/mosaics in a spec2d file.

    Parameters
    ----------
    file : :obj:`str`, `Path`_
        spec2d file.

    Returns
    -------
    :obj:`list`
        Detector/mosaic names (e.g., ``['DET01', 'MSC02']``).
    """
    return fits.getval(file, 'HIERARCH ALLSPEC2D_DETS').split(',')


def parse_detname(file, det=None):
    """
    Return the name of a detector in a spec2d file.

    Parameters
    ----------
    file : :obj:`str`, `Path`_
        spec2d file.
    det : :obj:`str`, :obj:`int`, optional
        Detector name (e.g., ``'DET01'`` or ``'MSC01'``) or number.  If None,
        use the first detector in the file.  If a number, the name is
        constructed assuming the reduction is for a single detector.

    Returns
    -------
    :obj:`str`
        Detector name.
    """
    if det is None:
        return get_detnames(file)[0]
    try:
        _det = int(det)
    except (TypeError, ValueError):
        return str(det)
    return DetectorContainer.get_name(_det)


def read_spec1d_traces(file, detname):
    """
    Read the object traces for one detector from a spec1d file.

    Only the traces and the few header keywords needed for display are read,
    which is much faster than reading the full :class:`~pypeit.specobjs.SpecObjs`
    object for files with many objects (e.g., fiber-fed IFUs).

    Parameters
    ----------
    file : :obj:`str`, `Path`_
        spec1d file.
    detname : :obj:`str`
        Name of the detector or mosaic.

    Returns
    -------
    :obj:`tuple`
        Four objects: the traces, shape :math:`(N_{\\rm spec}, N_{\\rm obj})`,
        the object names, and the flags for objects extracted using the
        slit-mask design and manually extracted objects.  All are None if
        there are no objects for this detector.
    """
    traces, names, maskdef_extract, manual_extract = [], [], [], []
    with io.fits_open(file) as hdu:
        for h in hdu[1:]:
            hdr = h.header
            if hdr.get('DMODCLS') != 'SpecObj' or hdr.get('DET') != detname:
                continue
            traces.append(h.data['TRACE_SPAT'])
            name = hdr.get('NAME', h.name)
            objname = hdr.get('MASKDEF_OBJNAME')
            names.append(f'{name}     OBJNAME:{objname}' if objname is not None else name)
            maskdef_extract.append(hdr.get('MASKDEF_EXTRACT') is True)
            manual_extract.append(hdr.get('HAND_EXTRACT_FLAG') is True)
    if len(traces) == 0:
        return None, None, None, None
    return (np.column_stack(traces), np.array(names), np.array(maskdef_extract),
            np.array(manual_extract))


def _cuts(image, gpm, lower=1.0, upper=4.0):
    """Cut levels from sigma-clipped statistics of the good pixels."""
    if not np.any(gpm):
        return None
    mean, med, sigma = sigma_clipped_stats(image[gpm], sigma_lower=5.0, sigma_upper=5.0)
    return (mean - lower*sigma, mean + upper*sigma)


def _old_datamodel_message(file, chk_version):
    """Error/warning message for a file with an incompatible datamodel."""
    try:
        # Try to get the pypeit version used to write this file
        file_pypeit_version = fits.getval(file, 'VERSPYP', 0)
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
    return (f'Your installed version of PypeIt ({__version__}) cannot be used to parse '
            f'{file}, which was reduced using version {file_pypeit_version}.  You '
            'are strongly encouraged to re-reduce your data using this (or, better yet, '
            'the most recent) version of PypeIt.  ' + addendum)


def _read_old_datamodel(file, detname):
    """
    Read the components needed for display directly from the FITS extensions.

    Returns the images, the image and model good-pixel masks, the bitmask
    image, the wavelength image, and the slit data (edges, mask, spat_id,
    maskdef_id, and slit/order IDs).
    """
    slit_left = slit_right = slit_mask = slit_spat_id = slit_mask_id = slit_slid_IDs = None
    with io.fits_open(file) as hdu:
        names = [h.name for h in hdu]
        if not any([detname in n for n in names]):
            raise PypeItError(f'Provided file has no extensions including {detname}.')
        for ext in ['SCIIMG', 'SKYMODEL', 'OBJMODEL', 'IVARMODEL']:
            if f'{detname}-{ext}' not in names:
                raise PypeItError(f'{file} missing extension {detname}-{ext}.')

        sciimg = hdu[f'{detname}-SCIIMG'].data
        skymodel = hdu[f'{detname}-SKYMODEL'].data
        objmodel = hdu[f'{detname}-OBJMODEL'].data
        ivarmodel = hdu[f'{detname}-IVARMODEL'].data
        _ext = f'{detname}-BPMMASK'
        bpmmask = hdu[_ext].data if _ext in names else None
        _ext = f'{detname}-WAVEIMG'
        waveimg = hdu[_ext].data if _ext in names else None

        img_gpm = np.ones(sciimg.shape, dtype=bool) if bpmmask is None else bpmmask == 0

        _ext = f'{detname}-SLITS'
        if _ext not in names:
            log.warning(f'{file} missing extension {_ext}; cannot show slit edges.')
        else:
            data = hdu[_ext].data
            columns = hdu[_ext].columns.names
            slit_spat_id = data['spat_id'] if 'spat_id' in columns else None
            slit_left = data['left_tweak'].T if 'left_tweak' in columns \
                            else (data['left_init'].T if 'left_init' in columns else None)
            slit_right = data['right_tweak'].T if 'right_tweak' in columns \
                            else (data['right_init'].T if 'right_init' in columns else None)
            slit_mask = data['mask'] if 'mask' in columns else None
            slit_mask_id = data['maskdef_id'] if 'maskdef_id' in columns else None
            hdr = hdu[f'{detname}-SCIIMG'].header
            if 'SCI_SPAT_FLEXURE' in hdr and slit_left is not None and slit_right is not None:
                sci_spat_flexure = float(hdr['SCI_SPAT_FLEXURE'])
                slit_left = slit_left + sci_spat_flexure
                slit_right = slit_right + sci_spat_flexure
                log.info(f'Offseting slits by {sci_spat_flexure} pixels.')
            pypeline = hdr['PYPELINE'] if 'PYPELINE' in hdr else None
            if pypeline in ['MultiSlit', 'SlicerIFU', 'Fiber']:
                slit_slid_IDs = slit_spat_id
            elif pypeline == 'Echelle':
                slit_slid_IDs = data['ech_order'] if 'ech_order' in columns else None

    return (sciimg, skymodel, objmodel, ivarmodel, img_gpm, img_gpm.copy(), bpmmask, waveimg,
            (slit_left, slit_right, slit_mask, slit_spat_id, slit_mask_id, slit_slid_IDs))


def load_spec2d_display(file, det=None, channels=None, prefix='', showmask=None,
                        removetrace=False, ignore_extract_mask=False, spat_id=None,
                        maskID=None, try_old=False):
    """
    Read a spec2d file and construct what is needed to display one detector.

    Parameters
    ----------
    file : :obj:`str`, `Path`_
        spec2d file.  The spec1d file, if it exists, is found by replacing
        ``spec2d`` with ``spec1d`` in the file name.
    det : :obj:`str`, :obj:`int`, optional
        Detector to show; see :func:`parse_detname`.
    channels : array-like, optional
        Indices of the channels to show; see :data:`CHANNEL_NAMES`.  If None,
        show all channels.
    prefix : :obj:`str`, optional
        Prefix for the channel names.
    showmask : :obj:`list`, optional
        If not None, include a channel with the mask.  If an empty list, the
        mask bit values are shown; otherwise the image flags the pixels with
        any of the listed mask flags.
    removetrace : :obj:`bool`, optional
        Do not overplot the object traces in the skysub, sky_resid, and resid
        channels.
    ignore_extract_mask : :obj:`bool`, optional
        Ignore the extraction mask when constructing the model images.
    spat_id : :obj:`int`, optional
        Only show the slit with this spatial ID.
    maskID : :obj:`int`, optional
        Only show the slit with this slit-mask ID.
    try_old : :obj:`bool`, optional
        Attempt to read files written with an older datamodel.

    Returns
    -------
    :class:`Spec2DDisplay`
        The images and overlays to display.
    """
    _file = Path(file).absolute()
    detname = parse_detname(_file, det)
    detnames = get_detnames(_file)
    chk_version = not try_old
    show_channels = list(range(len(CHANNEL_NAMES))) if channels is None \
                        else [int(c) for c in channels]
    try:
        pyp_spec = fits.getval(_file, 'PYP_SPEC', 0)
    except (KeyError, OSError):
        pyp_spec = None

    # Try to read the Spec2DObj using the current datamodel, but allowing
    # for the datamodel version to be different
    try:
        spec2DObj = spec2dobj.Spec2DObj.from_file(_file, detname, chk_version=chk_version)
    except (PypeItDataModelError, PypeItBitMaskError):
        message = _old_datamodel_message(_file, chk_version)
        if chk_version:
            raise PypeItError(message)
        log.warning(message)
        spec2DObj = None

    spec1d_file = None
    if spec2DObj is None:
        sciimg, skymodel, objmodel, ivarmodel, img_gpm, model_gpm, bpmmask, waveimg, slits \
                = _read_old_datamodel(_file, detname)
        slit_left, slit_right, slit_mask, slit_spat_id, slit_mask_id, slit_slid_IDs = slits
        # Ignore any spec1d file
    else:
        # Use the parsed SpecObjs object
        sciimg = spec2DObj.sciimg
        skymodel = spec2DObj.skymodel
        objmodel = spec2DObj.objmodel
        ivarmodel = spec2DObj.ivarmodel
        bpmmask = spec2DObj.bpmmask
        waveimg = spec2DObj.waveimg

        img_gpm = spec2DObj.select_flag(invert=True)
        if not np.any(img_gpm):
            log.warning('The full science image is masked!')

        model_gpm = img_gpm.copy()
        if ignore_extract_mask:
            model_gpm |= spec2DObj.select_flag(flag='EXTRACT')

        if spec2DObj.sci_spat_flexure is not None:
            log.info(f'Offseting slits by {spec2DObj.sci_spat_flexure}')
        slit_left, slit_right, slit_mask \
                = spec2DObj.slits.select_edges(flexure=spec2DObj.sci_spat_flexure)
        slit_spat_id = spec2DObj.slits.spat_id
        slit_mask_id = spec2DObj.slits.maskdef_id
        slit_slid_IDs = spec2DObj.slits.slitord_id

        # Object traces from spec1d file
        spec1d_file = _file.with_name(_file.name.replace('spec2d', 'spec1d'))
        if spec1d_file.suffix == '.gz':
            spec1d_file = spec1d_file.with_suffix('')
        if not spec1d_file.is_file():
            log.warning(f'Could not find spec1d file: {spec1d_file}\n'
                        'No objects were extracted.')
            spec1d_file = None

    out = Spec2DDisplay(file=_file, detname=detname, detnames=detnames, pyp_spec=pyp_spec,
                        waveimg=waveimg, prefix=prefix, removetrace=removetrace,
                        showmask=showmask, sciimg=sciimg, skymodel=skymodel,
                        objmodel=objmodel, ivarmodel=ivarmodel, img_gpm=img_gpm,
                        model_gpm=model_gpm, bpmmask=bpmmask)

    # Select the slits to show
    if slit_mask is not None and slit_left is not None and slit_right is not None:
        # TODO: This may be too restrictive, i.e. ignore BADFLTCALIB??
        slit_gpm = slit_mask == 0
        if slit_spat_id is not None and spat_id is not None:
            slit_gpm &= slit_spat_id == spat_id
        # TODO: Should this be 'elif'?  Why not just 'if'?  I.e., are
        # spat_id and maskdef_id mutually exclusive?
        elif slit_mask_id is not None and maskID is not None:
            slit_gpm &= slit_mask_id == maskID
        out.left = slit_left[:, slit_gpm]
        out.right = slit_right[:, slit_gpm]
        out.slit_ids = None if slit_slid_IDs is None else slit_slid_IDs[slit_gpm]
        out.maskdef_ids = None if slit_mask_id is None else slit_mask_id[slit_gpm]

    # Object traces for this detector
    if spec1d_file is not None:
        try:
            out.traces, out.trace_names, out.maskdef_extract, out.manual_extract \
                    = read_spec1d_traces(spec1d_file, detname)
        except (KeyError, ValueError, TypeError):
            # Fall back to reading the full SpecObjs object
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=False)
            in_det = np.where(sobjs.DET == detname)[0]
            if len(in_det) > 0:
                out.traces = np.column_stack([sobjs[k]['TRACE_SPAT'] for k in in_det])
                names = []
                for k in in_det:
                    objname = sobjs[k].MASKDEF_OBJNAME
                    names.append(f'{sobjs[k].NAME}     OBJNAME:{objname}'
                                 if objname is not None else sobjs[k].NAME)
                out.trace_names = np.array(names)
                out.maskdef_extract = np.array([sobjs[k].MASKDEF_EXTRACT is True
                                                for k in in_det])
                out.manual_extract = np.array([sobjs[k].hand_extract_flag is True
                                               for k in in_det])
        if out.traces is None:
            log.warning('spec1d file found, but no objects were extracted for this detector.')

    # The channels to show
    for key in (['mask'] if showmask is not None else []) + show_channels:
        chan = out.make_channel(key)
        if chan is not None:
            out.channels.append(chan)
    return out
