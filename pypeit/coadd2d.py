"""
Module for performing two-dimensional coaddition of spectra.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import dataclasses
from pathlib import Path

from IPython import embed

from astropy.table import Table, vstack
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from scipy import ndimage

from pypeit import log
from pypeit import PypeItError
from pypeit import extraction
from pypeit import find_objects
from pypeit import slittrace
from pypeit import specobjs
from pypeit import utils
from pypeit.images import pypeitimage
from pypeit.core import coadd
from pypeit.core import findobj_skymask
from pypeit.core import parse 
from pypeit.core.wavecal import wvutils
from pypeit.core.moment import moment1d
from pypeit.manual_extract import ManualExtractionObj
from pypeit.spec2dobj import Spec2DObj


@dataclasses.dataclass
class CoAdd2dStack:
    """ CoAdd 2D Stack dataclass

    Attributes
    ----------
    specobjs_list : list
        List of :class:`~pypeit.specobjs.SpecObjs` objects, one per exposure,
        containing the 1D extracted objects associated with the loaded
        ``spec2d`` inputs for the selected detector or mosaic. This list may be
        empty if matching ``spec1d`` files are not available.
    slits_list : list
        List of :class:`~pypeit.slittrace.SlitTraceSet` objects, one per
        exposure, containing the slit or order definitions for the input
        frames.
    slitmask_stack : list
        List of two-dimensional slit-ID images, one per exposure, produced from
        the corresponding :attr:`slits_list` entry using
        :func:`~pypeit.slittrace.SlitTraceSet.slit_img`. Each image labels the
        detector pixels by slit or order, including any recorded spatial
        flexure correction.
    sciimg_stack : list
        List of two-dimensional science images, one per exposure. These are the
        primary input images to be rectified and coadded.
    sciivar_stack : list
        List of two-dimensional inverse-variance images, one per exposure,
        aligned with :attr:`sciimg_stack`. If the exposure times differ among
        the input frames, these images are rescaled consistently with the
        science images before coaddition.
    skymodel_stack : list
        List of two-dimensional sky-model images, one per exposure, aligned
        with :attr:`sciimg_stack`. These are used when constructing sky-
        subtracted images and are also rescaled if the exposure times are
        homogenized.
    mask_stack : list
        List of bad-pixel mask images, one per exposure, containing the raw
        mask values from the input :class:`~pypeit.spec2dobj.Spec2DObj`
        instances.
    waveimg_stack : list
        List of two-dimensional wavelength images, one per exposure, giving the
        wavelength associated with each detector pixel in the corresponding
        science image.
    exptime_stack : list
        List of exposure times, one per input exposure, read from the primary
        headers of the input ``spec2d`` data.
    exptime_coadd : float
        Reference exposure time adopted for the coadd. This is taken to be the
        median input exposure time, using the higher of the two middle values
        for even-length stacks. If the individual exposure times differ by more
        than the allowed tolerance, the science, sky, and inverse-variance
        images are rescaled to this effective exposure time before coaddition.
    redux_path : :obj:`~pathlib.Path`
        Reduction path associated with the current coaddition context. At
        present this is set to the current working directory when the stack is
        constructed.
    detectors : list
        List of detector or detector-mosaic metadata objects, one per exposure,
        copied from the loaded :class:`~pypeit.spec2dobj.Spec2DObj` instances.
        These are later used for bookkeeping and for propagating detector
        properties into the coadded products.
    spectrograph : str
        Name of the spectrograph associated with the stack.
    pypeline : str
        Name of the PypeIt reduction mode associated with the stack, e.g.
        multislit or echelle.
    maskdef_designtab_list : list
        List of slitmask-design tables, one per exposure, taken from the input
        ``spec2d`` objects. These are used when propagating mask-design
        metadata into the coadded products when such information is available.
    spat_flexure_list : list
        List of spatial-flexure corrections, one per exposure, used when
        constructing slit masks and when relating slit geometry in the input
        frames to the coadded frame.
    """
    specobjs_list:list
    slits_list:list
    slitmask_stack:list
    sciimg_stack:list
    sciivar_stack:list
    skymodel_stack:list
    mask_stack:list
    waveimg_stack:list
    exptime_stack:list
    exptime_coadd:float
    redux_path:Path
    detectors:list
    spectrograph:str
    pypeline:str
    maskdef_designtab_list:list
    spat_flexure_list:list


#TODO We should decide which parameters go in through the parset 
# and which parameters are passed in to the method as arguments
class CoAdd2D:

    """
    Main driver for two-dimensional spectral coaddition.

    This class coordinates loading a set of reduced ``spec2d`` products, selecting
    the slits or orders to combine, determining relative offsets and weights among
    the input exposures, rectifying the input images onto a common coordinate grid,
    and reducing the resulting pseudo-image into final coadded products.

    The base class provides the shared machinery used by both multislit and
    echelle coadds. Subclasses specialize the handling of reference objects,
    offset determination, weighting, and wavelength-grid construction for their
    respective reduction modes.

    Notes
    -----
    The overall workflow is:

        #. Load the per-exposure data products into a :class:`~pypeit.core.coadd2d.CoAdd2dStack`.
        #. Determine the subset of slits or orders that are valid for coaddition.
        #. Compute relative offsets among the exposures.
        #. Compute the weights used during coaddition.
        #. Rectify and coadd the input images slit-by-slit or order-by-order.
        #. Build a pseudo-image from the coadd outputs.
        #. Run the extraction and bookkeeping steps needed to write the coadded products to disk.

    See Also
    --------
    MultiSlitCoAdd2D
        Multislit and longslit implementation.
    EchelleCoAdd2D
        Echelle implementation.
    """
    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec2dfiles, spectrograph, par, det=1,
                     only_slits=None, exclude_slits=None,
                     sn_smooth_npix=None, bkg_redux=False, find_negative=False, show=False,
                     show_peaks=False, debug_offsets=False, debug=False):
        """
        Instantiate the appropriate :class:`CoAdd2D` subclass for a spectrograph.

        The subclass is selected by matching the spectrograph ``pypeline`` name to a
        child class of :class:`CoAdd2D`.

        Parameters
        ----------
        spec2dfiles : list
            List of input ``spec2d`` files or already-instantiated
            :class:`~pypeit.spec2dobj.Spec2DObj` objects.
        spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            Spectrograph object describing the instrument and reduction mode.
        par : :class:`~pypeit.par.parset.ParSet`
            Parameter set controlling the coaddition.
        det : int or tuple, optional
            Detector number or detector mosaic identifier to process.
        only_slits : list, optional
            List of slit or order identifiers to include.
        exclude_slits : list, optional
            List of slit or order identifiers to exclude.
        sn_smooth_npix : int, optional
            Number of pixels used when smoothing S/N estimates for weight
            determination.
        bkg_redux : bool, optional
            Flag indicating that the science images have already been background
            subtracted.
        find_negative : bool, optional
            If True, search for and mask negative traces produced by image
            differencing.
        show : bool, optional
            If True, display intermediate results.
        show_peaks : bool, optional
            If True, show QA plots for peak finding.
        debug_offsets : bool, optional
            If True, show QA plots related to offset determination.
        debug : bool, optional
            If True, enable additional debug output.

        Returns
        -------
        :class:`CoAdd2D`
            Instance of the subclass appropriate for the requested reduction mode.
        """

        return next(c for c in cls.__subclasses__()
                    if c.__name__ == (spectrograph.pypeline + 'CoAdd2D'))(
                        spec2dfiles, spectrograph, par, det=det, 
                        only_slits=only_slits, exclude_slits=exclude_slits,
                        sn_smooth_npix=sn_smooth_npix, bkg_redux=bkg_redux, find_negative=find_negative,
                        show=show, show_peaks=show_peaks, debug_offsets=debug_offsets, debug=debug)

    def __init__(self, spec2d, spectrograph, par, det=1, 
                 only_slits=None, exclude_slits=None, 
                 sn_smooth_npix=None, bkg_redux=False, find_negative=False, show=False,
                 show_peaks=False, debug_offsets=False, debug=False):
        """
        Initialize a two-dimensional coadd driver.

        Parameters
        ----------
        spec2d : list
            List of input ``spec2d`` files or
            :class:`~pypeit.spec2dobj.Spec2DObj` objects to coadd.
        spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            Spectrograph object describing the instrument and reduction mode.
        par : :class:`~pypeit.par.parset.ParSet`
            Parameter set controlling the coaddition.
        det : int or tuple, optional
            One-indexed detector number or detector-mosaic identifier to process.
        only_slits : list, optional
            List of slit or order identifiers to include in the coadd.
        exclude_slits : list, optional
            List of slit or order identifiers to exclude from the coadd.
        sn_smooth_npix : int, optional
            Number of pixels used when smoothing S/N estimates for weight
            determination. If None, a default based on the image size is used.
            TODO: for truncated echelle orders we should be doing something more
            intelligent.
        bkg_redux : bool, optional
            If True, the science images have already been background subtracted.
        find_negative : bool, optional
            If True, search for and mask negative traces produced by differenced
            science frames.
        show : bool, optional
            If True, display intermediate results in Ginga.
        show_peaks : bool, optional
            If True, show QA plots for object-finding peaks.
        debug_offsets : bool, optional
            If True, display QA related to automatic offset determination.
        debug : bool, optional
            If True, enable additional debug output.

        Raises
        ------
        PypeItError
            Raised if the input stack is inconsistent, e.g. if the exposures do not
            contain the same number of slits or have incompatible binning.
        """
        # Use Cases:
        # offsets
        #    1) offsets = 'auto' -- auto compute offsets from brightest object (if exists)
        #    2) offsets not 'auto' (i.e. a list) - use them
        #    -------------- only for Multislit --------------
        #    3) offsets = 'maskdef_offsets' - use `maskdef_offset` saved in SlitTraceSet
        #    4) offsets = 'header' - use the dither offsets recorded in the header
        # ===============================================================================
        # weights
        #    1) weights = 'auto' -- if brightest object exists auto compute weights,
        #                           otherwise use uniform weights
        #    2) weights = 'uniform' -- use uniform weights
        #    3) weights is a list - use them

        self.spec2d = spec2d
        self.spectrograph = spectrograph
        self.par = par

        # This can be a single integer for a single detector or a tuple for
        # multiple detectors placed in a mosaic.
        self.det = det

        # This is the string name of the detector or mosaic used when saving the
        # processed data to PypeIt's main output files
        self.detname = self.spectrograph.get_det_name(self.det)

        self.bkg_redux = bkg_redux
        self.find_negative = find_negative
        self.show = show
        self.show_peaks = show_peaks
        self.debug_offsets = debug_offsets
        self.debug = debug
        self.offsets = None
        self.coadd2d_stack = None
        self.pseudo_dict = None

        # Brightest object attributes used for both MultislitCoAdd2D and EchelleCoAdd2D
        # Array with shape = (nexp,) containing spat_pixpos_id (MultiSlit) or
        # ech_fracpos_id (Echelle) of the brightest object in each exposure
        self.obj_id_bri = None
        # Array with shape = (nexp,) containing the S/N of the brightest object in each exposure
        self.snr_bar_bri = None

        # This is a list of length self.nexp that is assigned by the compute_weights method
        self.use_weights = None
        self.wave_grid = None
        self.good_slits = None
        self.maskdef_offset = None

        # Load the stack_dict
        self.coadd2d_stack = self.load_coadd2d_stacks(self.spec2d)
        self.pypeline = self.spectrograph.pypeline

        self.nexp = len(self.spec2d)

        # Check that there are the same number of slits on every exposure
        nslits_list = [slits.nslits for slits in self.coadd2d_stack.slits_list]
        if not len(set(nslits_list)) == 1:
            raise PypeItError(
                'Not all of your files have the same number of slits. Check your inputs'
            )
        # This is the number of slits of the single (un-coadded) frames
        self.nslits_single = nslits_list[0]

        # Check that nspec is the same for all the exposures
        self.nspec_array = np.array([slits.nspec for slits in self.coadd2d_stack.slits_list])
        self.nspec_max = self.nspec_array.max()

        # Check that binning is the same for all the exposures
        binspec_list = [slits.binspec for slits in self.coadd2d_stack.slits_list]
        binspat_list = [slits.binspat for slits in self.coadd2d_stack.slits_list]
        if not len(set(binspec_list)) == 1:
            raise PypeItError(
                'Not all of your files have the same spectral binning. Check your inputs'
            )
        if not len(set(binspat_list)) == 1:
            raise PypeItError(
                'Not all of your files have the same spatial binning. Check your inputs'
            )
        self.binning = np.array([self.coadd2d_stack.slits_list[0].binspec,
                                 self.coadd2d_stack.slits_list[0].binspat])

        self.spat_ids = self.coadd2d_stack.slits_list[0].spat_id

        # If smoothing is not input, smooth by 10% of the maximum spectral dimension
        self.sn_smooth_npix = sn_smooth_npix if sn_smooth_npix is not None else 0.1*self.nspec_max

        # coadded frame parameters

        # get slit index that indicates which slits are good for coadding
        self.good_slits = self.good_slitindx(only_slits=only_slits, exclude_slits=exclude_slits)
        # get the number of slits that are going to be coadded
        self.nslits_coadded = self.good_slits.size

        # effective exposure time
        self.exptime_coadd = self.coadd2d_stack.exptime_coadd
        
        # define the wavelength grid for the 2d coadd
        self.wave_grid, self.wave_grid_mid, self.dsamp = self.get_wave_grid()
        
        
        # Handle the reference object
        self.handle_reference_obj()

        # get self.use_weights
        self.compute_weights()

        # get self.offsets
        self.compute_offsets()




    @staticmethod
    def default_par(spectrograph, inp_cfg=None, det=None, only_slits=None, exclude_slits=None):
        """
        Construct the default parameter set for two-dimensional coaddition.

        Parameters
        ----------
        spectrograph : str
            PypeIt spectrograph name.
        inp_cfg : dict, optional
            Existing configuration dictionary to update.
        det : list, str, or tuple, optional
            Detector or mosaic selection to place in the returned configuration.
        only_slits : list or str, optional
            Slit or order identifiers to include in the coadd.
        exclude_slits : list or str, optional
            Slit or order identifiers to exclude from the coadd.

        Returns
        -------
        dict
            Dictionary with the default coadd2d configuration entries.
        """
        cfg = dict(rdx=dict(spectrograph=spectrograph))
        if inp_cfg is not None:
            cfg = utils.recursive_update(cfg, dict(inp_cfg))
        if only_slits is not None and det is not None:
            log.warning('only_slits and det are mutually exclusive. Ignoring det.')
            _det = None
        else:
            _det = det

        if det is not None:
            cfg['rdx']['detnum'] = _det

        if only_slits is not None and exclude_slits is not None:
            log.warning('only_slits and exclude_slits are mutually exclusive. Ignoring exclude_slits.')
            _exclude_slits = None
        else:
            _exclude_slits = exclude_slits

        if only_slits is not None:
            utils.add_sub_dict(cfg, 'coadd2d')
            cfg['coadd2d']['only_slits'] = only_slits
        if _exclude_slits is not None:
            utils.add_sub_dict(cfg, 'coadd2d')
            cfg['coadd2d']['exclude_slits'] = _exclude_slits
        # TODO: Heliocentric for coadd2d needs to be thought through. Currently
        # turning it off.
        utils.add_sub_dict(cfg, 'calibrations')
        utils.add_sub_dict(cfg['calibrations'], 'wavelengths')
        cfg['calibrations']['wavelengths']['refframe'] = 'observed'
        # TODO: Flexure correction for coadd2d needs to be thought through.
        # Currently turning it off.
        utils.add_sub_dict(cfg, 'flexure')
        cfg['flexure']['spec_method'] = 'skip'
        # TODO: This is currently the default for 2d coadds, but we need a way
        # to toggle it on/off
        utils.add_sub_dict(cfg, 'reduce')
        utils.add_sub_dict(cfg['reduce'], 'findobj')
        cfg['reduce']['findobj']['skip_skysub'] = True

        return cfg

    @staticmethod
    def default_basename(spec2d_files):
        """
        Construct the default basename for a coadded ``spec2d`` product.

        Parameters
        ----------
        spec2d_files : list
            List of input ``spec2d`` filenames.

        Returns
        -------
        str
            Basename used when writing the coadded output products.
        """
        # Get the output basename
        frsthdr = fits.getheader(spec2d_files[0])
        lasthdr = fits.getheader(spec2d_files[-1])
        if 'FILENAME' not in frsthdr:
            raise PypeItError(f'Missing FILENAME keyword in {spec2d_files[0]}.  Set the basename '
                        'using the command-line option.')
        if 'FILENAME' not in lasthdr:
            raise PypeItError(f'Missing FILENAME keyword in {spec2d_files[-1]}.  Set the basename '
                        'using the command-line option.')
        if 'TARGET' not in frsthdr:
            raise PypeItError(f'Missing TARGET keyword in {spec2d_files[0]}.  Set the basename '
                        'using the command-line option.')
        return f"{frsthdr['FILENAME'].split('.fits')[0]}-" \
                f"{lasthdr['FILENAME'].split('.fits')[0]}-{frsthdr['TARGET'].replace(' ','')}"

    @staticmethod
    def output_paths(spec2d_files, par, coadd_dir=None):
        """
        Construct and create the science and QA output directories for coadd2d.

        Parameters
        ----------
        spec2d_files : list
            List of input ``spec2d`` filenames. The parent reduction directory is
            inferred from these paths if ``coadd_dir`` is not provided.
        par : :class:`~pypeit.par.pypeitpar.PypeItPar`
            Full parameter set. The routine uses and may update the values in
            ``par['rdx']['scidir']`` and ``par['rdx']['qadir']``.
        coadd_dir : str, optional
            Root directory for the coadd2d output. If None, the parent reduction
            directory inferred from ``spec2d_files`` is used.

        Returns
        -------
        sci_output_dir : str
            Path to the science output directory
        qa_output_dir : str
            Path to the QA output directory.

        Notes
        -----
        The required directories are created if they do not already exist.
        """
        # Science output directory
        if coadd_dir is not None:
            pypeit_scidir = Path(coadd_dir).absolute() / 'Science'
        else:
            pypeit_scidir = Path(spec2d_files[0]).parent
        coadd_scidir = pypeit_scidir.parent / f"{par['rdx']['scidir']}_coadd"
        if not coadd_scidir.exists():
            coadd_scidir.mkdir(parents=True)
        # QA directory
        par['rdx']['qadir'] += '_coadd'
        qa_path = pypeit_scidir.parent / par['rdx']['qadir'] / 'PNGs'
        if not qa_path.exists():
            qa_path.mkdir(parents=True)
        return str(coadd_scidir), str(qa_path)

    def good_slitindx(self, only_slits=None, exclude_slits=None):
        """
        Determine which slits or orders are eligible for two-dimensional coaddition.

        A slit or order is considered good only if it is not masked by the common
        bitmask across the full set of uncoadded exposures and it satisfies the
        user-supplied inclusion or exclusion criteria.

        Parameters
        ----------
        only_slits : list, optional
            Slit or order identifiers to include. If provided, only these identifiers
            are considered good.
        exclude_slits : list, optional
            Slit or order identifiers to exclude. Ignored if ``only_slits`` is
            provided.

        Returns
        -------
        numpy.ndarray
            One-dimensional array of indices selecting the good slits or orders in the
            uncoadded stack.
        """

        if exclude_slits is not None and only_slits is not None:
            log.warning('Both `only_slits` and `exclude_slits` are provided. They are mutually exclusive. '
                      'Using `only_slits` and ignoring `exclude_slits`')
            _exclude_slits = None
        else:
            _exclude_slits = exclude_slits

        # This creates a unified bpm common to all frames
        slits0 = self.coadd2d_stack.slits_list[0]
        # bpm for the first frame

        reduce_bpm = slits0.bitmask.flagged(slits0.mask,
                                            and_not=slits0.bitmask.exclude_for_reducing)
        for i in range(1, self.nexp):
            # update bpm with the info from the other frames
            slits = self.coadd2d_stack.slits_list[i]
            reduce_bpm |= slits.bitmask.flagged(slits.mask,
                                                and_not=slits.bitmask.exclude_for_reducing)
        # these are the good slit index according to the bpm mask
        good_slitindx = np.where(np.logical_not(reduce_bpm))[0]

        # If we want to coadd all the good slits
        if only_slits is None and _exclude_slits is None:
            return good_slitindx

        # If instead we want to coadd only a selected (by the user) number of slits
        if only_slits is not None:
            # these are the `slitord_id` of the slits that we want to coadd
            _only_slits = np.atleast_1d(only_slits)
            # create an array of slit index that are selected by the user and are also good slits
            good_onlyslits = np.array([], dtype=int)
            log.info('Coadding only the following slits:')
            for islit in _only_slits:
                if islit not in slits0.slitord_id[good_slitindx]:
                    # Warnings for the slits that are selected by the user but NOT good slits
                    log.warning('Slit {} cannot be coadd because masked'.format(islit))
                else:
                    log.info(f'Slit {islit}')
                    indx = np.where(slits0.slitord_id[good_slitindx] == islit)[0]
                    good_onlyslits = np.append(good_onlyslits, good_slitindx[indx])
            return good_onlyslits

        # if we want to exclude some slits (selected by the user) from coadding
        # these are the `slitord_id` of the slits that we want to exclude
        _exclude_slits = np.atleast_1d(_exclude_slits)
        # create an array of slit index that are excluded by the user
        exclude_slitindx = np.array([], dtype=int)
        log.info('Excluding the following slits:')
        for islit in _exclude_slits:
            if islit in slits0.slitord_id[good_slitindx]:
                log.info(f'Slit {islit}')
                exclude_slitindx = np.append(exclude_slitindx,
                                             np.where(slits0.slitord_id[good_slitindx] == islit)[0][0])
        # these are the good slit index excluding the slits that are selected by the user
        return np.delete(good_slitindx, exclude_slitindx)

    def optimal_weights(self, uniq_obj_id, order=None, weight_method='auto'):
        """
        Compute the optimal exposure weights for a set of reference objects.

        The routine extracts the relevant one-dimensional spectra for the supplied
        reference-object identifiers and passes them to
        :func:`pypeit.core.coadd.sn_weights` to determine the weighting of each input
        exposure.

        Parameters
        ----------
        uniq_obj_id : numpy.ndarray
            Array of unique object identifiers, one per exposure, used to select the
            reference object for the weight calculation.
        order : int, optional
            Echelle order to use when computing the weights. Ignored for multislit
            reductions.
        weight_method : {'auto', 'constant', 'uniform', 'wave_dependent', 'relative', 'ivar'}, optional
            Weighting algorithm passed to :func:`pypeit.core.coadd.sn_weights`.

        Returns
        -------
        rms_sn : numpy.ndarray
            Array of root-mean-square S/N value for each input spectra. Shape = (nexp,)
        weights : list
            List of  len(nexp) containing the signal-to-noise squared weights to be
            applied to the spectra. This output is aligned with the vector (or
            vectors) provided in waves which is read in by this routine, i.e. it is a
            list of arrays of type `numpy.ndarray`_  with the same shape as those in waves.

        Notes
        -----
        This method is used by both the multislit and echelle subclasses when the
        weights are determined from the brightest detected reference object.
        """

        # Grab the traces, flux, wavelength and noise for this uniq_obj_id.
        waves, fluxes, ivars, gpms = [], [], [], []

        for iexp, sobjs in enumerate(self.coadd2d_stack.specobjs_list):
            ithis = sobjs.slitorder_uniq_id_indices(uniq_obj_id[iexp], order=order)
            if not np.any(ithis):
                raise PypeItError(f'Object {uniq_obj_id[iexp]} provided not valid. Optimal weights cannot be determined.')
            order_str = f' on slit/order {order}' if order is not None else ''
            # check if OPT_COUNTS is available
            if sobjs[ithis][0].has_opt_ext() and np.any(sobjs[ithis][0].OPT_MASK):
                wave_iexp, flux_iexp, ivar_iexp, gpm_iexp = sobjs[ithis][0].get_opt_ext()
                waves.append(wave_iexp)
                fluxes.append(flux_iexp)
                ivars.append(ivar_iexp)
                gpms.append(gpm_iexp)
            # check if BOX_COUNTS is available
            elif sobjs[ithis][0].has_box_ext() and np.any(sobjs[ithis][0].BOX_MASK):
                wave_iexp, flux_iexp, ivar_iexp, gpm_iexp = sobjs[ithis][0].get_box_ext()
                waves.append(wave_iexp)
                fluxes.append(flux_iexp)
                ivars.append(ivar_iexp)
                gpms.append(gpm_iexp)

                log.warning(
                    f'Optimal extraction not available for object {uniq_obj_id[iexp]} '
                    f'{order_str} in file {iexp}. Using box extraction.'
                )
            else:
                raise PypeItError(
                    f'Optimal weights cannot be determined because flux not available for object '
                    f'= {uniq_obj_id[iexp]} {order_str} in file {iexp}. '
                )

        # TODO For now just use the zero as the reference for the wavelengths? Perhaps we should be rebinning the data though?
        rms_sn, weights = coadd.sn_weights(fluxes, ivars, gpms, sn_smooth_npix=self.sn_smooth_npix, weight_method=weight_method)
        return rms_sn, weights

    def coadd(self, interp_dspat=True):
        """
        Perform the two-dimensional coaddition.

        This routine loops over the selected good slits or orders, rectifies the
        input images for each one onto a common coordinate system, and combines the
        resulting rebinned images using the previously determined offsets and weights.

        Parameters
        ----------
        interp_dspat : bool, optional
            If True, interpolate the spatial sampling when constructing the coadd
            grids.

        Returns
        -------
        list
            List of per-slit or per-order coadd products. Each element contains the
            data needed by :meth:`create_pseudo_image` to build the final pseudo-image.

        Notes
        -----
        The exact contents of the returned list depend on the pypeline-specific
        coaddition implementation, but each entry represents one coadded slit or
        order.
        """

        coadd_list = []
        for slit_idx in self.good_slits:
            _slitord_id = self.coadd2d_stack.slits_list[0].slitord_id
            log.info(f'Performing 2D coadd for slit/order {_slitord_id[slit_idx]} ({slit_idx + 1}/{self.nslits_single})')

            # mask identifying the current slit in each exposure
            thismask_stack = [np.abs(slitmask - self.spat_ids[slit_idx]) <= self.par['coadd2d']['spat_toler']
                              for slitmask in self.coadd2d_stack.slitmask_stack]

            # check if the slit is found in every exposure
            if not np.all([np.any(thismask) for thismask in thismask_stack]):
                log.warning(
                    f'Slit/order {_slitord_id[slit_idx]} was not found in every file. 2D coadd '
                    'cannot be performed on this slit. Try increasing the parameter spat_toler'
                )
                continue

            # reference trace
            ref_trace_stack = self.reference_trace_stack(slit_idx, offsets=self.offsets, uniq_obj_id=self.obj_id_bri)

            # Perform the 2d coadd
            # NOTE: mask_stack is a gpm, and this is called inmask_stack in
            # compute_coadd2d, and outmask in coadd_dict is also a gpm
            mask_stack = [mask == 0 for mask in self.coadd2d_stack.mask_stack]
            coadd_dict = coadd.compute_coadd2d(ref_trace_stack, self.coadd2d_stack.sciimg_stack,
                                               self.coadd2d_stack.sciivar_stack,
                                               self.coadd2d_stack.skymodel_stack,
                                               mask_stack,
                                               thismask_stack,
                                               self.coadd2d_stack.waveimg_stack,
                                               self.wave_grid, self.par['coadd2d']['spat_samp_fact'],
                                               maskdef_dict=self.get_maskdef_dict(slit_idx, ref_trace_stack),
                                               weights=self._get_weights(indx=slit_idx), interp_dspat=interp_dspat)
            coadd_list.append(coadd_dict)

        if len(coadd_list) == 0:
            raise PypeItError(
                'All the slits were missing in one or more files. 2D coadd cannot be performed.'
            )

        return coadd_list

    def create_pseudo_image(self, coadd_list):
        """
        Assemble the per-slit or per-order coadds into a pseudo-image.

        Parameters
        ----------
        coadd_list : list
            List of coadded slit or order products produced by :meth:`coadd`.

        Returns
        -------
        dict
            Dictionary containing the pseudo-image and associated metadata needed for
            downstream reduction and output writing.

        Notes
        -----
        The pseudo-image is the rectified, stacked image representation consumed by
        the later extraction steps.

        .. todo::

            THIS UNDOCUMENTED CODE PROBABLY SHOULD GENERATE AND RETURN
            STANDARD PYPEIT OBJCTS INSTEAD OF SOME UNDEFINED DICT
        """

        # Check that self.nslit is equal to len(coadd_list)
        if self.nslits_coadded != len(coadd_list):
            raise PypeItError('Wrong number of slits for the 2d coadded frame')

        nspec_vec = np.zeros(self.nslits_coadded,dtype=int)
        nspat_vec = np.zeros(self.nslits_coadded,dtype=int)
        for islit, cdict in enumerate(coadd_list):
            nspec_vec[islit]=cdict['nspec']
            nspat_vec[islit]=cdict['nspat']

        # Determine the size of the pseudo image
        nspat_pad = 10
        nspec_pseudo = nspec_vec.max()
        nspat_pseudo = int(np.sum(nspat_vec) + (self.nslits_coadded + 1)*nspat_pad)  # Cast for SlitTraceSet
        spec_vec_pseudo = np.arange(nspec_pseudo)
        shape_pseudo = (nspec_pseudo, nspat_pseudo)
        imgminsky_pseudo = np.zeros(shape_pseudo)
        sciivar_pseudo = np.zeros(shape_pseudo)
        waveimg_pseudo = np.zeros(shape_pseudo)
        waveimg_mid_pseudo = np.zeros(shape_pseudo)
        tilts_pseudo = np.zeros(shape_pseudo)
        spat_img_pseudo = np.zeros(shape_pseudo)
        nused_pseudo = np.zeros(shape_pseudo, dtype=int)
        inmask_pseudo = np.zeros(shape_pseudo, dtype=bool)
        wave_mid = np.zeros((nspec_pseudo, self.nslits_coadded))
        wave_mask = np.zeros((nspec_pseudo, self.nslits_coadded),dtype=bool)
        wave_min = np.zeros((nspec_pseudo, self.nslits_coadded))
        wave_max = np.zeros((nspec_pseudo, self.nslits_coadded))
        dspat_mid = np.zeros((nspat_pseudo, self.nslits_coadded))

        spat_left = nspat_pad
        slit_left = np.zeros((nspec_pseudo, self.nslits_coadded))
        slit_righ = np.zeros((nspec_pseudo, self.nslits_coadded))
        spec_min1 = np.zeros(self.nslits_coadded)
        spec_max1 = np.zeros(self.nslits_coadded)

        # maskdef info
        all_maskdef_ids = np.array([cc['maskdef_id'] for cc in coadd_list])
        if None not in all_maskdef_ids:
            maskdef_id = np.zeros(self.nslits_coadded, dtype=int)
            maskdef_objpos = np.zeros(self.nslits_coadded)
            maskdef_slitcen = np.zeros((nspec_pseudo, self.nslits_coadded))
            maskdef_designtab = Table()
        else:
            maskdef_id = None
            maskdef_objpos = None
            maskdef_slitcen = None
            maskdef_designtab = None

        nspec_grid = self.wave_grid_mid.size
        for islit, coadd_dict in enumerate(coadd_list):
            spat_righ = spat_left + nspat_vec[islit]
            ispec = slice(0,nspec_vec[islit])
            ispat = slice(spat_left,spat_righ)
            imgminsky_pseudo[ispec, ispat] = coadd_dict['imgminsky']
            sciivar_pseudo[ispec, ispat] = coadd_dict['sciivar']
            waveimg_pseudo[ispec, ispat] = coadd_dict['waveimg']
            # NOTE: inmask is a gpm
            inmask_pseudo[ispec, ispat] = coadd_dict['outmask']
            image_temp = (coadd_dict['dspat'] - coadd_dict['dspat_mid'][0] + spat_left) #*coadd_dict['outmask']
            # spat_img_pseudo is the sub-pixel image position on the rebinned pseudo image
            spat_img_pseudo[ispec, ispat] = image_temp
            nused_pseudo[ispec, ispat] = coadd_dict['nused']
            wave_min[ispec, islit] = coadd_dict['wave_min']
            wave_max[ispec, islit] = coadd_dict['wave_max']
            wave_mid[ispec, islit] = coadd_dict['wave_mid']
            # waveimg_mid_pseudo image containing the bin centers that the data was rebinned onto
            waveimg_mid_pseudo[ispec, ispat] = np.repeat(wave_mid[ispec, islit][:, np.newaxis], nspat_vec[islit], axis=1)
            # Patch locations where the waveimg is zero with the midpoints of the grid. This prevents discontinuities
            # in the wavelength image. This means howver that the 2d wavelength image has wavelengths with
            # two different meanings, i.e. where unmasked they are averaged rebinned wavelengths, but where masked
            # it is the original grid.
            # TODO THink about whether we should just use the fixed grid wavelengths throughout as the waveimg rather than
            # have this hybrid defintion.
            waveimg_pseudo[ispec, ispat][np.logical_not(inmask_pseudo[ispec, ispat])] = \
                waveimg_mid_pseudo[ispec, ispat][np.logical_not(inmask_pseudo[ispec, ispat])]
            wave_mask[ispec, islit] = True
            tilts_pseudo[ispec, ispat] = (waveimg_pseudo[ispec, ispat] - coadd_dict['wave_min'][0])/(coadd_dict['wave_max'][-1] - coadd_dict['wave_min'][0])

            # Fill in the rest of the wave_mid with the corresponding points in the wave_grid
            #wave_this = wave_mid[wave_mask[:,islit], islit]
            #ind_upper = np.argmin(np.abs(self.wave_grid_mid - wave_this.max())) + 1
            #if nspec_vec[islit] != nspec_pseudo:
            #    wave_mid[nspec_vec[islit]:, islit] = self.wave_grid_mid[ind_upper:ind_upper + (nspec_pseudo-nspec_vec[islit])]


            dspat_mid[ispat, islit] = coadd_dict['dspat_mid']
            slit_left[:,islit] = np.full(nspec_pseudo, spat_left)
            slit_righ[:,islit] = np.full(nspec_pseudo, spat_righ)
            spec_max1[islit] = nspec_vec[islit]-1
            spat_left = spat_righ + nspat_pad

            # maskdef info
            if None not in all_maskdef_ids:
                maskdef_id[islit] = coadd_dict['maskdef_id']
                maskdef_objpos[islit] = coadd_dict['maskdef_objpos']
                maskdef_slitcen[:, islit] = np.full(nspec_pseudo, coadd_dict['maskdef_slitcen'] + slit_left[:,islit])
                if coadd_dict['maskdef_designtab'] is not None:
                    maskdef_designtab = vstack([maskdef_designtab, coadd_dict['maskdef_designtab']])

        slits_pseudo \
                = slittrace.SlitTraceSet(slit_left, slit_righ, self.pypeline, detname=self.detname,
                                         nspat=nspat_pseudo, PYP_SPEC=self.spectrograph.name,
                                         specmin=spec_min1, specmax=spec_max1,
                                         maskdef_id=maskdef_id, maskdef_objpos=maskdef_objpos,
                                         maskdef_offset=0., maskdef_slitcen=maskdef_slitcen,
                                         maskdef_designtab=maskdef_designtab)

        # change value of spat_id in maskdef_designtab
        # needs to be done here because spat_id is computed in slittrace
        if maskdef_designtab is not None:
            slits_pseudo.maskdef_designtab['SPAT_ID'] = slits_pseudo.spat_id

        # assign ech_order if exist
        slits_pseudo.ech_order = self.coadd2d_stack.slits_list[0].ech_order[self.good_slits] \
            if self.coadd2d_stack.slits_list[0].ech_order is not None else None
        slitmask_pseudo = slits_pseudo.slit_img()
        # This is a kludge to deal with cases where bad wavelengths result in large regions where the slit is poorly sampled,
        # which wreaks havoc on the local sky-subtraction
        min_slit_frac = 0.70
        spec_min = np.zeros(self.nslits_coadded)
        spec_max = np.zeros(self.nslits_coadded)
        for islit in range(self.nslits_coadded):
            spat_id = slits_pseudo.spat_id[islit]
            slit_width = np.sum(inmask_pseudo & (slitmask_pseudo == spat_id), axis=1)
            slit_width_img = np.outer(slit_width, np.ones(nspat_pseudo))
            med_slit_width = np.median(slit_width_img[slitmask_pseudo == spat_id])
            # TODO -- need inline docs
            nspec_eff = np.sum(slit_width > min_slit_frac*med_slit_width)
            nsmooth = int(np.fmax(np.ceil(nspec_eff*0.02),10))
            slit_width_sm = ndimage.filters.median_filter(slit_width, size=nsmooth, mode='reflect')
            igood = (slit_width_sm > min_slit_frac*med_slit_width)
            # TODO -- need inline docs
            spec_min[islit] = spec_vec_pseudo[igood].min()
            spec_max[islit] = spec_vec_pseudo[igood].max()
            bad_pix = (slit_width_img < min_slit_frac*med_slit_width) & (slitmask_pseudo == spat_id)
            inmask_pseudo[bad_pix] = False

        # Update slits_pseudo
        slits_pseudo.specmin = spec_min
        slits_pseudo.specmax = spec_max

        return dict(nspec=nspec_pseudo, nspat=nspat_pseudo, imgminsky=imgminsky_pseudo,
                    sciivar=sciivar_pseudo, inmask=inmask_pseudo, tilts=tilts_pseudo,
                    waveimg=waveimg_pseudo, waveimg_mid=waveimg_mid_pseudo, spat_img=spat_img_pseudo, slits=slits_pseudo,
                    wave_mask=wave_mask, wave_mid=wave_mid, wave_min=wave_min, wave_max=wave_max)

    def reduce(self, pseudo_dict, show=False, clear_ginga=True, show_peaks=False, show_skysub_fit=False, basename=None):
        """
        Reduce a coadded pseudo-image into final science products.

        Parameters
        ----------
        pseudo_dict : dict
            Dictionary describing the pseudo-image and associated metadata, typically
            returned by :meth:`create_pseudo_image`.
        show : bool, optional
            If True, display intermediate reduction results.
        clear_ginga : bool, optional
            If True, clear the Ginga display before showing new content.
        show_peaks : bool, optional
            If True, show QA plots related to object-finding peaks.
        show_skysub_fit : bool, optional
            If True, show QA for the sky-subtraction fit.
        basename : str, optional
            Basename to use when writing the output products.

        Returns
        -------
        tuple
            The reduced data products generated from the pseudo-image.

        Notes
        -----
        This method performs the coadd-specific analogue of the standard PypeIt
        extraction sequence for a single frame.
        """

        show = self.show if show is None else show
        show_peaks = self.show_peaks if show_peaks is None else show_peaks
        # NOTE: inmask is a gpm
        sciImage = pypeitimage.PypeItImage(pseudo_dict['imgminsky'], ivar=pseudo_dict['sciivar'],
                                           bpm=np.logical_not(pseudo_dict['inmask']))
        sciImage.detector = self.coadd2d_stack.detectors[0].copy()
        # update platescale in the detector object to reflect the resampling done in coadd2d.
        # This is done to be able to propagate the correct spatial sampling to FindObjects and Extract
        sciImage.detector.platescale *= self.par['coadd2d']['spat_samp_fact']
        # if this is a mosaic, we need to update the detector object of each detector in the mosaic
        if "detectors" in sciImage.detector.keys():
            for det in sciImage.detector['detectors']:
                det.platescale *= self.par['coadd2d']['spat_samp_fact']

        slitmask_pseudo = pseudo_dict['slits'].slit_img()
        sciImage.build_mask(slitmask=slitmask_pseudo)

        # Make changes to parset specific to 2d coadds
        parcopy = self.par.copy()
        # Enforce low order traces since we are rectified
        parcopy['reduce']['findobj']['trace_npoly'] = int(np.clip(parcopy['reduce']['findobj']['trace_npoly'],None,3))
        # Manual extraction.
        manual_obj = None
        if self.par['coadd2d']['manual'] is not None and len(self.par['coadd2d']['manual']) > 0:
            manual_obj = ManualExtractionObj.by_fitstbl_input('None', self.par['coadd2d']['manual'], self.spectrograph)
        # Get bpm mask. There should not be any masked slits because we excluded those already
        # before the coadd, but we need to pass a bpm to FindObjects and Extract
        slits = pseudo_dict['slits']

        # Initiate FindObjects object
        objFind = find_objects.FindObjects.get_instance(sciImage, pseudo_dict['slits'], self.spectrograph, parcopy,
                                                        'science_coadd2d', tilts=pseudo_dict['tilts'],
                                                        bkg_redux=self.bkg_redux, manual=manual_obj,
                                                        find_negative=self.find_negative, basename=basename,
                                                        clear_ginga=clear_ginga, show=show)
        if show:
            gpm = sciImage.select_flag(invert=True)
            objFind.show('image', image=pseudo_dict['imgminsky']*gpm.astype(float), chname='imgminsky', slits=True)

        global_sky_pseudo, sobjs_obj = objFind.run(show_peaks=show or show_peaks, show_skysub_fit=show_skysub_fit)

        # maskdef stuff
        if parcopy['reduce']['slitmask']['assign_obj'] and slits.maskdef_designtab is not None:
            # Get pixel scale, binned and resampled (if requested), i.e., pixel scale of the pseudo image
            resampled_pixscale = parse.parse_binning(sciImage.detector.binning)[1]*sciImage.detector.platescale

            # Assign slitmask design information to detected objects
            slits.assign_maskinfo(sobjs_obj, resampled_pixscale, None, TOLER=parcopy['reduce']['slitmask']['obj_toler'])

            if parcopy['reduce']['slitmask']['extract_missing_objs'] is True:
                # Set the FWHM for the extraction of missing objects
                fwhm = slits.get_maskdef_extract_fwhm(sobjs_obj, resampled_pixscale,
                                                      parcopy['reduce']['slitmask']['missing_objs_fwhm'],
                                                      parcopy['reduce']['findobj']['find_fwhm'])
                # Assign undetected objects
                sobjs_obj = slits.mask_add_missing_obj(sobjs_obj, None, fwhm,
                                                       parcopy['reduce']['slitmask']['missing_objs_boxcar_rad']/resampled_pixscale)

        # Initiate Extract object
        exTract = extraction.Extract.get_instance(sciImage, pseudo_dict['slits'], sobjs_obj, self.spectrograph, parcopy,
                                                  'science_coadd2d', global_sky=None, tilts=pseudo_dict['tilts'],
                                                  waveimg=pseudo_dict['waveimg'], bkg_redux=self.bkg_redux,
                                                  basename=basename, show=show)

        skymodel_pseudo, _, objmodel_pseudo, ivarmodel_pseudo, outmask_pseudo, sobjs, _, _, _ = exTract.run(
            model_noise=False, spat_pix=pseudo_dict['spat_img'])

        # Add the rest to the pseudo_dict
        pseudo_dict['skymodel'] = skymodel_pseudo
        pseudo_dict['objmodel'] = objmodel_pseudo
        pseudo_dict['ivarmodel'] = ivarmodel_pseudo
        pseudo_dict['outmask'] = outmask_pseudo
        pseudo_dict['sobjs'] = sobjs
        self.pseudo_dict=pseudo_dict

        return pseudo_dict['imgminsky'], pseudo_dict['sciivar'], skymodel_pseudo, \
               objmodel_pseudo, ivarmodel_pseudo, outmask_pseudo, sobjs, sciImage.detector, slits, \
               pseudo_dict['tilts'], pseudo_dict['waveimg']

    @staticmethod
    def offsets_report(offsets, pixscale, offsets_method):
        """
        Log a summary of the offsets used for coaddition.

        Parameters
        ----------
        offsets : array-like
            Offsets, in pixels, applied to the input exposures.
        pixscale : float
            Pixel scale used to convert the offsets to arcseconds for reporting.
        offsets_method : str
            Description of how the offsets were determined.
        """

        if offsets_method is not None and offsets is not None:
            msg_string = '\n---------------------------------------------------------------------------------\n'
            msg_string += f' Summary of offsets from {offsets_method}     \n'
            msg_string += '---------------------------------------------------------------------------------\n'
            msg_string += '          file#      offset (pixels)    offset (arcsec)\n'
            for iexp, off in enumerate(offsets):
                msg_string += f'            {iexp:2d}            {off:6.2f}              {off*pixscale:6.3f}\n'
            msg_string += '---------------------------------------------------------------------------------'
            log.info(msg_string)

    def offset_slit_cen(self, slitid, offsets):
        """
        Return the slit center traces shifted by the exposure offsets.

        Parameters
        ----------
        slitid : int
            Slit or order identifier.
        offsets : array-like
            Offsets, in pixels, to apply to the slit centers for each exposure.

        Returns
        -------
        numpy.ndarray
            Stack of shifted slit-center traces with one column per exposure.
        """
        return [slits.center[:,slitid] - offsets[iexp]
                    for iexp, slits in enumerate(self.coadd2d_stack.slits_list)]

    def get_wave_grid(self):
        """
        Construct the wavelength grid used for the two-dimensional coadd.

        Returns
        -------
        wave_grid : numpy.ndarray
            New wavelength grid, not masked
        wave_grid_mid : numpy.ndarray
            New wavelength grid evaluated at the centers of the wavelength
            bins, that is this grid is simply offset from wave_grid by
            dsamp/2.0, in either linear space or log10 depending on whether
            linear or (log10 or velocity) was requested.  For iref or
            concatenate the linear wavelength sampling will be calculated.
        dsamp : float
            The pixel sampling for wavelength grid created.

        Notes
        -----
        The exact grid construction depends on the reduction mode and is controlled
        by :meth:`wave_method`.
        """
        nobjs_tot = int(np.array([len(spec) for spec in self.coadd2d_stack.specobjs_list]).sum())
        # TODO: Do we need this flag since we can determine whether or not we have specobjs from nobjs_tot?
        #  This all seems a bit hacky
        if self.par['coadd2d']['use_slits4wvgrid'] or nobjs_tot==0:
            nslits_tot = np.sum([slits.nslits for slits in self.coadd2d_stack.slits_list])
            waves, gpms = [], []
            box_radius = 3.
            #indx = 0
            # Loop on the exposures
            for iexp, (waveimg, slitmask, slits) in enumerate(zip(self.coadd2d_stack.waveimg_stack,
                                                self.coadd2d_stack.slitmask_stack,
                                                self.coadd2d_stack.slits_list)):
                slits_left, slits_righ, _ = slits.select_edges()
                row = np.arange(slits_left.shape[0])
                # Loop on the slits
                for kk, spat_id in enumerate(slits.spat_id):
                    mask = slitmask == spat_id
                    # Create apertures at 5%, 50%, and 95% of the slit width to cover full range of wavelengths
                    # on this slit
                    trace_spat = slits_left[:, kk][:,np.newaxis] +  np.outer((slits_righ[:,kk] - slits_left[:,kk]),[0.05,0.5,0.95])
                    box_denom = moment1d(waveimg * mask > 0.0, trace_spat, 2 * box_radius, row=row)[0]
                    wave_box = moment1d(waveimg * mask, trace_spat, 2 * box_radius,
                                    row=row)[0] / (box_denom + (box_denom == 0.0))
                    gpm_box = box_denom > 0.
                    waves += [wave for (wave, gpm) in zip(wave_box.T, gpm_box.T) if np.any(gpm)]
                    gpms += [(wave > 0.) & gpm for (wave, gpm) in zip(wave_box.T, gpm_box.T) if np.any(gpm)]

        else:
            waves, gpms = [], []
            for iexp, spec_this in enumerate(self.coadd2d_stack.specobjs_list):
                for spec in spec_this:
                    # NOTE: BOX extraction usage needed for quicklook
                    good_opt_ext = spec.has_opt_ext() and np.any(spec.OPT_MASK)
                    good_box_ext = spec.has_box_ext() and np.any(spec.BOX_MASK)
                    if good_opt_ext or good_box_ext:
                        waves.append(spec.OPT_WAVE if good_opt_ext else spec.BOX_WAVE)
                        gpms.append(spec.OPT_MASK if good_opt_ext else spec.BOX_MASK)
                        # TODO -- OPT_MASK is likely to become a bpm with int values
                        #gpm[:self.nspec_array[iexp], indx] = spec.OPT_MASK
                        #indx += 1

        return wvutils.get_wave_grid(waves=waves, gpms=gpms, wave_method=self.wave_method(),
                                                                spec_samp_fact=self.par['coadd2d']['spec_samp_fact'])

    def load_coadd2d_stacks(self, spec2d:list, chk_version:bool=False) -> CoAdd2dStack:
        """
        Load the input ``spec2d`` products into a :class:`CoAdd2dStack`.

        Parameters
        ----------
        spec2d : list
            List of input ``spec2d`` filenames or
            :class:`~pypeit.spec2dobj.Spec2DObj` objects.
        chk_version : bool, optional
            If True, check the version compatibility of the loaded files.

        Returns
        -------
        :class:`CoAdd2dStack`
            Dataclass containing the stacked science images, masks, coordinate
            images, slit definitions, exposure times, and related metadata needed for
            coaddition.

        Notes
        -----
        If the exposure times differ among the inputs beyond the allowed tolerance,
        the science, sky, and inverse-variance images are rescaled to a common
        effective exposure time before being stored in the returned stack.
        """
        # Grab the files
        #head2d_list = []

        # Image stacks
        sciimg_stack = []
        waveimg_stack = []
        skymodel_stack = []
        sciivar_stack = []
        mask_stack = []
        slitmask_stack = []
        exptime_stack = []
        #tilts_stack = []
        # Object stacks
        specobjs_list = []
        slits_list = []
        nfiles =len(spec2d)
        detectors_list = []
        maskdef_designtab_list = []
        spat_flexure_list = []
        for f in spec2d:
            if isinstance(f, Spec2DObj):
                # If spec2d is a list of objects
                s2dobj = f
            else:
                # If spec2d is a list of files, option to also use spec1ds
                s2dobj = Spec2DObj.from_file(f, self.detname, chk_version=chk_version)
                spec1d_file = f.replace('spec2d', 'spec1d')
                if Path(spec1d_file).is_file():
                    sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=chk_version)
                    this_det = sobjs.DET == self.detname
                    specobjs_list.append(sobjs[this_det])
            # TODO the code should run without a spec1d file, but we need to implement that
            slits_list.append(s2dobj.slits)
            detectors_list.append(s2dobj.detector)
            maskdef_designtab_list.append(s2dobj.maskdef_designtab)
            spat_flexure_list.append(s2dobj.sci_spat_flexure)

            sciimg_stack.append(s2dobj.sciimg)
            exptime_stack.append(s2dobj.head0['EXPTIME'])
            waveimg_stack.append(s2dobj.waveimg)
            skymodel_stack.append(s2dobj.skymodel)
            sciivar_stack.append(s2dobj.ivarmodel)
            mask_stack.append(s2dobj.bpmmask.mask)
            slitmask_stack.append(s2dobj.slits.slit_img(flexure=s2dobj.sci_spat_flexure))

        # check if exptime is consistent for all images
        exptime_coadd = np.percentile(exptime_stack, 50., method='higher')
        isclose_exptime = np.isclose(exptime_stack, exptime_coadd, atol=1.)
        if not np.all(isclose_exptime):
            log.warning('Exposure time is not consistent (within 1 sec) for all frames being coadded! '
                      f'Scaling each image by the median exposure time ({exptime_coadd} s) before coadding.')
            exp_scale = exptime_coadd / exptime_stack
            for iexp in range(nfiles):
                if not isclose_exptime[iexp]:
                    sciimg_stack[iexp] *= exp_scale[iexp]
                    skymodel_stack[iexp] *= exp_scale[iexp]
                    sciivar_stack[iexp] /= exp_scale[iexp]**2

        return CoAdd2dStack(specobjs_list=specobjs_list,
                    slits_list=slits_list,
                    slitmask_stack=slitmask_stack,
                    sciimg_stack=sciimg_stack,
                    sciivar_stack=sciivar_stack,
                    skymodel_stack=skymodel_stack, 
                    mask_stack=mask_stack,
                    waveimg_stack=waveimg_stack,
                    exptime_stack=exptime_stack,
                    exptime_coadd=exptime_coadd,
                    redux_path=Path.cwd(),
                    detectors=detectors_list,
                    spectrograph=self.spectrograph.name,
                    pypeline=self.spectrograph.pypeline,
                    maskdef_designtab_list=maskdef_designtab_list,
                    spat_flexure_list=spat_flexure_list)
    #                    tilts_stack=tilts_stack, waveimg_stack=waveimg_stack,

    def check_input(self, input, type):
        """
        Normalize a coadd configuration input to an internal representation.

        Parameters
        ----------
        input : object
            User-supplied value to validate and normalize.
        type : str
            Expected input type or semantic category.

        Returns
        -------
        object
            Validated and normalized value.

        Raises
        ------
        PypeItError
            Raised if the input does not conform to the expected type or format.
        """
        if type != 'weights' and type != 'offsets':
            raise PypeItError('Unrecognized type for check_input')
        if isinstance(input, (list, np.ndarray)):
            if len(input) != self.nexp:
                raise PypeItError(
                    f'If {type} are input it must be a list/array with same number of elements '
                    'as files.'
                )
            return np.atleast_1d(input).tolist() if type == 'weights' else np.atleast_1d(input)
        raise PypeItError(f'Unrecognized format for {type}')

    def compute_offsets(self):
        """
        Determine the relative offsets among the input exposures.

        Notes
        -----
        This base-class method is intended to be overridden by subclasses that define
        how offsets are determined for specific reduction modes.
        """
        log.info('Get Offsets')
        # binned pixel scale of the frames to be coadded
        pixscale = parse.parse_binning(self.coadd2d_stack.detectors[0].binning)[1]*self.coadd2d_stack.detectors[0].platescale
        # 1) offsets are provided in the header of the spec2d files
        if self.par['coadd2d']['offsets'] == 'header':
            log.info('Using offsets from header')
            dithoffs = [self.spectrograph.get_meta_value(f, 'dithoff') for f in self.spec2d]
            if None in dithoffs:
                raise PypeItError('Dither offsets keyword not found for one or more spec2d files. '
                           'Choose another option for `offsets`')
            dithoffs_pix = - np.array(dithoffs) / pixscale
            self.offsets = dithoffs_pix[0] - dithoffs_pix
            self.offsets_report(self.offsets, pixscale, 'header keyword')

        elif self.obj_id_bri is None and self.par['coadd2d']['offsets'] == 'auto':
            raise PypeItError('Offsets cannot be computed because no unique reference object '
                       'with the highest S/N was found. To continue, provide offsets in `Coadd2DPar`')

        # 2) a list of offsets is provided by the user (no matter if we have a bright object or not)
        elif isinstance(self.par['coadd2d']['offsets'], (list, np.ndarray)):
            log.info('Using user input offsets')
            # use them
            self.offsets = self.check_input(self.par['coadd2d']['offsets'], 'offsets')
            self.offsets_report(self.offsets, pixscale, 'user input')

        # 3) parset `offsets` is = 'maskdef_offsets' (no matter if we have a bright object or not)
        elif self.par['coadd2d']['offsets'] == 'maskdef_offsets':
            self.maskdef_offset = np.array([slits.maskdef_offset for slits in self.coadd2d_stack.slits_list])
            # Check if maskdef_offset is actually recoded in the SlitTraceSet
            if np.any(self.maskdef_offset == None):
                raise PypeItError(
                    'maskdef_offsets are not recoded in the SlitTraceSet for one or more files. '
                    'They cannot be used.'
                )
            # the offsets computed during the main reduction (`run_pypeit`) are used
            log.info('Determining offsets using maskdef_offset recoded in SlitTraceSet')
            self.offsets = self.maskdef_offset[0] - self.maskdef_offset
            self.offsets_report(self.offsets, pixscale, 'maskdef_offset')

        # 4) parset `offsets` = 'auto' but we have a bright object
        elif self.par['coadd2d']['offsets'] == 'auto' and self.obj_id_bri is not None:
            # see child method
            pass
        else:
            raise PypeItError('Invalid value for `offsets`')

    def compute_weights(self):
        """
        Determine the exposure weights used in the coadd.

        Notes
        -----
        This base-class method is intended to be overridden by subclasses that define
        how weights are determined for specific reduction modes.
        """
        log.info('Get Weights')

        # 1) User input weight
        if isinstance(self.par['coadd2d']['weights'], (list, np.ndarray)):
            # use those inputs
            self.use_weights = self.check_input(self.par['coadd2d']['weights'], 'weights')
            log.info('Using user input weights')

        # 2) No bright object and parset `weights` is 'auto' or 'uniform',
        # or Yes bright object but the user wants to still use uniform weights
        elif ((self.obj_id_bri is None) and (self.par['coadd2d']['weights'] in ['auto', 'uniform'])) or \
                ((self.obj_id_bri is not None) and (self.par['coadd2d']['weights'] == 'uniform')):
            if self.par['coadd2d']['weights'] == 'auto':
                # TODO maybe better behavior here would be to crash out to force the user to change the weight method explicitly
                # to 'uniform'. What I don't like here is that we are using uniform weights even though the user requested 'auto'
                # and they might miss the warning. Its debatable though.
                
                # warn if the user had put `auto` in the parset
                log.warning('Weights cannot be computed because no unique reference object '
                          'with the highest S/N was found. Using uniform weights instead.')
            elif self.par['coadd2d']['weights'] == 'uniform':
                log.info('Using uniform weights')
            # use uniform weights
            self.use_weights = (np.ones(self.nexp) / float(self.nexp)).tolist()

        # 3) Bright object exists and parset `weights` is equal to 'auto'
        elif (self.obj_id_bri is not None) and (self.par['coadd2d']['weights'] == 'auto'):
            # see child method
            pass
        else:
            raise PypeItError('Invalid value for `weights`')

    def _get_weights(self, indx=None):
        """
        Return the weights appropriate for a given slit or order.

        Parameters
        ----------
        indx : int, optional
            Slit or order index for which to return the weights.

        Returns
        -------
        list
            List of weights, one per exposure.
        """
        return self.use_weights

    @staticmethod
    def unpack_specobj(spec, spatord_id=None):
        """
        Extract trace and S/N information from a :class:`~pypeit.specobjs.SpecObj`.

        Parameters
        ----------
        spec : :class:`~pypeit.specobjs.SpecObj`
            Object to unpack.
        spatord_id : int, optional
            Slit spatial identifier or echelle order identifier used to select or
            validate the object.

        Returns
        -------
        flux : numpy.ndarray
            Flux array from the :class:`~pypeit.specobjs.SpecObj`
        ivar : numpy.ndarray
            Inverse variance array from the :class:`~pypeit.specobjs.SpecObj`
        gpm : numpy.ndarray
            Good pixel mask array from the :class:`~pypeit.specobjs.SpecObj`

        Notes
        -----
        This helper provides a common interface for the object-handling code used by
        the multislit and echelle subclasses.
        """
        # Get the slit/order ID if not provided
        if spatord_id is None:
            spatord_id = spec.ECH_ORDER if spec.ECH_ORDER is not None else spec.SLITID

        # get OBJID, which is different for Echelle and MultiSlit
        objid = spec.ECH_FRACPOS_ID if spec.ECH_FRACPOS_ID is not None else spec.SPAT_PIXPOS_ID

        # check if OPT_COUNTS is available
        if spec.has_opt_ext() and np.any(spec.OPT_MASK):
            _, flux, ivar, gpm = spec.get_opt_ext()
        # check if BOX_COUNTS is available
        elif spec.has_box_ext() and np.any(spec.BOX_MASK):
            _, flux, ivar, gpm = spec.get_box_ext()
            log.warning(f'Optimal extraction not available for obj_id {objid} '
                      f'in slit/order {spatord_id}. Using box extraction.')
        else:
            log.warning(f'Optimal and Boxcar extraction not available for obj_id {objid} in slit/order {spatord_id}.')
            _, flux, ivar, gpm = None, None, None, None

        return flux, ivar, gpm

    def get_brightest_obj(self, specobjs_list, spat_ids):
        """
        Identify the brightest reference object in each exposure.

        Parameters
        ----------
        specobjs_list : list
            List of :class:`~pypeit.specobjs.SpecObjs` containers, one per exposure.
        spat_ids : numpy.ndarray
            Slit or order identifiers used to constrain the object selection.

        Returns
        -------
        tuple
            Information describing the selected brightest objects.

        Notes
        -----
        This base-class method is intended to be overridden by subclasses because the
        object identifiers differ between multislit and echelle reductions.
        """
        raise PypeItError('The get_brightest_obj() method should be overloaded by the child class.')
    
    def handle_reference_obj(self):
        """
        Interpret the user-supplied reference-object selection.

        Notes
        -----
        This base-class method is overridden by subclasses.
        """
        
        raise PypeItError('The handle_reference_obj() method should be overloaded by the child class.')


    def reference_trace_stack(self, slitid, offsets=None, uniq_obj_id=None):
        """
        Construct the stack of reference traces used for rectification and coaddition.

        Parameters
        ----------
        slitid : int
            Slit or order identifier.
        offsets : array-like, optional
            Offsets applied to the reference traces for each exposure.
        uniq_obj_id : array-like, optional
            Object identifiers selecting the reference trace in each exposure.

        Returns
        -------
        numpy.ndarray
            Stack of reference traces with one column per exposure.

        Notes
        -----
        This base-class method is overridden by subclasses.
        """
        raise PypeItError('The reference_trace_stack() method should be overloaded by the child class.')


    def get_maskdef_dict(self, slit_idx, ref_trace_stack):
        """
        Collect mask-design metadata for a slit in the coadd.

        Parameters
        ----------
        slit_idx : int
            Index of the slit in the uncoadded frames.
        ref_trace_stack : numpy.ndarray
            Stack of reference traces about which the images are rectified and
            coadded. The shape is ``(nspec, nimgs)``.

        Returns
        -------
        dict
            Dictionary containing the mask-design metadata associated with the slit,
            including identifiers, object positions, slit centers, and design-table
            information.

        Notes
        -----
        This base-class method is overridden by subclasses.
        """

        return dict(maskdef_id=None, maskdef_objpos=None, maskdef_slitcen=None, maskdef_designtab=None)

    def wave_method(self):
        """
        Return the wavelength-grid construction method used for coadd2d.

        Returns
        -------
        str
            Name of the wavelength-grid method.

        Notes
        -----
        This base-class method is overridden by subclasses.
        """

        raise PypeItError('The wave_method() method should be overloaded by the child class.')

# Multislit can coadd with:
# 1) input offsets or if offsets is None, it will find the brightest trace and compute them
# 2) specified weights, or if weights is None and auto_weights=True, it will compute weights using the brightest object

class MultiSlitCoAdd2D(CoAdd2D):
    """
    Two-dimensional coaddition driver for multislit and longslit reductions.

    This subclass implements the multislit-specific handling of reference
    objects, offsets, weights, and wavelength-grid construction on top of the
    shared :class:`CoAdd2D` framework.

    Notes
    -----
    Multislit coadds can use either user-supplied offsets or offsets derived from
    a reference object or mask-design information. Weights may be user supplied,
    uniform, or determined automatically from the brightest detected object.
    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, 
                 only_slits=None, exclude_slits=None, sn_smooth_npix=None,
                 bkg_redux=False, find_negative=False, show=False, show_peaks=False, debug_offsets=False, debug=False):
        """
        Initialize a multislit or longslit two-dimensional coadd driver.

        Parameters
        ----------
        spec2d_files : list
            List of input ``spec2d`` files or
            :class:`~pypeit.spec2dobj.Spec2DObj` objects.
        spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            Spectrograph object.
        par : :class:`~pypeit.par.parset.ParSet`
            Parameter set controlling the coaddition.
        det : int or tuple, optional
            Detector or detector-mosaic identifier.
        only_slits : list, optional
            Slits to include.
        exclude_slits : list, optional
            Slits to exclude.
        sn_smooth_npix : int, optional
            S/N smoothing scale used for weight estimation.
        bkg_redux : bool, optional
            If True, the input science frames have already been background
            subtracted.
        find_negative : bool, optional
            If True, search for negative traces.
        show : bool, optional
            If True, display intermediate results.
        show_peaks : bool, optional
            If True, show QA plots for object-finding peaks.
        debug_offsets : bool, optional
            If True, show QA related to offset determination.
        debug : bool, optional
            If True, enable additional debug output.
        """

        # Attributes specifically used by MultislitCoAdd2D
        # This is an integer, which is the spatial slit id of the slit with the brightest object.
        # Used for both offsets (if offsets='auto') and weights (if weights='auto').
        # Can be user specified if user_obj_ids is provided
        self.spatid_bri = None

        # This will be an array of the object spatial pixel positions used for auto weights in each exposure
        self.spat_pixpos_id_weights = None

        super().__init__(spec2d_files, spectrograph, det=det, 
                         only_slits=only_slits, exclude_slits=exclude_slits,
                         sn_smooth_npix=sn_smooth_npix, bkg_redux=bkg_redux, find_negative=find_negative, par=par,
                         show=show, show_peaks=show_peaks, debug_offsets=debug_offsets,
                         debug=debug)


    def handle_reference_obj(self):
        """
        Interpret the multislit reference-object syntax.

        Notes
        -----
        This method parses the user-requested reference object used to compute
        offsets and weights in multislit coadds.
        """
        # Check if 1) specobjs exist AND 2) either of 'offsets' or 'weights' are 'auto'; otherwise return
        if not self.coadd2d_stack.specobjs_list or 'auto' not in [self.par['coadd2d']['offsets'], self.par['coadd2d']['weights']]:
            return

        # If no `user_obj_ids` are passed in, find the brightest object in the stack and obtain the relevant information
        if self.par['coadd2d']['user_obj_ids'] is None:
            self.spat_pixpos_id_bri, self.spat_pixpos_bri, self.spatid_bri, self.snr_bar_bri = self.get_brightest_obj(self.coadd2d_stack.specobjs_list, self.spat_ids)
            self.obj_id_bri = self.spat_pixpos_id_bri        
            return
        
        # The user passed in user_obj_ids that we will use these for the brighest object to 
        # be optionally used for weights.
        if self.par['coadd2d']['weights'] != 'auto':
            raise PypeItError('Parameter `user_obj_ids` can only be used if weights are set to `auto`.')
        if len(self.par['coadd2d']['user_obj_ids']) != self.nexp:
            raise PypeItError(
                'Parameter `user_obj_ids` must have the same number of elements as files.'
            )

        user_obj_exist = np.zeros(self.nexp, dtype=bool)
        # Get the flux, ivar, gpm, and spatial pixel position of the user object
        fluxes, ivars, gpms, spatids, spat_pixpos = [], [], [], [], []
        # Place the input objects into `spat_pixpos_id_bri`
        self.spat_pixpos_id_bri = np.array(self.par['coadd2d']['user_obj_ids'])
        # Loop over specobjs
        for i, sobjs in enumerate(self.coadd2d_stack.specobjs_list):
            # Get the index of the user-requested object in this slit
            user_idx = sobjs.slitorder_uniq_id_indices(self.spat_pixpos_id_bri[i])
            if np.any(user_idx):
                this_sobj = sobjs[user_idx][0]
                flux_iobj, ivar_iobj, gpm_iobj = self.unpack_specobj(this_sobj)
                if flux_iobj is not None and ivar_iobj is not None and gpm_iobj is not None:
                    fluxes.append(flux_iobj)
                    ivars.append(ivar_iobj)
                    gpms.append(gpm_iobj)
                    spat_pixpos.append(this_sobj.SPAT_PIXPOS)
                    spatids.append(this_sobj.SLITID)
                    user_obj_exist[i] = True

        # Check that the user object exists in all the exposures
        if not np.all(user_obj_exist):
            raise PypeItError(
                    'Not all of the spat_pixpos_ids provided through `user_obj_ids` exist '
                    'in all of the files.'
            )

        # Check that all spatids are within the spat_toler of each other
        if not np.all(np.abs(spatids - np.mean(spatids[0])) <= self.par['coadd2d']['spat_toler']):
            raise PypeItError('Not all spatial IDs are within spat_toler of each other')
        self.spatid_bri = int(np.rint(np.mean(spatids)))
        self.spat_pixpos_bri = np.array(spat_pixpos)
        self.snr_bar_bri, _ = coadd.calc_snr(fluxes, ivars, gpms)                
        self.obj_id_bri = self.spat_pixpos_id_bri                


    # TODO When we run multislit, we actually compute the rebinned images twice. Once here to compute the offsets
    # and another time to weighted_combine the images in compute2d. This could be sped up
    # TODO The reason we rebin the images for the purposes of computing the offsets is to deal with combining
    # data that are dithered in the spectral direction. In these situations you cannot register the two dithered
    # reference objects into the same frame without first rebinning them onto the same grid.
    def compute_offsets(self):
        """
        Determine the exposure offsets for a multislit coadd.

        The offsets are defined relative to the first exposure and may be supplied by
        the user, inferred from slitmask-design information, or computed from the
        brightest detected reference object.
        """
        super().compute_offsets()

        # If no bright object ID or the offsets are not 'auto', return
        if self.obj_id_bri is None or self.par['coadd2d']['offsets'] != 'auto':
            return

        # Set boolean
        use_obj_ids = self.par['coadd2d']['user_obj_ids'] is not None

        # Compute offsets using the bright object / user-supplied object IDs
        if use_obj_ids:
            offsets_method = f'user object on slitid = {self.spatid_bri}'
        else:
            offsets_method = f'brightest object found on slit: {self.spatid_bri} with avg SNR={np.mean(self.snr_bar_bri):5.2f}'

        log.info(f'Determining offsets using {offsets_method}')
        thismask_stack = [np.abs(slitmask - self.spatid_bri) <= self.par['coadd2d']['spat_toler'] for slitmask in self.coadd2d_stack.slitmask_stack]

        slitidx_bri = np.where(np.abs(self.spat_ids - self.spatid_bri) <= self.par['coadd2d']['spat_toler'])[0][0]
        # TODO Need to think abbout whether we have multiple tslits_dict for each exposure or a single one
        trace_stack_bri = [slits.center[:, slitidx_bri]
                                for slits in self.coadd2d_stack.slits_list]

        # Determine the wavelength grid that we will use for the current slit/order
        ## TODO: Should the spatial and spectral samp_facts here match those of the final coadded data, or she would
        ## compute offsets at full resolution??
        wave_bins = coadd.get_wave_bins(thismask_stack, self.coadd2d_stack.waveimg_stack, self.wave_grid)
        dspat_bins, dspat_stack = coadd.get_spat_bins(thismask_stack, trace_stack_bri)

        sci_list = [[sciimg - skymodel for sciimg, skymodel in zip(self.coadd2d_stack.sciimg_stack, self.coadd2d_stack.skymodel_stack)]]
        var_list = [[utils.inverse(sciivar) for sciivar in self.coadd2d_stack.sciivar_stack]]
        if use_obj_ids:
            ny, _ = self.coadd2d_stack.sciimg_stack[0].shape
            obj_list = [[(spat, ny/2) for spat in self.spat_pixpos_bri]]
        else:
            obj_list = None

        log.info('Rebinning Images')
        mask_stack = [mask == 0 for mask in self.coadd2d_stack.mask_stack]
        sci_list_rebin, var_list_rebin, norm_rebin_stack, _, obj_list_rebin = coadd.rebin2d(
            wave_bins, 
            dspat_bins, 
            self.coadd2d_stack.waveimg_stack, 
            dspat_stack, 
            thismask_stack, 
            mask_stack, 
            sci_list, 
            var_list, 
            obj_list,
        )

        # Build up the masks
        thismask = np.ones_like(sci_list_rebin[0][0,:,:],dtype=bool)
        nspec_pseudo, nspat_pseudo = thismask.shape
        slit_left = np.full(nspec_pseudo, 0.0)
        slit_righ = np.full(nspec_pseudo, nspat_pseudo)
        inmask = norm_rebin_stack > 0
        traces_rect = np.zeros((nspec_pseudo, self.nexp))
        user_obj_dspats = []


        # Loop over exposures
        for iexp in range(self.nexp):

            if use_obj_ids:
                spat, spec = obj_list_rebin[0][iexp]
                detnum = self.coadd2d_stack.detectors[iexp].det
                # Treat the object finding as a manual process
                manual_obj = ManualExtractionObj.by_fitstbl_input(
                    self.spec2d[iexp],
                    f"{detnum}:{spat}:{spec}:{self.par['reduce']['findobj']['find_fwhm']}",
                    self.spectrograph)
                manual_extract_dict = manual_obj.dict_for_objfind(self.spectrograph.get_det_name(detnum))
            else:
                # Otherwise, indicate auto-find
                manual_extract_dict = None

            # Perform the extraction
            sobjs_exp = findobj_skymask.objs_in_slit(
                sci_list_rebin[0][iexp,:,:], 
                utils.inverse(var_list_rebin[0][iexp,:,:]), 
                thismask, 
                slit_left, 
                slit_righ,
                inmask=inmask[iexp,:,:],
                fwhm=self.par['reduce']['findobj']['find_fwhm'],
                trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
                maxshift=self.par['reduce']['findobj']['trace_maxshift'],
                maxdev=self.par['reduce']['findobj']['trace_maxdev'],
                numiterfit=self.par['reduce']['findobj']['find_numiterfit'],
                ncoeff=self.par['reduce']['findobj']['trace_npoly'],
                snr_thresh=self.par['reduce']['findobj']['snr_thresh'],
                nperslit=1 if self.par['coadd2d']['user_obj_ids'] is None else None, 
                find_min_max=self.par['reduce']['findobj']['find_min_max'],
                spec_min_max=self.par['reduce']['findobj']['trace_min_max'],
                hand_extract_dict=manual_extract_dict,
                show_trace=self.debug_offsets, 
                show_peaks=self.debug_offsets)

            # Perform QA on the extracted objects
            if len(sobjs_exp) == 0:
                raise PypeItError(
                    f'No objects found in the rebinned image for file {iexp} (used to compute '
                    'the offsets).  Check `FindObjPar` parameters and try to adjust '
                    '`snr_thresh`.'
                )

            # Add the rectified traces to the stack
            if use_obj_ids:
                # find the spectrum of the user object and the corresponding trace
                idx_orig = self.coadd2d_stack.specobjs_list[iexp].slitorder_uniq_id_indices(
                    self.par['coadd2d']['user_obj_ids'][iexp])
                trace_orig = self.coadd2d_stack.specobjs_list[iexp][idx_orig].TRACE_SPAT
                # find the slit of the user object and its left edge
                slitidx_orig = \
                    self.coadd2d_stack.slits_list[iexp].spat_id == self.coadd2d_stack.specobjs_list[iexp][idx_orig].SLITID
                left_edge_orig = self.coadd2d_stack.slits_list[iexp].select_edges(
                    flexure=self.coadd2d_stack.spat_flexure_list[iexp])[0][:, slitidx_orig][:,0]
                # Compute the mean median offset betweeh the original trace and the left edge of the slit
                dist_to_left = np.median(trace_orig - left_edge_orig)
                # Identify the trace in the sobjs_exp from the rebinned image that is closest to the original trace taking this offset into account
                dspat_exp_orig = np.abs(np.median(sobjs_exp.TRACE_SPAT - dist_to_left, axis=1))
                dspat_ex_orig_min = dspat_exp_orig.min()
                if dspat_ex_orig_min < self.par['coadd2d']['spat_toler']:
                    traces_rect[:, iexp] = sobjs_exp[np.argmin(dspat_exp_orig)].TRACE_SPAT
                    user_obj_dspats.append(dspat_ex_orig_min)
                else:
                    raise PypeItError(f'Could not identify an object in the rebinned image corresponding '
                                f'to the trace for the user object {self.par["coadd2d"]["user_obj_ids"][iexp]} '
                                f'in file {iexp+1} within the specified spatial '
                                f'tolerance ={self.par["coadd2d"]["spat_toler"]}')
            else: 
                traces_rect[:, iexp] = sobjs_exp.TRACE_SPAT

        # After looping through all exposures, announce the offsets for user-specified objects
        if use_obj_ids:
            log.info(f'The median distance between the original traces and those in the '
                        f'rebinned image for the user_obj_ids is {np.median(user_obj_dspats):.2f} pixels')

        # Now deterimine the offsets. Arbitrarily set the zeroth trace to the reference
        med_traces_rect = np.median(traces_rect,axis=0)
        self.offsets = med_traces_rect[0] - med_traces_rect
        # TODO create a QA with this
        if self.debug_offsets:
            for iexp in range(self.nexp):
                plt.plot(traces_rect[:, iexp], linestyle='--', label='original trace')
                plt.plot(traces_rect[:, iexp] + self.offsets[iexp], label='shifted traces')
                plt.legend()
            plt.show()

        # Binned pixel scale of the frames to be coadded
        pixscale = parse.parse_binning(self.coadd2d_stack.detectors[0].binning)[1]*self.coadd2d_stack.detectors[0].platescale
        self.offsets_report(self.offsets, pixscale, offsets_method)


    def compute_weights(self):
        """
        Determine the exposure weights for a multislit coadd.

        This method sets the internal ``use_weights`` attribute using either the
        user-supplied weights, uniform weights, or weights computed from the
        brightest detected reference object.
        """

        super().compute_weights()

        # If no bright object ID or the weights are not 'auto', return
        if self.obj_id_bri is None or self.par['coadd2d']['weights'] != 'auto':
            return

        # compute weights using bright object
        # TODO add a parset for weight_method in optimal_weights. The default is currently 'auto'
        _, self.use_weights = self.optimal_weights(self.obj_id_bri)
        if self.par['coadd2d']['user_obj_ids'] is not None:
            log.info(f'Weights computed using a unique reference object in slit={self.spatid_bri} provided by the user')
        else:
            log.info(f'Weights computed using a unique reference object in slit={self.spatid_bri} with the highest S/N')
        self.snr_report(self.spatid_bri, self.spat_pixpos_bri, self.snr_bar_bri)


    def get_brightest_obj(self, specobjs_list, slit_spat_ids):
        """
        Identify the brightest reference object in each multislit exposure.

        Parameters
        ----------
        specobjs_list : list
            List of :class:`~pypeit.specobjs.SpecObjs` containers, one per exposure.
        slit_spat_ids : numpy.ndarray
            Array of slit spatial identifiers used to constrain the search.

        Returns
        -------
        tuple
            Tuple containing:

            - the reference-object ``SPAT_PIXPOS`` identifiers for each exposure,
            - the corresponding spatial positions,
            - the slit ``SPAT_ID`` containing the highest-S/N object, and
            - the average S/N values for the selected object in each exposure.
        """
        log.info('Finding brightest object')
        nexp = len(specobjs_list)
        nslits = slit_spat_ids.size

        slit_snr_max = np.zeros((nslits, nexp), dtype=float)
        bpm = np.ones(slit_snr_max.shape, dtype=bool)
        spat_pixpos_id_max = np.zeros((nslits, nexp), dtype=int)
        spat_pixpos_max = np.zeros((nslits, nexp), dtype=float)
        # Loop over each exposure, slit, find the brightest object on that slit for every exposure
        for iexp, sobjs in enumerate(specobjs_list):
            log.info(f"Working on file {iexp}")
            for islit, spat_id in enumerate(slit_spat_ids):
                if len(sobjs) == 0:
                    continue
                ithis = np.abs(sobjs.SLITID - spat_id) <= self.par['coadd2d']['spat_toler']
                if np.any(ithis):
                    spat_pixpos_id_this = sobjs[ithis].SPAT_PIXPOS_ID
                    spat_pixpos_this = sobjs[ithis].SPAT_PIXPOS
                    fluxes, ivars, gpms = [], [], []
                    for spec in sobjs[ithis]:
                        flux_iobj, ivar_iobj, gpm_iobj = self.unpack_specobj(spec, spatord_id=spat_id)
                        if flux_iobj is not None and ivar_iobj is not None and gpm_iobj is not None:
                            fluxes.append(flux_iobj)
                            ivars.append(ivar_iobj)
                            gpms.append(gpm_iobj)

                    # if there are objects on this slit left, we can proceed with computing rms_sn
                    if len(fluxes) > 0:
                        rms_sn, _ = coadd.calc_snr(fluxes, ivars, gpms)
                        imax = np.argmax(rms_sn)
                        slit_snr_max[islit, iexp] = rms_sn[imax]
                        spat_pixpos_id_max[islit, iexp] = spat_pixpos_id_this[imax]
                        spat_pixpos_max[islit, iexp] = spat_pixpos_this[imax]
                        bpm[islit, iexp] = False

        # If a slit has bpm = True for some exposures and not for others, set bpm = True for all exposures
        # Find the rows where any of the bpm values are True
        bpm_true_idx = np.array([np.any(b) for b in bpm])
        if np.any(bpm_true_idx):
            # Flag all exposures in those rows
            bpm[bpm_true_idx, :] = True

        # Find the highest snr object among all the slits
        if np.all(bpm):
            log.warning(
                'You do not appear to have a unique reference object that was traced as the '
                'highest S/N ratio on the same slit of every file. Try increasing the parameter '
                '`spat_toler`.'
            )
            return None, None, None, None
        else:
            # mask the bpm
            slit_snr_max_masked = np.ma.array(slit_snr_max, mask=bpm)
            slit_snr = np.mean(slit_snr_max_masked, axis=1)
            slitid = np.argmax(slit_snr)
            snr_bar_mean = slit_snr[slitid]
            snr_bar = slit_snr_max[slitid, :]
            spat_pix_pos_id = spat_pixpos_id_max[slitid, :]
            spat_pixpos = spat_pixpos_max[slitid, :]

        return spat_pix_pos_id, spat_pixpos, slit_spat_ids[slitid], snr_bar

    def snr_report(self, slitid, spat_pixpos, snr_bar):
        """
        Log an S/N report for the multislit reference object.

        Parameters
        ----------
        slitid : int
            ``SPAT_ID`` of the slit containing the reference object.
        spat_pixpos : numpy.ndarray
            Spatial positions of the reference object, one per exposure.
        snr_bar : numpy.ndarray
            Average S/N values of the reference object, one per exposure.
        """

        # Print out a report on the SNR
        msg_string = '\n-------------------------------------\n'
        msg_string += '  Summary for highest S/N object\n'
        msg_string += f'      found on slitid = {slitid}            \n'
        msg_string += '-------------------------------------\n'
        msg_string += '      file#   spat_pixpos     S/N\n'
        msg_string += '-------------------------------------\n'
        for iexp, (spat,snr) in enumerate(zip(spat_pixpos, snr_bar)):
            msg_string += f'       {iexp:2d}      {spat:7.1f}      {snr:5.2f}\n'
        msg_string += '-------------------------------------'
        log.info(msg_string)


    # TODO add an option here to actually use the reference trace for cases where they are on the same slit and it is
    # single slit???
    def reference_trace_stack(self, slitid, offsets=None, uniq_obj_id=None):
        """
        Construct the reference-trace stack for a multislit slit.

        Parameters
        ----------
        slitid : int
            Slit identifier for which to build the reference traces.
        offsets : list or numpy.ndarray, optional
            Offsets to apply to the slit-center traces, one per exposure.
        uniq_obj_id : list or numpy.ndarray, optional
            Object identifiers. This argument is accepted for interface
            compatibility but is not used in multislit mode.

        Returns
        -------
        list
            List of reference traces, one per exposure, for the requested slit.
        """

        return self.offset_slit_cen(slitid, offsets)

    def get_maskdef_dict(self, slit_idx, ref_trace_stack):
        """
        Collect mask-design metadata for a multislit slit in the coadd.

        Parameters
        ----------
        slit_idx : int
            Index of the slit in the uncoadded frames.
        ref_trace_stack : numpy.ndarray
            Stack of reference traces about which the images are rectified and
            coadded. The shape is ``(nspec, nimgs)``.

        Returns
        -------
        dict
            Dictionary containing mask-design metadata for the slit.
        """

        # maskdef info
        if self.par['calibrations']['slitedges']['use_maskdesign'] and \
                self.coadd2d_stack.slits_list[0].maskdef_id is not None and \
                self.coadd2d_stack.slits_list[0].maskdef_objpos is not None and \
                self.coadd2d_stack.maskdef_designtab_list[0] is not None and \
                self.par['coadd2d']['offsets'] == 'maskdef_offsets':
            # maskdef_designtab info for only this slit
            this_idx = self.coadd2d_stack.maskdef_designtab_list[0]['SPAT_ID'] == self.spat_ids[slit_idx]
            this_maskdef_designtab = self.coadd2d_stack.maskdef_designtab_list[0][this_idx]
            # remove columns that are irrelevant in the coadd2d frames
            this_maskdef_designtab.remove_columns(['TRACEID', 'TRACESROW', 'TRACELPIX', 'TRACERPIX',
                                                   'SLITLMASKDEF', 'SLITRMASKDEF'])
            this_maskdef_designtab.meta['MASKRMSL'] = 0.
            this_maskdef_designtab.meta['MASKRMSR'] = 0.

            # maskdef_id for this slit
            imaskdef_id = self.coadd2d_stack.slits_list[0].maskdef_id[slit_idx]

            # maskdef_slitcen (slit center along the spectral direction) and
            # maskdef_objpos (expected position of the target, as distance from left slit edge) for this slit

            # These are the binned maskdef_slitcen positions w.r.t. the center of the slit in ref_trace_stack
            slit_cen_dspat_vec = np.zeros(self.nexp)
            # These are the binned maskdef_objpos positions w.r.t. the center of the slit in ref_trace_stack
            objpos_dspat_vec = np.zeros(self.nexp)
            for iexp in range(self.nexp):
                # get maskdef_slitcen
                mslitcen_pixpos = self.coadd2d_stack.slits_list[iexp].maskdef_slitcen
                if mslitcen_pixpos.ndim < 2:
                    mslitcen_pixpos = mslitcen_pixpos[:, None]
                maskdef_slitcen_pixpos = mslitcen_pixpos[self.nspec_array[0]//2, slit_idx] + self.maskdef_offset[iexp]

                # get maskdef_objpos
                # find left edge
                slits_left, _, _ = \
                    self.coadd2d_stack.slits_list[iexp].select_edges(flexure=self.coadd2d_stack.spat_flexure_list[iexp])
                # targeted object spat pix
                maskdef_obj_pixpos = \
                    self.coadd2d_stack.slits_list[iexp].maskdef_objpos[slit_idx] + self.maskdef_offset[iexp] \
                    + slits_left[slits_left.shape[0]//2, slit_idx]

                # change reference system
                ref_trace = ref_trace_stack[iexp]
                nspec_this = ref_trace.shape[0]
                slit_cen_dspat_vec[iexp] = (maskdef_slitcen_pixpos - ref_trace[nspec_this // 2]) / self.par['coadd2d']['spat_samp_fact']

                objpos_dspat_vec[iexp] = (maskdef_obj_pixpos - ref_trace[nspec_this // 2]) / self.par['coadd2d']['spat_samp_fact']

            imaskdef_slitcen_dspat = np.mean(slit_cen_dspat_vec)
            imaskdef_objpos_dspat = np.mean(objpos_dspat_vec)

        else:
            this_maskdef_designtab = None
            imaskdef_id = None
            imaskdef_slitcen_dspat = None
            imaskdef_objpos_dspat = None

        return dict(maskdef_id=imaskdef_id, maskdef_objpos=imaskdef_objpos_dspat,
                    maskdef_slitcen=imaskdef_slitcen_dspat, maskdef_designtab=this_maskdef_designtab)

    def wave_method(self):
        """
        Return the wavelength-grid method used for multislit coadd2d.

        Returns
        -------
        str
            Name of the wavelength-grid method.
        """
        return self.par['coadd2d']['wave_method'] if self.par['coadd2d']['wave_method'] is not None else 'linear'


# Echelle can either stack with:
# 1) input offsets or if offsets is None, it will find the objid of brightest trace and stack all orders relative to the trace of this object.
# 2) specified weights, or if weights is None and auto_weights=True,
#    it will use wavelength dependent weights determined from the spectrum of the brightest objects objid on each order

class EchelleCoAdd2D(CoAdd2D):
    """
    Two-dimensional coaddition driver for echelle reductions.

    This subclass implements the echelle-specific handling of reference objects,
    offsets, weights, and wavelength-grid construction on top of the shared
    :class:`CoAdd2D` framework.

    Notes
    -----
    Echelle coadds can be aligned either by user-supplied offsets or by the trace
    of a selected reference object. Weights may be user supplied, uniform, or
    derived order-by-order from the brightest reference object.
    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, 
                 only_slits=None, exclude_slits=None, sn_smooth_npix=None,
                 bkg_redux=False, find_negative=False, show=False, show_peaks=False, debug_offsets=False, debug=False):
        """
        Initialize an echelle two-dimensional coadd driver.

        Parameters
        ----------
        spec2d_files : list
            List of input ``spec2d`` files or
            :class:`~pypeit.spec2dobj.Spec2DObj` objects.
        spectrograph : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            Spectrograph object.
        par : :class:`~pypeit.par.parset.ParSet`
            Parameter set controlling the coaddition.
        det : int or tuple, optional
            Detector or detector-mosaic identifier.
        only_slits : list, optional
            Orders to include.
        exclude_slits : list, optional
            Orders to exclude.
        sn_smooth_npix : int, optional
            S/N smoothing scale used for weight estimation.
        bkg_redux : bool, optional
            If True, the input science frames have already been background
            subtracted.
        find_negative : bool, optional
            If True, search for negative traces.
        show : bool, optional
            If True, display intermediate results.
        show_peaks : bool, optional
            If True, show QA plots for object-finding peaks.
        debug_offsets : bool, optional
            If True, show QA related to offset determination.
        debug : bool, optional
            If True, enable additional debug output.
        """
        super().__init__(spec2d_files, spectrograph, det=det, 
                         only_slits=only_slits, exclude_slits=exclude_slits,
                         sn_smooth_npix=sn_smooth_npix, bkg_redux=bkg_redux, find_negative=find_negative,
                         par=par, show=show, show_peaks=show_peaks, debug_offsets=debug_offsets,
                         debug=debug)


    def handle_reference_obj(self):
        """
        Interpret the echelle reference-object syntax.

        Notes
        -----
        This method parses the user-requested reference object used to compute
        offsets and weights in echelle coadds.
        """

        # If a user-input object to compute offsets and weights is provided, check if it exists and get the needed info
        if len(self.coadd2d_stack.specobjs_list) > 0 and self.par['coadd2d']['user_obj_ids'] is not None:
            if len(self.par['coadd2d']['user_obj_ids']) != self.nexp:
                raise PypeItError(
                    f'Parameter `user_obj_ids` {self.par["coadd2d"]["user_obj_ids"]} must have '
                    f'the same number of elements as files {self.nexp}.'
                )
            else:
                # does it exists?
                user_obj_exist = np.zeros((self.nexp,self.nslits_single), dtype=bool)
                orders= self.coadd2d_stack.slits_list[0].slitord_id
                for iexp, sobjs in enumerate(self.coadd2d_stack.specobjs_list):
                    for iord, ord in enumerate(orders):
                        # check if the object exists in this exposure
                        ind = sobjs.slitorder_uniq_id_indices(self.par['coadd2d']['user_obj_ids'][iexp], order=ord)
                        if (len(ind) == 0) or (not np.any(ind)):
                            raise PypeItError(
                                'Object with user_obj_id '
                                f'{self.par["coadd2d"]["user_obj_ids"][iexp]} does not exist in '
                                f'file {iexp+1} for order {ord}.'
                            )
                        flux, ivar, mask = self.unpack_specobj(sobjs[ind][0])
                        if flux is not None and ivar is not None and mask is not None:
                                user_obj_exist[iexp, iord] = True

                            
                if not np.all(user_obj_exist):
                    raise PypeItError(
                        'Object provided through `user_obj_ids` does not exist in all the files.'
                    )

                # get the needed info about the user object
                self.obj_id_bri = np.array(self.par['coadd2d']['user_obj_ids'])

        elif len(self.coadd2d_stack.specobjs_list) > 0 and (self.par['coadd2d']['offsets'] == 'auto' or self.par['coadd2d']['weights'] == 'auto'):
            self.obj_id_bri, self.snr_bar_bri = \
                self.get_brightest_obj(self.coadd2d_stack.specobjs_list, self.coadd2d_stack.slits_list[0].slitord_id)
        
        

    def compute_offsets(self):
        """
        Determine the exposure offsets for an echelle coadd.

        The offsets are defined relative to the first exposure and may be either
        user supplied or derived from the trace of the selected reference object.
        """
        super().compute_offsets()

        # adjustment for echelle to case 2): a list of offsets is provided by the user
        if isinstance(self.offsets, (list, np.ndarray)):
            self.obj_id_bri = None

        # adjustment for echelle to case 4) parset `offsets` = 'auto' but we have a bright object
        elif self.par['coadd2d']['offsets'] == 'auto' and self.obj_id_bri is not None:
            # offsets are not determined, but the bright object is used to construct
            # a reference trace (this is done in coadd using method `reference_trace_stack`)
            self.offsets = None
            if self.par['coadd2d']['user_obj_ids'] is not None:
                log.info('Reference trace about which 2d coadd is performed is computed using user object')
            else:
                log.info('Reference trace about which 2d coadd is performed is computed using the brightest object')

    def compute_weights(self):
        """
        Determine the exposure weights for an echelle coadd.

        This method sets the internal ``use_weights`` attribute using either the
        user-supplied weights, uniform weights, or wavelength-dependent weights
        computed from the brightest detected reference object.
        """
        super().compute_weights()

        # adjustment for echelle to case 3) Bright object exists and parset `weights` is equal to 'auto'
        if (self.obj_id_bri is not None) and (self.par['coadd2d']['weights'] == 'auto'):
            # computing a list of weights for all the slitord_ids that we than parse in coadd
            slitord_ids = self.coadd2d_stack.slits_list[0].slitord_id
            self.use_weights = []
            for order in slitord_ids:
                _, iweights = self.optimal_weights(self.obj_id_bri, order=order)
                self.use_weights.append(iweights)
            if self.par['coadd2d']['user_obj_ids'] is not None:
                log.info('Weights computed using a unique reference object provided by the user')
                # TODO: implement something here to print out the snr_report
            else:
                log.info('Weights computed using a unique reference object with the highest S/N')
                self.snr_report(self.snr_bar_bri)

    def _get_weights(self, indx=None):
        """
        Determine the exposure weights for an echelle coadd.

        This method sets the internal ``use_weights`` attribute using either the
        user-supplied weights, uniform weights, or wavelength-dependent weights
        computed from the brightest detected reference object.
        """

        # if this is echelle data and the parset 'weights' is set to 'auto',
        # then the weights are computed per order, i.e., every order has a
        # different set of weights in each exposure (len(self.use_weights[indx]) = nexp)
        if self.par['coadd2d']['weights'] == 'auto' and indx is None:
            raise PypeItError('The index of the slit/order must be provided when using auto weights for Echelle data.')
        return self.use_weights[indx] if self.par['coadd2d']['weights'] == 'auto' else super()._get_weights()


    def get_brightest_obj(self, specobjs_list, orders):
        """
        Identify the brightest reference object in each echelle exposure.

        Parameters
        ----------
        specobjs_list : list
            List of :class:`~pypeit.specobjs.SpecObjs` containers, one per exposure.
        orders : numpy.ndarray
            Array of order identifiers over which to search for the brightest object.

        Returns
        -------
        tuple
            Tuple containing:

            - the reference-object echelle fractional-position identifiers for each
              exposure, and
            - the average S/N values of the selected object in each exposure.

        """
        log.info('Finding brightest object')
        nexp = len(specobjs_list)

        fracpos_id = np.zeros(nexp, dtype=int)
        snr_bar = np.zeros(nexp)
        for iexp, sobjs in enumerate(specobjs_list):
            log.info("Working on file {}".format(iexp))
            uni_fracpos_id = np.unique(sobjs.ECH_FRACPOS_ID)
            nobjs = len(uni_fracpos_id)
            order_snr = np.zeros((orders.size, nobjs), dtype=float)
            bpm = np.ones((orders.size, nobjs), dtype=bool)
            for iord, ord in enumerate(orders):
                for iobj in range(nobjs):
                    ind = sobjs.slitorder_uniq_id_indices(uni_fracpos_id[iobj], order=ord)
                    flux, ivar, mask = self.unpack_specobj(sobjs[ind][0], spatord_id=sobjs[ind][0].ECH_ORDER)

                    if flux is not None and ivar is not None and mask is not None:
                        rms_sn, _ = coadd.calc_snr([flux], [ivar], [mask])
                        order_snr[iord, iobj] = rms_sn[0]
                        bpm[iord, iobj] = False

            # If there are orders that have bpm = True for some objs and not for others, set bpm = True for all objs
            # Find the rows where any of the bpm values are True
            bpm_true_idx = np.array([np.any(b) for b in bpm])
            if np.any(bpm_true_idx):
                # Flag all objs in those rows
                bpm[bpm_true_idx, :] = True

            # Compute the average SNR and find the brightest object
            if not np.all(bpm):
                # mask the bpm
                order_snr_masked = np.ma.array(order_snr, mask=bpm)
                snr_bar_vec = np.mean(order_snr_masked, axis=0)
                fracpos_id[iexp] = uni_fracpos_id[snr_bar_vec.argmax()]
                snr_bar[iexp] = snr_bar_vec[snr_bar_vec.argmax()]
        if 0 in snr_bar:
            log.warning(
                'You do not appear to have a unique reference object that was traced as the '
                'highest S/N ratio for every file.'
            )
            return None, None
        return fracpos_id, snr_bar

    def snr_report(self, snr_bar):
        """
        Log an S/N report for the echelle reference object.

        Parameters
        ----------
        snr_bar : numpy.ndarray
            Average S/N values of the selected reference object, one per exposure.
        """

        # Print out a report on the SNR
        msg_string = '\n-------------------------------------'
        msg_string += '\n  Summary for highest S/N object'
        msg_string += '\n-------------------------------------'
        msg_string += '\n          file#        S/N'
        for iexp, snr in enumerate(snr_bar):
            msg_string += f'\n            {iexp:d}         {snr:5.2f}'
        msg_string += '\n-------------------------------------'
        log.info(msg_string)


    def reference_trace_stack(self, slitid, offsets=None, uniq_obj_id=None):
        """
        Construct the reference-trace stack for an echelle order.

        Parameters
        ----------
        slitid : int
            Order identifier for which to build the reference traces.
        offsets : list or numpy.ndarray, optional
            Offsets to apply to the slit-center traces, one per exposure.
        uniq_obj_id : list or numpy.ndarray, optional
            Object identifiers selecting the reference trace in each exposure. For
            echelle reductions, these are ``ECH_FRACPOS_ID`` values.

        Returns
        -------
        list
            List of reference traces, one per exposure, for the requested order.

        Raises
        ------
        PypeItError
            Raised if neither or both of ``offsets`` and ``uniq_obj_id`` are
            provided.

        Notes
        -----
        Two alignment modes are supported:

        - offset-based alignment about the slit or order center, and
        - object-based alignment about a selected reference trace.
        """

        # check inputs
        if offsets is not None and uniq_obj_id is not None:
            raise PypeItError('You can only input offsets or an uniq_obj_id, but not both')
        if offsets is None and uniq_obj_id is None:
            raise PypeItError('You must input either offsets or a uniq_obj_id to determine the stack of '
                       'reference traces')

        # if offset is provided, we stack about the center of the slit
        if isinstance(offsets, (list, np.ndarray)):
            return self.offset_slit_cen(slitid, offsets)

        # if uniq_obj_id is provided, we stack about the trace of the object
        orders = self.coadd2d_stack.slits_list[0].slitord_id
        specobjs_list = self.coadd2d_stack.specobjs_list
        ref_trace_stack = []
        for iexp, sobjs in enumerate(specobjs_list):
            ithis = sobjs.slitorder_uniq_id_indices(uniq_obj_id[iexp], order=orders[slitid])
            ref_trace_stack.append(sobjs[ithis][0].TRACE_SPAT)
        return ref_trace_stack

    def wave_method(self):
        """
        Return the wavelength-grid method used for echelle coadd2d.

        Returns
        -------
        str
            Name of the wavelength-grid method.
        """
        return 'log10' if self.par['coadd2d']['wave_method'] is None else self.par['coadd2d']['wave_method']