"""
Module for guiding 1D Wavelength Calibration

.. include:: ../include/links.rst

"""
import inspect
import json

import numpy as np
from matplotlib import pyplot as plt

from linetools.utils import jsonify
from astropy.table import Table
from astropy.io import fits

from pypeit import msgs
from pypeit.core import arc, qa
from pypeit.core import fitting
from pypeit.core import parse
from pypeit.core.wavecal import autoid, wv_fitting, wvutils
from pypeit.core.gui.identify import Identify
from pypeit import datamodel
from pypeit import calibframe
from pypeit.core.wavecal import echelle


from IPython import embed

class WaveCalib(calibframe.CalibFrame):
    """
    Calibration frame containing the wavelength calibration.

    All of the items in the datamodel are required for instantiation, although
    they can be None (but shouldn't be)

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_wavecalib.rst

    """
    version = '1.1.2'

    # Calibration frame attributes
    calib_type = 'WaveCalib'
    calib_file_format = 'fits'

    # NOTE:
    #   - Internals are identical to the base class
    #   - Datamodel already contains CalibFrame base elements, so no need to
    #     include it here.

    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'wv_fits': dict(otype=np.ndarray, atype=wv_fitting.WaveFit,
                                 descr='WaveFit to each 1D wavelength solution'),
                 #'wv_fit2d': dict(otype=fitting.PypeItFit,
                 #                 descr='2D wavelength solution (echelle)'),
                 'wv_fit2d': dict(otype=np.ndarray, atype=fitting.PypeItFit,
                                  descr='2D wavelength solution(s) (echelle).  If there is more '
                                        'than one, they must be aligned to the separate detectors '
                                        'analyzed'),
                 'fwhm_map': dict(otype=np.ndarray, atype=fitting.PypeItFit,
                                  descr='A fit that determines the spectral FWHM at every location of every slit'),
                 'det_img': dict(otype=np.ndarray, atype=np.integer,
                                  descr='Detector image which indicates which pixel in the mosaic '
                                        'corresponds to which detector; used occasionally by '
                                        'echelle.  Currently only saved if ech_separate_2d=True'),
                 'arc_spectra': dict(otype=np.ndarray, atype=np.floating,
                                     descr='2D array: 1D extracted spectra, slit by slit '
                                           '(nspec, nslits)'),
                 'nslits': dict(otype=int,
                                descr='Total number of slits.  This can include masked slits'),
                 'spat_ids': dict(otype=np.ndarray, atype=np.integer, 
                                  descr='Slit spat_ids. Named distinctly from that in WaveFit '),
                 'ech_orders': dict(otype=np.ndarray, atype=np.integer,
                                   descr='Echelle order ID numbers.  Defined only for echelle.'),
                 'strpar': dict(otype=str, descr='Parameters as a string'),
                 'lamps': dict(otype=str,
                               descr='List of arc lamps used for the wavelength calibration')}

    def __init__(self, wv_fits=None, fwhm_map=None, nslits=None, spat_ids=None, ech_orders=None,
                 PYP_SPEC=None, strpar=None, wv_fit2d=None, arc_spectra=None, lamps=None,
                 det_img=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _bundle(self):
        """
        Override base class function to write one HDU per image.  Any extras are
        in the HDU header of the primary image.

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension.
        """
        _d = []

        # Spat_ID are always first
        if self.spat_ids is None:
            msgs.error('Cannot write WaveCalib without spat_ids!')
        _d.append(dict(spat_ids=self.spat_ids))
        # Echelle orders
        if self.ech_orders is not None:
            _d.append(dict(ech_orders=self.ech_orders))

        # Rest of the datamodel
        for key in self.keys():
            if key in ['spat_ids', 'ech_orders']:
                continue
            # Skip None
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray and \
                key not in ['wv_fits', 'wv_fit2d', 'fwhm_map']:
                _d.append({key: self[key]})
            # TODO: Can we put all the WAVEFIT and PYPEITFIT at the end of the
            # list of HDUs?  This would mean ARC_SPECTRA is always in the same
            # extension number, regardless of the number of slits.
            elif key == 'wv_fits':
                for ss, wv_fit in enumerate(self[key]):
                    # TODO: Are we writing empty extensions if any of the
                    # elements of self[key] are None?  If so, is this required
                    # behavior?  Why?
                    # Naming
                    # TODO: Shouldn't this name match the dkey below?
                    #   Oddly enough it is coded correctly below
                    dkey = 'WAVEFIT-{}'.format(self.spat_ids[ss])
                    # Generate a dummy?
                    if wv_fit is None:
                        echorder = self.ech_orders[ss] if self.ech_orders is not None else None
                        kwv_fit = wv_fitting.WaveFit(self.spat_ids[ss], ech_order=echorder)
                    else:
                        kwv_fit = wv_fit
                    # This is required to deal with a single HDU WaveFit() bundle
                    if kwv_fit.pypeitfit is None:
                        dkey = 'SPAT_ID-{}_WAVEFIT'.format(self.spat_ids[ss])
                    # Save
                    _d.append({dkey: kwv_fit})
            elif key == 'wv_fit2d':
                for ss, wv_fit2d in enumerate(self[key]):
                    dkey = f'WAVE2DFIT-{ss}'
                    _d.append({dkey: wv_fit2d})
            elif key == 'fwhm_map':
                for ss, fwhm_fit in enumerate(self[key]):
                    dkey = 'SPAT_ID-{}_FWHMFIT'.format(self.spat_ids[ss])
                    # Generate a dummy?
                    if fwhm_fit is None:
                        _fwhm_fit = fitting.PypeItFit()
                    else:
                        _fwhm_fit = fwhm_fit
                    # Save
                    _d.append({dkey: _fwhm_fit})
            else: # Add to header of the spat_id image
                _d[0][key] = self[key]
        # Return
        return _d

    @classmethod
    def from_hdu(cls, hdu, chk_version=True, **kwargs):
        """
        Instantiate the object from an HDU extension.

        This overrides the base-class method. Overriding this method is
        preferrable to overriding the ``_parse`` method because it makes it
        easier to deal with the :class:`~pypeit.datamodel.DataContainer` nesting
        of this object.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            chk_version (:obj:`bool`, optional):
                If True, raise an error if the datamodel version or
                type check failed. If False, throw a warning only.
            kwargs (:obj:`dict`, optional):
                Used for consistency with base class. Ignored.
        """
        # Run the default parser to get most of the data. This won't parse the
        # extensions with the WAVEFIT, WAVE2DFIT, or FWHMFIT results.
        d, version_passed, type_passed, parsed_hdus = super()._parse(hdu)
        # Check
        cls._check_parsed(version_passed, type_passed, chk_version=chk_version)

        # Get the list of all extensions
        ext = [h.name for h in hdu] if isinstance(hdu, fits.HDUList) else [hdu.name]

        # Get the SPAT_IDs
        if 'SPAT_IDS' in parsed_hdus:
            # Use the ones parsed above
            spat_ids = d['spat_ids']
        else:
            # This line parses all the spat_ids from the extension names,
            # filters out any None values from the list, and gets the unique set
            # of integers
            spat_ids = np.unique(list(filter(None.__ne__,
                            [wv_fitting.WaveFit.parse_spatid_from_hduext(e) for e in ext])))

        # Parse all the WAVEFIT extensions
        wave_fits = []
        for spat_id in spat_ids:
            _ext = wv_fitting.WaveFit.hduext_prefix_from_spatid(spat_id)+'WAVEFIT'
            if _ext not in ext:
                continue
            # TODO: I (KBW) don't think we should be writing empty HDUs
            if len(hdu[_ext].data) == 0:
                wave_fits += [wv_fitting.WaveFit(hdu[_ext].header['SPAT_ID'],
                                                 ech_order=hdu[_ext].header.get('ECH_ORDER'))]
            else:
                wave_fits += [wv_fitting.WaveFit.from_hdu(hdu, spat_id=spat_id,
                                                          chk_version=chk_version)]
        if len(wave_fits) > 0:
            d['wv_fits'] = np.asarray(wave_fits)
                
        # Parse all the WAVE2DFIT extensions
        # TODO: It would be good to have the WAVE2DFIT extensions follow the
        # same naming convention as the WAVEFIT extensions...
        wave2d_fits = [fitting.PypeItFit.from_hdu(hdu[e], chk_version=chk_version)
                            for e in ext if 'WAVE2DFIT' in e]
        if len(wave2d_fits) > 0:
            d['wv_fit2d'] = np.asarray(wave2d_fits)

        # Parse all the FWHMFIT extensions
        fwhm_fits = []
        for _ext in ext:
            if 'FWHMFIT' not in _ext:
                continue
            # TODO: I (KBW) don't think we should be writing empty HDUs
            fwhm_fits += [fitting.PypeItFit() if len(hdu[_ext].data) == 0 \
                            else fitting.PypeItFit.from_hdu(hdu[_ext], chk_version=chk_version)]
        if len(fwhm_fits) > 0:
            d['fwhm_map'] = np.asarray(fwhm_fits)

        # Instantiate the object
        self = cls.from_dict(d=d)
        # This is a CalibFrame, so parse the relevant keys for the naming system
        self.calib_keys_from_header(hdu[parsed_hdus[0]].header)
        # Return the constructed object
        return self

    @property
    def par(self):
        return json.loads(self.strpar)

    def chk_synced(self, slits):
        """
        Confirm the slits in WaveCalib are aligned to that in SlitTraceSet

        Barfs if not

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`):

        """
        if not np.array_equal(self.spat_ids, slits.spat_id):
            msgs.error('Your wavelength solutions are out of sync with your slits.  Remove '
                       'Calibrations and restart from scratch.')

    def build_fwhmimg(self, tilts, slits, initial=False, spat_flexure=None):
        """
        Generates an image of the instrument spectral FWHM (units=pixels) at every pixel on the detector.

        Args:
            tilts (`numpy.ndarray`_):
                Image holding tilts
            slits (:class:`pypeit.slittrace.SlitTraceSet`):
                Properties of the slits
            initial (bool, optional):
                If True, the initial slit locations will be used. Otherwise, the tweaked edges will be used.
            spat_flexure (float, optional):
                Spatial flexure correction in pixels.

        Returns:
            `numpy.ndarray`_: The spectral FWHM image.
        """
        # Check spatial flexure type
        if (spat_flexure is not None) and (not isinstance(spat_flexure, float)):
            msgs.error("Spatial flexure must be None or float")
        # Generate the slit mask and slit edges - pad slitmask by 1 for edge effects
        slitmask = slits.slit_img(pad=1, initial=initial, flexure=spat_flexure)
        slits_left, slits_right, _ = slits.select_edges(initial=initial, flexure=spat_flexure)
        # Build a map of the spectral FWHM
        fwhmimg = np.zeros(tilts.shape)
        for sl, spat_id in enumerate(slits.spat_id):
            this_mask = slitmask == spat_id
            spec, spat = np.where(this_mask)
            spat_loc = (spat - slits_left[spec, sl]) / (slits_right[spec, sl] - slits_left[spec, sl])
            fwhmimg[this_mask] = self.fwhm_map[sl].eval(spec, spat_loc)
        return fwhmimg

    def build_waveimg(self, tilts, slits, spat_flexure=None, spec_flexure=None):
        """
        Main algorithm to build the wavelength image

        Only applied to good slits, which means any non-flagged or flagged
         in the exclude_for_reducing list

        Args:
            tilts (`numpy.ndarray`_):
                Image holding tilts
            slits (:class:`pypeit.slittrace.SlitTraceSet`):
                Properties of the slits
            spat_flexure (float, optional):
                Spatial flexure correction in pixels.
            spec_flexure (float, `numpy.ndarray`_, optional):
                Spectral flexure correction in pixels. If a float,
                the same spectral flexure correction will be applied
                to all slits. If a numpy array, the length of the
                array should be the same as the number of slits. The
                value of each element is the spectral shift in pixels
                to be applied to each slit.

        Returns:
            `numpy.ndarray`_: The wavelength image.
        """
        # Check spatial flexure type
        if (spat_flexure is not None) and (not isinstance(spat_flexure, float)):
            msgs.error("Spatial flexure must be None or float")
        # Check spectral flexure type
        if spec_flexure is None: spec_flex = np.zeros(slits.nslits)
        elif isinstance(spec_flexure, float): spec_flex = spec_flexure*np.ones(slits.nslits)
        elif isinstance(spec_flexure, np.ndarray):
            spec_flex = spec_flexure.copy()
            assert(spec_flexure.size == slits.nslits)
        spec_flex /= (slits.nspec - 1)

        # Setup
        #ok_slits = slits.mask == 0
#        bpm = slits.mask.astype(bool)
#        bpm &= np.logical_not(slits.bitmask.flagged(slits.mask, flag=slits.bitmask.exclude_for_reducing))
        bpm = slits.bitmask.flagged(slits.mask, and_not=slits.bitmask.exclude_for_reducing)
        ok_slits = np.logical_not(bpm)
        #
        image = np.zeros_like(tilts)
        slitmask = slits.slit_img(flexure=spat_flexure, exclude_flag=slits.bitmask.exclude_for_reducing)

        # Separate detectors for the 2D solutions?
        if self.par['ech_separate_2d']:
            # Error checking
            if self.det_img is None:
                msgs.error("This WaveCalib object was not generated with ech_separate_2d=True")
            # Grab slit_img
            slit_img = slits.slit_img()
        
        # Unpack some 2-d fit parameters if this is echelle
        for islit in np.where(ok_slits)[0]:
            slit_spat = slits.spat_id[islit]
            thismask = (slitmask == slit_spat)
            if not np.any(thismask):
                msgs.error("Something failed in wavelengths or masking..")
            if self.par['echelle'] and self.par['ech_2dfit']:
                # evaluate solution --
                if self.par['ech_separate_2d']:
                    ordr_det = slits.det_of_slit(
                        slit_spat, self.det_img,
                        slit_img=slit_img)
                    # There are ways for this to go sour..
                    #  if the seperate solutions are not aligned with the detectors
                    #  or if one reruns with a different number of detectors
                    #  without regeneating
                    #  But that would be bad practice
                    idx_fit2d = ordr_det-1  
                else:
                    idx_fit2d = 0
                image[thismask] = self.wv_fit2d[idx_fit2d].eval(
                    tilts[thismask] + spec_flex[islit], 
                    x2=np.full_like(tilts[thismask], 
                                    slits.ech_order[islit]))
                image[thismask] /= slits.ech_order[islit]
            else:
                iwv_fits = self.wv_fits[islit]
                image[thismask] = iwv_fits.pypeitfit.eval(
                    tilts[thismask] + spec_flex[islit])
        # Return
        return image

    def wave_diagnostics(self, print_diag=False):
        """
        Create a table with wavecalib diagnostics

        Args:
            print_diag (:obj:`bool`, optional):
                If True, the diagnostic table is printed to screen
        Returns:
            `astropy.table.Table`_: wavecalib diagnostics table

        """
        # wavelength range of calibrated arc spectra
        minWave = np.array([0 if wvfit.wave_soln is None else wvfit.wave_soln[0] for wvfit in self.wv_fits])
        maxWave = np.array([0 if wvfit.wave_soln is None else wvfit.wave_soln[-1] for wvfit in self.wv_fits])

        # wavelength range of fitted ID'd lines
        lines_wmin = np.array([0 if wvfit is None or wvfit.pypeitfit is None else
                              wvfit.wave_fit[wvfit.pypeitfit.gpm == 1][0] for wvfit in self.wv_fits])
        lines_wmax = np.array([0 if wvfit is None or wvfit.pypeitfit is None else
                              wvfit.wave_fit[wvfit.pypeitfit.gpm == 1][-1] for wvfit in self.wv_fits])

        # wavelength coverage of fitted ID'd lines
        lines_waverange = lines_wmax - lines_wmin
        spec_waverange = maxWave - minWave
        lines_cov = [0 if spec_waverange[i] == 0 else
                     lines_waverange[i] / spec_waverange[i] * 100 for i in range(self.wv_fits.size)]

        # Generate a table
        diag = Table()
        # Slit number
        diag['N.'] = np.arange(self.wv_fits.size)
        diag['N.'].format = 'd'
        # spat_id or order number
        diag['SpatOrderID'] = [wvfit.spat_id for wvfit in self.wv_fits] if self.ech_orders is None \
            else [wvfit.ech_order for wvfit in self.wv_fits]
        # Central wave, delta wave
        diag['minWave'] = minWave
        diag['minWave'].format = '0.1f'
        diag['Wave_cen'] = [0 if wvfit.cen_wave is None else wvfit.cen_wave for wvfit in self.wv_fits]
        diag['Wave_cen'].format = '0.1f'
        diag['maxWave'] = maxWave
        diag['maxWave'].format = '0.1f'
        diag['dWave'] = [0 if wvfit.cen_disp is None else wvfit.cen_disp for wvfit in self.wv_fits]
        diag['dWave'].format = '0.3f'
        # Number of good lines
        diag['Nlin'] = [0 if wvfit.pypeitfit is None else np.sum(wvfit.pypeitfit.gpm) for wvfit in self.wv_fits]
        diag['IDs_Wave_range'] = ['{:9.3f} - {:9.3f}'.format(lines_wmin[i], lines_wmax[i]) for i in range(self.wv_fits.size)]
        diag['IDs_Wave_cov(%)'] = lines_cov
        diag['IDs_Wave_cov(%)'].format = '0.1f'
        # FWHM
        diag['measured_fwhm'] = [0. if wvfit.fwhm is None else wvfit.fwhm for wvfit in self.wv_fits]
        diag['measured_fwhm'].format = '0.1f'
        # RMS
        diag['RMS'] = [0 if wvfit.rms is None else wvfit.rms for wvfit in self.wv_fits]
        diag['RMS'].format = '0.3f'
        if print_diag:
            # Print it
            print('')
            diag.pprint_all()
        return diag


class BuildWaveCalib:
    """
    Class to guide wavelength calibration

    Args:
        msarc (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Arc image, created by the ArcImage class
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit edges
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        par (:class:`~pypeit.par.pypeitpar.WavelengthSolutionPar`):
            The parameters used for the wavelength solution.
            Uses ``['calibrations']['wavelengths']``.
        meta_dict (dict, optional):
            Dictionary containing meta information required for wavelength
            calibration. Specifically for non-fixed format echelles this dict
            must contain the following keys:

               - ``'echangle'``:  the echelle angle
               - ``'xdangle'``: the cross-disperser angle
               - ``'dispmame'``: the disperser name

        det (int, optional):
            Detector number
        msbpm (`numpy.ndarray`_, optional):
            Bad pixel mask image
        qa_path (str, optional):
            For QA

    Attributes:
        steps (list):
            List of the processing steps performed
        wv_calib (dict):
            Primary output.  Keys 0, 1, 2, 3 are solution for individual
            previously slit steps
        arccen (`numpy.ndarray`_):
            (nwave, nslit) Extracted arc(s) down the center of the slit(s)
        maskslits (`numpy.ndarray`_):
            Slits to ignore because they were not extracted. WARNING:
            Outside of this Class, it is best to regenerate the mask
            using  make_maskslits()
        gpm (`numpy.ndarray`_):
            Good pixel mask
            Eventually, we might attach this to self.msarc although that would then
            require that we write it to disk with self.msarc.image
        nonlinear_counts (float):
            Specifies saturation level for the arc lines
        wvc_bpm (`numpy.ndarray`_):
            Mask for slits attempted to have a wv_calib solution
    """

    # TODO: Is this used anywhere?
    frametype = 'wv_calib'

    def __init__(self, msarc, slits, spectrograph, par, lamps,
                 meta_dict=None, det=1, qa_path=None, msbpm=None):

        # TODO: This should be a stop-gap to avoid instantiation of this with
        # any Nones.
        if None in [msarc, slits, spectrograph, par, lamps]:
            msgs.error('CODING ERROR: Cannot instantiate BuildWaveCalib with Nones.')

        # Required parameters
        self.msarc = msarc
        self.binspectral = parse.parse_binning(self.msarc.detector.binning)[0]
        self.slits = slits
        self.spectrograph = spectrograph
        self.par = par
        self.lamps = lamps
        self.meta_dict = meta_dict

        # Optional parameters
        self.bpm = self.msarc.select_flag(flag='BPM') if msbpm is None else msbpm.astype(bool)
        if self.bpm.shape != self.msarc.shape:
            msgs.error('Bad-pixel mask is not the same shape as the arc image.')
        self.qa_path = qa_path
        self.det = det

        # Attributes
        self.steps = []     # steps executed
        self.wv_calib = {}  # main output
        self.arccen = None  # central arc spectrum

        # Get the non-linear count level
        # TODO: This is currently hacked to deal with Mosaics
        try:
            self.nonlinear_counts = self.msarc.detector.nonlinear_counts()
        except:
            self.nonlinear_counts = 1e10

        # --------------------------------------------------------------
        # TODO: Build another base class that does these things for both
        # WaveTilts and WaveCalib?

        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with arcimages that
        # have a different binning then the trace images used to defined
        # the slits
        if self.slits is not None and self.msarc is not None:
            # Redo?
            if self.par['redo_slits'] is not None:
                if self.par['echelle'] and self.slits.ech_order is not None:
                    idx = np.in1d(self.slits.ech_order, self.par['redo_slits'])
                    # Turn off mask
                    self.slits.mask[idx] = self.slits.bitmask.turn_off(
                            self.slits.mask[idx], 'BADWVCALIB')
                else:
                    idx = np.in1d(self.slits.spat_id, self.par['redo_slits'])
                    self.slits.mask[idx] = self.slits.bitmask.turn_off(
                            self.slits.mask[idx], 'BADWVCALIB')

            # Load up slits
            # TODO -- Allow for flexure
            slits_left, slits_right, mask = self.slits.select_edges(initial=True, flexure=None)  # Grabs all, init slits + flexure
            self.orders = self.slits.ech_order  # Can be None
#            self.spat_coo = self.slits.spatial_coordinates()  # All slits, even masked
            # Internal mask for failed wv_calib analysis
            # TODO -- Allow for an option to re-attempt those previously flagged as BADWVCALIB?
            self.wvc_bpm = np.invert(mask == 0)
            ## We want to keep the 'BOXSLIT', which mask value is 2. But we don't want to keep 'BOXSLIT'
            ## with other bad flag (for which the mask value would be > 2)
            #self.wvc_bpm = mask > 2
            self.wvc_bpm_init = self.wvc_bpm.copy()
            # Slitmask -- Grabs only unmasked, initial slits
            #self.slitmask_science = self.slits.slit_img(initial=True, flexure=None, exclude_flag=['BOXSLIT'])
            self.slitmask_science = self.slits.slit_img(initial=True, flexure=None)
            # Resize
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.msarc.image.shape
            # slitcen is padded to include slits that may be masked, for convenience in coding downstream
            self.slits_left = arc.resize_slits2arc(self.shape_arc, self.shape_science, slits_left)
            self.slits_right = arc.resize_slits2arc(self.shape_arc, self.shape_science, slits_right)
            self.slitcen = (self.slits_left+self.slits_right)/2
            #self.slitcen = arc.resize_slits2arc(self.shape_arc, self.shape_science, (self.slits_left+self.slits_right)/2)
            self.slitmask = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            # Mask
            # TODO: The bpm defined above is already a boolean and cannot be None.
            gpm = np.logical_not(self.bpm)
            self.gpm = arc.resize_mask2arc(self.shape_arc, gpm)
            # We want even the saturated lines in full_template for the cross-correlation
            #   They will be excised in the detect_lines() method on the extracted arc
            if self.par['method'] != 'full_template':
                self.gpm &= self.msarc.image < self.nonlinear_counts

        else:
            self.orders = None
            self.wvc_bpm = None
            self.wvc_bpm_init = None
            self.slitmask_science = None
            self.shape_science = None
            self.shape_arc = None
            self.slitcen = None
            self.slits_left = None
            self.slits_right = None
            self.slitmask = None
            self.gpm = None

    def build_wv_calib(self, arccen, method, skip_QA=False,
                       prev_wvcalib=None):
        """
        Main routine to generate the wavelength solutions in a loop over slits
        Wrapper to arc.simple_calib or arc.calib_with_arclines

        self.maskslits is updated for slits that fail

        Args:
            method : str
              'simple' -- arc.simple_calib
              'arclines' -- arc.calib_with_arclines
              'holy-grail' -- wavecal.autoid.HolyGrail
              'reidentify' -- wavecal.auotid.ArchiveReid
              'identify' -- wavecal.identify.Identify
              'full_template' -- wavecal.auotid.full_template
            skip_QA (bool, optional)
            prev_wvcalib (WaveCalib, optional):
                Previous wavelength calibration

        Returns:
            dict:  self.wv_calib
        """
        # Obtain a list of good slits
        ok_mask_idx = np.where(np.invert(self.wvc_bpm))[0]

        # print to screen the slit widths if maskdef_designtab is available
        if self.slits.maskdef_designtab is not None:
            msgs.info("Slit widths (arcsec): {}".format(np.round(self.slits.maskdef_designtab['SLITWID'].data, 2)))

        # Generate a map of the instrumental spectral FWHM
        # TODO nsample should be a parameter
        fwhm_map = autoid.map_fwhm(self.msarc.image, self.gpm, self.slits_left, self.slits_right, self.slitmask,
                                   nsample=10, slit_bpm=self.wvc_bpm, specord=self.par['fwhm_spec_order'],
                                   spatord=self.par['fwhm_spat_order'])
        # Calculate the typical spectral FWHM down the centre of the slit
        measured_fwhms = np.zeros(arccen.shape[1], dtype=object)
        for islit in range(arccen.shape[1]):
            if islit not in ok_mask_idx:
                continue
            # Measure the spectral FWHM (in pixels) at the midpoint of the slit
            # (i.e. the midpoint in both the spectral and spatial directions)
            measured_fwhms[islit] = fwhm_map[islit].eval(self.msarc.image.shape[0]//2, 0.5)

        # Save for redo's
        self.measured_fwhms = measured_fwhms

        # Obtain calibration for all slits
        if method == 'holy-grail':
            # Sometimes works, sometimes fails
            arcfitter = autoid.HolyGrail(arccen, self.lamps, par=self.par, 
                                         ok_mask=ok_mask_idx,
                                         measured_fwhms=self.measured_fwhms,
                                         nonlinear_counts=self.nonlinear_counts,
                                         spectrograph=self.spectrograph.name)
            patt_dict, final_fit = arcfitter.get_results()
        elif method == 'identify':
            raise NotImplementedError('method = identify not yet implemented')
            final_fit = {}
            # Manually identify lines
            msgs.info("Initializing the wavelength calibration tool")
            #embed(header='line 222 wavecalib.py')
            for slit_idx in ok_mask_idx:
                arcfitter = Identify.initialise(arccen, self.lamps, self.slits, slit=slit_idx, par=self.par)
                final_fit[str(slit_idx)] = arcfitter.get_results()
                arcfitter.store_solution(final_fit[str(slit_idx)], "", self.binspectral,
                                         specname=self.spectrograph.name,
                                         gratname="UNKNOWN", dispangl="UNKNOWN")
        elif method == 'reidentify':
            # Now preferred
            # Slit positions
            arcfitter = autoid.ArchiveReid(
                arccen, self.lamps, self.par,
                ech_fixed_format=self.spectrograph.ech_fixed_format, 
                ok_mask=ok_mask_idx,
                measured_fwhms=self.measured_fwhms,
                orders=self.orders,
                nonlinear_counts=self.nonlinear_counts)
            patt_dict, final_fit = arcfitter.get_results()

            # Grab arxiv for redo later?
            if self.par['echelle']: 
                # Hold for later usage
                self.wave_soln_arxiv, self.arcspec_arxiv = arcfitter.get_arxiv(self.orders)
                self.arccen = arccen
        elif method == 'full_template':
            # Now preferred
            if self.binspectral is None:
                msgs.error("You must specify binspectral for the full_template method!")
            final_fit, order_vec = autoid.full_template(arccen, self.lamps, self.par, ok_mask_idx, self.det,
                                             self.binspectral, slit_ids=self.slits.slitord_id,
                                             measured_fwhms=self.measured_fwhms,
                                             nonlinear_counts=self.nonlinear_counts,
                                             nsnippet=self.par['nsnippet'], 
                                             x_percentile=self.par['cc_percent_ceil'])

            # Grab arxiv for redo later?
            if self.par['echelle']: 
                # Hold for later usage
                self.slits.ech_order = order_vec[:self.slits.nslits]
                self.arccen = arccen
        elif self.par['method'] == 'echelle':
            # Echelle calibration files
            angle_fits_file, composite_arc_file = self.spectrograph.get_echelle_angle_files()

            # Identify the echelle orders
            msgs.info("Finding the echelle orders")
            order_vec, wave_soln_arxiv, arcspec_arxiv = echelle.identify_ech_orders(
                    arccen, self.meta_dict['echangle'],
                    self.meta_dict['xdangle'],
                    self.meta_dict['dispname'],
                    angle_fits_file,
                    composite_arc_file,
                    pad=self.par['echelle_pad'],
                    cc_percent_ceil = self.par['cc_percent_ceil'], debug=False)
            # Put the order numbers in the slit object
            self.slits.ech_order = order_vec
            msgs.info(f"The observation covers the following orders: {order_vec}")

            patt_dict, final_fit = autoid.echelle_wvcalib(
                arccen, order_vec, arcspec_arxiv, wave_soln_arxiv,
                self.lamps, self.par, ok_mask=ok_mask_idx,
                measured_fwhms=self.measured_fwhms,
                nonlinear_counts=self.nonlinear_counts,
                debug_all=False,
                redo_slits=np.atleast_1d(self.par['redo_slits']) if self.par['redo_slits'] is not None else None)

            # Save as internals in case we need to redo
            self.wave_soln_arxiv = wave_soln_arxiv
            self.arcspec_arxiv = arcspec_arxiv
            self.arccen = arccen

        else:
            msgs.error('Unrecognized wavelength calibration method: {:}'.format(method))

        # Build the DataContainer
        if self.par['redo_slits'] is not None:
            # If we are only redoing slits, we start from the
            #  previous wv_calib and update only the (good) redone slits
            self.wv_calib = prev_wvcalib
            # Update/reset items
            self.wv_calib.arc_spectra = arccen
            # Save the new fits (if they meet tolerance)
            for key in final_fit.keys():
                idx = int(key)
                # get FWHM for this slit
                fwhm = autoid.set_fwhm(self.par, measured_fwhm=self.measured_fwhms[idx], verbose=True)
                # get rms threshold for this slit
                wave_rms_thresh = round(self.par['rms_thresh_frac_fwhm'] * fwhm, 3)
                if final_fit[key]['rms'] < wave_rms_thresh:
                    self.wv_calib.wv_fits[idx] = final_fit[key]
                    self.wv_calib.wv_fits[idx].spat_id = self.slits.spat_id[idx]
                    self.wv_calib.wv_fits[idx].ech_order = self.slits.ech_order[idx] if self.slits.ech_order is not None else None
                    self.wv_calib.wv_fits[idx].fwhm = self.measured_fwhms[idx]
        else: # Generate the DataContainer from scratch
            # Loop on WaveFit items
            tmp = []
            for idx in range(self.slits.nslits):
                echorder = self.slits.ech_order[idx] if self.slits.ech_order is not None else None
                item = final_fit.pop(str(idx))
                if item is None:  # Add an empty WaveFit
                    tmp.append(wv_fitting.WaveFit(self.slits.spat_id[idx], ech_order=echorder))
                else:
                    # This is for I/O naming
                    item.spat_id = self.slits.spat_id[idx]
                    item.ech_order = echorder
                    # add measured fwhm
                    item['fwhm'] = measured_fwhms[idx]
                    tmp.append(item)
            self.wv_calib = WaveCalib(wv_fits=np.asarray(tmp),
                                      fwhm_map=fwhm_map,
                                      arc_spectra=arccen,
                                      nslits=self.slits.nslits,
                                      spat_ids=self.slits.spat_id,
                                      ech_orders=self.slits.ech_order,
                                      PYP_SPEC=self.spectrograph.name,
                                      lamps=','.join(self.lamps))
        # Inherit the calibration frame naming from self.msarc
        # TODO: Should throw an error here if these calibration frame naming
        # elements are not defined by self.msarc...
        self.wv_calib.copy_calib_internals(self.msarc)

        # Update mask
        self.update_wvmask()

        # QA
        if not skip_QA:
            ok_mask_idx = np.where(np.invert(self.wvc_bpm))[0]
            for slit_idx in ok_mask_idx:
                msgs.info(f"Preparing wavelength calibration QA for slit {slit_idx+1}/{self.slits.nslits}")
                # Obtain the output QA name for the wavelength solution
                outfile = qa.set_qa_filename(
                    self.wv_calib.calib_key, 'arc_fit_qa', 
                    slit=self.slits.slitord_id[slit_idx],
                    out_dir=self.qa_path)
                # Save the wavelength solution fits
                autoid.arc_fit_qa(self.wv_calib.wv_fits[slit_idx], 
                                  title=f'Arc Fit QA for slit/order: {self.slits.slitord_id[slit_idx]}',
                                  outfile=outfile, log=self.par['qa_log'])

                # Obtain the output QA name for the spectral resolution map
                outfile_fwhm = qa.set_qa_filename(self.wv_calib.calib_key, 'arc_fwhm_qa',
                                                  slit=self.slits.slitord_id[slit_idx],
                                                  out_dir=self.qa_path)
                # Save the wavelength solution fits
                autoid.arc_fwhm_qa(self.wv_calib.fwhm_map[slit_idx],
                                   self.slits.slitord_id[slit_idx], self.slits.slitord_txt,
                                   outfile=outfile_fwhm)

        # Return
        self.steps.append(inspect.stack()[0][3])
        return self.wv_calib

    def redo_echelle_orders(self, bad_orders:np.ndarray, dets:np.ndarray, order_dets:np.ndarray,
                            bad_orders_maxfrac:float=0.1, frac_rms_thresh:float=1.1):
        """ Attempt to redo the wavelength calibration for a set 
        of bad echelle orders

        Args:
            bad_orders (np.ndarray): Array of bad order numbers
            dets (np.ndarray): detectors of the spectrograph
                Multiple numbers for mosaic (typically)
            order_dets (np.ndarray): Orders on each detector
            bad_orders_maxfrac (float): Maximum fraction of bad orders
                in a detector for attempting a refit
            frac_rms_thresh (float): Fractional change in the RMS threshold
                for accepting a refit

        Returns:
            bool: True if any of the echelle orders were successfully redone
        """

        # Make this outside the for loop..
        #bad_orders = self.slits.ech_order[np.where(bad_rms)[0]]
        fixed = False

        for idet in range(len(dets)):
            in_det = np.in1d(bad_orders, order_dets[idet])
            if not np.any(in_det):
                continue
            msgs.info(f"Attempting to refit bad orders in detector={dets[idet]}")
            # Are there few enough?
            max_bad = int(len(order_dets[idet])*bad_orders_maxfrac)
            if np.sum(in_det) > max_bad:
                msgs.warn(f"Too many bad orders in detector={dets[idet]} to attempt a refit.")
                continue
            # Loop
            for order in bad_orders[in_det]:
                iord = np.where(self.slits.ech_order == order)[0][0]
                # Predict the wavelengths
                nspec = self.arccen.shape[0]
                spec_vec_norm = np.linspace(0,1,nspec)
                wv_order_mod = self.wv_calib.wv_fit2d[idet].eval(spec_vec_norm, 
                                    x2=np.ones_like(spec_vec_norm)*order)/order
                # get FWHM for this order
                fwhm = autoid.set_fwhm(self.par, measured_fwhm=self.measured_fwhms[iord], verbose=True)
                # get rms threshold for this order
                wave_rms_thresh = round(self.par['rms_thresh_frac_fwhm'] * fwhm, 3)
                # Link me
                tcent, spec_cont_sub, patt_dict_slit, tot_llist = autoid.match_to_arxiv(
                    self.lamps, self.arccen[:,iord], wv_order_mod, 
                    self.arcspec_arxiv[:, iord],  self.wave_soln_arxiv[:,iord],
                    self.par['nreid_min'], 
                match_toler=self.par['match_toler'], 
                nonlinear_counts=self.nonlinear_counts, 
                sigdetect=wvutils.parse_param(self.par, 'sigdetect', iord),
                fwhm=fwhm)

                if not patt_dict_slit['acceptable']:
                    msgs.warn(f"Order {order} is still not acceptable after attempt to reidentify.")
                    continue

                # Fit me -- RMS may be too high again
                n_final = wvutils.parse_param(self.par, 'n_final', iord)
                # TODO - Make this cheaper
                final_fit = wv_fitting.fit_slit(
                    self.arccen[:,iord], patt_dict_slit, tcent, tot_llist,
                    match_toler=self.par['match_toler'], 
                    func=self.par['func'], 
                    n_first=self.par['n_first'],
                    sigrej_first=self.par['sigrej_first'],
                    n_final=n_final, 
                    sigrej_final=self.par['sigrej_final'])
                msgs.info(f"New RMS for redo of order={order}: {final_fit['rms']}")

                # Keep?
                if final_fit['rms'] < frac_rms_thresh*wave_rms_thresh:
                    msgs.info('Updating wavelength solution.')
                    # TODO -- This is repeated from build_wv_calib()
                    #  Would be nice to consolidate
                    # QA
                    outfile = qa.set_qa_filename(
                        self.wv_calib.calib_key, 'arc_fit_qa', 
                        slit=order,
                        out_dir=self.qa_path)
                    autoid.arc_fit_qa(final_fit,
                                    title=f'Arc Fit QA for slit/order: {order}',
                                    outfile=outfile, log=self.par['qa_log'])
                    # This is for I/O naming
                    final_fit.spat_id = self.slits.spat_id[iord]
                    final_fit.ech_order = self.slits.ech_order[iord] if self.slits.ech_order is not None else None
                    final_fit.fwhm = self.measured_fwhms[iord]
                    # Save the wavelength solution fits
                    self.wv_calib.wv_fits[iord] = final_fit
                    self.wvc_bpm[iord] = False
                    fixed = True
                else:
                    msgs.warn(f'New RMS is too high (>{frac_rms_thresh}xRMS threshold). '
                              f'Not updating wavelength solution.')
        #
        return fixed

    def echelle_2dfit(self, wv_calib, debug=False, skip_QA=False):
        """
        Fit a two-dimensional wavelength solution for echelle data.

        Primarily a wrapper for :func:`pypeit.core.arc.fit2darc`,
        using data unpacked from the ``wv_calib`` dictionary.

        Parameters
        ----------
        wv_calib : :class:`pypeit.wavecalib.WaveCalib`
            Wavelength calibration object
        debug : :obj:`bool`, optional
            Show debugging info
        skip_QA : :obj:`bool`, optional
            Flag to skip construction of the nominal QA plots.

        Returns
        -------
        fit2ds : list of :class:`pypeit.fitting.PypeItFit`
            Contains information from 2-d fit.  Frequently a list of 1 fit.  The
            main exception is for a mosaic when one sets
            ``echelle_separate_2d=True``.
        dets : list
            List of integers for the detector numbers.
        save_order_dets: list
            List of integer lists providing list of the orders.
        """
        if self.spectrograph.pypeline != 'Echelle':
            msgs.error('Cannot execute echelle_2dfit for a non-echelle spectrograph.')

        msgs.info('Fitting 2-d wavelength solution for echelle....')

        # Obtain a list of good slits
        ok_mask_idx = np.where(np.logical_not(self.wvc_bpm))[0]
        ok_mask_order = self.slits.slitord_id[ok_mask_idx]
        nspec = self.msarc.image.shape[0]

        # Prep
        if self.par['ech_separate_2d']:
            slit_img = self.slits.slit_img()
            # Grab the detectors in the mosaice (1-based indexing)
            dets = np.unique(self.msarc.det_img)
            dets = dets[dets > 0]
        else:
            # The value here is irrelevant
            dets = [1]

        # Loop on detectors
        fit2ds = []
        save_order_dets = []
        for idet in dets:
            order_in_dets = []
            msgs.info('Fitting detector {:d}'.format(idet))
            # Init
            all_wave = np.array([], dtype=float)
            all_pixel = np.array([],dtype=float)
            all_order = np.array([],dtype=float)

            # Loop to grab the good orders
            for ii in range(wv_calib.nslits):
                iorder = self.slits.ech_order[ii]

                # Separate detector analysis?
                if self.par['ech_separate_2d']:
                    spat_id = wv_calib.spat_ids[ii]
                    # What is the most common detector for this order?
                    ordr_det = self.slits.det_of_slit(
                        spat_id, self.msarc.det_img,
                        slit_img=slit_img)
                    # Correct detector?
                    if ordr_det != idet:
                        continue

                # Need to record this whether or not it is ok
                order_in_dets.append(iorder)                                                        

                # Is it ok?
                if iorder not in ok_mask_order:
                    continue

                # Slurp
                mask_now = wv_calib.wv_fits[ii].pypeitfit.bool_gpm
                all_wave = np.append(all_wave, wv_calib.wv_fits[ii]['wave_fit'][mask_now])
                all_pixel = np.append(all_pixel, wv_calib.wv_fits[ii]['pixel_fit'][mask_now])
                all_order = np.append(all_order, np.full_like(wv_calib.wv_fits[ii]['pixel_fit'][mask_now],
                                                            float(iorder)))

            # Fit
            if len(all_order) < 2:
                msgs.warn(f"Fewer than 2 orders to fit for detector {idet}.  Skipping")
                save_order_dets.append([])
                # Add a dummy fit
                fit2ds.append(fitting.PypeItFit())
                continue

            fit2d = arc.fit2darc(all_wave, all_pixel, all_order, nspec,
                                    nspec_coeff=self.par['ech_nspec_coeff'],
                                    norder_coeff=self.par['ech_norder_coeff'],
                                    sigrej=self.par['ech_sigrej'], debug=debug)
            # Save
            save_order_dets.append(order_in_dets)
            fit2ds.append(fit2d)
            self.steps.append(inspect.stack()[0][3])

            # QA
            if not skip_QA:
                if wv_calib.calib_key is None:
                    msgs.warn('WaveCalib object provided does not have a defined calibration '
                              'key.  The QA files will not include this key in the file name, '
                              'meaning that existing QA files may be overwritten.')
                    calib_key = '' 
                else:
                    calib_key = wv_calib.calib_key
                # Separate detectors?
                if self.par['ech_separate_2d']:
                    det_str = f'_{idet}'
                else:
                    det_str = ''
                # Global QA
                outfile_global = qa.set_qa_filename(
                    calib_key+det_str, 'arc_fit2d_global_qa',
                    out_dir=self.qa_path)
                arc.fit2darc_global_qa(fit2d, nspec, outfile=outfile_global)
                # Order QA
                outfile_orders = qa.set_qa_filename(
                    calib_key+det_str, 'arc_fit2d_orders_qa',
                    out_dir=self.qa_path)
                arc.fit2darc_orders_qa(fit2d, nspec, outfile=outfile_orders)

        return fit2ds, dets, save_order_dets


    # TODO: JFH this method is identical to the code in wavetilts.
    # SHould we make it a separate function?
    def extract_arcs(self, slitIDs=None):
        """
        Extract the arcs down each slit/order

        Wrapper to arc.get_censpec()

        Parameters
        ----------
        slitIDs : :obj:`list`, optional
            A list of the slit IDs to extract (if None, all slits will be
            extracted)

        Returns
        -------
        arccen : `numpy.ndarray`_
            arc spectrum for all slits, shape is (nspec, nslit):
        wvc_bpm : `numpy.ndarray`_
            boolean array containing a mask indicating which slits are good. True
            = masked (bad).
        """
        # Do it on the slits not masked in self.slitmask
        arccen, arccen_bpm, arc_maskslit = arc.get_censpec(
            self.slitcen, self.slitmask, self.msarc.image,
            gpm=self.gpm, slit_bpm=self.wvc_bpm,
            slitIDs=slitIDs)
        # Step
        self.steps.append(inspect.stack()[0][3])

        # Update the mask
        self.wvc_bpm |= arc_maskslit

        return arccen, self.wvc_bpm


    def update_wvmask(self):
        """
        (re)Generate the mask for wv_calib based on its contents
        This is the safest way to go...

        Args:
            nslit (int): Number of slits/orders

        Returns:
            `numpy.ndarray`_: self.wvc_bpm, boolean array -- True = masked, i.e. do not use

        """
        # Update mask based on wv_calib
        for kk, fit in enumerate(self.wv_calib.wv_fits):
            if fit is None or fit.pypeitfit is None:
                self.wvc_bpm[kk] = True


    def run(self, skip_QA=False, debug=False,
            prev_wvcalib=None):
        """
        Main method for wavelength calibration

        Code flow:
          1. Extract 1D arc spectra down the center of each unmasked slit/order
          2. Generate the 1D wavelength fits
          3. If echelle, perform 2D fit(s).
          4. Deal with masking
          5. Return a WaveCalib object

        Args:
            skip_QA (bool, optional): Skip QA?
            prev_wvcalib (WaveCalib, optional):
                Previous wavelength calibration object (from disk, typically)

        Returns:
            WaveCalib:  wavelength calibration object

        """
        ###############
        # Extract an arc down each slit
        self.arccen, self.wvc_bpm = self.extract_arcs()

        # Fill up the calibrations and generate QA
        self.wv_calib = self.build_wv_calib(
            self.arccen, self.par['method'],
            skip_QA=skip_QA,
            prev_wvcalib=prev_wvcalib)

        # Fit 2D?
        if self.par['echelle'] and self.par['ech_2dfit']:
            # Assess the fits
            rms = np.array([999. if wvfit.rms is None else wvfit.rms for wvfit in self.wv_calib.wv_fits])
            # get used FWHM
            fwhm = self.par['fwhm'] if self.measured_fwhms is None or self.par['fwhm_fromlines'] is False \
                else self.measured_fwhms.astype(float)
            # get rms threshold for all orders
            wave_rms_thresh = np.round(self.par['rms_thresh_frac_fwhm'] * fwhm, 3)
            bad_rms = rms > wave_rms_thresh
            if np.any(bad_rms):
                self.wvc_bpm[bad_rms] = True
                msgs.warn("Masking one or more bad orders (RMS)")
            # Fit
            fit2ds, dets, order_dets = self.echelle_2dfit(
                self.wv_calib, skip_QA = skip_QA, debug=debug)

            # Save
            self.wv_calib.wv_fit2d = np.array(fit2ds)
            # Save det_img?
            if self.par['ech_separate_2d']:
                self.wv_calib.det_img = self.msarc.det_img.copy()

            # Try a second attempt with 1D, if needed
            if np.any(bad_rms):
                bad_orders = self.slits.ech_order[np.where(bad_rms)[0]]
                any_fixed = self.redo_echelle_orders(bad_orders, dets, order_dets,
                                                     bad_orders_maxfrac=self.par['bad_orders_maxfrac'],
                                                     frac_rms_thresh=self.par['frac_rms_thresh'])

                # Do another full 2D?
                if any_fixed:
                    fit2ds, _, _ = self.echelle_2dfit(self.wv_calib, skip_QA = skip_QA, debug=debug)
                    # Save
                    self.wv_calib.wv_fit2d = np.array(fit2ds)

            # Check that we have at least one good 2D fit
            if not np.any([fit2d.success for fit2d in self.wv_calib.wv_fit2d]):
                msgs.error("No successful 2D Wavelength fits.  Cannot proceed.")

        # Deal with mask
        self.update_wvmask()

        # Any masked during this analysis?
        wv_masked = np.where(np.invert(self.wvc_bpm_init) & self.wvc_bpm)[0]
        if len(wv_masked) > 0:
            self.slits.mask[wv_masked] = self.slits.bitmask.turn_on(
                    self.slits.mask[wv_masked], 'BADWVCALIB')

        # Pack up
        sv_par = self.par.data.copy()
        j_par = jsonify(sv_par)
        self.wv_calib['strpar'] = json.dumps(j_par)#, sort_keys=True, indent=4, separators=(',', ': '))

        return self.wv_calib



    def show(self, item, slit=None):
        """
        Show one of the class internals

        Args:
            item (str):
              'spec' -- Show the fitted points and solution;  requires slit
              'fit' -- Show fit QA; requires slit
            slit (int, optional):

        Returns:

        """
        if item == 'spec':
            # spec
            spec = self.wv_calib[str(slit)]['spec']
            # tcent
            tcent = self.wv_calib[str(slit)]['tcent']
            yt = np.zeros_like(tcent)
            for jj,t in enumerate(tcent):
                it = int(np.round(t))
                yt[jj] = np.max(spec[it-1:it+1])
            # Plot
            plt.clf()
            ax=plt.gca()
            ax.plot(spec, drawstyle='steps-mid')
            ax.scatter(tcent, yt, color='red', marker='*')
            ax.set_xlabel('Pixel')
            ax.set_ylabel('Counts')
            plt.show()
        elif item == 'fit':
            autoid.arc_fit_qa(self.wv_calib[str(slit)])

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt

