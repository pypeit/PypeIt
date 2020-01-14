"""
Implements the flat-field class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
import numpy as np

from IPython import embed

from pypeit import msgs
from pypeit import ginga
from pypeit import masterframe

from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import flat
from pypeit.core import save
from pypeit.core import load
from pypeit.core import pixels
from pypeit.core import procimg


class FlatField(calibrationimage.CalibrationImage, masterframe.MasterFrame):
    """
    Builds pixel-level flat-field and the illumination flat-field.

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the instrument used to
            take the observations.  See usage by
            :class:`pypeit.processimages.ProcessImages` base class.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the flat frames.
        files (:obj:`list`, optional):
            The list of files to process.  Can be an empty list.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames
        msbias (`numpy.ndarray`_, :obj:`str`, optional):
            Either an image with the bias to be subtracted, or a string
            providing the method to use for bias correction.
        msbpm (`numpy.ndarray`_, optional):
            Bad pixel mask image
        flatpar (:class:`pypeit.par.pypeitpar.FlatFieldPar`, optional):
            User-level parameters for constructing the flat-field
            corrections.  If None, the default parameters are used.
        slits (:class:`pypeit.edgetrace.SlitTraceSet`):
            The current slit traces.
        tilts_dict (:obj:`dict`, optional):
            The current wavelength tilt traces; see
            :class:`pypeit.wavetilts.WaveTilts`.
        reuse_masters (:obj:`bool`, optional):
            Reuse already created master files from disk.

    Attributes:
        rawflatimg (PypeItImage):
        mspixelflat (ndarray):
            Normalized flat
        msillumflat (ndarray):
            Illumination flat
    """

    # Frame type is a class attribute
    frametype = 'pixelflat'
    master_type = 'Flat'

    @classmethod
    def from_master_file(cls, master_file, par=None):
        """
        Instantiate the class from a master file

        Args:
            master_file (str):
            par (:class:`pypeit.par.pypeitpar.PypeItPar`, optional):
                Full par set

        Returns:
            :class:`pypeit.flatfield.FlatField`:
                With the flat images loaded up

        """
        # Spectrograph
        spectrograph, extras = masterframe.items_from_master_file(master_file)
        head0 = extras[0]
        # Par
        if par is None:
            par = spectrograph.default_pypeit_par()
        # Master info
        master_dir = head0['MSTRDIR']
        master_key = head0['MSTRKEY']
        # Instantiate
        slf = cls(spectrograph, par['calibrations']['pixelflatframe'], master_dir=master_dir, master_key=master_key,
                  reuse_masters=True)
        # Load
        rawflatimg, slf.mspixelflat, slf.msillumflat = slf.load(ifile=master_file)
        # Convert rawflatimg to a PypeItImage
        slf.rawflatimg = pypeitimage.PypeItImage(rawflatimg)
        # Return
        return slf

    def __init__(self, spectrograph, par, files=None, det=1, master_key=None,
                 master_dir=None, reuse_masters=False, flatpar=None, msbias=None, msbpm=None,
                 slits=None, tilts_dict=None):

        # Image processing parameters
        if not isinstance(par, pypeitpar.FrameGroupPar):
            msgs.error('Must provide a FrameGroupPar instance as the parameters for FlatField.')
        self.par = par

        # Instantiate the base classes
        #   - Basic processing of the raw images
        calibrationimage.CalibrationImage.__init__(self, spectrograph, det, self.par['process'], files=files)
        #   - Construction and interface as a master frame
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

        # FieldFlattening parameters
        self.flatpar = pypeitpar.FlatFieldPar() if flatpar is None else flatpar

        # Input master frame data
        self.msbias = msbias
        self.msbpm = msbpm
        self.slits = slits
        self.tilts_dict = tilts_dict

        # Attributes unique to this Object
        self.rawflatimg = None      # Un-normalized pixel flat as a PypeItImage
        self.mspixelflat = None     # Normalized pixel flat
        self.msillumflat = None     # Slit illumination flat
        self.flat_model = None      # Model flat

        # Completed steps
        self.steps = []

        # Child-specific Internals
#        self.extrap_slit = None
#        self.msblaze = None
#        self.blazeext = None
#        self.slit_profiles = None
#        self.ntckx = None
#        self.ntcky = None

    @property
    def nslits(self):
        """
        Return the number of slits.  Pulled directly from :attr:`slits`, if it exists.
        """
        return 0 if self.slits is None else self.slits.nslits

    def build_pixflat(self, trim=True, force=False):
        """
        Generate the flat image.

        Args:
            trim (:obj:`bool`, optional):
                Trim the image down to just the data section.
            force (:obj:`bool`, optional):
                Force the flat to be reconstructed if it already exists

        Returns:
            pypeitimage.PypeItImage:  The image with the unnormalized pixel-flat data.
        """
        if self.rawflatimg is None or force:
            # Process steps
            self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
            ## JFH We never need untrimmed images. Why is this even an option?
            if trim:
                self.process_steps += ['trim']
            self.process_steps += ['apply_gain']
            self.process_steps += ['orient']
            # Turning this on leads to substantial edge-tracing problems when last tested
            #     JXP November 22, 2019
            #if self.par['cr_reject']:
            #    self.process_steps += ['crmask']
            self.steps.append(inspect.stack()[0][3])
            # Do it
            self.rawflatimg = super(FlatField, self).build_image(bias=self.msbias,
                                                                 bpm=self.msbpm,
                                                                 ignore_saturation=True)
        return self.rawflatimg

    # TODO Need to add functionality to use a different frame for the ilumination flat, e.g. a sky flat
    def run(self, debug=False, show=False, maskslits=None):
        """
        Generate normalized pixel and illumination flats

        Code flow::
           1.  Generate the pixelflat image (if necessary)
           2.  Prepare b-spline knot spacing
           3.  Loop on slits/orders
               a. Calculate the slit profile
               b. Normalize
               c. Save

        Args:
            debug (:obj:`bool`, optional):
                Run in debug mode.
            show (:obj:`bool`, optional):
                Show the results in the ginga viewer.
            maskslits (np.ndarray, optional):
               Array specifying whether a slit is good.
               True = bad

        Returns:
            `numpy.ndarray`_: Two arrays are returned, the normalized
            pixel flat data and the slit illumination correction data.
        """
        # Mask
        if maskslits is None:
            maskslits = np.zeros(self.nslits, dtype=bool)

        # Build the pixel flat (as needed)
        self.build_pixflat()

        # Prep tck (sets self.ntckx, self.ntcky)
        #self._prep_tck()

        # Setup
        self.mspixelflat = np.ones_like(self.rawflatimg.image)
        self.msillumflat = np.ones_like(self.rawflatimg.image)
        self.flat_model = np.zeros_like(self.rawflatimg.image)
#        self.slitmask = pixels.tslits2mask(self.tslits_dict)
        self.slitmask = self.slits.slit_img()

        final_tilts = np.zeros_like(self.rawflatimg.image)

        # If we are tweaking slits allocate the new aray to hold tweaked slit boundaries
        if self.flatpar['tweak_slits']:
            self.slits.init_tweaks()
#            self.tslits_dict['slit_left_tweak'] = np.zeros_like(self.tslits_dict['slit_left'])
#            self.tslits_dict['slit_righ_tweak'] = np.zeros_like(self.tslits_dict['slit_righ'])

        inmask = np.ones_like(self.rawflatimg.image, dtype=bool) if self.msbpm is None \
                    else np.invert(self.msbpm)
        nonlinear_counts = self.spectrograph.nonlinear_counts(det=self.det)

        # Loop on slits
        for slit in range(self.nslits):
            # Is this a good slit??
            if maskslits[slit]:
                msgs.info('Skipping bad slit: {}'.format(slit))
                continue

            msgs.info('Computing flat field image for slit: {:d}/{:d}'.format(slit,self.nslits-1))

            # Fit flats for a single slit
            this_tilts_dict = {'tilts':self.tilts_dict['tilts'],
                               'coeffs':self.tilts_dict['coeffs'][:,:,slit].copy(),
                               'slitcen':self.tilts_dict['slitcen'][:,slit].copy(),
                               'func2d':self.tilts_dict['func2d']}

            pixelflat, illumflat, flat_model, tilts_out, thismask_out, slit_left_out, \
                    slit_righ_out \
                            = flat.fit_flat(self.rawflatimg.image, this_tilts_dict, self.tslits_dict,
                                           slit, inmask=inmask, nonlinear_counts=nonlinear_counts,
                                           spec_samp_fine=self.flatpar['spec_samp_fine'],
                                           spec_samp_coarse=self.flatpar['spec_samp_coarse'],
                                           spat_samp=self.flatpar['spat_samp'],
                                           tweak_slits=self.flatpar['tweak_slits'],
                                           tweak_slits_thresh=self.flatpar['tweak_slits_thresh'],
                                           tweak_slits_maxfrac=self.flatpar['tweak_slits_maxfrac'],
                                           debug=debug)

            self.mspixelflat[thismask_out] = pixelflat[thismask_out]
            self.msillumflat[thismask_out] = illumflat[thismask_out]
            self.flat_model[thismask_out] = flat_model[thismask_out]

            # Did we tweak slit boundaries? If so, update the tslits_dict and the tilts_dict
            if self.flatpar['tweak_slits']:
                # TODO: Why do we need slit_left, slit_left_orig, and
                # slit_left_tweak? Shouldn't we only need two of these?
                
                self.slits.left[:,slit] = slit_left_out
                self.slits.left_tweak[:,slit] = slit_left_out

                self.slits.right[:,slit] = slit_righ_out
                self.slits.right_tweak[:,slit] = slit_righ_out
#                self.tslits_dict['slit_left'][:,slit] = slit_left_out
#                self.tslits_dict['slit_righ'][:,slit] = slit_righ_out
#                self.tslits_dict['slit_left_tweak'][:,slit] = slit_left_out
#                self.tslits_dict['slit_righ_tweak'][:,slit] = slit_righ_out
                final_tilts[thismask_out] = tilts_out[thismask_out]

        # If we tweaked the slits update the tilts_dict
        if self.flatpar['tweak_slits']:
            self.tilts_dict['tilts'] = final_tilts

        if show:
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show(slits=True, wcs_match = True)

        # If illumination flat fielding is turned off, set the illumflat to be None.
        if not self.flatpar['illumflatten']:
            msgs.warn('No illumination flat will be applied to your data (illumflatten=False).')
            self.msillumflat = None

        # Return
        return self.mspixelflat, self.msillumflat

    def show(self, slits=True, wcs_match=True):
        """
        Show all of the flat field products

        Args:
            slits (bool, optional):
            wcs_match (bool, optional):

        """
        viewer, ch = ginga.show_image(self.mspixelflat, chname='pixeflat', cuts=(0.9, 1.1),
                                      wcs_match=wcs_match, clear=True)
        viewer, ch = ginga.show_image(self.msillumflat, chname='illumflat', cuts=(0.9, 1.1),
                                      wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.rawflatimg, chname='flat', wcs_match=wcs_match)
        viewer, ch = ginga.show_image(self.flat_model, chname='flat_model', wcs_match=wcs_match)

        if slits and self.slits is not None:
            ginga.show_slits(viewer, ch, self.slits.left, self.slits.right, self.slits.id)

    def save(self, outfile=None, overwrite=True):
        """
        Save the flat-field master data to a FITS file

        Extensions are:
            RAWFLAT
            PIXELFLAT
            ILLUMFLAT
            MODEL

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Setup the items
        hdr = self.build_master_header(steps=self.steps, raw_files=self.file_list)
        data = [self.rawflatimg.image, self.mspixelflat, self.msillumflat, self.flat_model]
        extnames = ['RAWFLAT', 'PIXELFLAT', 'ILLUMFLAT', 'MODEL']

        # Save to a multi-extension FITS
        save.write_fits(hdr, data, _outfile, extnames=extnames)
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load the flat-field data from a save master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.

        Returns:
            tuple: Returns three `numpy.ndarray`_ objects with the raw
            flat-field image, the normalized pixel flat, and the
            illumination flat.
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return None, None, None
        # Load
        ext = ['RAWFLAT', 'PIXELFLAT', 'ILLUMFLAT', 'MODEL']
        self.rawflatimg, self.mspixelflat, self.msillumflat, self.flat_model, head0 = load.load_multiext_fits(master_file, ext)
        # Return
        return self.rawflatimg, self.mspixelflat, self.msillumflat

    # flat is self.rawflatimg.image
    # tilts_dict is self.tilts_dict
    # tslits_dict_in is self.slits
    # inmask = np.invert(self.msbpm)

            nonlinear_counts = self.spectrograph.nonlinear_counts(det=self.det)

                            = flat.fit_flat(
                                self.rawflatimg.image, this_tilts_dict, self.tslits_dict,
                                           slit, inmask=inmask, nonlinear_counts=nonlinear_counts,
                                           spec_samp_fine=self.flatpar['spec_samp_fine'],
                                           spec_samp_coarse=self.flatpar['spec_samp_coarse'],
                                           spat_samp=self.flatpar['spat_samp'],
                                           tweak_slits=self.flatpar['tweak_slits'],
                                           tweak_slits_thresh=self.flatpar['tweak_slits_thresh'],
                                           tweak_slits_maxfrac=self.flatpar['tweak_slits_maxfrac'],
                                           debug=debug)


            pixelflat, illumflat, flat_model, tilts_out, thismask_out, slit_left_out, \
                    slit_righ_out \




    def fit(self, debug=False, pad=5.):
        """
        Compute pixelflat and illumination flat from a flat field image.

        Parameters
        ----------
        flat :  float ndarray, shape (nspec, nspat)
            Flat field image in units of electrons.
        tilts_dict: dict
              Dictionary containing wavelength tilts image and other
              information indicating how wavelengths move across the slit


        tslits_dict: dict
              Dictionary with information on the slit boundaries
        slit: int
              Slit currently being considered
        inmask: boolean ndarray, shape (nspec, nspat), default inmask = None, optional
          Input mask for pixels not to be included in sky subtraction fits. True = Good (not masked), False = Bad (masked)

        spec_samp_fine: float, default = 1.2, optional
          bspline break point spacing in units of pixels for spectral fit to flat field blaze function.

        spec_samp_coarse: float, default = 50.0, optional
          bspline break point spacing in units of pixels for 2-d bspline-polynomial fit to flat field image residuals.
          This should be a large number unless you are trying to fit a sky flat with lots of features.

        spat_samp: float, default = 5.0, optional
          Spatial sampling for spatial slit illumination function. This is the width of the median filter in pixels used to
          determine the slit illumination function, and thus sets the minimum scale on which the illumination function will
          have features.

        trim_edg: tuple of floats  (left_edge, right_edge), default (3,3), optional
          indicates how many pixels to trim from left and right slit edges for creating the edgemask, which is used to mask
          the edges from the initial (fine) spectroscopic fit to the blaze function.

        pad: int, default = 5, optional
          Padding window used to create expanded slitmask images used for re-determining slit boundaries. Tilts are also
          computed using this expanded slitmask in cases the slit boundaries need to be moved outward.

        npoly: int, default = None, optional
          Order of polynomial for 2-d bspline-polynomial fit to flat field image residuals. The code determines the order of
          these polynomials to each slit automatically depending on the slit width, which is why the default is None.
          Do not attempt to set this paramter unless you know what you are doing.


        tweak_slits: bool, default = True, optional
          Slit edges will be tweaked such the left and right bounadaries intersect the location where the illumination
          function falls below tweak_slits_thresh (see below) of its maximum value near the center (moving out from the center)

        tweak_slits_thresh: float, default = 0.93, optional
          If tweak_slits is True, this sets the illumination function threshold used to tweak the slits

        tweak_slits_maxfrac: float, default = 0.10, optional
          Maximum fractinoal amount (of slit width) allowed for each trimming the left and right slit boundaries, i.e. the
          default is 10% which means slits would shrink by at most 20% (10% on each side)

        debug: bool, default = False, optional
          Show plots useful for debugging. This will block further execution of the code until the plot windows are closed.

        Returns
        -------
        pixeflat:   ndarray with same shape as flat
          Pixelflat gives pixel-to-pixel variations of detector response. Values are centered about unity.

        illumflat:  ndarray with same shape as flat
          Illumination flat gives variations of the slit illumination function across the spatial direction of the detect.
          Values are centered about unity. The slit illumination function is computed by dividing out the spectral response and
          collapsing out the spectral direction.

        flat_model:  ndarray with same shape as flat
          Full 2-d model image of the input flat image in units of electrons.  The pixelflat is defined to be flat/flat_model.

        tilts: ndarray with same shape as flat
          Tilts image fit for this slit evaluated using the new slit boundaries

        thismask_out: ndarray with same shape as flat, bool
           Boolean mask indicating which pixels are on the slit now with the new slit boundaries

        slit_left_out: ndarray with shape (nspec,)
           Tweaked left slit bounadries

        slit_righ_out: ndarray with shape (nspec,)
           Tweaked right slit bounadries

        Notes
        -----
    
        Revision History
            - 11-Mar-2005  First version written by Scott Burles.
            - 2005-2018    Improved by J. F. Hennawi and J. X. Prochaska
            - 3-Sep-2018 Ported to python by J. F. Hennawi and significantly improved
        """

        spec_samp_fine = self.flatpar['spec_samp_fine']
        spat_samp = self.flatpar['spat_samp']

        flat, tilts_dict, tslits_dict_in, slit, inmask=None, spec_samp_fine=1.2,
                 spec_samp_coarse=50.0, spat_samp=5.0, npoly=None, trim_edg=(3.0,3.0), pad=5.0,
                 tweak_slits=True, tweak_slits_thresh=0.93, tweak_slits_maxfrac=0.10,
                 nonlinear_counts=1e10, debug=False):

        nspec, nspat = self.rawflatimg.image.shape
        flat = self.rawflatimg.image
        gpm = np.ones_like(flat, dtype=bool) if self.msbpm is None else np.invert(self.msbpm)
        nonlinear_counts = self.spectrograph.nonlinear_counts(det=self.det)

        median_slit_width = np.median(self.right - self.left, axis=0)

        # Create both padded and unpadded slit ID images
        slitid_img = self.slits.slit_img()
        padded_slitid_img = self.slits.slit_img(pad=pad)

        for slit in range(self.slits.nslits):

            # Find the pixels on the slit and within a padded region around each slit
            onslit = slitid_img == slit

            # Check for saturation of the flat. If there are not enough
            # pixels do not attempt a fit, and continue to the next
            # slit.  TODO: set the threshold to a parameter?
            good_frac = np.sum(onslit & (flat < nonlinear_counts))/np.sum(onslit)
            if good_frac < 0.5:
                msgs.warn('Only {:4.2f}'.format(100*good_frac)
                          + '% of the pixels on this slit are not saturated.' + msgs.newline()
                          + 'Consider raising nonlinear_counts={:5.3f}'.format(nonlinear_counts) +
                          msgs.newline() + 'Not attempting to flat field slit {:d}'.format(slit))
                continue

            # Demand at least 10 pixels per row (on average) per degree
            # of the polynomial
            if npoly is None:
                # Approximate number of pixels sampling each spatial pixel
                # for this (original) slit.
                npercol = np.fmax(np.floor(np.sum(onslit)/nspec),1.0)
                npoly  = max(1, int(np.ceil(npercol/10.)))
                msgs.info('Fitting flat-field ')

            # TODO: Warn the user if npoly is provided but higher than
            # the nominal calculation if it is not provided?

            # TODO: Make npoly a parameter in fitpar?

            # Create an image with the spatial coordinates relative to the left edge
            coo_img, offslit_trimmed = self.slits.spatial_coordinate_image(slitids=slit, full=True)

            # Find pixels on the padded and trimmed slit coordinates
            onslit_padded = padded_slitid_img == slit
            onslit_trimmed = np.invert(offslit_trimmed)

            # ----------------------------------------------------------
            # Collapse the slit spatially and fit the spectral function

            # Create the tilts image for this slit
            # TODO: Is the copy needed?
            tilts = tracewave.fit2tilts(flat.shape, self.tilts_dict['coeffs'][:,:,slit].copy(),
                                        self.tilts_dict['func2d'])
            # Convert the tilt image to an image with the spectral pixel index
            spec_coo = tilts * (nspec-1)

            # Get the log of the flat flux for fitting
            flat_log = np.log(np.fmax(flat, 1.0))
            gpm_log = (flat > 1.0) & gpm
            ivar_log = gpm_log.astype(float)/0.5**2 # set errors to just be 0.5 in the log

            # Only include the trimmed set of pixels in the flat-field
            # fit along the spectral direction.
            fit_gpm = onslit_trimmed & gpm # & (flat < nonlinear_counts)
            nfit_spec = np.sum(fit_gpm)
            ntot_spec = np.sum(onslit)
            msgs.info('Spectral fit of flatfield for {0}/{1} '.format(nfit_spec, ntot_spec)
                      + ' pixels.')
            if nfit_spec/ntot_spec < 0.5:
                # TODO: Shouldn't this raise an exception or continue to the next slit instead?
                msgs.warn('Flatfield fit is to only {:4.2f}'.format(100*nfit_spec/ntot_spec)
                          + '% of the pixels on this slit.' + msgs.newline()
                          + '          Something appears to be wrong here')

            # TODO: This essentially requires that spec_img[fit_gpm]
            # always returns a vector, right?
            srt = np.argsort(spec_coo[fit_gpm])
            spec_coo_fit = spec_coo[fit_gpm][srt]
            flat_log_fit = flat_log[fit_gpm][srt]
            ivar_log_fit = ivar_log[fit_gpm][srt]
            gpm_log_fit = gpm_log[fit_gpm][srt]

            # rejection threshold for spectral fit in log(image)
            logrej = 0.5

            # Fit the spectral direction of the blaze.
            # TODO: Figure out how to deal with the fits going crazy at
            # the edges of the chip in spec direction
            spec_set_fine, outmask_spec, specfit, _, exit_status \
                    = utils.bspline_profile(spec_coo_fit, flat_log_fit, ivar_log_fit,
                                            np.ones_like(spec_fit), inmask=gpm_log_fit, nord=4,
                                            upper=logrej, lower=logrej,
                                            kwargs_bspline={'bkspace': spec_samp_fine},
                                            kwargs_reject={'groupbadpix': True, 'maxrej': 5})

            # TODO: Check exit status?

            # Debugging/checking spectral fit
            if debug:
                utils.bspline_profile_qa(spec_coo_fit, flat_log_fit, spec_set_fine, outmask_spec,
                                         specfit, xlabel='Spectral Pixel',
                                         ylabel='log(flat counts)',
                                         title='Spectral Fit for slit={:d}'.format(slit))

            # Evaluate and save
            spec_model = np.ones_like(flat)
            spec_model[onslit_padded], _ = np.exp(spec_set_fine.value(spec_coo[onslit_padded]))
            norm_spec = np.ones_like(flat)
            norm_spec[onslit_padded] = flat[onslit_padded] / np.fmax(spec_model[onslit_padded], 1.0)

            # ----------------------------------------------------------


            # ----------------------------------------------------------
            # To fit the spatial response, first normalize out the
            # spectral response, and then collapse the slit spectrally.
            # This assumes that the spatial response is, to first
            # order, independent of the spectral position.

            # Determine maximum counts in median filtered flat
            # spectrum. Only fit pixels > 0.1 of this maximum.
            specfit_interp = interp1d(spec_coo_fit, specfit, kind='linear', bounds_error=False,
                                      fill_value=-np.inf)
            spec_sm = utils.fast_running_median(np.exp(specfit_interp(np.arange(nspec))),
                                                np.fmax(np.ceil(0.10*nspec).astype(int),10))
            fit_spat = padded_slitmask & gpm & (spec_model > 1.0) \
                        & (spec_model > 0.1*np.amax(spec_sm)) & (norm_spec > 0.0) \
                        & (norm_spec < 1.7) #& (flat < nonlinear_counts)

            nfit_spat = np.sum(fit_spat)
            ntot_spat = np.sum(padded_slitmask)
            msgs.info('Spatial fit to flatfield for {:}'.format(nfit_spec) + ' pixels')
            if nfit_spat/ntot_spat < 0.5:
                # TODO: Shouldn't this raise an exception or continue to the next slit instead?
                msgs.warn('Spatial flatfield fit is to only {:4.2f}'.format(100*spat_frac) 
                          + '% of the pixels on this slit.' + msgs.newline()
                          + '              Something appears to be wrong here')

            # Sort the pixels by their spatial position
            srt = np.argsort(coo_img[fit_spat])
            spat_coo_fit = coo_img[fit_spat][srt]
            norm_spec_fit = norm_spec[fit_spat][srt]

            # Assume the density of samples in any given spatial
            # coordinate is roughly the same at all spatial positions.
            # Calculate the fraction of the slit width for the median
            # filter as set by the ``spat_samp`` parameter.
            med_width = int(np.ceil(nfit_spat * spat_samp / median_slit_width[:,slit]))
            # Median filter the data
            norm_flat_spat = utils.fast_running_median(norm_spec_fit, med_width)
            # Gaussian filter the data with a kernel that is 1/20th of
            # the median-filter width (or at least 0.5 pixels where
            # here a "pixel" is just the index of the data to fit)
            norm_flat_spat = scipy.ndimage.filters.gaussian_filter1d(norm_flat_spat,
                                                                     np.fmax(med_width/20.0, 0.5),
                                                                     mode='nearest')

            # Make sure that the normalized and filtered flat is finite!
            if np.any(np.invert(np.isfinite(norm_flat_spat))):
                msgs.error('Inifinities in slit illumination function computation!')

            # Determine the breakpoint spacing from the sampling of the
            # spatial coordinates. Use breakpoints at a spacing of a
            # 1/10th of a pixel, but do not allow a bsp smaller than
            # the typical sampling. Use the bspline class to determine
            # the breakpoints:
            spat_bspl = pydl.bspline(spat_coo_fit,nord=4,
                                     bkspace=np.fmax(1.0/median_slit_width[:,slit]/10.0,
                                                     1.2*np.median(np.diff(spat_coo_fit)))
            spat_set, outmask_spat, spatfit, _, exit_status \
                    = utils.bspline_profile(spat_coo_fit, norm_flat_spat,
                                            np.ones_like(norm_flat_spat),
                                            np.ones_like(norm_flat_spat), nord=4, upper=5.0,
                                            lower=5.0, fullbkpt=spat_bspl.breakpoints)

            # Evaluate the illumination profile
            illumflat = np.ones_like(flat)
            illumflat[onslit_padded], _ = spat_set.value(coo_img[onslit_padded])

            # Construct the spectrally and spatially normalized flat
            norm_spec_spat = np.ones_like(flat)
            norm_spec_spat[onslit_padded] = flat[onslit_padded] \
                                                / np.fmax(spec_model[onslit_padded], 1.0) \
                                                / np.fmax(illumflat[onslit_padded], 0.01)

        if tweak_slits:
            slit_left_out, slit_righ_out, tweak_dict \
                    = tweak_slit_edges(self.slits.left[:,slit], self.slits.right[:,slit], spat_coo_fit,
                                       norm_flat_spat, tweak_slits_thresh, tweak_slits_maxfrac)
            # Recreate all the quantities we need based on the tweaked slits
            tslits_dict_out = copy.deepcopy(tslits_dict_in)
            tslits_dict_out['slit_left'][:,slit] = slit_left_out
            tslits_dict_out['slit_righ'][:,slit] = slit_righ_out
            slitmask_out = pixels.tslits2mask(tslits_dict_out)
            thismask_out = (slitmask_out == slit)
            ximg_out, edgmask_out = pixels.ximg_and_edgemask(slit_left_out, slit_righ_out,
                                                             thismask_out, trim_edg=trim_edg)
            # Note that nothing changes with the tilts, since these were
            # already extrapolated across the whole image.
        else:
            # Generate the edgemask using the original slit boundaries and
            # thismask_in
            slit_left_out = np.copy(self.slits.left[:,slit])
            slit_righ_out = np.copy(self.slits.right[:,slit])
            thismask_out = onslit
            ximg_out = coo_img

        # Add an approximate pixel axis at the top
        if debug:
            plt.clf()
            ax = plt.gca()
            ax.plot(spat_coo_fit, norm_spec_fit, color='k', marker='o', markersize=0.4, mfc='k',
                    fillstyle='full',linestyle='None', label = 'all pixels')
            #ax.plot(spat_coo_fit[~imed], norm_spec_fit[~imed], color='darkred', marker='+',markersize=4.0,
            #        mfc='red', fillstyle='full', linestyle='None', label = 'masked')
            #ax.plot(spat_coo_fit[imed], normfit[imed], color='orange', label = 'median spatial profile')
            ax.plot(spat_coo_fit, spatfit, color='cornflowerblue',
                    label='final slit illumination function')
            ymin = np.fmax(0.8 * spatfit.min(), 0.5)
            ymax = 1.2*spatfit.max()
            ax.set_ylim((np.fmax(0.8 * spatfit.min(), 0.5), 1.2 * spatfit.max()))
            ax.set_xlim(spat_coo_fit.min(), spat_coo_fit.max())
            plt.vlines(0.0, ymin, ymax, color='lightgreen', linestyle=':', linewidth=2.0,
                       label='original left edge', zorder=8)
            plt.vlines(1.0, ymin, ymax, color='red', linestyle=':', linewidth=2.0,
                       label='original right edge', zorder=9)
            if tweak_slits:
                if tweak_dict['tweak_left']:
                    label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) \
                                + ' % of max of left illumprofile'
                    plt.hlines(tweak_slits_thresh*tweak_dict['norm_max_left'], spat_coo_fit.min(), 0.5,
                               color='lightgreen', linewidth=3.0, label=label, zorder=10)
                    plt.vlines(tweak_dict['xleft'], ymin, ymax, color='lightgreen', linestyle='--',
                               linewidth=3.0, label='tweaked left edge', zorder=11)
                if tweak_dict['tweak_righ']:
                    label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) \
                                + ' % of max of right illumprofile'
                    plt.hlines(tweak_slits_thresh * tweak_dict['norm_max_righ'], 0.5, spat_coo_fit.max(),
                               color='red', linewidth=3.0, label=label, zorder=10)
                    plt.vlines(tweak_dict['xrigh'], ymin, ymax, color='red', linestyle='--',
                               linewidth=3.0, label='tweaked right edge', zorder=20)
            plt.legend()
            plt.xlabel('Normalized Slit Position')
            plt.ylabel('Normflat Spatial Profile')
            plt.title('Illumination Function Fit for slit={:d}'.format(slit))
            plt.show()

        msgs.info('Performing illumination + scattembedered light flat field fit')

        # Flat field pixels for fitting spectral direction
        isrt_spec = np.argsort(spec_coo[thismask_out])
        pix_twod = spec_coo[thismask_out][isrt_spec]
        ximg_twod = ximg_out[thismask_out][isrt_spec]
        norm_twod = norm_spec_spat[thismask_out][isrt_spec]

        fitmask = inmask[thismask_out][isrt_spec] & (np.abs(norm_twod - 1.0) < 0.30)
        # Here we ignore the formal photon counting errors and simply
        # assume that a typical error per pixel. This guess is somewhat
        # aribtrary. We then set the rejection threshold with sigrej_illum
        var_value = 0.01
        norm_twod_ivar = fitmask.astype(float)/(var_value**2)
        sigrej_illum = 4.0

        poly_basis = pydl.fpoly(2.0*ximg_twod - 1.0, npoly).T

        # Perform the full 2d fit now
        twod_set, outmask_twod, twodfit, _ , exit_status \
                = utils.bspline_profile(pix_twod, norm_twod, norm_twod_ivar, poly_basis,
                                        inmask=fitmask, nord=4, upper=sigrej_illum, lower=sigrej_illum,
                                        kwargs_bspline={'bkspace': spec_samp_coarse},
                                        kwargs_reject={'groupbadpix': True, 'maxrej': 10})

        if debug:
            resid = (norm_twod - twodfit)
            badpix = np.invert(outmask_twod) & fitmask
            goodpix = outmask_twod & fitmask
            plt.clf()
            ax = plt.gca()
            ax.plot(pix_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k',
                    fillstyle='full', linestyle='None', label='good points')
            ax.plot(pix_twod[badpix], resid[badpix], color='red', marker='+', markersize=0.5,
                    mfc='red', fillstyle='full', linestyle='None', label='masked')
            plt.hlines(sigrej_illum*var_value, pix_twod.min(), pix_twod.max(), color='lawngreen',
                       linestyle='--', label='rejection thresholds', zorder=10, linewidth=2.0)
            plt.hlines(-sigrej_illum*var_value,pix_twod.min(),pix_twod.max(), color='lawngreen',
                       linestyle='--', zorder=10, linewidth=2.0)
            ax.set_ylim((-0.05,0.05))
            ax.set_xlim((pix_twod.min(), pix_twod.max()))
            plt.legend()
            plt.xlabel('Spectral Pixel')
            plt.ylabel('Residuals from pixelflat 2-d fit')
            plt.title('Spectral Residuals for slit={:d}'.format(slit))
            plt.show()

            plt.clf()
            ax = plt.gca()
            ax.plot(ximg_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k',
                    fillstyle='full', linestyle='None', label='good points')
            ax.plot(ximg_twod[badpix], resid[badpix], color='red', marker='+', markersize=0.5,
                    mfc='red', fillstyle='full', linestyle='None', label='masked')
            plt.hlines(sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen',
                       linestyle='--', label='rejection thresholds', zorder=10, linewidth=2.0)
            plt.hlines(-sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen',
                       linestyle='--', zorder=10,linewidth=2.0)
            ax.set_ylim((-0.05, 0.05))
            ax.set_xlim(-0.02, 1.02)
            plt.legend()
            plt.xlabel('Normalized Slit Position')
            plt.ylabel('Residuals from pixelflat 2-d fit')
            plt.title('Spatial Residuals for slit={:d}'.format(slit))
            plt.show()

        # Evaluate and save
        twod_model = np.ones_like(flat)
        twod_this = np.zeros_like(twodfit)
        twod_this[isrt_spec] = twodfit
        twod_model[thismask_out] = twod_this

        # Compute all the final output images output
        pixelflat = np.ones_like(flat)
        flat_model = np.ones_like(flat)
        flat_model[thismask_out] = twod_model[thismask_out]*np.fmax(illumflat[thismask_out],0.05) \
                                        * np.fmax(spec_model[thismask_out],1.0)
        pixelflat[thismask_out] = flat[thismask_out]/flat_model[thismask_out]

        # TODO: Add some code here to treat the edges and places where fits
        # go bad?
        # Set the pixelflat to 1.0 wherever the flat was nonlinear
        pixelflat[flat >= nonlinear_counts] = 1.0

        # Do not apply pixelflat field corrections that are greater than
        # 100% to avoid creating edge effects, etc.
        # TODO: Should we do the same for the illumflat??
        #pixelflat = np.fmax(np.fmin(pixelflat, 2.0), 0.5)
        pixelflat = np.clip(pixelflat, 0.5, 2.0)

        return pixelflat, illumflat, flat_model, tilts, thismask_out, slit_left_out, slit_righ_out


