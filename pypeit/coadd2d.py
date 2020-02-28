"""
Module for performing two-dimensional coaddition of spectra.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import copy

from IPython import embed

import numpy as np
import scipy
from matplotlib import pyplot as plt

from astropy.io import fits

from pypeit import msgs
from pypeit.masterframe import MasterFrame
from pypeit.waveimage import WaveImage
from pypeit.wavetilts import WaveTilts
from pypeit import specobjs
from pypeit import edgetrace
from pypeit import reduce
from pypeit.core import extract
from pypeit.core import coadd, pixels
from pypeit.core import parse
from pypeit.images import scienceimage
from pypeit.spectrographs import util
from pypeit import calibrations


class CoAdd2D(object):

    """
    Main routine to run the extraction for 2d coadds.

    Algorithm steps are as follows:
        - Fill this in.

    This performs 2d coadd specific tasks, and then also performs some
    of the tasks analogous to the pypeit.extract_one method. Docs coming
    soon....
    """


    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec2dfiles, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                     ir_redux=False, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        """
        Instantiate the CoAdd2d subclass appropriate for the provided spectrograph.

        The class must be subclassed this class CoAdd2d.

        Parameters
        ----------
            See the documenation for the __init__ function below

        Returns
        -------
            :class:`CoAdd2d`: One of the subclasses with :class:`CoAdd2d` as its
            base.
        """

        return next(c for c in cls.__subclasses__() if c.__name__ == (spectrograph.pypeline + 'CoAdd2D'))(
            spec2dfiles, spectrograph, par, det=det, offsets=offsets, weights=weights, sn_smooth_npix=sn_smooth_npix,
            ir_redux=ir_redux, show=show, show_peaks=show_peaks, debug_offsets=debug_offsets, debug=debug, **kwargs_wave)

    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        """

        Parameters
        ----------
            spec2d_files (list):
               List of spec2d files
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
                The instrument used to collect the data to be reduced.
            par (:class:`pypeit.par.parset.ParSet`):
                Parset object
            det (int): optional
                Detector to reduce
            offsets (ndarray): default=None
                Spatial offsets to be applied to each image before coadding. For the default mode of None, images
                are registered automatically using the trace of the brightest object. Input offsets are not yet supported.
            weights (str, list or ndarray):
                Mode for the weights used to coadd images. Options are 'auto' (default), 'uniform', or list/array of
                weights with shape = (nexp,) can be input and will be applied to the image. Note 'auto' is not allowed
                if offsets are input, and if set this will cause an exception.
            sn_smooth_npix (int): optional, default=None
                Number of pixels to median filter by when computing S/N used to decide how to scale and weight spectra. If
                set to None, the code will simply take 10% of the image size in the spectral direction.
                TODO: for truncated echelle orders we should be doing something more intelligent.
            ir_redux (bool): optional, default=False
                Is this an near-IR reduction, True=yes. This parameter is passed to pypeit.reduce for determining the
                reduction steps.
            show (bool): optional, default=False
                Show results to ginga
            show_peaks (bool): optional, default=False
                Show the QA for object finding algorithm peak finding to the screen.
            debug_offset (bool): optional, default=False
                Show QA for debugging the automatic determination of offsets to the screen.
            debug (bool): optional, default=False
                Show QA for debugging.
            **kwargs_wave
                Keyword arguments pass to `pypeit.core.coadd.get_wvae_grid` which determine how the wavelength grid
                is created for the 2d coadding.
        """

        ## Use Cases:
        #  1) offsets is None -- auto compute offsets from brightest object, so then default to auto_weights=True
        #  2) offsets not None, weights = None (uniform weighting) or weights is not None (input weights)
        #  3) offsets not None, auto_weights=True (Do not support)
        if offsets is not None and 'auto' in weights:
            msgs.error("Automatic weights cannot be computed for input offsets. "
                       "Set weights='uniform' or input an array of weights with shape (nexp,)")
        self.spec2d_files = spec2d_files
        self.spectrograph = spectrograph
        self.par = par
        self.det = det
        self.offsets = offsets
        self.weights = weights
        self.ir_redux = ir_redux
        self.show = show
        self.show_peaks = show_peaks
        self.debug_offsets = debug_offsets
        self.debug = debug
        self.stack_dict = None
        self.psuedo_dict = None

        self.objid_bri = None
        self.slitid_bri  = None
        self.snr_bar_bri = None


        # Load the stack_dict
        self.stack_dict = self.load_coadd2d_stacks(self.spec2d_files)
        self.pypeline = self.spectrograph.pypeline

        # Check that there are the same number of slits on every exposure
        nslits_list = []
        for tslits_dict in self.stack_dict['tslits_dict_list']:
            nspec, nslits_now = tslits_dict['slit_left'].shape
            nslits_list.append(nslits_now)
        if not len(set(nslits_list))==1:
            msgs.error('Not all of your exposures have the same number of slits. Check your inputs')
        self.nslits = nslits_list[0]
        self.nexp = len(self.stack_dict['specobjs_list'])
        self.nspec = nspec
        self.binning = np.array([self.stack_dict['tslits_dict_list'][0]['binspectral'],
                                 self.stack_dict['tslits_dict_list'][0]['binspatial']])

        # If smoothing is not input, smooth by 10% of the spectral dimension
        self.sn_smooth_npix = sn_smooth_npix if sn_smooth_npix is not None else 0.1*self.nspec

    def optimal_weights(self, slitorderid, objid, const_weights=False):
        """
        Determine optimal weights for 2d coadds. This script grabs the information from SpecObjs list for the
        object with specified slitid and objid and passes to coadd.sn_weights to determine the optimal weights for
        each exposure.

        Args:
            slitorderid (int):
               The slit or order id that has the brightest object whose S/N will be used to determine the weight for each frame.
            objid (np.ndarray):
               Array of object indices with  shape = (nexp,) of the brightest object whose S/N will be used to determine the weight for each frame.
            const_weights (bool):
               Use constant weights for coadding the exposures. Default=False

        Returns:
            rms_sn : ndarray, shape = (len(specobjs_list),)
                Root mean square S/N value for each input spectra
            weights : ndarray, shape (len(specobjs_list),)
                Weights to be applied to the spectra. These are signal-to-noise squared weights.
        """

        nexp = len(self.stack_dict['specobjs_list'])
        nspec = self.stack_dict['specobjs_list'][0][0].TRACE_SPAT.shape[0]
        # Grab the traces, flux, wavelength and noise for this slit and objid.
        flux_stack = np.zeros((nspec, nexp), dtype=float)
        ivar_stack = np.zeros((nspec, nexp), dtype=float)
        wave_stack = np.zeros((nspec, nexp), dtype=float)
        mask_stack = np.zeros((nspec, nexp), dtype=bool)

        for iexp, sobjs in enumerate(self.stack_dict['specobjs_list']):
            ithis = sobjs.slitorder_objid_indices(slitorderid, objid[iexp])
            flux_stack[:, iexp] = sobjs[ithis].OPT_COUNTS
            ivar_stack[:, iexp] = sobjs[ithis].OPT_COUNTS_IVAR
            wave_stack[:, iexp] = sobjs[ithis].OPT_WAVE
            mask_stack[:, iexp] = sobjs[ithis].OPT_MASK

        # TODO For now just use the zero as the reference for the wavelengths? Perhaps we should be rebinning the data though?
        rms_sn, weights = coadd.sn_weights(wave_stack, flux_stack, ivar_stack, mask_stack, self.sn_smooth_npix,
                                           const_weights=const_weights)
        return rms_sn, weights.T


    def coadd(self, only_slits=None):

        only_slits = [only_slits] if (only_slits is not None and
                                      isinstance(only_slits, (int, np.int, np.int64, np.int32))) else only_slits
        good_slits = np.arange(self.nslits) if only_slits is None else only_slits

        coadd_list = []
        for islit in good_slits:
            msgs.info('Performing 2d coadd for slit: {:d}/{:d}'.format(islit, self.nslits - 1))
            ref_trace_stack = self.reference_trace_stack(islit, offsets=self.offsets, objid=self.objid_bri)
            thismask_stack = self.stack_dict['slitmask_stack'] == islit
            # TODO Can we get rid of this one line simply making the weights returned by parse_weights an
            # (nslit, nexp) array?
            # This one line deals with the different weighting strategies between MultiSlit echelle. Otherwise, we
            # would need to copy this method twice in the subclasses
            if 'auto_echelle' in self.use_weights:
                rms_sn, weights = self.optimal_weights(islit, self.objid_bri)
            else:
                weights = self.use_weights
            # Perform the 2d coadd
            coadd_dict = coadd.compute_coadd2d(ref_trace_stack, self.stack_dict['sciimg_stack'],
                                           self.stack_dict['sciivar_stack'],
                                           self.stack_dict['skymodel_stack'], self.stack_dict['mask_stack'] == 0,
                                           self.stack_dict['tilts_stack'], thismask_stack,
                                           self.stack_dict['waveimg_stack'],
                                           self.wave_grid, weights=weights)
            coadd_list.append(coadd_dict)

        return coadd_list


    def create_psuedo_image(self, coadd_list):
        """ THIS UNDOCUMENTED CODE PROBABLY SHOULD GENERATE AND RETURN
        STANDARD PYPEIT OBJCTS INSTEAD OF SOME UNDEFINED DICT"""



        nspec_vec = np.zeros(self.nslits,dtype=int)
        nspat_vec = np.zeros(self.nslits,dtype=int)
        for islit, cdict in enumerate(coadd_list):
            nspec_vec[islit]=cdict['nspec']
            nspat_vec[islit]=cdict['nspat']

        # Determine the size of the psuedo image
        nspat_pad = 10
        nspec_psuedo = nspec_vec.max()
        nspat_psuedo = np.sum(nspat_vec) + (self.nslits + 1)*nspat_pad
        spec_vec_psuedo = np.arange(nspec_psuedo)
        shape_psuedo = (nspec_psuedo, nspat_psuedo)
        imgminsky_psuedo = np.zeros(shape_psuedo)
        sciivar_psuedo = np.zeros(shape_psuedo)
        waveimg_psuedo = np.zeros(shape_psuedo)
        tilts_psuedo = np.zeros(shape_psuedo)
        spat_img_psuedo = np.zeros(shape_psuedo)
        nused_psuedo = np.zeros(shape_psuedo, dtype=int)
        inmask_psuedo = np.zeros(shape_psuedo, dtype=bool)
        wave_mid = np.zeros((nspec_psuedo, self.nslits))
        wave_mask = np.zeros((nspec_psuedo, self.nslits),dtype=bool)
        wave_min = np.zeros((nspec_psuedo, self.nslits))
        wave_max = np.zeros((nspec_psuedo, self.nslits))
        dspat_mid = np.zeros((nspat_psuedo, self.nslits))

        spat_left = nspat_pad
        slit_left = np.zeros((nspec_psuedo, self.nslits))
        slit_righ = np.zeros((nspec_psuedo, self.nslits))
        spec_min1 = np.zeros(self.nslits)
        spec_max1 = np.zeros(self.nslits)

        nspec_grid = self.wave_grid_mid.size
        for islit, coadd_dict in enumerate(coadd_list):
            spat_righ = spat_left + nspat_vec[islit]
            ispec = slice(0,nspec_vec[islit])
            ispat = slice(spat_left,spat_righ)
            imgminsky_psuedo[ispec, ispat] = coadd_dict['imgminsky']
            sciivar_psuedo[ispec, ispat] = coadd_dict['sciivar']
            waveimg_psuedo[ispec, ispat] = coadd_dict['waveimg']
            tilts_psuedo[ispec, ispat] = coadd_dict['tilts']
            # spat_img_psuedo is the sub-pixel image position on the rebinned psuedo image
            inmask_psuedo[ispec, ispat] = coadd_dict['outmask']
            image_temp = (coadd_dict['dspat'] -  coadd_dict['dspat_mid'][0] + spat_left)*coadd_dict['outmask']
            spat_img_psuedo[ispec, ispat] = image_temp
            nused_psuedo[ispec, ispat] = coadd_dict['nused']
            wave_min[ispec, islit] = coadd_dict['wave_min']
            wave_max[ispec, islit] = coadd_dict['wave_max']
            wave_mid[ispec, islit] = coadd_dict['wave_mid']
            wave_mask[ispec, islit] = True
            # Fill in the rest of the wave_mid with the corresponding points in the wave_grid
            #wave_this = wave_mid[wave_mask[:,islit], islit]
            #ind_upper = np.argmin(np.abs(self.wave_grid_mid - wave_this.max())) + 1
            #if nspec_vec[islit] != nspec_psuedo:
            #    wave_mid[nspec_vec[islit]:, islit] = self.wave_grid_mid[ind_upper:ind_upper + (nspec_psuedo-nspec_vec[islit])]


            dspat_mid[ispat, islit] = coadd_dict['dspat_mid']
            slit_left[:,islit] = np.full(nspec_psuedo, spat_left)
            slit_righ[:,islit] = np.full(nspec_psuedo, spat_righ)
            spec_max1[islit] = nspec_vec[islit]-1
            spat_left = spat_righ + nspat_pad

        slitcen = (slit_left + slit_righ)/2.0
        tslits_dict_psuedo = dict(slit_left=slit_left, slit_righ=slit_righ, slitcen=slitcen,
                                  nspec=nspec_psuedo, nspat=nspat_psuedo, pad=0,
                                  nslits = self.nslits, binspectral=1, binspatial=1, spectrograph=self.spectrograph.spectrograph,
                                  spec_min=spec_min1, spec_max=spec_max1,
                                  maskslits=np.zeros(slit_left.shape[1], dtype=np.bool))

        slitmask_psuedo = pixels.tslits2mask(tslits_dict_psuedo)
        # This is a kludge to deal with cases where bad wavelengths result in large regions where the slit is poorly sampled,
        # which wreaks havoc on the local sky-subtraction
        min_slit_frac = 0.70
        spec_min = np.zeros(self.nslits)
        spec_max = np.zeros(self.nslits)
        for islit in range(self.nslits):
            slit_width = np.sum(inmask_psuedo*(slitmask_psuedo == islit),axis=1)
            slit_width_img = np.outer(slit_width, np.ones(nspat_psuedo))
            med_slit_width = np.median(slit_width_img[slitmask_psuedo == islit])
            nspec_eff = np.sum(slit_width > min_slit_frac*med_slit_width)
            nsmooth = int(np.fmax(np.ceil(nspec_eff*0.02),10))
            slit_width_sm = scipy.ndimage.filters.median_filter(slit_width, size=nsmooth, mode='reflect')
            igood = (slit_width_sm > min_slit_frac*med_slit_width)
            spec_min[islit] = spec_vec_psuedo[igood].min()
            spec_max[islit] = spec_vec_psuedo[igood].max()
            bad_pix = (slit_width_img < min_slit_frac*med_slit_width) & (slitmask_psuedo == islit)
            inmask_psuedo[bad_pix] = False

        # Update with tslits_dict_psuedo
        tslits_dict_psuedo['spec_min'] = spec_min
        tslits_dict_psuedo['spec_max'] = spec_max

        psuedo_dict = dict(nspec=nspec_psuedo, nspat=nspat_psuedo, imgminsky=imgminsky_psuedo, sciivar=sciivar_psuedo,
                           inmask=inmask_psuedo, tilts=tilts_psuedo,
                           waveimg=waveimg_psuedo, spat_img = spat_img_psuedo,
                           tslits_dict=tslits_dict_psuedo,
                           wave_mask=wave_mask, wave_mid=wave_mid, wave_min=wave_min, wave_max=wave_max)

        return psuedo_dict

    def reduce(self, psuedo_dict, show=None, show_peaks=None):

        show = self.show if show is None else show
        show_peaks = self.show_peaks if show_peaks is None else show_peaks

        # Generate a ScienceImage
        sciImage = scienceimage.ScienceImage(self.spectrograph, self.det,
                                                      self.par['scienceframe']['process'],
                                                      psuedo_dict['imgminsky'],
                                                      psuedo_dict['sciivar'],
                                                      np.zeros_like(psuedo_dict['inmask']),  # Dummy bpm
                                                      rn2img=np.zeros_like(psuedo_dict['inmask']),  # Dummy rn2img
                                                      crmask=np.invert(psuedo_dict['inmask']))
        slitmask_psuedo = pixels.tslits2mask(psuedo_dict['tslits_dict'])
        sciImage.build_mask(slitmask=slitmask_psuedo)

        # Make changes to parset specific to 2d coadds
        parcopy = copy.deepcopy(self.par)
        parcopy['scienceimage']['findobj']['trace_npoly'] = 3        # Low order traces since we are rectified
        #parcopy['scienceimage']['find_extrap_npoly'] = 1  # Use low order for trace extrapolation
        # Instantiate Calibrations class
        caliBrate = calibrations.MultiSlitCalibrations(None, parcopy['calibrations'], self.spectrograph)
        caliBrate.tslits_dict = psuedo_dict['tslits_dict']
        caliBrate.tilts_dict = dict(tilts=psuedo_dict['tilts'])
        caliBrate.mswave = psuedo_dict['waveimg']
        #
        # redux = reduce.instantiate_me(sciImage, self.spectrograph, psuedo_dict['tslits_dict'], parcopy, psuedo_dict['tilts'],
        redux=reduce.instantiate_me(sciImage, self.spectrograph, parcopy, caliBrate,
                                    ir_redux=self.ir_redux, objtype='science_coadd2d',
                                    det=self.det, binning=self.binning, show=show)

        if show:
            redux.show('image', image=psuedo_dict['imgminsky']*(sciImage.mask == 0), chname = 'imgminsky', slits=True, clear=True)
        # Object finding
        sobjs_obj, nobj, skymask_init = redux.find_objects(sciImage.image, show_peaks=show_peaks)
        # Local sky-subtraction
        global_sky_psuedo = np.zeros_like(psuedo_dict['imgminsky']) # No global sky for co-adds since we go straight to local
        skymodel_psuedo, objmodel_psuedo, ivarmodel_psuedo, outmask_psuedo, sobjs = redux.local_skysub_extract(
            caliBrate.mswave, global_sky_psuedo, sobjs_obj, spat_pix=psuedo_dict['spat_img'], model_noise=False,
            show_profile=show, show=show)

        if self.ir_redux:
            sobjs.purge_neg()

        # TODO: Removed this, but I'm not sure that's what you want...
#        # Add the information about the fixed wavelength grid to the sobjs
#        for spec in sobjs:
#            idx = spec.slit_orderindx
#            # Fill
#            spec.BOX_WAVE_GRID_MASK, spec.OPT_WAVE_GRID_MASK = [psuedo_dict['wave_mask'][:,idx]]*2
#            spec.BOX_WAVE_GRID, spec.OPT_WAVE_GRID = [psuedo_dict['wave_mid'][:,idx]]*2
#            spec.BOX_WAVE_GRID_MIN, spec.OPT_WAVE_GRID_MIN = [psuedo_dict['wave_min'][:,idx]]*2
#            spec.BOX_WAVE_GRID_MAX, spec.OPT_WAVE_GRID_MAX = [psuedo_dict['wave_max'][:,idx]]*2

        # Add the rest to the psuedo_dict
        psuedo_dict['skymodel'] = skymodel_psuedo
        psuedo_dict['objmodel'] = objmodel_psuedo
        psuedo_dict['ivarmodel'] = ivarmodel_psuedo
        psuedo_dict['outmask'] = outmask_psuedo
        psuedo_dict['sobjs'] = sobjs
        self.psuedo_dict=psuedo_dict

        return psuedo_dict['imgminsky'], psuedo_dict['sciivar'], skymodel_psuedo, objmodel_psuedo, ivarmodel_psuedo, outmask_psuedo, sobjs


    def save_masters(self, master_dir):

        # Write out the psuedo master files to disk
        master_key_dict = self.stack_dict['master_key_dict']

        # TODO: These saving operations are a temporary kludge
        waveImage = WaveImage(None, None, None, self.spectrograph,  # spectrograph is needed for header
                              None, None, master_key=master_key_dict['arc'],
                              master_dir=master_dir)
        waveImage.save(image=self.psuedo_dict['waveimg'])

        edges = edgetrace.EdgeTraceSet.from_tslits_dict(self.psuedo_dict['tslits_dict'],
                                                        master_key_dict['trace'], master_dir)
        edges.save()

    def snr_report(self, snr_bar, slitid=None):

        # Print out a report on the SNR
        msg_string = msgs.newline() + '-------------------------------------'
        msg_string += msgs.newline() + '  Summary for highest S/N object'
        if slitid is not None:
            msg_string += msgs.newline() + '      found on slitid = {:d}            '.format(slitid)
        msg_string += msgs.newline() + '-------------------------------------'
        msg_string += msgs.newline() + '           exp#        S/N'
        for iexp, snr in enumerate(snr_bar):
            msg_string += msgs.newline() + '            {:d}         {:5.2f}'.format(iexp, snr)

        msg_string += msgs.newline() + '-------------------------------------'
        msgs.info(msg_string)

    def get_good_slits(self, only_slits):

        only_slits = [only_slits] if (only_slits is not None and
                                        isinstance(only_slits, (int, np.int, np.int64, np.int32))) else only_slits
        good_slits = np.arange(self.nslits) if only_slits is None else only_slits
        return good_slits

    def offset_slit_cen(self, slitid, offsets):

        nexp = len(offsets)
        tslits_dict_list = self.stack_dict['tslits_dict_list']
        nspec, nslits = tslits_dict_list[0]['slit_left'].shape
        ref_trace_stack = np.zeros((nspec, nexp))
        for iexp, tslits_dict in enumerate(tslits_dict_list):
            ref_trace_stack[:, iexp] = (tslits_dict['slit_left'][:, slitid] +
                                        tslits_dict['slit_righ'][:, slitid])/2.0 - offsets[iexp]
        return ref_trace_stack

    def get_wave_grid(self, **kwargs_wave):
        """
        Routine to create a wavelength grid for 2d coadds using all of the wavelengths of the extracted objects. Calls
        coadd1d.get_wave_grid.

        Args:
            **kwargs_wave (dict):
                Optional argumments for coadd1d.get_wve_grid function

        Returns:
            tuple: Returns the following:
                - wave_grid (np.ndarray): New wavelength grid, not
                  masked
                - wave_grid_mid (np.ndarray): New wavelength grid
                  evaluated at the centers of the wavelength bins, that
                  is this grid is simply offset from wave_grid by
                  dsamp/2.0, in either linear space or log10 depending
                  on whether linear or (log10 or velocity) was
                  requested.  For iref or concatenate the linear
                  wavelength sampling will be calculated.
                - dsamp (float): The pixel sampling for wavelength grid
                  created.
        """

        nobjs_tot = np.array([len(spec) for spec in self.stack_dict['specobjs_list']]).sum()
        waves = np.zeros((self.nspec, nobjs_tot))
        masks = np.zeros_like(waves, dtype=bool)
        indx = 0
        for spec_this in self.stack_dict['specobjs_list']:
            for spec in spec_this:
                waves[:, indx] = spec.OPT_WAVE
                masks[:, indx] = spec.OPT_MASK
                indx += 1

        wave_grid, wave_grid_mid, dsamp = coadd.get_wave_grid(waves, masks=masks, **kwargs_wave)

        return wave_grid, wave_grid_mid, dsamp

    def load_coadd2d_stacks(self, spec2d_files):
        """
        Routine to read in required images for 2d coadds given a list of spec2d files.

        Args:
            spec2d_files: list
               List of spec2d filenames
            det: int
               detector in question

        Returns:
            dict: Dictionary containing all the images and keys required
            for perfomring 2d coadds.
        """

        # Get the detector string
        sdet = parse.get_dnum(self.det, prefix=False)

        # Get the master dir

        redux_path = os.getcwd()

        # Grab the files
        head2d_list = []
        tracefiles = []
        waveimgfiles = []
        tiltfiles = []
        spec1d_files = []
        for f in spec2d_files:
            head = fits.getheader(f)
            if os.path.exists(head['PYPMFDIR']):
                master_path = head['PYPMFDIR']
            else:
                master_dir = os.path.basename(head['PYPMFDIR'])
                master_path = os.path.join(os.path.split(os.path.split(f)[0])[0], master_dir)

            trace_key = '{0}_{1:02d}'.format(head['TRACMKEY'], self.det)
            wave_key = '{0}_{1:02d}'.format(head['ARCMKEY'], self.det)

            head2d_list.append(head)
            spec1d_files.append(f.replace('spec2d', 'spec1d'))
            tracefiles.append(os.path.join(master_path,
                            '{0}.gz'.format(MasterFrame.construct_file_name('Edges', trace_key))))
#                                           MasterFrame.construct_file_name('Trace', trace_key)))
            waveimgfiles.append(os.path.join(master_path,
                                             MasterFrame.construct_file_name('Wave', wave_key)))
            tiltfiles.append(os.path.join(master_path,
                                          MasterFrame.construct_file_name('Tilts', wave_key)))

        nfiles = len(spec2d_files)

        specobjs_list = []
        head1d_list = []
        tslits_dict_list = []
        # TODO Sort this out with the correct detector extensions etc.
        # Read in the image stacks
        waveimgfile, tiltfile, tracefile = None, None, None
        for ifile in range(nfiles):
            # Load up the calibs, if needed
            if waveimgfiles[ifile] != waveimgfile:
                waveimg = WaveImage.from_master_file(waveimgfiles[ifile]).image
            if tiltfile != tiltfiles[ifile]:
                tilts = WaveTilts.from_master_file(tiltfiles[ifile]).tilts_dict
            # Save
            waveimgfile = waveimgfiles[ifile]
            tiltfile = tiltfiles[ifile]
            #
            hdu = fits.open(spec2d_files[ifile])
            # One detector, sky sub for now
            names = [hdu[i].name for i in range(len(hdu))]
            # science image
            try:
                exten = names.index('DET{:s}-PROCESSED'.format(sdet))
            except:  # Backwards compatability
                coadd.det_error_msg(exten, sdet)
            sciimg = hdu[exten].data
            # skymodel
            try:
                exten = names.index('DET{:s}-SKY'.format(sdet))
            except:  # Backwards compatability
                coadd.det_error_msg(exten, sdet)
            skymodel = hdu[exten].data
            # Inverse variance model
            try:
                exten = names.index('DET{:s}-IVARMODEL'.format(sdet))
            except ValueError:  # Backwards compatability
                coadd.det_error_msg(exten, sdet)
            sciivar = hdu[exten].data
            # Mask
            try:
                exten = names.index('DET{:s}-MASK'.format(sdet))
            except ValueError:  # Backwards compatability
                coadd.det_error_msg(exten, sdet)
            mask = hdu[exten].data
            if ifile == 0:
                # the two shapes accomodate the possibility that waveimg and tilts are binned differently
                shape_wave = (nfiles, waveimg.shape[0], waveimg.shape[1])
                shape_sci = (nfiles, sciimg.shape[0], sciimg.shape[1])
                waveimg_stack = np.zeros(shape_wave, dtype=float)
                tilts_stack = np.zeros(shape_wave, dtype=float)
                sciimg_stack = np.zeros(shape_sci, dtype=float)
                skymodel_stack = np.zeros(shape_sci, dtype=float)
                sciivar_stack = np.zeros(shape_sci, dtype=float)
                mask_stack = np.zeros(shape_sci, dtype=float)
                slitmask_stack = np.zeros(shape_sci, dtype=float)

            # Slit Traces and slitmask
            if tracefile != tracefiles[ifile]:
                tslits_dict \
                    = edgetrace.EdgeTraceSet.from_file(tracefiles[ifile]).convert_to_tslits_dict()
            tracefile = tracefiles[ifile]
            #
            tslits_dict_list.append(tslits_dict)
            slitmask = pixels.tslits2mask(tslits_dict)
            slitmask_stack[ifile, :, :] = slitmask
            waveimg_stack[ifile, :, :] = waveimg
            tilts_stack[ifile, :, :] = tilts['tilts']
            sciimg_stack[ifile, :, :] = sciimg
            sciivar_stack[ifile, :, :] = sciivar
            mask_stack[ifile, :, :] = mask
            skymodel_stack[ifile, :, :] = skymodel

            # Specobjs
            if os.path.isfile(spec1d_files[ifile]):
                sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_files[ifile])
                head1d_list.append(sobjs.header)
                this_det = sobjs.DET == self.det
                specobjs_list.append(sobjs[this_det])

        # slitmask_stack = np.einsum('i,jk->ijk', np.ones(nfiles), slitmask)

        # Fill the master key dict
        head2d = head2d_list[0]
        master_key_dict = {}
        master_key_dict['frame'] = head2d['FRAMMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['bpm'] = head2d['BPMMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['bias'] = head2d['BIASMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['arc'] = head2d['ARCMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['trace'] = head2d['TRACMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['flat'] = head2d['FLATMKEY'] + '_{:02d}'.format(self.det)

        # TODO In the future get this stuff from the headers once data model finalized
        spectrograph = util.load_spectrograph(tslits_dict['spectrograph'])

        stack_dict = dict(specobjs_list=specobjs_list, tslits_dict_list=tslits_dict_list,
                          slitmask_stack=slitmask_stack,
                          sciimg_stack=sciimg_stack, sciivar_stack=sciivar_stack,
                          skymodel_stack=skymodel_stack, mask_stack=mask_stack,
                          tilts_stack=tilts_stack, waveimg_stack=waveimg_stack,
                          head1d_list=head1d_list, head2d_list=head2d_list,
                          redux_path=redux_path,
                          master_key_dict=master_key_dict,
                          spectrograph=spectrograph.spectrograph,
                          pypeline=spectrograph.pypeline)

        return stack_dict

# Multislit can coadd with:
# 1) input offsets or if offsets is None, it will find the brightest trace and compute them
# 2) specified weights, or if weights is None and auto_weights=True, it will compute weights using the brightest object

# Echelle can either stack with:
# 1) input offsets or if offsets is None, it will find the objid of brightest trace and stack all orders relative to the trace of this object.
# 2) specified weights, or if weights is None and auto_weights=True,
#    it will use wavelength dependent weights determined from the spectrum of the brightest objects objid on each order

class MultiSlitCoAdd2D(CoAdd2D):
    """
    Child of Coadd2d for Multislit and Longslit reductions. For documentation see CoAdd2d parent class above.

        # Multislit can coadd with:
        # 1) input offsets or if offsets is None, it will find the brightest trace and compute them
        # 2) specified weights, or if weights is None and auto_weights=True, it will compute weights using the brightest object


    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        super(MultiSlitCoAdd2D, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                        sn_smooth_npix=sn_smooth_npix, ir_redux=ir_redux, par=par,
                                        show=show, show_peaks=show_peaks, debug_offsets=debug_offsets,
                                        debug=debug, **kwargs_wave)


        ## Use Cases:
        #  1) offsets is None -- auto compute offsets from brightest object, so then default to auto_weights=True
        #  2) offsets not None, weights = None (uniform weighting) or weights is not None (input weights)
        #  3) offsets not None, auto_weights=True (Do not support)

        # Default wave_method for Multislit is linear
        kwargs_wave['wave_method'] = 'linear' if 'wave_method' not in kwargs_wave else kwargs_wave['wave_method']
        self.wave_grid, self.wave_grid_mid, self.dsamp = self.get_wave_grid(**kwargs_wave)

        if offsets is None:
            self.objid_bri, self.slitid_bri, self.snr_bar_bri, self.offsets = self.compute_offsets()

        self.use_weights = self.parse_weights(weights)

    def parse_weights(self, weights):

        if 'auto' in weights:
            rms_sn, use_weights = self.optimal_weights(self.slitid_bri, self.objid_bri, const_weights=True)
            return use_weights
        elif 'uniform' in weights:
            return 'uniform'
        elif isinstance(weights, (list, np.ndarray)):
            if len(weights) != self.nexp:
                msgs.error('If weights are input it must be a list/array with same number of elements as exposures')
            return weights
        else:
            msgs.error('Unrecognized format for weights')

    # TODO When we run multislit, we actually compute the rebinned images twice. Once here to compute the offsets
    # and another time to weighted_combine the images in compute2d. This could be sped up
    def compute_offsets(self):

        objid_bri, slitid_bri, snr_bar_bri = self.get_brightest_obj(self.stack_dict['specobjs_list'], self.nslits)
        msgs.info('Determining offsets using brightest object on slit: {:d} with avg SNR={:5.2f}'.format(slitid_bri,np.mean(snr_bar_bri)))
        thismask_stack = self.stack_dict['slitmask_stack'] == slitid_bri
        trace_stack_bri = np.zeros((self.nspec, self.nexp))
        # TODO Need to think abbout whether we have multiple tslits_dict for each exposure or a single one
        for iexp in range(self.nexp):
            trace_stack_bri[:,iexp] = (self.stack_dict['tslits_dict_list'][iexp]['slit_left'][:,slitid_bri] +
                                       self.stack_dict['tslits_dict_list'][iexp]['slit_righ'][:,slitid_bri])/2.0
        # Determine the wavelength grid that we will use for the current slit/order
        wave_bins = coadd.get_wave_bins(thismask_stack, self.stack_dict['waveimg_stack'], self.wave_grid)
        dspat_bins, dspat_stack = coadd.get_spat_bins(thismask_stack, trace_stack_bri)

        sci_list = [self.stack_dict['sciimg_stack'] - self.stack_dict['skymodel_stack']]
        var_list = []

        msgs.info('Rebinning Images')
        sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack = rebin2d(
            wave_bins, dspat_bins, self.stack_dict['waveimg_stack'], dspat_stack, thismask_stack,
            (self.stack_dict['mask_stack'] == 0), sci_list, var_list)
        thismask = np.ones_like(sci_list_rebin[0][0,:,:],dtype=bool)
        nspec_psuedo, nspat_psuedo = thismask.shape
        slit_left = np.full(nspec_psuedo, 0.0)
        slit_righ = np.full(nspec_psuedo, nspat_psuedo)
        inmask = norm_rebin_stack > 0
        traces_rect = np.zeros((nspec_psuedo, self.nexp))
        sobjs = specobjs.SpecObjs()
        #specobj_dict = {'setup': 'unknown', 'slitid': 999, 'orderindx': 999, 'det': self.det, 'objtype': 'unknown',
        #                'pypeline': 'MultiSLit' + '_coadd_2d'}
        for iexp in range(self.nexp):
            sobjs_exp, _ = extract.objfind(sci_list_rebin[0][iexp,:,:], thismask, slit_left, slit_righ,
                                           inmask=inmask[iexp,:,:], ir_redux=self.ir_redux,
                                           fwhm=self.par['scienceimage']['findobj']['find_fwhm'],
                                           trim_edg=self.par['scienceimage']['findobj']['find_trim_edge'],
                                           npoly_cont=self.par['scienceimage']['findobj']['find_npoly_cont'],
                                           maxdev=self.par['scienceimage']['findobj']['find_maxdev'],
                                           ncoeff=3, sig_thresh=self.par['scienceimage']['findobj']['sig_thresh'], nperslit=1,
                                           show_trace=self.debug_offsets, show_peaks=self.debug_offsets)
            sobjs.add_sobj(sobjs_exp)
            traces_rect[:, iexp] = sobjs_exp.TRACE_SPAT
        # Now deterimine the offsets. Arbitrarily set the zeroth trace to the reference
        med_traces_rect = np.median(traces_rect,axis=0)
        offsets = med_traces_rect[0] - med_traces_rect
        # Print out a report on the offsets
        msg_string = msgs.newline()  + '---------------------------------------------'
        msg_string += msgs.newline() + ' Summary of offsets for highest S/N object   '
        msg_string += msgs.newline() + '         found on slitid = {:d}              '.format(slitid_bri)
        msg_string += msgs.newline() + '---------------------------------------------'
        msg_string += msgs.newline() + '           exp#      offset                  '
        for iexp, off in enumerate(offsets):
            msg_string += msgs.newline() + '            {:d}        {:5.2f}'.format(iexp, off)

        msg_string += msgs.newline() + '-----------------------------------------------'
        msgs.info(msg_string)
        if self.debug_offsets:
            for iexp in range(self.nexp):
                plt.plot(traces_rect[:, iexp], linestyle='--', label='original trace')
                plt.plot(traces_rect[:, iexp] + offsets[iexp], label='shifted traces')
                plt.legend()
            plt.show()

        return objid_bri, slitid_bri, snr_bar_bri, offsets

    def get_brightest_obj(self, specobjs_list, nslits):

        """
        Utility routine to find the brightest object in each exposure given a specobjs_list for MultiSlit reductions.

        Args:
            specobjs_list: list
               List of SpecObjs objects.
            echelle: bool, default=True, optional

        Returns:
            tuple: Returns the following:
                - objid: ndarray, int, shape (len(specobjs_list),):
                  Array of object ids representing the brightest object
                  in each exposure
                - slitid (int): Slit that highest S/N ratio object is on
                  (only for pypeline=MultiSlit)
                - snr_bar: ndarray, float, shape (len(list),): Average
                  S/N over all the orders for this object
        """
        nexp = len(specobjs_list)
        nspec = specobjs_list[0][0].TRACE_SPAT.shape[0]

        slit_snr_max = np.full((nslits, nexp), -np.inf)
        objid_max = np.zeros((nslits, nexp), dtype=int)
        # Loop over each exposure, slit, find the brighest object on that slit for every exposure
        for iexp, sobjs in enumerate(specobjs_list):
            for islit in range(nslits):
                ithis = sobjs.SLITID == islit
                nobj_slit = np.sum(ithis)
                if np.any(ithis):
                    objid_this = sobjs[ithis].OBJID
                    flux = np.zeros((nspec, nobj_slit))
                    ivar = np.zeros((nspec, nobj_slit))
                    wave = np.zeros((nspec, nobj_slit))
                    mask = np.zeros((nspec, nobj_slit), dtype=bool)
                    for iobj, spec in enumerate(sobjs[ithis]):
                        flux[:, iobj] = spec.OPT_COUNTS
                        ivar[:, iobj] = spec.OPT_COUNTS_IVAR
                        wave[:, iobj] = spec.OPT_WAVE
                        mask[:, iobj] = spec.OPT_MASK
                    rms_sn, weights = coadd.sn_weights(wave, flux, ivar, mask, None, const_weights=True)
                    imax = np.argmax(rms_sn)
                    slit_snr_max[islit, iexp] = rms_sn[imax]
                    objid_max[islit, iexp] = objid_this[imax]
        # Find the highest snr object among all the slits
        slit_snr = np.mean(slit_snr_max, axis=1)
        slitid = slit_snr.argmax()
        snr_bar_mean = slit_snr[slitid]
        snr_bar = slit_snr_max[slitid, :]
        objid = objid_max[slitid, :]
        if (snr_bar_mean == -np.inf):
            msgs.error('You do not appear to have a unique reference object that was traced as the highest S/N '
                       'ratio on the same slit of every exposure')

        self.snr_report(snr_bar, slitid=slitid)

        return objid, slitid, snr_bar

    # TODO add an option here to actually use the reference trace for cases where they are on the same slit and it is
    # single slit???
    def reference_trace_stack(self, slitid, offsets=None, objid=None):

        return self.offset_slit_cen(slitid, offsets)


class EchelleCoAdd2D(CoAdd2D):
    """
    Coadd Echelle reductions.
    
    For documentation see :class:`CoAdd2D`.

    Echelle can either stack with:

        - input ``offsets`` or if ``offsets`` is None, it will find
          the ``objid`` of brightest trace and stack all orders
          relative to the trace of this object.

        - specified ``weights``, or if ``weights`` is None and
          ``auto_weights`` is True, it will use wavelength dependent
          weights determined from the spectrum of the brightest
          objects ``objid`` on each order

    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        super(EchelleCoAdd2D, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                      sn_smooth_npix=sn_smooth_npix, ir_redux=ir_redux, par=par,
                                      show=show, show_peaks=show_peaks, debug_offsets=debug_offsets, debug=debug,
                                      **kwargs_wave)

        # Default wave_method for Echelle is log10
        kwargs_wave['wave_method'] = 'log10' if 'wave_method' not in kwargs_wave else kwargs_wave['wave_method']
        self.wave_grid, self.wave_grid_mid, self.dsamp = self.get_wave_grid(**kwargs_wave)

        self.objid_bri = None
        self.slitid_bri  = None
        self.snr_bar_bri = None
        if offsets is None:
            self.objid_bri, self.slitid_bri, self.snr_bar_bri = self.get_brightest_obj(self.stack_dict['specobjs_list'], self.nslits)

        self.use_weights = self.parse_weights(weights)

    def parse_weights(self, weights):

        if 'auto' in weights:
            return 'auto_echelle'
        elif 'uniform' in weights:
            return 'uniform'
        elif isinstance(weights, (list, np.ndarray)):
            if len(weights) != self.nexp:
                msgs.error('If weights are input it must be a list/array with same number of elements as exposures')
            return weights
        else:
            msgs.error('Unrecognized format for weights')

    def get_brightest_obj(self, specobjs_list, nslits):
        """
        Utility routine to find the brightest object in each exposure given a specobjs_list for Echelle reductions.

        Args:
            specobjs_list: list
               List of SpecObjs objects.
            echelle: bool, default=True, optional

        Returns:
            tuple: Returns the following:
                - objid: ndarray, int, shape (len(specobjs_list),):
                  Array of object ids representing the brightest object
                  in each exposure
                - snr_bar: ndarray, float, shape (len(list),): Average
                  S/N over all the orders for this object
        """
        nexp = len(specobjs_list)

        objid = np.zeros(nexp, dtype=int)
        snr_bar = np.zeros(nexp)
        # norders = specobjs_list[0].ech_orderindx.max() + 1
        for iexp, sobjs in enumerate(specobjs_list):
            uni_objid = np.unique(sobjs.ECH_OBJID)
            nobjs = len(uni_objid)
            order_snr = np.zeros((nslits, nobjs))
            for iord in range(nslits):
                for iobj in range(nobjs):
                    ind = (sobjs.ECH_ORDERINDX == iord) & (sobjs.ECH_OBJID == uni_objid[iobj])
                    flux = sobjs[ind][0].OPT_COUNTS
                    ivar = sobjs[ind][0].OPT_COUNTS_IVAR
                    wave = sobjs[ind][0].OPT_WAVE
                    mask = sobjs[ind][0].OPT_MASK
                    rms_sn, weights = coadd.sn_weights(wave, flux, ivar, mask, self.sn_smooth_npix, const_weights=True)
                    order_snr[iord, iobj] = rms_sn

            # Compute the average SNR and find the brightest object
            snr_bar_vec = np.mean(order_snr, axis=0)
            objid[iexp] = uni_objid[snr_bar_vec.argmax()]
            snr_bar[iexp] = snr_bar_vec[snr_bar_vec.argmax()]

        self.snr_report(snr_bar)

        return objid, None, snr_bar

    def reference_trace_stack(self, slitid, offsets=None, objid=None):
        """
        Utility function for determining the reference trace about which 2d coadds are performed.
        There are two modes of operation to determine the reference trace for the 2d coadd of a given slit/order:

         1) offsets: we stack about the center of the slit for the slit in question with the input offsets added
         2) ojbid: we stack about the trace ofa reference object for this slit given for each exposure by the input objid

        Either offsets or objid must be provided, but the code will raise an exception if both are provided.

        Args:
            slitid (int):
               The slit or order that we are currently considering
            stack_dict (dict):
               Dictionary containing all the images and keys required for perfomring 2d coadds.
            offsets (list or np.ndarray):
               An array of offsets with the same dimensionality as the nexp, the numer of images being coadded.
            objid: (list or np.ndarray):
               An array of objids with the same dimensionality as the nexp, the number of images being coadded.

        Returns:
            ref_trace_stack

            ref_trace_stack (np.ndarray):
                An array with shape (nspec, nexp) containing the reference trace for each of the nexp exposures.

        """

        if offsets is not None and objid is not None:
            msgs.errror('You can only input offsets or an objid, but not both')
        nexp = len(offsets) if offsets is not None else len(objid)
        if offsets is not None:
            return self.offset_slit_cen(slitid, offsets)
        elif objid is not None:
            specobjs_list = self.stack_dict['specobjs_list']
            nspec = specobjs_list[0][0].TRACE_SPAT.shape[0]
            # Grab the traces, flux, wavelength and noise for this slit and objid.
            ref_trace_stack = np.zeros((nspec, nexp), dtype=float)
            for iexp, sobjs in enumerate(specobjs_list):
                ithis = (sobjs.ECH_ORDERINDX == slitid) & (sobjs.ECH_OBJID == objid[iexp])
                ref_trace_stack[:, iexp] = sobjs[ithis].TRACE_SPAT
            return ref_trace_stack
        else:
            msgs.error('You must input either offsets or an objid to determine the stack of reference traces')
            return None

