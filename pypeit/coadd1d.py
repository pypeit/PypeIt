
"""
Coadding module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy import stats

from pypeit.spectrographs.util import load_spectrograph
from pypeit.onespec import OneSpec
from pypeit.orderstack import OrderStack
from pypeit import utils
from pypeit import sensfunc
from pypeit import specobjs
from pypeit import msgs
from pypeit.core import coadd, flux_calib
from pypeit.history import History



class CoAdd1D:

    @classmethod
    def get_instance(cls, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None,
                     setup_id=None, debug=False, show=False, chk_version=True):
        """
        Superclass factory method which generates the subclass instance. See
        :class:`CoAdd1D` instantiation for argument descriptions.
        """
        pypeline = fits.getheader(spec1dfiles[0])['PYPELINE'] + 'CoAdd1D'
        return next(c for c in utils.all_subclasses(CoAdd1D) if c.__name__ == pypeline)(
            spec1dfiles, objids, spectrograph=spectrograph, par=par, sensfuncfile=sensfuncfile,
            setup_id=setup_id, debug=debug, show=show, chk_version=chk_version)

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None,
                 setup_id=None, debug=False, show=False, chk_version=True):
        """

        Args:
            spec1dfiles (list):
                List of strings which are the spec1dfiles
            objids (list):
                List of strings which are the objids for the object in each
                spec1d file that you want to coadd.
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`, optional):
                Spectrograph object.
            par (:class:`pypeit.par.pypeitpar.Coadd1DPar`, optional):
                PypeIt parameter set object for Coadd1D
            sensfuncile (str or list of strings, optional):
                File or list of files holding the sensitivity function. This is
                required for echelle coadds only.
            setup_id (str or list of strings, optional):
                A string or list of strings identifiying the setup IDs to coadd.
                This is only used for echelle coadds where a loop over the
                different echelle setups is performed.  If None, it will be
                assumed that all the input files, objids, and sensfuncfiles
                correspond to the same setup.
            debug (bool, optional)
                Debug. Default = False
            show (bool, optional):
                Debug. Default = True
            chk_version (:obj:`bool`, optional):
                When reading in existing files written by PypeIt, perform strict
                version checking to ensure a valid file.  If False, the code
                will try to keep going, but this may lead to faults and quiet
                failures.  User beware!
        """
        # Instantiate attributes
        self.spec1dfiles = spec1dfiles
        self.objids = objids

        # Optional
        if spectrograph is not None:
            self.spectrograph = spectrograph
        else:
            header = fits.getheader(spec1dfiles[0])
            self.spectrograph = load_spectrograph(header['PYP_SPEC'])
        if par is None:
            self.par = self.spectrograph.default_pypeit_par()['coadd1d']
        else:
            self.par = par
        #
        self.debug = debug
        self.show = show
        self.chk_version = chk_version
        self.nexp = len(self.spec1dfiles) # Number of exposures
        self.coaddfile = None
        self.gpm_exp = np.ones(self.nexp, dtype=bool).tolist()  # list of bool indicating the exposures that have been coadded

    def run(self):
        """
        Runs the coadding
        """

        # Coadd the data
        # if there are multiple orders/slits, the stacks will be extracted from the coadd1d object, otherwise it'll be None
        self.wave_grid_mid, self.wave_coadd, self.flux_coadd, self.ivar_coadd, self.gpm_coadd, self.order_stacks = self.coadd()

        # Scale to a filter magnitude?
        if self.par['filter'] != 'none':
            scale = flux_calib.scale_in_filter(self.wave_coadd, self.flux_coadd, self.gpm_coadd, self.par)
            self.flux_coadd *= scale
            self.ivar_coadd = self.ivar_coadd / scale**2



    def load(self):
        """
        Load the arrays we need for performing coadds. Dummy method overloaded by children.
        """
        msgs.error('This method is undefined in the base classes and should only be called by the subclasses')


    def save(self, coaddfile, telluric=None, obj_model=None, overwrite=True):
        """
        Generate a :class:`~pypeit.onespec.OneSpec` object and write it to disk.

        Args:
            coaddfile (str):
                File to output coadded spectrum to.
            telluric (`numpy.ndarray`_, optional):
                Telluric model.
            obj_model (str, optional):
                Name of the object model
            overwrite (bool, optional):
               Overwrite existing file?
        """
        self.coaddfile = coaddfile
        # Generate the spectrum container object
        onespec = OneSpec(wave=self.wave_coadd, wave_grid_mid=self.wave_grid_mid, flux=self.flux_coadd,
                          PYP_SPEC=self.spectrograph.name, ivar=self.ivar_coadd,
                          sigma=np.sqrt(utils.inverse(self.ivar_coadd)),
                          mask=self.gpm_coadd.astype(int),
                          ext_mode=self.par['ex_value'], fluxed=self.par['flux_value'])

        # TODO This is a hack, not sure how to merge the headers at present
        onespec.head0 = self.headers[0]

        # Add history entries for coadding.
        history = History()
        history.add_coadd1d(self.spec1dfiles, self.objids, gpm_exp=self.gpm_exp)

        # Add on others
        if telluric is not None:
            onespec.telluric = telluric
        if obj_model is not None:
            onespec.obj_model = obj_model

        #save the order stacks from the echelle reduction
        if self.order_stacks is not None and (fits.getheader(self.spec1dfiles[0])['PYPELINE'] == 'Echelle'):
            for setup_num, setup_val in enumerate(self.unique_setups):
                if len(self.unique_setups) > 1:
                    wave_stack = np.array(self.order_stacks[0][setup_num])
                    flux_stack = np.array(self.order_stacks[1][setup_num])
                    ivar_stack = np.array(self.order_stacks[2][setup_num])
                    mask_stack = np.array(self.order_stacks[3][setup_num]).astype(int)
                    sigma_stack = np.sqrt(utils.inverse(ivar_stack))
                else:
                    wave_stack = np.array(self.order_stacks[0])
                    flux_stack = np.array(self.order_stacks[1])
                    ivar_stack = np.array(self.order_stacks[2])
                    mask_stack = np.array(self.order_stacks[3]).astype(int)
                    sigma_stack = np.sqrt(utils.inverse(ivar_stack))
                orderstack = OrderStack(wave_stack, 
                                        flux_stack, 
                                        ivar_stack=ivar_stack, 
                                        mask_stack=mask_stack, 
                                        sigma_stack=sigma_stack,
                                        PYP_SPEC=self.spectrograph.name, 
                                        ext_mode=self.par['ex_value'], fluxed=self.par['flux_value'],
                                        setup_name = setup_val)
                orderstack.head0 = self.headers[setup_num]
                orderstack.to_file(coaddfile.split('.fits')[0] + '_orderstack' + setup_val + '.fits', history=history, overwrite=overwrite)
        # Write
        onespec.to_file(coaddfile, history=history, overwrite=overwrite)

    def coadd(self):
        """
        Dummy method overloaded by sub-classes
        """
        raise NotImplementedError(f'coadding function not defined for {self.__class__.__name__}!')


class MultiSlitCoAdd1D(CoAdd1D):
    """
    Child of CoAdd1d for Multislit and Longslit reductions.
    """

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None, setup_id=None, 
                 debug=False, show=False, chk_version=False):
        """
        See :class:`CoAdd1D` instantiation for argument descriptions.
        """
        super().__init__(spec1dfiles, objids, spectrograph=spectrograph, par=par, sensfuncfile=sensfuncfile,
                         setup_id=setup_id, debug=debug, show=show, chk_version=chk_version)


    def load(self):
        """
        Load the arrays we need for performing coadds.

        Returns
        -------
        waves : list of float `numpy.ndarray`_
            List of wavelength arrays. The length of the list is nexp. The
            arrays can have different shapes.
        fluxes : list of float `numpy.ndarray`_
            List of flux arrays. The arrays can have different shapes, but all
            are aligned with what is in waves.
        ivars : list of float `numpy.ndarray`_
            List of inverse variance arrays. The arrays can have different
            shapes, but all are aligned with what is in waves.
        gpms : list of bool `numpy.ndarray`_
            List of good pixel mask variance arrays. The arrays can have
            different shapes, but all are aligned with what is in waves.
        headers : list of header objects
            List of headers of length nexp
        """
        waves, fluxes, ivars, gpms, headers = [], [], [], [], []
        for iexp in range(self.nexp):
            sobjs = specobjs.SpecObjs.from_fitsfile(self.spec1dfiles[iexp],
                                                    chk_version=self.chk_version)
            indx = sobjs.name_indices(self.objids[iexp])
            if not np.any(indx):
                msgs.error(
                    "No matching objects for {:s}.  Odds are you input the wrong OBJID".format(self.objids[iexp]))
            if np.sum(indx) > 1:
                msgs.error("Error in spec1d file for exposure {:d}: "
                           "More than one object was identified with the OBJID={:s} in file={:s}".format(
                    iexp, self.objids[iexp], self.spec1dfiles[iexp]))
            wave_iexp, flux_iexp, ivar_iexp, gpm_iexp, _, _, _, header = \
                sobjs[indx].unpack_object(ret_flam=self.par['flux_value'], extract_type=self.par['ex_value'])
            waves.append(wave_iexp)
            fluxes.append(flux_iexp)
            ivars.append(ivar_iexp)
            gpms.append(gpm_iexp)
            header_out = header.copy()
            if 'RA' in sobjs[indx][0].keys() and 'DEC' in sobjs[indx][0].keys():
                header_out['RA_OBJ'] = sobjs[indx][0]['RA']
                header_out['DEC_OBJ'] = sobjs[indx][0]['DEC']
            headers.append(header_out)

        return waves, fluxes, ivars, gpms, headers

    def check_exposures(self):
        """
        Check if there are bad exposures.
        Exposures with flux masked everywhere are always removed.
        Exposures that are considered bad based on their S/N compared to the
        average S/N among all the exposures, are removed only if self.par['sigrej_exp'] is set.
        The attributes self.waves, self.fluxes, self.ivars, self.gpms need to be defined.

        Returns
        -------
        gpm_exp: list of bool
            List of boolean that indicates which exposures
            have been coadded. The length of the list is nexp.
        _waves : list of float `numpy.ndarray`_
            Updated list of wavelength arrays.
        _fluxes : list of float `numpy.ndarray`_
            Updated list of flux arrays.
        _ivars : list of float `numpy.ndarray`_
            Updated list of inverse variance arrays.
        _gpms : list of bool `numpy.ndarray`_
            Updated list of good pixel mask variance arrays.
        """

        # initialize the exposures lists
        _waves = [wave for wave in self.waves]
        _fluxes = [flux for flux in self.fluxes]
        _ivars = [ivar for ivar in self.ivars]
        _gpms = [gpm for gpm in self.gpms]
        _spec1dfiles = [spec1dfile for spec1dfile in self.spec1dfiles]
        _objids = [objid for objid in self.objids]

        # good exposures index
        goodindx_exp = np.arange(self.nexp)

        # check if there are exposures that are completely masked out, i.e., gpms = False for all spectral pixels
        masked_exps = [np.all(np.logical_not(gpm)) for gpm in _gpms]
        if np.any(masked_exps):
            msgs.warn(f'The following exposure(s) is/are completely masked out. It/They will not be coadded.')
            [msgs.warn(f"Exposure {i}: {fname.split('/')[-1]}  {obj}")
             for i, (fname, obj, masked_exp) in enumerate(zip(_spec1dfiles, _objids, masked_exps)) if masked_exp]
            # remove masked out exposure
            _waves = [wave for (wave, masked_exp) in zip(_waves, masked_exps) if not masked_exp]
            _fluxes = [flux for (flux, masked_exp) in zip(_fluxes, masked_exps) if not masked_exp]
            _ivars = [ivar for (ivar, masked_exp) in zip(_ivars, masked_exps) if not masked_exp]
            _gpms = [gpm for (gpm, masked_exp) in zip(_gpms, masked_exps) if not masked_exp]
            _spec1dfiles = [spec1dfile for (spec1dfile, masked_exp) in zip(_spec1dfiles, masked_exps) if not masked_exp]
            _objids = [objid for (objid, masked_exp) in zip(_objids, masked_exps) if not masked_exp]
            # update good exposures index
            goodindx_exp = goodindx_exp[np.logical_not(masked_exps)]

        # check if there is still at least 1 exposure left
        if len(_fluxes) < 1:
            msgs.error('At least 1 unmasked exposures are required for coadding.')

        # check if there is any bad exposure by comparing the rms_sn with the median rms_sn among all exposures
        if len(_fluxes) > 2:
            # Evaluate the rms_sn
            rms_sn, _ = coadd.calc_snr(_fluxes, _ivars, _gpms)
            # some stats
            mean, med, sigma = stats.sigma_clipped_stats(rms_sn, sigma_lower=2., sigma_upper=2.)
            _sigrej = self.par['sigrej_exp'] if self.par['sigrej_exp'] is not None else 10.0
            # we set thresh_value to never be less than 0.2
            thresh_value = round(0.2 + med + _sigrej * sigma, 2)
            bad_exps = rms_sn > thresh_value
            if np.any(bad_exps):
                warn_msg = f'The following exposure(s) has/have S/N > {thresh_value:.2f} ' \
                           f'({_sigrej} sigma above the median S/N in the stack).'
                if self.par['sigrej_exp'] is not None:
                        warn_msg += ' It/They WILL NOT BE COADDED.'
                msgs.warn(warn_msg)
                [msgs.warn(f"Exposure {i}: {fname.split('/')[-1]}  {obj}")
                 for i, (fname, obj, bad_exp) in enumerate(zip(_spec1dfiles, _objids, bad_exps)) if bad_exp]
                if self.par['sigrej_exp'] is not None:
                    # remove bad exposure
                    _waves = [wave for (wave, bad_exp) in zip(_waves, bad_exps) if not bad_exp]
                    _fluxes = [flux for (flux, bad_exp) in zip(_fluxes, bad_exps) if not bad_exp]
                    _ivars = [ivar for (ivar, bad_exp) in zip(_ivars, bad_exps) if not bad_exp]
                    _gpms = [gpm for (gpm, bad_exp) in zip(_gpms, bad_exps) if not bad_exp]
                    _spec1dfiles = [spec1dfile for (spec1dfile, bad_exp) in zip(_spec1dfiles, bad_exps) if not bad_exp]
                    _objids = [objid for (objid, bad_exp) in zip(_objids, bad_exps) if not bad_exp]
                    # update good exposures index
                    goodindx_exp = goodindx_exp[np.logical_not(bad_exps)]

        # gpm for the exposures, i.e., which exposures have been coadded
        gpm_exp = np.zeros(self.nexp, dtype=bool)
        gpm_exp[goodindx_exp] = True

        return gpm_exp.tolist(), _waves, _fluxes, _ivars, _gpms

    def coadd(self):
        """
        Perform coadd for for Multi/Longslit data using multi_combspec

        Returns:
            tuple: see objects returned by
            :func:`~pypeit.core.coadd.multi_combspec`.
        """
        # Load the data
        self.waves, self.fluxes, self.ivars, self.gpms, self.headers = self.load()
        # check if there are bad exposures and remove them
        self.gpm_exp, _waves, _fluxes, _ivars, _gpms = self.check_exposures()
        # Perform the coadd
        wave_grid_mid, wave_coadd, flux_coadd, ivar_coadd, gpm_coadd = \
            coadd.multi_combspec(_waves, _fluxes, _ivars, _gpms,
            sn_smooth_npix=self.par['sn_smooth_npix'], wave_method=self.par['wave_method'],
            dv=self.par['dv'], dwave=self.par['dwave'], dloglam=self.par['dloglam'],
            wave_grid_min=self.par['wave_grid_min'], wave_grid_max=self.par['wave_grid_max'],
            spec_samp_fact=self.par['spec_samp_fact'], ref_percentile=self.par['ref_percentile'],
            maxiter_scale=self.par['maxiter_scale'], sigrej_scale=self.par['sigrej_scale'],
            scale_method=self.par['scale_method'], sn_min_medscale=self.par['sn_min_medscale'],
            sn_min_polyscale=self.par['sn_min_polyscale'], weight_method = self.par['weight_method'],
            maxiter_reject=self.par['maxiter_reject'], lower=self.par['lower'], upper=self.par['upper'],
            maxrej=self.par['maxrej'], sn_clip=self.par['sn_clip'], debug=self.debug, show=self.show)

        # return
        return wave_grid_mid, wave_coadd, flux_coadd, ivar_coadd, gpm_coadd, None



class EchelleCoAdd1D(CoAdd1D):
    """
    Child of CoAdd1d for Echelle reductions.
    """

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None,
                 setup_id=None, debug=False, show=False, chk_version=True):
        """
        See :class:`CoAdd1D` instantiation for argument descriptions.
        """
        super().__init__(spec1dfiles, objids, spectrograph=spectrograph, par=par, sensfuncfile=sensfuncfile,
                         setup_id=setup_id, debug=debug, show=show,
                         chk_version=chk_version)

        if sensfuncfile is None:
            msgs.error('sensfuncfile is a required argument for echelle coadding')

        self.sensfuncfile = self.nexp * [sensfuncfile] if isinstance(sensfuncfile, str) else sensfuncfile
        nsens = len(self.sensfuncfile)
        if nsens == 1:
            self.sensfuncfile = self.nexp * [self.sensfuncfile[0]]
            nsens = self.nexp
        if nsens != self.nexp:
            msgs.error('Must enter either one sensfunc file for all exposures or one sensfunc file for '
                       f'each exposure.  Entered {nsens} files for {self.nexp} exposures.')

        if setup_id is None:
            self.setup_id = self.nexp*['A']
        else:
            self.setup_id = self.nexp*[setup_id] if isinstance(setup_id, str) else setup_id
        nsetup = len(self.setup_id)
        if nsetup == 1:
            self.setup_id = self.nexp * [self.setup_id[0]]
            nsetup = self.nexp
        if nsetup != self.nexp:
            msgs.error('Must enter either a single setup_id for all exposures or one setup_id for '
                       f'each exposure.  Entered {nsetup} files for {self.nexp} exposures.')


        self.unique_setups = np.unique(self.setup_id).tolist()
        self.nsetups_unique = len(self.unique_setups)


    def coadd(self):
        """
        Perform coadd for echelle data using ech_combspec

        Returns
        -------
        wave_grid_mid : `numpy.ndarray`_
            Wavelength grid (in Angstrom) evaluated at the bin centers,
            uniformly-spaced either in lambda or log10-lambda/velocity. See
            core.wavecal.wvutils.py for more.  shape=(ngrid,)
        wave_coadd : `numpy.ndarray`_
            Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the
            wavelengths of each spectrum that contriuted to a bin in the input
            wave_grid wavelength grid. It thus has ngrid elements, whereas
            wave_grid has ngrid+1 elements to specify the ngrid total number
            of bins. Note that wave_giant_stack is NOT simply the wave_grid
            bin centers, since it computes the weighted average;
        flux_coadd : `numpy.ndarray`_
            Final stacked spectrum on wave_stack wavelength grid.  shape=(ngrid,)
        ivar_coadd : `numpy.ndarray`_
            Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to
            weighting and masking. shape=(ngrid,)
        gpm_coadd : `numpy.ndarray`_
            Mask for stacked spectrum on wave_stack wavelength grid. True=Good. shape=(ngrid,)
        order_stacks_output : `numpy.ndarray`_, None
            Stacked orders.  None if not echelle data.  shape=(norders, ngrid)
        """

        # Load the data
        self.waves, self.fluxes, self.ivars, self.gpms, self.weights_sens, self.headers = self.load()
        wave_grid_mid, (wave_coadd, flux_coadd, ivar_coadd, gpm_coadd),  order_stacks \
                = coadd.ech_combspec(self.waves, self.fluxes, self.ivars, self.gpms, self.weights_sens,
                                     setup_ids=self.unique_setups,
                                     nbests=self.par['nbests'],
                                     wave_method=self.par['wave_method'],
                                     dv=self.par['dv'], dwave=self.par['dwave'], dloglam=self.par['dloglam'],
                                     wave_grid_min=self.par['wave_grid_min'],
                                     wave_grid_max=self.par['wave_grid_max'],
                                     spec_samp_fact=self.par['spec_samp_fact'],
                                     ref_percentile=self.par['ref_percentile'],
                                     maxiter_scale=self.par['maxiter_scale'], 
                                     sigrej_scale=self.par['sigrej_scale'],
                                     scale_method=self.par['scale_method'],
                                     sn_min_medscale=self.par['sn_min_medscale'],
                                     sn_min_polyscale=self.par['sn_min_polyscale'],
                                     maxiter_reject=self.par['maxiter_reject'],
                                     lower=self.par['lower'], upper=self.par['upper'],
                                     maxrej=self.par['maxrej'], sn_clip=self.par['sn_clip'],
                                     debug=self.debug, show=self.show, show_exp=self.show)
        
        
        return wave_grid_mid, wave_coadd, flux_coadd, ivar_coadd, gpm_coadd, order_stacks

    def load_ech_arrays(self, spec1dfiles, objids, sensfuncfiles):
        """
        Load the arrays we need for performing coadds for a single setup.

        Args:
            spec1dfiles (list):
                List of spec1d files for this setup.
            objids (list):
                List of objids. This is aligned with spec1dfiles
            sensfuncfile (list):
                List of sensfuncfiles. This is aligned with spec1dfiles and objids

        Returns:
            tuple: waves, fluxes, ivars, gpms, header. Each array has shape =
            (nspec, norders, nexp)

        """
        nexp = len(spec1dfiles)
        for iexp in range(nexp):
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfiles[iexp], chk_version=self.chk_version)
            indx = sobjs.name_indices(objids[iexp])
            if not np.any(indx):
                msgs.error("No matching objects for {:s}.  Odds are you input the wrong OBJID".format(objids[iexp]))
            wave_iexp, flux_iexp, ivar_iexp, gpm_iexp, _, _, _, header = \
                    sobjs[indx].unpack_object(ret_flam=self.par['flux_value'], extract_type=self.par['ex_value'])
            # This np.atleast2d hack deals with the situation where we are wave_iexp is actually Multislit data, i.e. we are treating
            # it like an echelle spectrograph with a single order. This usage case arises when we want to use the
            # echelle coadding code to combine echelle and multislit data
            if wave_iexp.ndim == 1:
                wave_iexp, flux_iexp, ivar_iexp, gpm_iexp = np.atleast_2d(wave_iexp).T, np.atleast_2d(flux_iexp).T, np.atleast_2d(ivar_iexp).T, np.atleast_2d(gpm_iexp).T
            weights_sens_iexp = sensfunc.SensFunc.sensfunc_weights(sensfuncfiles[iexp], wave_iexp,
                                                                   debug=self.debug,
                                                                   chk_version=self.chk_version)
            # Allocate arrays on first iteration
            # TODO :: We should refactor to use a list of numpy arrays, instead of a 2D numpy array.
            if iexp == 0:
                waves = np.zeros(wave_iexp.shape + (nexp,))
                fluxes = np.zeros_like(waves)
                ivars = np.zeros_like(waves)
                weights_sens = np.zeros_like(waves)
                gpms = np.zeros_like(waves, dtype=bool)
                header_out = header
                if 'RA' in sobjs[indx][0].keys() and 'DEC' in sobjs[indx][0].keys():
                    header_out['RA_OBJ']  = sobjs[indx][0]['RA']
                    header_out['DEC_OBJ'] = sobjs[indx][0]['DEC']

            # Store the information
            waves[...,iexp], fluxes[...,iexp], ivars[..., iexp], gpms[...,iexp], weights_sens[...,iexp] \
                = wave_iexp, flux_iexp, ivar_iexp, gpm_iexp, weights_sens_iexp


        return waves, fluxes, ivars, gpms, weights_sens, header_out


    def load(self):
        """
        Load the arrays we need for performing echelle coadds.

        Returns
        -------
        waves : list
            List of arrays with the wavelength arrays for each setup. The length
            of the list equals the number of unique setups and each arrays in
            the list has shape = (nspec, norders, nexp)
        fluxes : list
            List of arrays with the flux arrays for each setup. The length of
            the list equals the number of unique setups and each arrays in the
            list has shape = (nspec, norders, nexp)
        ivars : list
            List of arrays with the ivar arrays for each setup. The length of
            the list equals the number of unique setups and each arrays in the
            list has shape = (nspec, norders, nexp)
        gpms : list
            List of arrays with the gpm arrays for each setup. The length of the
            list equals the number of unique setups and each arrays in the list
            has shape = (nspec, norders, nexp)
        weights_sens : list
            List of arrays with the sensfunc weights for each setup. The length
            of the list equals the number of unique setups and each arrays in
            the list has shape = (nspec, norders, nexp)
        headers : list
            List of headers for each setup. The length of the list is the number
            of unique setups.

        """

        _setup = np.asarray(self.setup_id)
        _sensfuncfiles = np.asarray(self.sensfuncfile)
        _spec1dfiles = np.asarray(self.spec1dfiles)
        _objids = np.asarray(self.objids)
        waves, fluxes, ivars, gpms, weights_sens, headers = [], [], [], [], [], []
        combined = [waves, fluxes, ivars, gpms, weights_sens, headers]
        for uniq_setup in self.unique_setups:
            setup_indx = _setup == uniq_setup
            loaded = self.load_ech_arrays(_spec1dfiles[setup_indx], _objids[setup_indx], _sensfuncfiles[setup_indx])
            for c, l in zip(combined, loaded):
                c.append(l)



        return waves, fluxes, ivars, gpms, weights_sens, headers


class SlicerIFUCoAdd1D(MultiSlitCoAdd1D):
    """
    Child of MultiSlitCoAdd1d for SlicerIFU reductions.
    """

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None, setup_id=None, debug=False, show=False, chk_version=False):
        """
        See :class:`CoAdd1D` instantiation for argument descriptions.
        """
        super().__init__(spec1dfiles, objids, spectrograph=spectrograph, par = par, sensfuncfile = sensfuncfile,
                         setup_id=setup_id, debug = debug, show = show, chk_version=chk_version)
