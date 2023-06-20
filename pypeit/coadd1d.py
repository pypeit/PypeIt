
"""
Coadding module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit.spectrographs.util import load_spectrograph
from pypeit.onespec import OneSpec
from pypeit import sensfunc
from pypeit import specobjs
from pypeit import msgs
from pypeit.core import coadd, flux_calib
from pypeit.history import History


class CoAdd1D:

    @classmethod
    def get_instance(cls, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None, setup_id=None,
                     debug=False, show=False):
        """
        Superclass factory method which generates the subclass instance. See __init__ docs for arguments.
        """
        pypeline = fits.getheader(spec1dfiles[0])['PYPELINE'] + 'CoAdd1D'
        return next(c for c in cls.__subclasses__() if c.__name__ == pypeline)(
            spec1dfiles, objids, spectrograph=spectrograph, par=par, sensfuncfile=sensfuncfile, setup_id=setup_id,
            debug=debug, show=show)

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None, setup_id=None,
                 debug=False, show=False):
        """

        Args:
            spec1dfiles (list):
               List of strings which are the spec1dfiles
            objids (list):
               List of strings which are the objids for the object in each spec1d file that you want to coadd
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`, optional):
            par (:class:`pypeit.par.pypeitpar.Coadd1DPar`, optional):
               Pypeit parameter set object for Coadd1D
            sensfuncile (str, optional):
               File holding the sensitivity function. This is required for echelle coadds only.
            debug (bool, optional)
               Debug. Default = False
            show (bool, optional):
               Debug. Default = True
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
        self.nexp = len(self.spec1dfiles) # Number of exposures
        self.coaddfile = None

    def run(self):
        """
        Runs the coadding
        """

        # Coadd the data
        self.wave_grid_mid, self.wave_coadd, self.flux_coadd, self.ivar_coadd, self.gpm_coadd = self.coadd()
        # Scale to a filter magnitude?
        if self.par['filter'] != 'none':
            scale = flux_calib.scale_in_filter(self.wave_coadd, self.flux_coadd, self.gpm_coadd, self.par)
            self.flux_coadd *= scale
            self.ivar_coadd = self.ivar_coadd / scale**2



    def load(self):
        """
        Load the arrays we need for performing coadds. Dummy method overloaded by children.

        Returns:
            None
        """

        return None

    def save(self, coaddfile, telluric=None, obj_model=None, overwrite=True):
        """
        Generate a :class:`OneSpec` object and write it to disk.

        Args:
            coaddfile (str):
               File to output coadded spectrum to.
            telluric (`numpy.ndarray`_):
            obj_model (str):
            overwrite (bool):
               Overwrite existing file?
        """
        self.coaddfile = coaddfile
        wave_gpm = self.wave_coadd > 1.0
        # Generate the spectrum container object
        onespec = OneSpec(wave=self.wave_coadd[wave_gpm], wave_grid_mid=self.wave_grid_mid[wave_gpm], flux=self.flux_coadd[wave_gpm],
                          PYP_SPEC=self.spectrograph.name, ivar=self.ivar_coadd[wave_gpm],
                          mask=self.gpm_coadd[wave_gpm].astype(int),
                          ext_mode=self.par['ex_value'], fluxed=self.par['flux_value'])

        # TODO This is a hack, not sure how to merge the headers at present
        onespec.head0 = self.headers[0]

        # Add history entries for coadding.
        history = History()
        history.add_coadd1d(self.spec1dfiles, self.objids)

        # Add on others
        if telluric is not None:
            onespec.telluric  = telluric[wave_gpm]
        if obj_model is not None:
            onespec.obj_model = obj_model[wave_gpm]
        # Write
        onespec.to_file(coaddfile, history=history, overwrite=overwrite)

    def coadd(self):
        """
        Dummy method overloaded by sub-classes

        Returns:
            :obj:`tuple`:  four items
              - wave
              - flux
              - ivar
              - gpm

        """
        return (None,)*4


class MultiSlitCoAdd1D(CoAdd1D):
    """
    Child of CoAdd1d for Multislit and Longslit reductions
    """

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None, setup_id=None, debug=False, show=False):
        """
        See `CoAdd1D` doc string
        """
        super().__init__(spec1dfiles, objids, spectrograph=spectrograph, par = par, sensfuncfile = sensfuncfile,
                         setup_id=setup_id, debug = debug, show = show)


    def load(self):
        """
        Load the arrays we need for performing coadds.

        Returns:
            tuple:
               - waves, fluxes, ivars, gpms, header
        """
        waves, fluxes, ivars, gpms, headers = [], [], [], [], []
        for iexp in range(self.nexp):
            sobjs = specobjs.SpecObjs.from_fitsfile(self.spec1dfiles[iexp], chk_version=self.par['chk_version'])
            indx = sobjs.name_indices(self.objids[iexp])
            if not np.any(indx):
                msgs.error(
                    "No matching objects for {:s}.  Odds are you input the wrong OBJID".format(self.objids[iexp]))
            wave_iexp, flux_iexp, ivar_iexp, gpm_iexp, trace_spec, trace_spat, meta_spec, header = \
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

    def coadd(self):
        """
        Perform coadd for for Multi/Longslit data using multi_combspec

        Returns:
            tuple
              - wave_grid_mid, wave, flux, ivar, gpm

        """
        # Load the data
        self.waves, self.fluxes, self.ivars, self.gpms, self.headers = self.load()
        # Perform and return the coadd
        return coadd.multi_combspec(
            self.waves, self.fluxes, self.ivars, self.gpms,
            sn_smooth_npix=self.par['sn_smooth_npix'], wave_method=self.par['wave_method'],
            dv=self.par['dv'], wave_grid_min=self.par['wave_grid_min'], wave_grid_max=self.par['wave_grid_max'],
            spec_samp_fact=self.par['spec_samp_fact'], ref_percentile=self.par['ref_percentile'],
            maxiter_scale=self.par['maxiter_scale'], sigrej_scale=self.par['sigrej_scale'],
            scale_method=self.par['scale_method'], sn_min_medscale=self.par['sn_min_medscale'],
            sn_min_polyscale=self.par['sn_min_polyscale'], maxiter_reject=self.par['maxiter_reject'],
            lower=self.par['lower'], upper=self.par['upper'], maxrej=self.par['maxrej'], sn_clip=self.par['sn_clip'],
            debug=self.debug, show=self.show)


        return waves, fluxes, ivars, gpms, header_out


class EchelleCoAdd1D(CoAdd1D):
    """
    Child of CoAdd1d for Echelle reductions
    """

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfuncfile=None, setup_id=None,
                 debug=False, show=False):
        """
        See `CoAdd1D` doc string

        """
        super().__init__(spec1dfiles, objids, spectrograph=spectrograph, par = par, sensfuncfile = sensfuncfile,
                         setup_id=setup_id, debug = debug, show = show)

        if sensfuncfile is None:
            msgs.error('sensfuncfile is a required argument for echelle coadding')
        else:
            nsens = len(sensfuncfile)
            if not ((nsens == 1) or (nsens == self.nexp)):
                msgs.error('Invalid length of sensfuncfile len(sensfuncfile)={:d}'.format(len(nsens)))
            elif nsens == 1:
                _sensfuncfile = sensfuncfile if isinstance(sensfuncfile, str) else sensfuncfile[0]
                self.sensfuncfile = self.nexp*[_sensfuncfile]
            elif nsens == self.nexp:
                self.sensfuncfile = sensfuncfile

        if setup_id is None:
            self.setup_id = self.nexp*['A']
        else:
            len_setup = len(setup_id)
            if not ((len_setup == 1) or (len_setup == self.nexp)):
                msgs.error('Invalid length of setup_id len(setup_id)={:d}'.format(len(len_setup)))
            _setup_id = list(setup_id) if isinstance(setup_id, str) else setup_id
            self.setup_id = _setup_id if len(_setup_id) == self.nexp else self.nexp*_setup_id

        self.unique_setups = np.unique(self.setup_id).tolist()
        self.nsetups = len(self.unique_setups)

        self.sensfuncfile = sensfuncfile
        self.setup_id= setup_id

    def coadd(self):
        """
        Perform coadd for echelle data using ech_combspec

        Returns:
            tuple
              - wave_grid_mid, wave, flux, ivar, gpm

        """

        # Load the data
        self.waves, self.fluxes, self.ivars, self.gpms, self.weights_sens, self.headers = self.load()
        wave_grid_mid, (wave_coadd, flux_coadd, ivar_coadd, gpm_coadd), order_stacks \
                = coadd.ech_combspec(self.waves, self.fluxes, self.ivars, self.gpms, self.weights_sens,
                                     setup_ids=self.unique_setups,
                                     nbests=self.par['nbests'],
                                     sn_smooth_npix=self.par['sn_smooth_npix'],
                                     wave_method=self.par['wave_method'],
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
                                     debug=self.debug, show=self.show)


        return wave_grid_mid, wave_coadd, flux_coadd, ivar_coadd, gpm_coadd


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
            tuple:
               - waves, fluxes, ivars, gpms, header. Each array has shape = (nspec, norders, nexp)

        """
        nexp = len(spec1dfiles)
        for iexp in range(nexp):
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfiles[iexp], chk_version=self.par['chk_version'])
            indx = sobjs.name_indices(objids[iexp])
            if not np.any(indx):
                msgs.error("No matching objects for {:s}.  Odds are you input the wrong OBJID".format(objids[iexp]))
            wave_iexp, flux_iexp, ivar_iexp, gpm_iexp, trace_spec, trace_spat, meta_spec, header = \
                    sobjs[indx].unpack_object(ret_flam=self.par['flux_value'], extract_type=self.par['ex_value'])
            weights_sens_iexp = sensfunc.SensFunc.sensfunc_weights(sensfuncfiles[iexp], wave_iexp, debug=self.debug)
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

        Returns:
            waves (list):
               List of arrays with the wavelength arrays for each setup. The length of the list
               equals the number of unique setups and each arrays in the list has shape = (nspec, norders, nexp)
            fluxes (list):
               List of arrays with the flux arrays for each setup. The length of the list
               equals the number of unique setups and each arrays in the list has shape = (nspec, norders, nexp)
            ivars (list):
               List of arrays with the ivar arrays for each setup. The length of the list
               equals the number of unique setups and each arrays in the list has shape = (nspec, norders, nexp)
            gpms (list):
               List of arrays with the gpm arrays for each setup. The length of the list
               equals the number of unique setups and each arrays in the list has shape = (nspec, norders, nexp)
            weights_sens (list):
               List of arrays with the sensfunc weights for each setup. The length of the list
               equals the number of unique setups and each arrays in the list has shape = (nspec, norders, nexp)
            headers (list):
               List of headers for each setup. The length of the list is the number of unique setups.

        """

        waves, fluxes, ivars, gpms, weights_sens, setup_ids, headers = [], [], [], [], [], [], []
        for uniq_setup in self.unique_setups:
            # TODO Is there a more python way to do this?
            setup_indx = np.array(self.setup_id)== uniq_setup
            sensfuncfiles = np.array(self.sensfuncfile)[setup_indx]
            spec1dfiles = np.array(self.spec1dfiles)[setup_indx]
            objids = np.array(self.objids)[setup_indx]
            wave, flux, ivar, gpm, weight_sens, header_out = self.load_ech_arrays(spec1dfiles, objids, sensfuncfiles)
            waves.append(wave)
            fluxes.append(flux)
            ivars.append(ivar)
            gpms.append(gpm)
            weights_sens.append(weight_sens)
            headers.append(header_out)


        return waves, fluxes, ivars, gpms, weights_sens, headers


