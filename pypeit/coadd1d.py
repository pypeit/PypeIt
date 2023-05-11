
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
    def get_instance(cls, spec1dfiles, objids, spectrograph=None, par=None, sensfile=None, debug=False, show=False):
        """
        Superclass factory method which generates the subclass instance. See __init__ docs for arguments.
        """
        pypeline = fits.getheader(spec1dfiles[0])['PYPELINE'] + 'CoAdd1D'
        return next(c for c in cls.__subclasses__() if c.__name__ == pypeline)(
            spec1dfiles, objids, spectrograph=spectrograph, par=par, sensfile=sensfile, debug=debug, show=show)

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfile=None, debug=False, show=False):
        """

        Args:
            spec1dfiles (list):
               List of strings which are the spec1dfiles
            objids (list):
               List of strings which are the objids for the object in each spec1d file that you want to coadd
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`, optional):
            par (:class:`pypeit.par.pypeitpar.Coadd1DPar`, optional):
               Pypeit parameter set object for Coadd1D
            sensfile (str, optional):
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
        self.sensfile = sensfile
        self.debug = debug
        self.show = show
        self.nexp = len(self.spec1dfiles) # Number of exposures
        self.coaddfile = None

    def run(self):
        """
        Runs the coadding
        """

        # Load the data
        self.waves, self.fluxes, self.ivars, self.gpms, self.header = self.load_arrays()
        # Coadd the data
        self.wave_grid_mid, self.wave_coadd, self.flux_coadd, self.ivar_coadd, self.gpm_coadd = self.coadd()
        # Scale to a filter magnitude?
        if self.par['filter'] != 'none':
            scale = flux_calib.scale_in_filter(self.wave_coadd, self.flux_coadd, self.gpm_coadd, self.par)
            self.flux_coadd *= scale
            self.ivar_coadd = self.ivar_coadd / scale**2

    def load_arrays(self):
        """
        Load the arrays we need for performing coadds.

        Returns:
            tuple:
               - waves, fluxes, ivars, gpms, header
        """
        for iexp in range(self.nexp):
            sobjs = specobjs.SpecObjs.from_fitsfile(self.spec1dfiles[iexp], chk_version=self.par['chk_version'])
            indx = sobjs.name_indices(self.objids[iexp])
            if not np.any(indx):
                msgs.error("No matching objects for {:s}.  Odds are you input the wrong OBJID".format(self.objids[iexp]))
            wave_iexp, flux_iexp, ivar_iexp, gpm_iexp, meta_spec, header = \
                    sobjs[indx].unpack_object(ret_flam=self.par['flux_value'], extract_type=self.par['ex_value'])
            # Allocate arrays on first iteration
            # TODO :: We should refactor to use a list of numpy arrays, instead of a 2D numpy array.
            if iexp == 0:
                waves = np.zeros(wave_iexp.shape + (self.nexp,))
                fluxes = np.zeros_like(waves)
                ivars = np.zeros_like(waves)
                gpms = np.zeros_like(waves, dtype=bool)
                header_out = header
                if 'RA' in sobjs[indx][0].keys() and 'DEC' in sobjs[indx][0].keys():
                    header_out['RA_OBJ']  = sobjs[indx][0]['RA']
                    header_out['DEC_OBJ'] = sobjs[indx][0]['DEC']
            # Check if the arrays need to be padded
            # TODO :: Remove the if/elif statement below once these 2D arrays have been converted to a list of 1D arrays
            if wave_iexp.shape[0] > waves.shape[0]:
                padv = [(0, wave_iexp.shape[0]-waves.shape[0]), (0, 0)]
                waves = np.pad(waves, padv, mode='constant', constant_values=(0, 0))
                fluxes = np.pad(fluxes, padv, mode='constant', constant_values=(0, 0))
                ivars = np.pad(ivars, padv, mode='constant', constant_values=(0, 1))
                gpms = np.pad(gpms, padv, mode='constant', constant_values=(False, False))
            elif wave_iexp.shape[0] < waves.shape[0]:
                padv = [0, waves.shape[0]-wave_iexp.shape[0]]
                wave_iexp = np.pad(wave_iexp, padv, mode='constant', constant_values=(0, 0))
                flux_iexp = np.pad(flux_iexp, padv, mode='constant', constant_values=(0, 0))
                ivar_iexp = np.pad(ivar_iexp, padv, mode='constant', constant_values=(0, 1))
                gpm_iexp = np.pad(gpm_iexp, padv, mode='constant', constant_values=(False, False))
            # Store the information
            waves[...,iexp], fluxes[...,iexp], ivars[..., iexp], gpms[...,iexp] \
                = wave_iexp, flux_iexp, ivar_iexp, gpm_iexp
        return waves, fluxes, ivars, gpms, header_out

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

        onespec.head0 = self.header

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

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfile=None, debug=False, show=False):
        """
        See `CoAdd1D` doc string
        """
        super().__init__(spec1dfiles, objids, spectrograph=spectrograph, par = par, sensfile = sensfile,
                         debug = debug, show = show)

    def coadd(self):
        """
        Perform coadd for for Multi/Longslit data using multi_combspec

        Returns:
            tuple
              - wave_grid_mid, wave, flux, ivar, gpm

        """
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





class EchelleCoAdd1D(CoAdd1D):
    """
    Child of CoAdd1d for Echelle reductions
    """

    def __init__(self, spec1dfiles, objids, spectrograph=None, par=None, sensfile=None, debug=False, show=False):
        """
        See `CoAdd1D` doc string

        """
        super().__init__(spec1dfiles, objids, spectrograph=spectrograph, par = par, sensfile = sensfile,
                         debug = debug, show = show)

    def coadd(self):
        """
        Perform coadd for for echelle data using ech_combspec

        Returns:
            tuple
              - wave_grid_mid, wave, flux, ivar, gpm

        """
        weights_sens = sensfunc.SensFunc.sensfunc_weights(self.sensfile, self.waves,
                                                          debug=self.debug)
        wave_grid_mid, (wave_coadd, flux_coadd, ivar_coadd, gpm_coadd), order_stacks \
                = coadd.ech_combspec(self.waves, self.fluxes, self.ivars, self.gpms, weights_sens,
                                     nbest=self.par['nbest'],
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


