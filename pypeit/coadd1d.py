
"""
Coadding module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import os

from IPython import embed

import numpy as np
from astropy.io import fits
from pypeit.spectrographs.util import load_spectrograph
from pypeit import utils
from pypeit import specobjs
from pypeit import msgs
from pypeit.core import coadd


class CoAdd1D(object):

    @classmethod
    def get_instance(cls, spec1dfiles, objids, par=None, sensfile=None, debug=False, show=False):
        """
        Superclass factory method which generates the subclass instance. See __init__ docs for arguments.
        """
        pypeline = fits.getheader(spec1dfiles[0])['PYPELINE'] + 'CoAdd1D'
        return next(c for c in cls.__subclasses__() if c.__name__ == pypeline)(
            spec1dfiles, objids, par=par, sensfile=sensfile, debug=debug, show=show)

    def __init__(self, spec1dfiles, objids, par=None, sensfile=None, debug=False, show=False):
        """

        Args:
            spec1dfiles (list):
               List of strings which are the spec1dfiles
            objids (list):
               List of strings which are the objids for the object in each spec1d file that you want to coadd
            par (parset):
               Pypeit parameter set object
            sensfile (str): optional
               File holding the sensitivity function. This is required for echelle coadds only.
            debug (bool): optional
               Debug. Default = False
            show (bool):
               Debug. Default = True
        """
        # Instantiate attributes
        self.spec1dfiles = spec1dfiles
        self.objids = objids
        self.sensfile = sensfile
        # Load the spectrograph only if par was not passed in to get default parset
        if par is None:
            header = fits.getheader(spec1dfiles[0])
            spectrograph = load_spectrograph(header['PYP_SPEC'])
            self.par = spectrograph.default_pypeit_par()['coadd1d']
        else:
            self.par = par
        self.debug = debug
        self.show = show
        self.nexp = len(self.spec1dfiles) # Number of exposures
        self.coaddfile = None

    def run(self):
        """
        Runs the coadding
        """

        # Load the data
        self.waves, self.fluxes, self.ivars, self.masks, self.header = self.load()
        # Coadd the data
        self.wave_coadd, self.flux_coadd, self.ivar_coadd, self.mask_coadd = self.coadd()

    def load(self):
        """
        Load the arrays we need for performing coadds.

        Returns:
            tuple:
               - waves, fluxes, ivars, masks, header
        """

        for iexp in range(self.nexp):
            sobjs = specobjs.SpecObjs.from_fitsfile(self.spec1dfiles[iexp])
            indx = sobjs.name_indices(self.objids[iexp])
            wave_iexp, flux_iexp, ivar_iexp, mask_iexp, meta_spec, header = \
                sobjs[indx].unpack_object(ret_flam=self.par['flux_value'])
            # Allocate arrays on first iteration
            if iexp == 0:
                waves = np.zeros(wave_iexp.shape + (self.nexp,))
                fluxes = np.zeros_like(waves)
                ivars = np.zeros_like(waves)
                masks = np.zeros_like(waves, dtype=bool)
                header_out = header

            waves[...,iexp], fluxes[...,iexp], ivars[..., iexp], masks[...,iexp] = wave_iexp, flux_iexp, ivar_iexp, mask_iexp

        return waves, fluxes, ivars, masks, header_out

    def save(self, coaddfile, telluric=None, obj_model=None, overwrite=True):
        """
        Routine to save 1d coadds to a fits file. This replaces save.save_coadd1d_to_fits

        Args:
            coaddfile (str):
               File to outuput coadded spectrum to.
            telluric (str):
               This is vestigial and should probably be removed.
            obj_model (str):
               This is vestigial and should probably be removed
            overwrite (bool):
               Overwrite existing file?

        """

        self.coaddfile = coaddfile
        ex_value = self.par['ex_value']
        # Estimate sigma from ivar
        sig = np.sqrt(utils.inverse(self.ivar_coadd))
        if (os.path.exists(self.coaddfile)) and (np.invert(overwrite)):
            hdulist = fits.open(self.coaddfile)
            msgs.info("Reading primary HDU from existing file: {:s}".format(self.coaddfile))
        else:
            msgs.info("Creating a new primary HDU.")
            prihdu = fits.PrimaryHDU()
            if self.header is None:
                msgs.warn('The primary header is none')
            else:
                prihdu.header = self.header
            hdulist = fits.HDUList([prihdu])

        wave_mask = self.wave_coadd > 1.0
        # Add Spectrum Table
        cols = []
        cols += [fits.Column(array=self.wave_coadd[wave_mask], name='{:}_WAVE'.format(ex_value), format='D')]
        cols += [fits.Column(array=self.flux_coadd[wave_mask], name='{:}_FLAM'.format(ex_value), format='D')]
        cols += [fits.Column(array=self.ivar_coadd[wave_mask], name='{:}_FLAM_IVAR'.format(ex_value), format='D')]
        cols += [fits.Column(array=sig[wave_mask], name='{:}_FLAM_SIG'.format(ex_value), format='D')]
        cols += [fits.Column(array=self.mask_coadd[wave_mask].astype(float), name='{:}_MASK'.format(ex_value), format='D')]
        if telluric is not None:
            cols += [fits.Column(array=telluric[wave_mask], name='TELLURIC', format='D')]
        if obj_model is not None:
            cols += [fits.Column(array=obj_model[wave_mask], name='OBJ_MODEL', format='D')]

        coldefs = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = 'OBJ0001-SPEC0001-{:}'.format(ex_value.capitalize())
        hdulist.append(tbhdu)

        if (os.path.exists(self.coaddfile)) and (np.invert(overwrite)):
            hdulist.writeto(self.coaddfile, overwrite=True)
            msgs.info("Appending 1D spectra to existing file {:s}".format(self.coaddfile))
        else:
            hdulist.writeto(self.coaddfile, overwrite=overwrite)
            msgs.info("Wrote 1D spectra to {:s}".format(self.coaddfile))



    def coadd(self):
        """
        Dummy method overloaded by sub-classes

        Returns:

        """
        return (None,)*4


class MultiSlitCoAdd1D(CoAdd1D):
    """
    Child of CoAdd1d for Multislit and Longslit reductions
    """

    def __init__(self, spec1dfiles, objids, par=None, sensfile=None, debug=False, show=False):
        """

            Args:
                spec1dfiles (list):
                   List of strings which are the spec1dfiles
                objids (list):
                   List of strings which are the objids for the object in each spec1d file that you want to coadd
                par (parset):
                   Pypeit parameter set object
                sensfile (str): optional
                   File holding the sensitivity function. This is required for echelle coadds only.
                debug (bool): optional
                   Debug. Default = False
                show (bool):
                   Debug. Default = True
        """

        super().__init__(spec1dfiles, objids, par=par, sensfile=sensfile, debug=debug, show=show)


    def coadd(self):
        """
        Perform coadd for for Multi/Longslit data using multi_combspec

        Returns:
            tuple
              - wave, flux, ivar, mask

        """
        return coadd.multi_combspec(
            self.waves, self.fluxes, self.ivars, self.masks,
            sn_smooth_npix=self.par['sn_smooth_npix'], wave_method=self.par['wave_method'],
            samp_fact=self.par['samp_fact'], ref_percentile=self.par['ref_percentile'],
            maxiter_scale=self.par['maxiter_scale'], sigrej_scale=self.par['sigrej_scale'],
            scale_method=self.par['scale_method'], sn_min_medscale=self.par['sn_min_medscale'],
            sn_min_polyscale=self.par['sn_min_polyscale'], maxiter_reject=self.par['maxiter_reject'],
            lower=self.par['lower'], upper=self.par['upper'], maxrej=self.par['maxrej'], sn_clip=self.par['sn_clip'],
            debug=self.debug, show=self.show)





class EchelleCoAdd1D(CoAdd1D):
    """
    Child of CoAdd1d for Echelle reductions
    """

    def __init__(self, spec1dfiles, objids, par=None, sensfile=None, debug=False, show=False):
        """

            Args:
                spec1dfiles (list):
                   List of strings which are the spec1dfiles
                objids (list):
                   List of strings which are the objids for the object in each spec1d file that you want to coadd
                par (parset):
                   Pypeit parameter set object
                sensfile (str): optional
                   File holding the sensitivity function. This is required for echelle coadds only.
                debug (bool): optional
                   Debug. Default = False
                show (bool):
                   Debug. Default = True
        """

        super().__init__(spec1dfiles, objids, par=par, sensfile=sensfile, debug=debug, show=show)

    def coadd(self):
        """
        Perform coadd for for echelle data using ech_combspec

        Returns:
            tuple
              - wave, flux, ivar, mask

        """
        (wave_coadd, flux_coadd, ivar_coadd, mask_coadd), order_stacks = coadd.ech_combspec(
            self.waves, self.fluxes, self.ivars, self.masks, self.sensfile,
            nbest=self.par['nbest'], sn_smooth_npix=self.par['sn_smooth_npix'], wave_method=self.par['wave_method'],
            samp_fact=self.par['samp_fact'], ref_percentile=self.par['ref_percentile'],
            maxiter_scale=self.par['maxiter_scale'], sigrej_scale=self.par['sigrej_scale'],
            scale_method=self.par['scale_method'], sn_min_medscale=self.par['sn_min_medscale'],
            sn_min_polyscale=self.par['sn_min_polyscale'], maxiter_reject=self.par['maxiter_reject'],
            lower=self.par['lower'], upper=self.par['upper'], maxrej=self.par['maxrej'], sn_clip=self.par['sn_clip'],
            debug = self.debug, show = self.show)

        return wave_coadd, flux_coadd, ivar_coadd, mask_coadd
