
"""
Coadding module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
import os

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit.spectrographs.util import load_spectrograph
from pypeit import specobjs
from pypeit import msgs
from pypeit.core import coadd, flux_calib
from pypeit import datamodel
from pypeit import io


class OneSpec(datamodel.DataContainer):
    """
    DataContainer to hold the products from :class:`pypeit.coadd1d.CoAdd1D`

    See the datamodel for argument descriptions

    Args:
        wave:
        flux:
        PYP_SPEC:

    Attributes:
        head0 (`astropy.io.fits.Header`):  Primary header
        spect_meta (:obj:`dict`): Parsed meta from the header
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Build from PYP_SPEC

    """
    version = '1.0.0'

    datamodel = {'wave': dict(otype=np.ndarray, atype=np.floating, descr='Wavelength array (Ang)'),
                 'flux': dict(otype=np.ndarray, atype=np.floating, descr='Flux array in units of counts/s or 10^-17 erg/s/cm^2/Ang'),
                 'ivar': dict(otype=np.ndarray, atype=np.floating, descr='Inverse variance array (matches units of flux)'),
                 'mask': dict(otype=np.ndarray, atype=np.integer, descr='Mask array (1=Good,0=Bad)'),
                 'telluric': dict(otype=np.ndarray, atype=np.floating, descr='Telluric model'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt: Spectrograph name'),
                 'obj_model': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Object model for tellurics'),
                 'ext_mode': dict(otype=str, descr='Extraction mode (options: BOX, OPT)'),
                 'fluxed': dict(otype=bool, descr='Boolean indicating if the spectrum is fluxed.'),
                 'spect_meta': dict(otype=dict, descr='header dict')}

    @classmethod
    def from_file(cls, ifile):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.from_file`
        to deal with the header

        Args:
            ifile (str):  Filename holding the object

        Returns:
            :class:`OneSpec`:

        """
        hdul = io.fits_open(ifile)
        slf = super(OneSpec, cls).from_hdu(hdul)

        # Internals
        slf.filename = ifile
        slf.head0 = hdul[0].header
        # Meta
        slf.spectrograph = load_spectrograph(slf.PYP_SPEC)
        slf.spect_meta = slf.spectrograph.parse_spec_header(slf.head0)
        #
        return slf


    def __init__(self, wave, flux, PYP_SPEC, ivar=None, mask=None, telluric=None,
                 obj_model=None, ext_mode=None, fluxed=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

    def _init_internals(self):
        self.head0 = None
        self.filename = None
        self.spectrograph = None
        self.spect_meta = None

    def to_file(self, ofile, primary_hdr=None, **kwargs):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.to_file`
        to deal with the header

        Args:
            ofile (:obj:`str`): Filename
            primary_hdr (`astropy.io.fits.Header`_, optional):
            **kwargs:  Passed to super.to_file()

        """
        if primary_hdr is None:
            primary_hdr = io.initialize_header(primary=True)
        # Build the header
        if self.head0 is not None and self.PYP_SPEC is not None:
            spectrograph = load_spectrograph(self.PYP_SPEC)
            subheader = spectrograph.subheader_for_spec(self.head0, self.head0)
        else:
            subheader = {}
        # Add em in
        for key in subheader:
            primary_hdr[key] = subheader[key]
        # Do it
        super(OneSpec, self).to_file(ofile, primary_hdr=primary_hdr,
                                     **kwargs)

class CoAdd1D(object):

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
        self.waves, self.fluxes, self.ivars, self.masks, self.header = self.load_arrays()
        # Coadd the data
        self.wave_coadd, self.flux_coadd, self.ivar_coadd, self.mask_coadd = self.coadd()
        # Scale to a filter magnitude?
        if self.par['filter'] != 'none':
            scale = flux_calib.scale_in_filter(self.wave_coadd, self.flux_coadd, self.mask_coadd, self.par)
            self.flux_coadd *= scale
            self.ivar_coadd = self.ivar_coadd / scale**2

    def load_arrays(self):
        """
        Load the arrays we need for performing coadds.

        Returns:
            tuple:
               - waves, fluxes, ivars, masks, header
        """

        for iexp in range(self.nexp):
            sobjs = specobjs.SpecObjs.from_fitsfile(self.spec1dfiles[iexp])
            indx = sobjs.name_indices(self.objids[iexp])
            if not np.any(indx):
                msgs.error("No matching objects for {:s}.  Odds are you input the wrong OBJID".format(self.objids[iexp]))
            wave_iexp, flux_iexp, ivar_iexp, mask_iexp, meta_spec, header = \
                    sobjs[indx].unpack_object(ret_flam=self.par['flux_value'], extract_type=self.par['ex_value'])
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
        wave_mask = self.wave_coadd > 1.0
        # Generate the DataContainer
        onespec = OneSpec(self.wave_coadd[wave_mask],
                          self.flux_coadd[wave_mask],
                          PYP_SPEC=self.spectrograph.name,
                          ivar=self.ivar_coadd[wave_mask],
                          mask=self.mask_coadd[wave_mask].astype(int),
                          ext_mode=self.par['ex_value'],
                          fluxed=self.par['flux_value'])
        onespec.head0 = fits.getheader(self.spec1dfiles[0])

        # Add on others
        if telluric is not None:
            onespec.telluric  = telluric[wave_mask]
        if obj_model is not None:
            onespec.obj_model = obj_model[wave_mask]
        # Write
        onespec.to_file(coaddfile, overwrite=overwrite)

    def coadd(self):
        """
        Dummy method overloaded by sub-classes

        Returns:
            :obj:`tuple`:  four items
              - wave
              - flux
              - ivar
              - mask

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
            debug = self.debug, show = self.show, extrap_sens=self.par['extrap_sens'])

        return wave_coadd, flux_coadd, ivar_coadd, mask_coadd
