"""
Provides a simple datamodel for a single spectrum.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect

import astropy
from IPython import embed
import numpy as np
from scipy.interpolate import interp1d

from pypeit import datamodel
from pypeit import io
from pypeit import log
from pypeit import PypeItError
from pypeit import utils
from pypeit.spectrographs.util import load_spectrograph


class OneSpec(datamodel.DataContainer):
    """
    DataContainer to hold single spectra, e.g., from
    :class:`~pypeit.coadd1d.CoAdd1D`.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_onespec.rst

    Args:
        wave
        wave_grid_mid
        flux
        PYP_SPEC
        ivar
        mask
        telluric
        obj_model
        ext_mode
        fluxed

    Attributes:
        head0 (`astropy.io.fits.Header`):
            Primary header
        spect_meta (:obj:`dict`):
            Parsed meta from the header
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Build from PYP_SPEC

    """
    version = '1.0.2'

    datamodel = {'wave': dict(otype=np.ndarray, atype=np.floating,
                              # TODO: The "weighted by pixel contributions" part
                              # should be better explained.
                              descr='Wavelength array (angstroms in vacuum), weighted by pixel '
                                    'contributions'),
                 'wave_grid_mid': dict(otype=np.ndarray, atype=np.floating,
                                       descr='Wavelength (angstroms in vacuum) evaluated at the '
                                             'bin centers of a grid that is uniformly spaced '
                                             'in either lambda or log10-lambda/velocity'),
                 'flux': dict(otype=np.ndarray, atype=np.floating,
                              descr='Flux array in units of counts/s or 10^-17 erg/s/cm^2/Ang; '
                                    'see ``fluxed``'),
                 'ivar': dict(otype=np.ndarray, atype=np.floating,
                              descr='Inverse variance array (matches units of flux)'),
                 'sigma': dict(otype=np.ndarray, atype=np.floating,
                              descr='One sigma noise array, equivalent to 1/sqrt(ivar) (matches units of flux)'),
                 'mask': dict(otype=np.ndarray, atype=np.integer,
                              descr='Mask array (1=Good,0=Bad)'),
                 'telluric': dict(otype=np.ndarray, atype=np.floating, descr='Telluric model'),
                 'PYP_SPEC': dict(otype=str, descr='``PypeIt`` spectrograph designation'),
                 'obj_model': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Object model for tellurics'),
                 'ext_mode': dict(otype=str, descr='Extraction mode (options: BOX, OPT)'),
                 'fluxed': dict(otype=bool, descr='Boolean indicating if the spectrum is fluxed.'),
                 # TODO: Needs a better description.  What's in the dictionary?
                 # Why isn't this dictionary expanded into its elements?  I.e.,
                 # shouldn't each element of this dictionary be a component of
                 # the datamodel?
                 'spect_meta': dict(otype=dict, descr='header dict')}

    internals = ['head0',
                 'filename',
                 'spectrograph',
                 'spect_meta',
                 'history']

    @classmethod
    def from_file(cls, ifile, verbose=True, chk_version=True, **kwargs):
        """
        Instantiate the object from an extension in the specified fits file.

        Over-load :func:`~pypeit.datamodel.DataContainer.from_file`
        to deal with the header
        
        Args:
            ifile (:obj:`str`, `Path`_):
                Fits file with the data to read
            verbose (:obj:`bool`, optional):
                Print informational messages (not currently used)
            chk_version (:obj:`bool`, optional):
                Passed to :func:`from_hdu`.
            kwargs (:obj:`dict`, optional):
                Arguments passed directly to :func:`from_hdu`.
        """
        with io.fits_open(ifile) as hdu:
            self = cls.from_hdu(hdu, chk_version=chk_version, **kwargs)
            self.filename = ifile
            self.head0 = hdu[0].header
            self.spectrograph = load_spectrograph(self.PYP_SPEC)
            self.spect_meta = self.spectrograph.parse_spec_header(self.head0)
        return self

    @classmethod
    def from_xspec_file(cls, ifile:str):
        hdul = io.fits_open(ifile)
        wave = hdul['WAVELENGTH'].data
        flux = hdul['FLUX'].data
        return cls(wave, None, flux, fluxed=False)

    @property
    def npix(self):
        """
        Number of pixels in the spectrum.

        Returns:
            :obj:`int`: Number of pixels in the spectrum.
        """
        return self.wave.size
    
    def __init__(self, wave, wave_grid_mid, flux, PYP_SPEC=None, ivar=None, sigma=None, mask=None, telluric=None,
                 obj_model=None, ext_mode=None, fluxed=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

    def _bundle(self):
        """
        Override the base class method simply to set the HDU extension name.
        """
        return super()._bundle(ext='SPECTRUM')

    def to_file(self, ofile, primary_hdr=None, history=None, **kwargs):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.to_file`
        to deal with the header

        Args:
            ofile (:obj:`str`): Filename
            primary_hdr (`astropy.io.fits.Header`_, optional):
            **kwargs:  Passed to super.to_file()

        """
        if primary_hdr is None:
            primary_hdr = io.initialize_header()
        # Build the header
        if self.head0 is not None and self.PYP_SPEC is not None:
            spectrograph = load_spectrograph(self.PYP_SPEC)
            subheader = spectrograph.subheader_for_spec(self.head0, self.head0,
                                                        extra_header_cards = ['RA_OBJ', 'DEC_OBJ'])
        else:
            subheader = {}
        # Add em in
        for key in subheader:
            primary_hdr[key] = subheader[key]

        # Add history
        if history is not None:
            history.write_to_header(primary_hdr)

        # Do it
        super().to_file(ofile, primary_hdr=primary_hdr, **kwargs)


    def rebin(self, new_wv, fill_value=0., grow_bad_sig=False):
        """ Rebin the spectrum to a new OneSpec object with the input array

        Uses simple linear interpolation.  The default (and only)
        option conserves counts (and flambda).

        WARNING: Do not trust either edge pixel of the new array.
        In fact the sig is set to 0 for each of these
        Also be aware that neighboring pixels are likely to be
        correlated in a manner that is not described by the error
        array.

        Parameters
        ----------
        new_wv : np.ndarray
            New wavelength array
        fill_value : float, optional
            Fill value at the edges
            Default = 0., but 'extrapolate' may be considered
        grow_bad_sig : bool, optional
            Allow sig<=0. values and grow them

        Returns
        -------
        newspec : OneSpec
        """
        # Save flux info to avoid unit issues
        flux = self.flux

        # Deal with nan
        badf = np.any([np.isnan(flux), np.isinf(flux)], axis=0)
        if np.sum(badf) > 0:
            log.warning("Ignoring pixels with NAN or INF in flux")
        gdf = ~badf
        flux = flux[gdf]

        # Check for bad pixels (not prepared for these)
        if self.sigma is not None:
            sig = self.sigma
            bad_sig = sig[gdf] <= 0.
            if np.sum(bad_sig) > 0:
                if not grow_bad_sig:
                    raise PypeItError(
                        'Data contains rejected pixels (sig=0). Use grow_bad_sig to proceed and '
                        'grow them.'
                    )
            bads = np.any([np.isnan(sig[gdf]), np.isinf(sig[gdf]**2)], axis=0)  # Latter is for way too large values
            bad_sig[bads] = True

        # Endpoints of original pixels
        npix = self.wave.size
        wvh = (self.wave + np.roll(self.wave, -1)) / 2.
        wvh[npix - 1] = self.wave[npix - 1] + \
                        (self.wave[npix - 1] - self.wave[npix - 2]) / 2.
        dwv = wvh - np.roll(wvh, 1)
        dwv[0] = 2 * (wvh[0] - self.wave[0])
        med_dwv = np.median(dwv)

        wvh = wvh[gdf]
        dwv = dwv[gdf]

        # Error
        if self.sigma is not None: 
            var = sig[gdf]**2
            var[bad_sig] = 0.
        else:
            var = np.ones_like(flux)

        # Cumulative Sum
        cumsum = np.cumsum(flux * dwv)
        cumvar = np.cumsum(var * dwv, dtype=np.float64)

        # Interpolate (loses the units)
        fcum = interp1d(wvh, cumsum, fill_value=fill_value, bounds_error=False)
        fvar = interp1d(wvh, cumvar, fill_value=0., bounds_error=False)

        # Endpoints of new pixels
        nnew = len(new_wv)
        nwvh = (new_wv + np.roll(new_wv, -1)) / 2.
        nwvh[nnew - 1] = new_wv[nnew - 1] + \
                        (new_wv[nnew - 1] - new_wv[nnew - 2]) / 2.
        # Pad starting point
        bwv = np.zeros(nnew + 1) 
        bwv[0] = new_wv[0] - (new_wv[1] - new_wv[0]) / 2.
        bwv[1:] = nwvh

        # Evaluate and put unit back
        newcum = fcum(bwv) 
        newvar = fvar(bwv) 

        # Rebinned flux, var, co
        new_fx = (np.roll(newcum, -1) - newcum)[:-1]
        new_var = (np.roll(newvar, -1) - newvar)[:-1]

        # Normalize (preserve counts and flambda)
        new_dwv = bwv - np.roll(bwv, 1)
        new_fx = new_fx / new_dwv[1:]
        # Preserve S/N (crudely)
        med_newdwv = np.median(new_dwv)
        new_var = new_var / (med_newdwv/med_dwv) / new_dwv[1:]

        # Return new spectrum
        if self.sigma is not None:
            # Create new_sig
            new_sig = np.zeros_like(new_var)
            gd = new_var > 0.
            new_sig[gd] = np.sqrt(new_var[gd])
            # Deal with bad pixels (grow_bad_sig should be True)
            bad = np.where(var <= 0.)[0]
            # Find nearby wavelengths in rebinned wavelength
            nearidxs = np.searchsorted(new_wv, self.wave[bad])
            # Pad arrays to enable vector operations
            pwv = np.concatenate([self.wave,[self.wave[-1]+dwv[-1]]])
            pndwv = np.concatenate([new_dwv,[new_dwv[-1]]])
            pnwv = np.concatenate([new_wv,[new_wv[-1]+new_dwv[-1]]])

            # Find distances between original bad wavelengths and nearby new ones
            ldiff = np.abs(new_wv[nearidxs-1]-pwv[bad]) - \
                    (pndwv[1:][nearidxs]+dwv[bad])/2
            rdiff = np.abs(pwv[bad]-pnwv[nearidxs]) - \
                    (pndwv[1:][nearidxs] + dwv[bad]) / 2
            # Set errors to 0; we have to mind the padding above
            new_sig[nearidxs[(ldiff<0)&(nearidxs<len(new_wv))]] = 0
            new_sig[nearidxs[(rdiff<0)&(nearidxs<len(new_wv))]] = 0

            # Zero out edge pixels -- not to be trusted
            igd = np.where(gd)[0]
            if len(igd) == 0:  # Should not get here!
                raise PypeItError('Not a single good pixel?!  Something went wrong...')
            new_sig[igd[0]] = 0.
            new_sig[igd[-1]] = 0.
        else:
            new_sig = None

        # Finish
        return OneSpec(new_wv, None, new_fx, sigma=new_sig)

    def gauss_smooth(self, fwhm):
        """
        Smooth the spectrum with a Gaussian kernal.

        .. warning::

            This does **not** propagate the errors!  The returned
            :class:`OneSpec` instance will *not* have any errors.  All other
            attributes are included in the returned
            :class:`~pypeit.onespec.OneSpec` instance; array attributes are
            copied.

        Parameters
        ----------
        fwhm : float
            FWHM of the Gaussian kernal in pixels.

        Returns
        -------
        :class:`~pypeit.onespec.OneSpec`
            The smoothed spectrum
        """
        return OneSpec(
            self.wave.copy(), None if self.wave_grid_mid is None else self.wave_grid_mid.copy(), 
            utils.convolve_psf(self.flux, fwhm), PYP_SPEC=self.PYP_SPEC, ivar=None, sigma=None,
            mask=None if self.mask is None else self.mask.copy(),
            telluric=None if self.telluric is None else self.telluric.copy(),
            obj_model=None if self.obj_model is None else self.obj_model.copy(),
            ext_mode=self.ext_mode, fluxed=self.fluxed
        )

