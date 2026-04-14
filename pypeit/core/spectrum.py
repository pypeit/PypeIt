"""
Temporary light-weight spectrum object.

Ideally this would be replaced by specutils.Spectrum

"""

from copy import deepcopy

from IPython import embed
import numpy as np

from pypeit import log
from pypeit import PypeItError
from pypeit import sampling
from pypeit import utils


class Spectrum:
    r"""
    A light-weight container class for a spectrum.

    The flux array can be either 1D or 2D.  If 2D, its shape is always assumed
    to be :math:`(N_{\rm pix},N_{\rm spec})`, where :math:`N_{\rm pix}` is the
    number of pixels per spectrum and :math:`N_{\rm spec}` is the number of
    spectra.  The wavelength array must always be provided as 1D; i.e., for a 2D
    flux array, the wavelength vector must be correct for each flux vector along
    the 1st axis of the 2D array.
    
    Parameters
    ----------
    wave : array-like
        Vacuum wavelengths in angstrom.  Must be 1D, and its length must match
        the first axis of ``flux``.
    flux : array-like
        Flux data at each wavelength.  Can be 2D or 1D; see class description.
    ivar : array-like, optional
        Inverse variance in the flux.  Shape must match ``flux``.  If None,
        assumes uncertainties are unknown.
    gpm : array-like, optional
        Boolean good-pixel mask.  Shape must match ``flux``.  If None, assumes
        all pixels are valid.
    meta : dict, optional
        Collection of relevant metadata.
    """
    def __init__(self, wave, flux, ivar=None, gpm=None, meta=None):

        # Throughout I use the copy method to ensure the original arrays are
        # copied
        self.flux = np.asarray(flux, dtype=float).copy()

        self.wave = np.asarray(wave, dtype=float).copy()
        if self.wave.ndim != 1:
            raise PypeItError('wavelength array must always be 1D in the spectrum object')
        if self.wave.size != self.flux.shape[0]:
            raise PypeItError('wavelength vector must match length of flux array')

        if ivar is None:
            self.ivar = None
        else:
            self.ivar = np.asarray(ivar, dtype=float).copy()
            if self.ivar.shape != self.flux.shape:
                raise PypeItError('Wavelength and inverse variance arrays do not have the same shape.')

        if gpm is None:
            self.gpm = np.ones(self.flux.shape, dtype=bool)
        else:
            self.gpm = np.asarray(gpm, dtype=bool).copy()
            if self.gpm.shape != self.flux.shape:
                raise PypeItError('Wavelength and good-pixel arrays do not have the same size.')

        self.meta = meta if meta is None else deepcopy(meta)

    @property
    def size(self):
        """
        The size of the flux array
        """
        return self.flux.size
    
    @property
    def shape(self):
        """
        The shape of the flux array
        """
        return self.flux.shape
    
    @property
    def ndim(self):
        """
        The dimensionality of the flux array
        """
        return self.flux.ndim
    
    def copy(self):
        """
        Return a deepcopy of the object.
        """
        _ivar = None if self.ivar is None else self.ivar.copy()
        _meta = None if self.meta is None else deepcopy(self.meta)
        return self.__class__(
            self.wave.copy(), self.flux.copy(), ivar=_ivar, gpm=self.gpm.copy(), meta=_meta
        )

    def multiply(self, a):
        """
        Multiply the spectrum by a scalar, vector, or another spectrum.

        *This modifies the spectral data in place.*  If uncertainties are
        available, they are propagated.  Any divisions by 0 result in an inverse
        variance of 0 and the good pixel mask is set to False.

        Parameters
        ----------
        a : scalar, array-like, :class:`pypeit.core.spectrum.Spectrum`
            Multiplicative factor.  If an array, its shape must match
            :attr:`flux`.  If a spectrum, the wavelength arrays of the two
            spectrum *must be identical*.
        """
        # Multiply by a scalar
        if isinstance(a, (int, np.integer, float, np.floating)):
            if float(a) == 0.:
                log.warning('Multiplicative factor is 0!')
            self.flux *= a
            if self.ivar is not None:
                if np.absolute(a) > 0:
                    self.ivar /= a**2
                else:
                    self.ivar *= 0.
            return

        # Workspace
        sqr_err_ratio = None
        a_gpm = None

        if isinstance(a, Spectrum):
            # Pull the necessary data out of the spectrum

            # Check the wavelength vectors
            # TODO: Loosen this; i.e., use isclose instead of array_equal?
            if not np.array_equal(a.wave, self.wave):
                raise PypeItError('To multiply two spectra, their wavelength vectors must be identical.')
            a_flux = a.flux
            a_gpm = a.gpm
            if a.ivar is not None:
                # Square of the ratio between the error and flux in a
                sqr_err_ratio = utils.inverse(a.flux**2 * a.ivar)
        else:
            # Convert the array-like object to a numpy array
            a_flux = np.asarray(a)

        # Check the input
        if a_flux.ndim > self.ndim:
            raise PypeItError(
                'Multiplication does not allow the dimensionality of the spectrum to change.  '
                f'The dimensionality of this spectrum is {self.ndim} and the multiplier is '
                f'{a_flux.ndim}.'
            )
        # Numpy broadcasting rules mean that the arithmetic operations performed
        # below should work, as long as the last a.ndim dimensions of a and this
        # spectrum match.
        if a_flux.shape != self.shape[:a_flux.ndim]:
            raise PypeItError(
                'Numpy will not be able to successfully broadcast arithmetic operations between '
                f'this spectrum, shape={self.shape}, and the multiplier, shape={a_flux.shape}.'
            )

        # Reshape the arrays, if necessary
        if a_flux.ndim != self.flux.ndim:
            a_flux = np.expand_dims(
                a_flux, tuple(np.arange(len(self.shape))[a_flux.ndim:].tolist())
            )
            if sqr_err_ratio is not None:
                sqr_err_ratio = np.expand_dims(
                    sqr_err_ratio, tuple(np.arange(len(self.shape))[sqr_err_ratio.ndim:].tolist())
                )
            if a_gpm is not None:
                a_gpm = np.expand_dims(
                    a_gpm, tuple(np.arange(len(self.shape))[a_gpm.ndim:].tolist())
                )

        # Add the error.  NOTE: This *must* be done before the multiplication by
        # a_flux below because that changes self.flux.
        if self.ivar is not None:
            # Square of the ratio between the error and flux in this spectrum
            if sqr_err_ratio is None:
                sqr_err_ratio = utils.inverse(self.flux**2 * self.ivar)
            else:
                sqr_err_ratio += utils.inverse(self.flux**2 * self.ivar)

        # Complete the multiplication
        self.flux *= a_flux
        if sqr_err_ratio is not None:
            # Propagate the error
            sqr_err = self.flux**2 * sqr_err_ratio
            self.ivar = utils.inverse(sqr_err)
            self.gpm[np.logical_not(self.ivar > 0)] = False
        if a_gpm is not None:
            # Propagate the good-pixel mask
            self.gpm &= a_gpm

    def inverse(self):
        """
        Replace the spectrum with its multiplicative inverse.

        *This modifies the spectrum in place.*  If uncertainties are available,
        they are propagated.  Any divisions by 0 result in an inverse variance
        of 0 and the good pixel mask is set to False.
        """
        if self.ivar is not None:
            self.ivar *= self.flux**4
            self.gpm[np.logical_not(self.ivar > 0)] = False
        self.flux = utils.inverse(self.flux)

    def to_magnitude(self, zeropoint=0.):
        r"""
        Convert the spectrum to magnitudes.

        For fluxes, :math:`f`, this returns

        .. math::

            m = -2.5 \log_{\rm 10} (f) + Z,

        where :math:`Z` is the provided zeropoint.

        *This modifies the spectrum in place.*  If uncertainties are available,
        they are propagated.  Any pixels with non-positive fluxes are masked.

        Parameters
        ----------
        zeropoint : float, optional
            The magnitude conversion zeropoint (see above)
        """
        if self.ivar is not None:
            self.ivar *= (self.flux * np.log(10) / 2.5)**2
        self.gpm[np.logical_not(self.flux > 0)] = False
        self.flux[np.logical_not(self.gpm)] = 0.
        self.flux[self.gpm] = -2.5 * np.log10(self.flux[self.gpm]) + zeropoint

    def resample(self, new_wave, pixel_fraction_threshold=0.8, conserve=False):
        r"""
        Resample the spectrum to a new wavelength array.

        If available, errors and masking are both propagated through the
        calculation.  This function is basically a wrapper for
        :class:`~pypeit.sampling.Resample`.

        Parameters
        ----------
        new_wave : array-like
            New wavelength vector for the spectrum
        pixel_fraction_threshold : float, optional
            The resampling calculates the fraction of each output pixel that has
            unmasked contributions from the original spectrum.  Fractions below
            this threshold will be masked in the output spectrum.
        conserve : bool, optional
            Conserve the flux in the resampled spectrum.  If the units of the
            spectrum are flux integrated over the pixel, this should typically
            be True; if the units are flux density (e.g., :math:`{\rm
            ergs/s/cm}^2{\rm /angstrom}`), this should typically be False.

        Returns
        -------
        :class:`~pypeit.core.spectrum.Spectrum`
            A new spectrum object with the resample data.
        """
        # Setup
        bpm = None if np.all(self.gpm) else np.logical_not(self.gpm).T
        if self.ivar is None:
            err = None
        else:
            err = np.zeros(self.ivar.shape, dtype=float)
            err[self.gpm] = np.sqrt(utils.inverse(self.ivar[self.gpm]))
            err = err.T

        # Resample accepts both 1D and 2D arrays, but it expects the 2D arrays
        # to have the spectra organized along the 2nd axis; i.e.,
        # (N_spec,N_pix) instead of (N_pix,N_spec).
        r = sampling.Resample(
            self.flux.T, e=err, mask=bpm, x=self.wave, newx=new_wave, conserve=conserve
        )
        ivar = None if err is None else utils.inverse(r.oute.T)**2
        return Spectrum(
            r.outx, r.outy.T, ivar=ivar, gpm=r.outf.T > pixel_fraction_threshold, meta=self.meta,
        )
