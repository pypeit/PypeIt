"""
Implements the objects used to construct sensitivity functions.

.. include:: ../include/links.rst
"""
import inspect

from IPython import embed

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import fits
from astropy import table

from pypeit import msgs
from pypeit import specobjs
from pypeit import utils
from pypeit.core import coadd
from pypeit.core import flux_calib
from pypeit.core import telluric
from pypeit.core import fitting
from pypeit.core.wavecal import wvutils
from pypeit.core import meta
from pypeit.spectrographs.util import load_spectrograph

from pypeit import datamodel

# TODO Add the data model up here as a standard thing using DataContainer.

# TODO Standard output location for sensfunc?

# TODO Add some QA plots, and plots to the screen if show is set.

class SensFunc(datamodel.DataContainer):
    r"""
    Base class for generating sensitivity functions from a standard-star
    spectrum.

    This class should not be instantated by itself; instead instantiate
    either :class:`UVISSensFunc` or :class:`IRSensFunc`, depending on the
    wavelength range of your data (UVIS for :math:`\lambda < 7000` angstrom,
    IR for :math:`\lambda > 7000` angstrom.)

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_sensfunc.rst

    Args:
        spec1dfile (:obj:`str`):
            PypeIt spec1d file for the standard file.
        sensfile (:obj:`str`):
            File name for the sensitivity function data.
        par (:class:`~pypeit.par.pypeitpar.SensFuncPar`, optional):
            The parameters required for the sensitivity function computation.
        debug (:obj:`bool`, optional):
            Run in debug mode, sending diagnostic information to the screen.
    """
    version = '1.0.1'
    """Datamodel version."""

    # TODO: Add this if we want to set the output float type for the np.ndarray
    # elements of the datamodel...
#    output_float_dtype = np.float32
#    """Regardless of datamodel, output floating-point data have this fixed bit size."""

    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'pypeline': dict(otype=str, descr='PypeIt pipeline reduction path'),
                 'spec1df': dict(otype=str,
                                 descr='PypeIt spec1D file used to for sensitivity function'),
                 'std_name': dict(otype=str, descr='Type of standard source'),
                 'std_cal': dict(otype=str,
                                 descr='File name (or shorthand) with the standard flux data'),
                 # TODO: Is it possible/useful to force the coordinates to always be floats
                 'std_ra': dict(otype=float, descr='RA of the standard source'),
                 'std_dec': dict(otype=float, descr='DEC of the standard source'),
                 'airmass': dict(otype=float, descr='Airmass of the observation'),
                 'exptime': dict(otype=float, descr='Exposure time'),
                 'telluric': dict(otype=telluric.Telluric,
                                  descr='Telluric model; see '
                                        ':class:`~pypeit.core.telluric.Telluric`'),
                 'sens': dict(otype=table.Table, descr='Table with the sensitivity function'),
                 'wave': dict(otype=np.ndarray, atype=float, descr='Wavelength vectors'),
                 'zeropoint': dict(otype=np.ndarray, atype=float,
                                   descr='Sensitivity function zeropoints'),
                 'throughput': dict(otype=np.ndarray, atype=float,
                                    descr='Spectrograph throughput measurements'),
                 'algorithm': dict(otype=str, descr='Algorithm used for the sensitivity calculation.')}
#                                    ,
#                 'wave_splice': dict(otype=np.ndarray, atype=float,
#                                     descr='Spliced-together wavelength vector'),
#                 'zeropoint_splice': dict(otype=np.ndarray, atype=float,
#                                          descr='Spliced-together sensitivity function zeropoints'),
#                 'throughput_splice': dict(otype=np.ndarray, atype=float,
#                                           descr='Spliced-together spectrograph throughput '
#                                                 'measurements')}
    """DataContainer datamodel."""

    internals = ['sensfile',
                 'spectrograph',
                 'par',
                 'qafile',
                 'thrufile',
                 'debug',
                 'wave_cnts',
                 'counts',
                 'counts_ivar',
                 'counts_mask',
                 'nspec_in',
                 'norderdet',
                 'wave_splice',
                 'zeropoint_splice',
                 'throughput_splice',
                 'steps',
                 'splice_multi_det',
                 'meta_spec',
                 'std_dict'
                ]

    _algorithm = None
    """Algorithm used for the sensitivity calculation."""

    @staticmethod
    def empty_sensfunc_table(norders, nspec, ncoeff=1):
        """
        Construct an empty `astropy.table.Table`_ for the sensitivity
        function.

        Args:
            norders (:obj:`int`):
                The number of slits/orders on the detector.
            nspec (:obj:`int`):
                The number of spectral pixels on the detector.

        Returns:
            `astropy.table.Table`_: Instance of the empty sensitivity
            function table.
        """
        return table.Table(data=[
            table.Column(name='SENS_WAVE', dtype=float, length=norders, shape=(nspec,),
                         description='Wavelength vector'),
            table.Column(name='SENS_COUNTS_PER_ANG', dtype=float, length=norders, shape=(nspec,),
                         description='Flux in counts per angstrom'),
            table.Column(name='SENS_ZEROPOINT', dtype=float, length=norders, shape=(nspec,),
                         description='Measured sensitivity zero-point data'),
            table.Column(name='SENS_ZEROPOINT_GPM', dtype=bool, length=norders, shape=(nspec,),
                         description='Good-pixel mask for the measured zero points'),
            table.Column(name='SENS_ZEROPOINT_FIT', dtype=float, length=norders, shape=(nspec,),
                         description='Best-fit smooth model to the zero points'),
            table.Column(name='SENS_ZEROPOINT_FIT_GPM', dtype=bool, length=norders, shape=(nspec,),
                         description='Good-pixel mask for the model zero points'),
            table.Column(name='SENS_COEFF', dtype=float, length=norders, shape=(ncoeff,),
                         description='Coefficients of smooth model fit to zero points'),
            table.Column(name='ECH_ORDERS', dtype=int, length=norders,
                         description='Echelle order for this specrum (echelle data only)'),
            table.Column(name='POLYORDER_VEC', dtype=int, length=norders,
                         description='Polynomial order for each slit/echelle (if applicable)'),
            table.Column(name='WAVE_MIN', dtype=float, length=norders,
                         description='Minimum wavelength included in the fit'),
            table.Column(name='WAVE_MAX', dtype=float, length=norders,
                         description='Maximum wavelength included in the fit')])

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec1dfile, sensfile, par, debug=False):
        """
        Instantiate the relevant subclass based on the algorithm provided in
        ``par``.
        """
        return next(c for c in cls.__subclasses__()
                    if c.__name__ == f"{par['algorithm']}SensFunc")(spec1dfile, sensfile, par=par,
                                                                    debug=debug)

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):

        # Instantiate as an empty DataContainer
        super().__init__()

        # Input and Output files
        self.spec1df = spec1dfile
        self.sensfile = sensfile

        # Spectrograph
        header = fits.getheader(self.spec1df)
        self.PYP_SPEC = header['PYP_SPEC']
        self.spectrograph = load_spectrograph(self.PYP_SPEC)
        self.pypeline = self.spectrograph.pypeline
        # TODO: This line is necessary until we figure out a way to instantiate
        # spectrograph objects with configuration specific information from
        # spec1d files.
        self.spectrograph.dispname = header['DISPNAME']

        # Get the algorithm parameters
        self.par = self.spectrograph.default_pypeit_par()['sensfunc'] if par is None else par
        # TODO: Check the type of the parameter object?

        # Set the algorithm in the datamodel
        self.algorithm = self.__class__._algorithm

        # QA and throughput plot filenames
        self.qafile = sensfile.replace('.fits', '') + '_QA.pdf'
        self.thrufile = sensfile.replace('.fits', '') + '_throughput.pdf'

        # Other
        self.debug = debug
        self.steps = []

        # Are we splicing together multiple detectors?
        self.splice_multi_det = True if self.par['multi_spec_det'] is not None else False

        # Read in the Standard star data
        sobjs_std = specobjs.SpecObjs.from_fitsfile(self.spec1df).get_std(
                            multi_spec_det=self.par['multi_spec_det'])

        if sobjs_std is None:
            msgs.error('There is a problem with your standard star spec1d file: {:s}'.format(self.spec1df))
        # Unpack standard
        wave, counts, counts_ivar, counts_mask, self.meta_spec, header = sobjs_std.unpack_object(ret_flam=False)

        # Perform any instrument tweaks
        wave_twk, counts_twk, counts_ivar_twk, counts_mask_twk \
                = self.spectrograph.tweak_standard(wave, counts, counts_ivar, counts_mask,
                                                   self.meta_spec)

        # Reshape to 2d arrays
        self.wave_cnts, self.counts, self.counts_ivar, self.counts_mask, self.nspec_in, \
            self.norderdet \
                = utils.spec_atleast_2d(wave_twk, counts_twk, counts_ivar_twk, counts_mask_twk)

        # If the user provided RA and DEC use those instead of what is in meta
        star_ra = self.meta_spec['RA'] if self.par['star_ra'] is None else self.par['star_ra']
        star_dec = self.meta_spec['DEC'] if self.par['star_dec'] is None else self.par['star_dec']
        # Convert to decimal deg, as needed
        star_ra, star_dec = meta.convert_radec(star_ra, star_dec)

        # Read in standard star dictionary
        self.std_dict = flux_calib.get_standard_spectrum(star_type=self.par['star_type'],
                                                         star_mag=self.par['star_mag'],
                                                         ra=star_ra, dec=star_dec)

    def _bundle(self):
        """
        Bundle the object for writing using
        :func:`~pypeit.datamodel.DataContainer.to_hdu`.
        """
        # All of the datamodel elements that can be written as a header are
        # written to the SENS extension.
        d = []
        for key in self.keys():
            if self[key] is None:
                continue
            if isinstance(self[key], datamodel.DataContainer) \
                    or isinstance(self[key], table.Table):
                d += [{key: self[key]}]
                continue
            if self.datamodel[key]['otype'] == np.ndarray:
                # TODO: We might want a general solution in
                # DataContainer that knows to convert (and revert)
                # boolean arrays to np.uint8 types. Does anyone know a
                # way around this? It's annoying that ImageHDUs cannot
                # have boolean type arrays.
                if issubclass(self[key].dtype.type, (bool, np.bool_)):
                    d += [{key: self[key].astype(np.uint8)}]
#                elif issubclass(self[key].dtype.type, (float, np.floating)):
#                    d += [{key: self[key].astype(self.output_float_dtype)}]
                else:
                    d += [{key: self[key]}]

        # If the spliced arrays have been set, use these to replace the
        # multislit/echelle wave, zeropoint, and throughput arrays
        if self.splice_multi_det:
            # TODO: I added this neurotic check, just to make sure...
            if self.wave_splice is None or self.zeropoint_splice is None \
                    or self.throughput_splice is None:
                msgs.error('CODING ERROR: Assumed if splice_multi_det is True, then the *_splice '
                           'arrays have all been defined.  Found a case where this is not true!')
            # Loop through this list of dictionaries
            for _d in d:
                if list(_d.keys())[0] not in ['wave', 'zeropoint', 'throughput']:
                    # Keep going if it's not one of the relevant attributes
                    continue
                # Replace with the spliced version.  Warning: This requires that
                # the relevant attribute names are paired as self.attr and
                # self.attr_splice
                attr = list(_d.keys())[0]
                _d[attr] = getattr(self, f'{attr}_splice')

        # Add all the non-array, non-DataContainer, non-Table elements to the
        # 'sens' extension. This is done after the first loop to just to make
        # sure that the order of extensions matches the order in the datamodel.
        for _d in d:
            if list(_d.keys()) != ['sens']:
                continue
            for key in self.keys():
                if self[key] is None or self.datamodel[key]['otype'] == np.ndarray \
                        or isinstance(self[key], datamodel.DataContainer) \
                        or isinstance(self[key], table.Table):
                    continue
                _d[key] = self[key]

        return d

    @classmethod
    def from_hdu(cls, hdu, hdu_prefix=None, chk_version=True):
        """
        Instantiate the object from an HDU extension.

        This overrides the base-class method, essentially just to handle the
        fact that the 'TELLURIC' extension is not called 'MODEL'.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            hdu_prefix (:obj:`str`, optional):
                Maintained for consistency with the base class but is
                not used by this method.
            chk_version (:obj:`bool`, optional):
                If True, raise an error if the datamodel version or
                type check failed. If False, throw a warning only.
        """
        # Run the default parser to get most of the data. This correctly parses
        # everything except for the Telluric.model data table.
        d, version_passed, type_passed, parsed_hdus = super()._parse(hdu, allow_subclasses=True)
        if not type_passed:
            msgs.error('The HDU(s) cannot be parsed by a {0} object!'.format(cls.__name__))
        if not version_passed:
            _f = msgs.error if chk_version else msgs.warn
            _f('Current version of {0} object in code (v{1})'.format(cls.__name__, cls.version)
               + ' does not match version used to write your HDU(s)!')

        # Load the telluric model, if it exists
        if 'TELLURIC' in [h.name for h in hdu]:
            # Instantiate the Telluric model from the header, if it exists
            # TODO: I don't like this, but this is the fastest work around
            d['telluric'] = telluric.Telluric.from_hdu(hdu['TELLURIC'], ext_pseudo='MODEL',
                                                       chk_version=chk_version)
        # Return the constructed object
        return super().from_dict(d=d)

    def compute_zeropoint(self):
        """
        Dummy method overloaded by subclasses
        """
        pass

    def run(self):
        """
        Execute the sensitivity function calculations.
        """
        # Compute the sensitivity function
        self.compute_zeropoint()

        # Extrapolate the zeropoint based on par['extrap_blu'], par['extrap_red']
        self.wave, self.zeropoint = self.extrapolate(samp_fact=self.par['samp_fact'])
        if self.splice_multi_det:
            self.wave_splice, self.zeropoint_splice = self.splice()

        # Compute the throughput
        self.throughput, self.throughput_splice = self.compute_throughput()

        # Write out QA and throughput plots
        self.write_QA()

    def eval_zeropoint(self, wave, iorddet):
        """
        Dummy method, overloaded by subclasses
        """
        pass

    def extrapolate(self, samp_fact=1.5):
        """
        Extrapolates the sensitivity function to cover an extra wavelength range
        set by the ``extrapl_blu`` and ``extrap_red`` parameters. This is
        important for making sure that the sensitivity function can be applied
        to data with slightly different wavelength coverage etc.

        Parameters
        ----------
        samp_fact : :obj:`float`
            Parameter governing the sampling of the wavelength grid used for the
            extrapolation.

        Returns
        -------
        wave_extrap : `numpy.ndarray`_
            Extrapolated wavelength array
        zeropoint_extrap : `numpy.ndarray`_
            Extrapolated sensitivity function
        """
        # Create a new set of oversampled and padded wavelength grids for the
        # extrapolation
        wave_extrap_min = self.sens['WAVE_MIN'].data * (1.0 - self.par['extrap_blu'])
        wave_extrap_max = self.sens['WAVE_MAX'].data * (1.0 + self.par['extrap_red'])
        nspec_extrap = 0

        # Find the maximum size of the wavewlength grids, since we want
        # everything to have the same
        for idet in range(self.norderdet):
            wave = self.wave_cnts if self.wave_cnts.ndim == 1 else self.wave_cnts[:, idet]
            dwave_data, dloglam_data, resln_guess, pix_per_sigma = wvutils.get_sampling(wave)
            nspec_now = np.ceil(samp_fact * (wave_extrap_max[idet] - wave_extrap_min[idet])
                                / dwave_data).astype(int)
            nspec_extrap = np.max([nspec_now, nspec_extrap])

        # Create the wavelength grid
        wave_extrap = np.outer(np.arange(nspec_extrap),
                               (wave_extrap_max - wave_extrap_min) / (nspec_extrap - 1)) \
                      + np.outer(np.ones(nspec_extrap), wave_extrap_min)
        zeropoint_extrap = np.zeros_like(wave_extrap)

        # Evaluate extrapolated zerpoint for all orders detectors
        for iorddet in range(self.norderdet):
            zeropoint_extrap[:, iorddet] = self.eval_zeropoint(wave_extrap[:,iorddet], iorddet)

        self.steps.append(inspect.stack()[0][3])
        return wave_extrap, zeropoint_extrap

    def splice(self):
        """
        Routine to splice together sensitivity functions into one global
        sensitivity function for spectrographs with multiple detectors extending
        across the wavelength direction.

        Returns
        -------
        wave_splice : `numpy.ndarray`_, shape is (nspec_splice, 1)
            wavelength array
        zeropoint_splice: `numpy.ndarray`_, shape is (nspec_splice, 1)
            zero-point array
        """

        msgs.info(f"Merging sensfunc for {self.norderdet} detectors {self.par['multi_spec_det']}")
        wave_splice_min = self.wave[self.wave > 1.0].min()
        wave_splice_max = self.wave[self.wave > 1.0].max()
        wave_splice_1d, _, _ = wvutils.get_wave_grid(waves=self.wave, wave_method='linear',
                                                     wave_grid_min=wave_splice_min,
                                                     wave_grid_max=wave_splice_max,
                                                     spec_samp_fact=1.0)
        zeropoint_splice_1d = np.zeros_like(wave_splice_1d)
        for idet in range(self.norderdet):
            wave_min = self.sens['WAVE_MIN'][idet]
            wave_max = self.sens['WAVE_MAX'][idet]
            if idet == 0:
                # If this is the bluest detector, extrapolate to wave_extrap_min
                wave_mask_min = wave_splice_min
                wave_mask_max = wave_max
            elif idet == (self.norderdet - 1):
                # If this is the reddest detector, extrapolate to wave_extrap_max
                wave_mask_min = wave_min
                wave_mask_max = wave_splice_max
            else:
                wave_mask_min = wave_min
                wave_mask_max = wave_max
            splice_wave_mask = (wave_splice_1d >= wave_mask_min) & (wave_splice_1d <= wave_mask_max)
            zeropoint_splice_1d[splice_wave_mask] \
                    = self.eval_zeropoint(wave_splice_1d[splice_wave_mask], idet)

        # Interpolate over gaps
        zeros = zeropoint_splice_1d == 0.
        if np.any(zeros):
            msgs.info("Interpolating over gaps (and extrapolating with fill_value=1, if need be)")
            interp_func = scipy.interpolate.interp1d(wave_splice_1d[np.invert(zeros)],
                                                     zeropoint_splice_1d[np.invert(zeros)],
                                                     kind='nearest', fill_value=0.,
                                                     bounds_error=False) #
            #kind='nearest', fill_value='extrapoloate', bounds_error=False)
            #  extrapolate fails for JXP, even on 1.4.1
            zero_values = interp_func(wave_splice_1d[zeros])
            zeropoint_splice_1d[zeros] = zero_values

        nspec_splice = wave_splice_1d.size
        self.steps.append(inspect.stack()[0][3])
        return wave_splice_1d.reshape(nspec_splice,1), zeropoint_splice_1d.reshape(nspec_splice,1)

    def compute_throughput(self):
        """
        Compute the spectroscopic throughput

        Returns
        -------
        throughput : `numpy.ndarray`_, :obj:`float`, shape is (nspec, norders)
            Throughput measurements

        throughput_splice : `numpy.ndarray`_, :obj:`float`, shape is (nspec_splice, norders)
            Throughput measurements for spliced spectra
        """

        # Set the throughput to be -1 in places where it is not defined.
        throughput = np.full_like(self.zeropoint, -1.0)
        for idet in range(self.wave.shape[1]):
            wave_gpm = (self.wave[:,idet] >= self.sens['WAVE_MIN'][idet]) \
                            & (self.wave[:,idet] <= self.sens['WAVE_MAX'][idet]) \
                            & (self.wave[:,idet] > 1.0)
            throughput[:,idet][wave_gpm] \
                    = flux_calib.zeropoint_to_throughput(self.wave[:,idet][wave_gpm],
                                                         self.zeropoint[:,idet][wave_gpm],
                                                         self.spectrograph.telescope.eff_aperture())
        if self.splice_multi_det:
            wave_gpm = (self.wave_splice >= np.amin(self.sens['WAVE_MIN'])) \
                            & (self.wave_splice <= np.amax(self.sens['WAVE_MAX'])) \
                            & (self.wave_splice > 1.0)
            throughput_splice = np.zeros_like(self.wave_splice)
            throughput_splice[wave_gpm] \
                    = flux_calib.zeropoint_to_throughput(self.wave_splice[wave_gpm],
                                                         self.zeropoint_splice[wave_gpm],
                                                         self.spectrograph.telescope.eff_aperture())
        else:
            throughput_splice = None

        return throughput, throughput_splice

    def write_QA(self):
        """
        Write out zeropoint QA files
        """
        utils.pyplot_rcparams()

        # Plot QA for zeropoint
        if 'Echelle' in self.spectrograph.pypeline:
            order_or_det = self.spectrograph.orders[np.arange(self.norderdet)]
            order_or_det_str = 'order'
        else:
            order_or_det = np.arange(self.norderdet) + 1
            order_or_det_str = 'det'

        spec_str = f' {self.spectrograph.name} {self.spectrograph.pypeline} ' \
                   f'{self.spectrograph.dispname} '
        zp_title = ['PypeIt Zeropoint QA for' + spec_str + order_or_det_str
                    + f'={order_or_det[idet]}' for idet in range(self.norderdet)]
        thru_title = [order_or_det_str + f'={order_or_det[idet]}'
                        for idet in range(self.norderdet)]

        is_odd = self.norderdet % 2 != 0
        npages = int(np.ceil(self.norderdet/2)) if is_odd else self.norderdet//2 + 1

        # TODO: PDF page logic is a bit complicated becauase we want to plot two
        # plots per page, but the number of pages depends on the number of
        # order/det. Consider just dumping out a set of plots or revamp once we
        # have a dashboard.
        with PdfPages(self.qafile) as pdf:
            for ipage in range(npages):
                figure, (ax1, ax2) = plt.subplots(2, figsize=(8.27, 11.69))
                if (2 * ipage) < self.norderdet:
                    flux_calib.zeropoint_qa_plot(self.sens['SENS_WAVE'][2*ipage],
                                                 self.sens['SENS_ZEROPOINT'][2*ipage],
                                                 self.sens['SENS_ZEROPOINT_GPM'][2*ipage],
                                                 self.sens['SENS_ZEROPOINT_FIT'][2*ipage],
                                                 self.sens['SENS_ZEROPOINT_FIT_GPM'][2*ipage],
                                                 title=zp_title[2*ipage], axis=ax1)
                if (2*ipage + 1) < self.norderdet:
                    flux_calib.zeropoint_qa_plot(self.sens['SENS_WAVE'][2*ipage+1],
                                                 self.sens['SENS_ZEROPOINT'][2*ipage+1],
                                                 self.sens['SENS_ZEROPOINT_GPM'][2*ipage+1],
                                                 self.sens['SENS_ZEROPOINT_FIT'][2*ipage+1],
                                                 self.sens['SENS_ZEROPOINT_FIT_GPM'][2*ipage+1],
                                                 title=zp_title[2*ipage+1], axis=ax2)

                if self.norderdet == 1:
                    # For single order/det just finish up after the first page
                    ax2.remove()
                    pdf.savefig()
                    plt.close('all')
                elif (self.norderdet > 1) & (ipage < npages-1):
                    # For multi order/det but not on the last page, finish up.
                    # No need to remove ax2 since there are always 2 plots per
                    # page except on the last page
                    pdf.savefig()
                    plt.close('all')
                else:
                    # For multi order/det but on the last page, add order/det
                    # summary plot to the final page Deal with even/odd page
                    # logic for axes
                    if is_odd:
                        # add order/det summary plot to axis 2 of current page
                        axis=ax2
                    else:
                        axis=ax1
                        ax2.remove()
                    for idet in range(self.norderdet):
                        # define the color
                        rr = (np.max(order_or_det) - order_or_det[idet]) \
                                / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
                        gg = 0.0
                        bb = (order_or_det[idet] - np.min(order_or_det)) \
                                / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
                        wave_gpm = self.sens['SENS_WAVE'][idet] > 1.0
                        axis.plot(self.sens['SENS_WAVE'][idet,wave_gpm],
                                  self.sens['SENS_ZEROPOINT_FIT'][idet,wave_gpm],
                                  color=(rr, gg, bb), linestyle='-', linewidth=2.5,
                                  label=thru_title[idet], zorder=5 * idet)

                    _wave_min = np.amin(self.sens['WAVE_MIN'])
                    _wave_max = np.amax(self.sens['WAVE_MAX'])

                    # If we are splicing, overplot the spliced zeropoint
                    if self.splice_multi_det:
                        wave_slice_gpm = (self.wave_splice >= _wave_min) \
                                            & (self.wave_splice <= _wave_max) \
                                            & (self.wave_splice > 1.0)
                        axis.plot(self.wave_splice[wave_slice_gpm].flatten(),
                                  self.zeropoint_splice[wave_slice_gpm].flatten(), color='black',
                                  linestyle='-', linewidth=2.5, label='Spliced Zeropoint',
                                  zorder=30, alpha=0.3)

                    wave_gpm = self.sens['SENS_WAVE'] > 1.0
                    axis.set_xlim((0.98 * _wave_min, 1.02 * _wave_max))
                    axis.set_ylim((0.95 * np.amin(self.sens['SENS_ZEROPOINT_FIT'][wave_gpm]),
                                   1.05 * np.amax(self.sens['SENS_ZEROPOINT_FIT'][wave_gpm])))
                    axis.legend(fontsize=14)
                    axis.set_xlabel('Wavelength (Angstroms)')
                    axis.set_ylabel('Zeropoint (AB mag)')
                    axis.set_title('PypeIt Zeropoints for' + spec_str, fontsize=12)

                    pdf.savefig()
                    plt.close('all')

        # Plot throughput curve(s) for all orders/det
        fig = plt.figure(figsize=(12,8))
        axis = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        for idet in range(self.wave.shape[1]):
            # define the color
            rr = (np.max(order_or_det) - order_or_det[idet]) \
                    / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
            gg = 0.0
            bb = (order_or_det[idet] - np.min(order_or_det)) \
                    / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
            gpm = (self.throughput[:, idet] >= 0.0)
            axis.plot(self.wave[gpm,idet], self.throughput[gpm,idet], color=(rr, gg, bb),
                      linestyle='-', linewidth=2.5, label=thru_title[idet], zorder=5*idet)
        if self.splice_multi_det:
            axis.plot(self.wave_splice[wave_slice_gpm].flatten(),
                      self.throughput_splice[wave_slice_gpm].flatten(), color='black',
                      linestyle='-', linewidth=2.5, label='Spliced Throughput', zorder=30,
                      alpha=0.3)

        axis.set_xlim((0.98*self.wave[self.throughput >=0.0].min(),
                       1.02*self.wave[self.throughput >=0.0].max()))
        axis.set_ylim((0.0, 1.05*self.throughput[self.throughput >=0.0].max()))
        axis.legend()
        axis.set_xlabel('Wavelength (Angstroms)')
        axis.set_ylabel('Throughput')
        axis.set_title('PypeIt Throughput for' + spec_str)
        fig.savefig(self.thrufile)

    @classmethod
    def sensfunc_weights(cls, sensfile, waves, debug=False, extrap_sens=True):
        """
        Get the weights based on the sensfunc

        Args:
            sensfile (str):
                the name of your fits format sensfile
            waves (ndarray): (nspec, norders, nexp) or (nspec, norders)
                wavelength grid for your output weights
            debug (bool): default=False
                show the weights QA

        Returns:
            ndarray: sensfunc weights evaluated on the input waves
            wavelength grid
        """
        sens = cls.from_file(sensfile)
    #    wave, zeropoint, meta_table, out_table, header_sens = sensfunc.SensFunc.load(sensfile)

        if waves.ndim == 2:
            nspec, norder = waves.shape
            nexp = 1
            waves_stack = np.reshape(waves, (nspec, norder, 1))
        elif waves.ndim == 3:
            nspec, norder, nexp = waves.shape
            waves_stack = waves
        else:
            msgs.error('Unrecognized dimensionality for waves')

        weights_stack = np.zeros_like(waves_stack)

        if norder != sens.zeropoint.shape[1]:
            msgs.error('The number of orders in {:} does not agree with your data. Wrong sensfile?'.format(sensfile))

        for iord in range(norder):
            for iexp in range(nexp):
                sensfunc_iord = flux_calib.get_sensfunc_factor(waves_stack[:,iord,iexp],
                                                               sens.wave[:,iord],
                                                               sens.zeropoint[:,iord], 1.0,
                                                               extrap_sens=extrap_sens)
                weights_stack[:,iord,iexp] = utils.inverse(sensfunc_iord)

        if debug:
            coadd.weights_qa(waves_stack, weights_stack, (waves_stack > 1.0), title='sensfunc_weights')

        if waves.ndim == 2:
            weights_stack = np.reshape(weights_stack, (nspec, norder))

        return weights_stack


# TODO Add a method which optionally merges sensfunc using the nsens > 1 logic


class IRSensFunc(SensFunc):
    r"""
    Determine a sensitivity functions from standard-star spectra. Should only
    be used with NIR spectra (:math:`\lambda > 7000` angstrom).

    Args:
        spec1dfile (:obj:`str`):
            PypeIt spec1d file for the standard file.
        sensfile (:obj:`str`):
            File name for the sensitivity function data.
        par (:class:`~pypeit.par.pypeitpar.SensFuncPar`, optional):
            The parameters required for the sensitivity function computation.
        debug (:obj:`bool`, optional):
            Run in debug mode.
    """

    _algorithm = 'IR'
    """Algorithm used for the sensitivity calculation."""

    def compute_zeropoint(self):
        """
        Calls routine to compute the sensitivity function.

        Returns
        -------
        TelObj : :class:`~pypeit.core.telluric.Telluric`
            Best-fitting telluric model
        """
        self.telluric = telluric.sensfunc_telluric(self.wave_cnts, self.counts, self.counts_ivar,
                                                   self.counts_mask, self.meta_spec['EXPTIME'],
                                                   self.meta_spec['AIRMASS'], self.std_dict,
                                                   self.par['IR']['telgridfile'],
                                                   polyorder=self.par['polyorder'],
                                                   ech_orders=self.meta_spec['ECH_ORDERS'],
                                                   resln_guess=self.par['IR']['resln_guess'],
                                                   resln_frac_bounds=self.par['IR']['resln_frac_bounds'],
                                                   sn_clip=self.par['IR']['sn_clip'],
                                                   mask_hydrogen_lines=self.par['mask_hydrogen_lines'],
                                                   maxiter=self.par['IR']['maxiter'],
                                                   lower=self.par['IR']['lower'],
                                                   upper=self.par['IR']['upper'],
                                                   delta_coeff_bounds=self.par['IR']['delta_coeff_bounds'],
                                                   minmax_coeff_bounds=self.par['IR']['minmax_coeff_bounds'],
                                                   tol=self.par['IR']['tol'],
                                                   popsize=self.par['IR']['popsize'],
                                                   recombination=self.par['IR']['recombination'],
                                                   polish=self.par['IR']['polish'],
                                                   disp=self.par['IR']['disp'], debug=self.debug,
                                                   debug_init=self.debug)

        # Copy the relevant metadata
        self.std_name = self.telluric.std_name
        self.std_cal = self.telluric.std_cal
        self.std_ra = self.telluric.std_ra
        self.std_dec = self.telluric.std_dec
        self.airmass = self.telluric.airmass
        self.exptime = self.telluric.exptime

        # Instantiate the main output data table
        self.sens = self.empty_sensfunc_table(self.telluric.norders, self.telluric.wave_grid.size,
                                              ncoeff=self.telluric.max_ntheta_obj)

        # For stupid reasons related to how astropy tables will let me store
        # this data I have to redundantly write out the wavelength grid for each
        # order.  In actuality a variable length wavelength grid is needed which
        # is a subset of the original fixed grid, but there is no way to write
        # this to disk.  Perhaps the better solution would be a list of of
        # astropy tables, one for each order, that would save some space, since
        # we could then write out only the subset used to the table column. I'm
        # going with this inefficient approach for now, since eventually we will
        # want to create a data container for these outputs. That said, coming
        # up with that standarized model is a bit tedious given the somewhat
        # heterogenous outputs of the two different fluxing algorithms.

        # Extract and construct the relevant data using the telluric model
        self.sens['SENS_COEFF'] = self.telluric.model['OBJ_THETA']
        self.sens['WAVE_MIN'] = self.telluric.model['WAVE_MIN']
        self.sens['WAVE_MAX'] = self.telluric.model['WAVE_MAX']
        self.sens['ECH_ORDERS'] = self.telluric.model['ECH_ORDERS']
        s = self.telluric.model['IND_LOWER']
        e = self.telluric.model['IND_UPPER']+1
        self.sens['POLYORDER_VEC'] = self.telluric.model['POLYORDER_VEC']
        for i in range(self.norderdet):
            # Compute and assign the zeropint_data from the input data and the
            # best-fit telluric model
            self.sens['SENS_WAVE'][i,s[i]:e[i]] = self.telluric.wave_grid[s[i]:e[i]]
            self.sens['SENS_ZEROPOINT_GPM'][i,s[i]:e[i]] = self.telluric.mask_arr[s[i]:e[i],i]
            self.sens['SENS_COUNTS_PER_ANG'][i,s[i]:e[i]] = self.telluric.flux_arr[s[i]:e[i],i]
            N_lam = self.sens['SENS_COUNTS_PER_ANG'][i,s[i]:e[i]] / self.exptime
            self.sens['SENS_ZEROPOINT'][i,s[i]:e[i]], _ \
                    = flux_calib.compute_zeropoint(self.sens['SENS_WAVE'][i,s[i]:e[i]], N_lam,
                                                   self.sens['SENS_ZEROPOINT_GPM'][i,s[i]:e[i]],
                                                   self.telluric.obj_dict_list[i]['flam_true'],
                                                   tellmodel=self.telluric.tellmodel_list[i])
            # TODO: func is always 'legendre' because that is what's set by
            # sensfunc_telluric
            self.sens['SENS_ZEROPOINT_FIT'][i,s[i]:e[i]] \
                    = fitting.evaluate_fit(self.sens['SENS_COEFF'][
                                                i,:self.sens['POLYORDER_VEC'][i]+2],
                                           self.telluric.func, self.sens['SENS_WAVE'][i,s[i]:e[i]],
                                           minx=self.sens['WAVE_MIN'][i],
                                           maxx=self.sens['WAVE_MAX'][i])
            self.sens['SENS_ZEROPOINT_FIT_GPM'][i,s[i]:e[i]] = self.telluric.outmask_list[i]

    def eval_zeropoint(self, wave, iorddet):
        """
        Evaluate the sensitivity function zero-points

        Parameters
        ----------
        wave : `numpy.ndarray`_, shape is (nspec)
            Wavelength array
        iorddet : :obj:`int`
            Order or detector (0-indexed)

        Returns
        -------
        zeropoint : `numpy.ndarray`_, shape is (nspec,)
        """
        return fitting.evaluate_fit(self.sens['SENS_COEFF'][
                                        iorddet,:self.telluric.model['POLYORDER_VEC'][iorddet]+2],
                                    self.telluric.func, wave, minx=self.sens['WAVE_MIN'][iorddet],
                                    maxx=self.sens['WAVE_MAX'][iorddet])


class UVISSensFunc(SensFunc):
    r"""
    Determine a sensitivity functions from standard-star spectra. Should only
    be used with UVIS spectra (:math:`\lambda < 7000` angstrom).

    Args:
        spec1dfile (:obj:`str`):
            PypeIt spec1d file for the standard file.
        sensfile (:obj:`str`):
            File name for the sensitivity function data.
        par (:class:`~pypeit.par.pypeitpar.SensFuncPar`, optional):
            The parameters required for the sensitivity function computation.
        debug (:obj:`bool`, optional):
            Run in debug mode.
    """

    _algorithm = 'UVIS'
    """Algorithm used for the sensitivity calculation."""

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):
        super().__init__(spec1dfile, sensfile, par=par, debug=debug)

        # Add some cards to the meta spec. These should maybe just be added
        # already in unpack object
        self.meta_spec['LATITUDE'] = self.spectrograph.telescope['latitude']
        self.meta_spec['LONGITUDE'] = self.spectrograph.telescope['longitude']

    def compute_zeropoint(self):
        """
        Calls routine to compute the sensitivity function.
        """
        meta_table, out_table = flux_calib.sensfunc(self.wave_cnts, self.counts, self.counts_ivar,
                                                    self.counts_mask, self.meta_spec['EXPTIME'],
                                                    self.meta_spec['AIRMASS'], self.std_dict,
                                                    self.meta_spec['LONGITUDE'],
                                                    self.meta_spec['LATITUDE'],
                                                    self.par['UVIS']['extinct_file'],
                                                    self.meta_spec['ECH_ORDERS'],
                                                    polyorder=self.par['polyorder'],
                                                    hydrogen_mask_wid=self.par['hydrogen_mask_wid'],
                                                    mask_hydrogen_lines=self.par['mask_hydrogen_lines'],
                                                    mask_helium_lines=self.par['mask_helium_lines'],
                                                    nresln=self.par['UVIS']['nresln'],
                                                    resolution=self.par['UVIS']['resolution'],
                                                    trans_thresh=self.par['UVIS']['trans_thresh'],
                                                    polycorrect=self.par['UVIS']['polycorrect'],
                                                    polyfunc=self.par['UVIS']['polyfunc'],
                                                    debug=self.debug)

        # Copy the relevant metadata
        self.std_name = meta_table['STD_NAME'][0]
        self.std_cal = meta_table['CAL_FILE'][0]
        self.std_ra = meta_table['STD_RA'][0]
        self.std_dec = meta_table['STD_DEC'][0]
        self.airmass = meta_table['AIRMASS'][0]
        self.exptime = meta_table['EXPTIME'][0]

        norder, nspec = out_table['SENS_ZEROPOINT'].shape

        # Instantiate the main output data table
        self.sens = self.empty_sensfunc_table(norder, nspec)

        # Copy the relevant data
        # NOTE: SENS_COEFF is empty!
        self.sens['SENS_WAVE'] = out_table['SENS_WAVE']
        self.sens['SENS_COUNTS_PER_ANG'] = out_table['SENS_COUNTS_PER_ANG']
        self.sens['SENS_ZEROPOINT'] = out_table['SENS_ZEROPOINT']
        self.sens['SENS_ZEROPOINT_GPM'] = out_table['SENS_ZEROPOINT_GPM']
        self.sens['SENS_ZEROPOINT_FIT'] = out_table['SENS_ZEROPOINT_FIT']
        self.sens['SENS_ZEROPOINT_FIT_GPM'] = out_table['SENS_ZEROPOINT_FIT_GPM']
        if self.meta_spec['ECH_ORDERS'] is not None:
            self.sens['ECH_ORDERS'] = self.meta_spec['ECH_ORDERS']
        self.sens['POLYORDER_VEC'] = np.full(norder, self.par['polyorder'])
        self.sens['WAVE_MIN'] = out_table['WAVE_MIN']
        self.sens['WAVE_MAX'] = out_table['WAVE_MAX']

    def eval_zeropoint(self, wave, iorddet):
        """
        Evaluate the sensitivity function zero-points

        Parameters
        ----------
        wave : `numpy.ndarray`_, shape is (nspec)
            Wavelength array
        iorddet : :obj:`int`
            Order or detector (0-indexed)

        Returns
        -------
        zeropoint : `numpy.ndarray`_, shape is (nspec,)
        """
        # This routine can extrapolate
        return scipy.interpolate.interp1d(self.sens['SENS_WAVE'][iorddet,:],
                                          self.sens['SENS_ZEROPOINT_FIT'][iorddet,:],
                                          bounds_error=False, fill_value='extrapolate')(wave)



