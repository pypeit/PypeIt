"""
Implements the objects used to construct sensitivity functions.

.. include:: ../include/links.rst
"""
import os
import inspect

from IPython import embed

import numpy as np
import scipy
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import fits
from astropy import table

from pypeit import msgs
from pypeit import specobjs
from pypeit import utils
from pypeit import io
from pypeit import datamodel
from pypeit.core import flux_calib
from pypeit.core import telluric
from pypeit.core import fitting
from pypeit.core import coadd
from pypeit.core.wavecal import wvutils
from pypeit.core import meta
from pypeit.spectrographs.util import load_spectrograph


# TODO Add the data model up here as a standard thing using DataContainer.

#TODO Should this be a master frame? I think not.
#TODO Standard output location for sensfunc?

# TODO Add some QA plots, and plots to the screen if show is set.

class SensFunc(datamodel.DataContainer):
    """
    Base class for generating sensitivity functions from a standard-star
    spectrum.

    This class should not be instantated by itself; instead instantiate
    either :class:`UVISSensFunc` or :class:`IRSensFunc`, depending on the
    wavelength range of your data (UVIS for :math:`\lambda < 7000` angstrom,
    IR for :math:`\lambda > 7000` angstrom.)

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

    version = '1.0.0'
    """Datamodel version."""

    datamodel = {'method': dict(otype=str, descr='PypeIt method used to construct the data'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'pypeline': dict(otype=str, descr='PypeIt pipeline reduction path'),
                 'spec1df': dict(otype=str,
                                 descr='PypeIt spec1D file used to construct sensitivity function'),
                 'stdname': dict(otype=str, descr='Type of standard source'),
                 'stdfile': dict(otype=str,
                                 descr='File name (or shorthand) with the standard flux data'),
                 'std_ra': dict(otype=float, descr='RA of the standard source'),
                 'std_dec': dict(otype=float, descr='DEC of the standard source'),
                 'airmass': dict(otype=float, descr='Airmass of the observation'),
                 'exptime': dict(otype=float, descr='Exposure time'),
                 'wave': dict(otype=np.ndarray, descr='Wavelength vectors'),
                 'zeropoint': dict(otype=np.ndarray, descr='Sensitivity function zero-points'),
                 'throughput': dict(otype=np.ndarray, descr='Estimated throughput'),
                 'telluric': dict(otype=telluric.TelluricData, descr='Telluric model'),
                 'sens': dict(otype=table.Table, descr='Table with the sensitivity function')}
    """DataContainer datamodel."""

    algorithm = None
    """Algorithm used for the sensitivity calculation."""

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec1dfile, sensfile, par=None, debug=False):
        return next(c for c in cls.__subclasses__()
                    if c.__name__ == '{0}SensFunc'.format(par['algorithm']))(spec1dfile, sensfile,
                                                                             par=par, debug=debug)

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):
        datamodel.DataContainer.__init__(self)

        # Arguments
        self.spec1df = spec1dfile
        self.sensfile = sensfile
        # Set spectrograph
        header = fits.getheader(self.spec1df)
        self.PYP_SPEC = header['PYP_SPEC']
        self.spectrograph = load_spectrograph(self.PYP_SPEC)
        self.spectrograph.dispname = header['DISPNAME']
        # Get the algorithm parameters
        self.par = self.spectrograph.default_pypeit_par()['sensfunc'] if par is None else par
        # TODO: Check the type of the parameter object?
        # QA and throughput plot filenames
        self.qafile = sensfile.replace('.fits', '') + '_QA.pdf'
        self.thrufile = sensfile.replace('.fits', '') + '_throughput.pdf'

        self.debug = debug
        # TODO: Is `steps` useful? Removing for now, but easy to bring back...
#        self.steps = []
        self.splice_multi_det = True if self.par['multi_spec_det'] is not None else False

        # Read in the Standard star data
        sobjs_std = specobjs.SpecObjs.from_fitsfile(self.spec1df).get_std(
                            multi_spec_det=self.par['multi_spec_det'])

        # Unpack standard
        wave, counts, counts_ivar, counts_mask, self.meta_spec, header \
                = sobjs_std.unpack_object(ret_flam=False)

        # Perform any instrument tweaks
        wave_twk, counts_twk, counts_ivar_twk, counts_mask_twk \
                = self.spectrograph.tweak_standard(wave, counts, counts_ivar, counts_mask,
                                                   self.meta_spec)

        # Reshape to 2d arrays
        self.wave, self.counts, self.counts_ivar, self.counts_mask, self.nspec, self.norders \
                = utils.spec_atleast_2d(wave_twk, counts_twk, counts_ivar_twk, counts_mask_twk)

        # If the user provided RA and DEC use those instead of what is in meta
        star_ra = self.meta_spec['RA'] if self.par['star_ra'] is None else self.par['star_ra']
        star_dec = self.meta_spec['DEC'] if self.par['star_dec'] is None else self.par['star_dec']
        # Convert to decimal deg, as need be
        star_ra, star_dec = meta.convert_radec(star_ra, star_dec)

        # Read in standard star dictionary
        self.std_dict = flux_calib.get_standard_spectrum(star_type=self.par['star_type'],
                                                         star_mag=self.par['star_mag'],
                                                         ra=star_ra, dec=star_dec)

    def _init_internals(self):
        """
        Add attributes not included in the data model.
        """
        self.norders = None
        self.nspec = None

        self.sensfile = None
        self.spectrograph = None
        self.par = None
        self.qafile = None
        self.thrufile = None
        self.debug = None
#        self.steps = None
        self.splice_multi_det = None
        self.meta_spec = None

        self.counts = None
        self.counts_ivar = None
        self.counts_mask = None
        self.std_dict = None

    @staticmethod
    def empty_sensfunc_table(norders, nspec, ncoeff=0):
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
            table.Column(name='WAVE', dtype=float, length=norders, shape=(nspec,),
                         description='Wavelength vector'),
            table.Column(name='COUNTS_PER_ANG', dtype=float, length=norders, shape=(nspec,),
                         description='Flux in counts per angstrom'),
            table.Column(name='ZEROPOINT', dtype=float, length=norders, shape=(nspec,),
                         description='Measured sensitivity zero-point data'),
            table.Column(name='ZEROPOINT_GPM', dtype=bool, length=norders, shape=(nspec,),
                         description='Good-pixel mask for the measured zero points'),
            table.Column(name='ZEROPOINT_FIT', dtype=float, length=norders, shape=(nspec,),
                         description='Best-fit smooth model to the zero points'),
            table.Column(name='ZEROPOINT_FIT_GPM', dtype=bool, length=norders, shape=(nspec,),
                         description='Good-pixel mask for the model zero points'),
            table.Column(name='COEFF', dtype=float, length=norders, shape=(ncoeff,),
                         description='Coefficients of smooth model fit to zero points'),
            table.Column(name='ECH_ORDER', dtype=int, length=norders,
                         description='Echelle order for this specrum (echelle data only)'),
            table.Column(name='WAVE_MIN', dtype=float, length=norders,
                         description='Minimum wavelength included in the fit'),
            table.Column(name='WAVE_MAX', dtype=float, length=norders,
                         description='Maximum wavelength included in the fit')])

    def compute_zeropoint(self):
        """
        Dummy method overloaded by subclasses
        """
        pass

    def eval_zeropoint(self, wave, iorddet):
        """
        Dummy method, overloaded by subclasses
        """
        pass

    def run(self):
        # Compute the sensitivity function
        self.compute_zeropoint()

        # Extrapolate the zeropoint based on par['extrap_blu'], par['extrap_red']
        self.wave_zp, self.zeropoint = self.extrapolate(samp_fact=self.par['samp_fact'])
        if self.splice_multi_det:
            self.wave_zp, self.zeropoint = self.splice()

        # Compute the throughput
        self.throughput = self.compute_throughput()

        # Write out QA and throughput plots
        self.write_QA()

        # If the zeropoint has just one order, or detectors were spliced, flatten the output
        # TODO -- Consider having self.splice() return a 2D array instead of 1D for multi_det
        if self.wave_zp.ndim == 2:
            if self.wave_zp.shape[1] == 1:
                self.wave_zp = self.wave_zp.flatten()
                self.zeropoint = self.zeropoint.flatten()
                self.throughput = self.throughput.flatten()

    def extrapolate(self, samp_fact=1.5):
        """
        Extrapolates the sensitivity function to cover an extra wavelength
        range set by the ``extrapl_blu`` and ``extrap_red`` parameters. This
        is important for making sure that the sensitivity function can be
        applied to data with slightly different wavelength coverage etc.

        Parameters
        ----------
        samp_fact: float
            Parameter governing the sampling of the wavelength grid used for
            the extrapolation.

        Returns
        -------
        wave_extrap: ndarray
            Extrapolated wavelength array
        zeropoint_extrap: ndarray
            Extrapolated sensfunc
        """

        # Create a new set of oversampled and padded wavelength grids for the
        # extrapolation
        wave_extrap_min = self.sens['WAVE_MIN'].data * (1.0 - self.par['extrap_blu'])
        wave_extrap_max = self.sens['WAVE_MAX'].data * (1.0 + self.par['extrap_red'])
        nspec_extrap = 0

        # Find the maximum size of the wavelength grids, since we want
        # everything to have the same wavelength range
        for idet in range(self.norders):
            wave = self.wave if self.wave.ndim == 1 else self.wave[:,idet]
            dwave_data, dloglam_data, resln_guess, pix_per_sigma = wvutils.get_sampling(wave)
            nspec_now = np.ceil(samp_fact * (wave_extrap_max[idet] - wave_extrap_min[idet])
                                 / dwave_data).astype(int)
            nspec_extrap = np.max([nspec_now, nspec_extrap])

        # Create the wavelength grid
        wave_extrap = np.outer(np.arange(nspec_extrap),
                               (wave_extrap_max - wave_extrap_min) / (nspec_extrap - 1)) \
                        + np.outer(np.ones(nspec_extrap), wave_extrap_min)
        zeropoint_extrap = np.zeros_like(wave_extrap)

        # Evaluate extrapolated zeropoint for all orders detectors
        for iorddet in range(self.norders):
            zeropoint_extrap[:,iorddet] = self.eval_zeropoint(wave_extrap[:,iorddet], iorddet)

#        self.steps.append(inspect.stack()[0][3])
        return wave_extrap, zeropoint_extrap

    def splice(self):
        """
        Routine to splice together sensitivity functions into one global
        sensitivity function for spectrographs with multiple detectors
        extending across the wavelength direction.

        Returns
        -------
        wave_splice: ndarray, shape (nspec_splice, 1)
        zeropoint_splice: ndarray, shape (nspec_splice, 1)
        """

        msgs.info('Merging sensfunc for {0} detectors {1}'.format(self.norders,
                                                                  self.par['multi_spec_det']))
        wave_splice_min = self.wave_zp[self.wave_zp > 1.0].min()
        wave_splice_max = self.wave_zp[self.wave_zp > 1.0].max()
        wave_splice_1d, _, _ = coadd.get_wave_grid(self.wave_zp, wave_method='linear',
                                                   wave_grid_min=wave_splice_min,
                                                   wave_grid_max=wave_splice_max,
                                                   spec_samp_fact=1.0)
        zeropoint_splice_1d = np.zeros_like(wave_splice_1d)
        for idet in range(self.norders):
            wave_min = self.sens['WAVE_MIN'][idet]
            wave_max = self.sens['WAVE_MAX'][idet]
            if idet == 0:
                # If this is the bluest detector, extrapolate to wave_extrap_min
                wave_mask_min = wave_splice_min
                wave_mask_max = wave_max
            elif idet == (self.norders - 1):
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
            indx = np.logical_not(zeros)
            interp_func = scipy.interpolate.interp1d(wave_splice_1d[indx],
                                                     zeropoint_splice_1d[indx], kind='nearest',
                                                     fill_value=0., bounds_error=False)
            zeropoint_splice_1d[zeros] = interp_func(wave_splice_1d[zeros])

        return wave_splice_1d.reshape(-1,1), zeropoint_splice_1d.reshape(-1,1)

    def compute_throughput(self):
        """
        Compute the spectroscopic throughput.

        Returns
        -------
        throughput : `numpy.ndarray`_
            Throughput measurement for each order/slit

        """
        # Set the throughput to be -1 in places where it is not defined.
        throughput = np.full_like(self.zeropoint, -1.0)
        for idet in range(self.wave_zp.shape[1]):
            wave_gpm = (self.wave_zp[:,idet] >= self.sens['WAVE_MIN'][idet]) \
                            & (self.wave_zp[:,idet] <= self.sens['WAVE_MAX'][idet]) \
                            & (self.wave_zp[:,idet] > 1.0)
            throughput[wave_gpm,idet] \
                    = flux_calib.zeropoint_to_throughput(self.wave_zp[wave_gpm,idet],
                                                         self.zeropoint[wave_gpm,idet],
                                                         self.spectrograph.telescope.eff_aperture())
        return throughput

    def write_QA(self):
        """
        Write out zeropoint QA files
        """
        utils.pyplot_rcparams()

        # Plot QA for zeropoint
        if self.spectrograph.pypeline == 'Echelle':
            order_or_det = self.spectrograph.orders[np.arange(self.norders)]
            order_or_det_str = 'order'
        else:
            order_or_det = np.arange(self.norders) + 1
            order_or_det_str = 'det'

        spec_str = f' {0} {1} {2} '.format(self.spectrograph.name, self.spectrograph.pypeline,
                                           self.spectrograph.dispname)
        zp_title = ['PypeIt Zeropoint QA for' + spec_str + order_or_det_str
                    + '={0}'.format(order_or_det[idet]) for idet in range(self.norders)]
        thru_title = [order_or_det_str + '={0}'.format(order_or_det[idet]) 
                        for idet in range(self.norders)]

        is_odd = self.norders % 2 != 0
        npages = int(np.ceil(self.norders/2)) if is_odd else self.norders//2 + 1

        # TODO: PDF page logic is a bit complicated because we want to plot
        # two plots per page, but the number of pages depends on the number of
        # order/det. Consider just dumping out a set of plots or revamp once we
        # have a dashboard
        with PdfPages(self.qafile) as pdf:
            for ipage in range(npages):
                figure, (ax1, ax2) = plt.subplots(2, figsize=(8.27, 11.69))
                for j, ax in zip([2*ipage, 2*ipage+1], [ax1,ax2]):
                    if j < self.norders:
                        flux_calib.zeropoint_qa_plot(self.sens['WAVE'][j],
                                                     self.sens['ZEROPOINT'][j],
                                                     self.sens['ZEROPOINT_GPM'][j],
                                                     self.sens['ZEROPOINT_FIT'][j],
                                                     self.sens['ZEROPOINT_FIT_GPM'][j],
                                                     title=zp_title[j], axis=ax)
                if self.norders == 1:
                    # For single order/det just finish up after the first page
                    ax2.remove()
                    pdf.savefig()
                    plt.close('all')
                elif (self.norders > 1) & (ipage < npages-1):
                    # For multi order/det but not on the last page, finish up.
                    # No need to remove ax2 since there are always 2 plots per
                    # page except on the last page
                    pdf.savefig()
                    plt.close('all')
                else:
                    # For multi order/det but on the last page, add order/det
                    # summary plot to the final page. Deal with even/odd page
                    # logic for axes
                    if is_odd:
                        # add order/det summary plot to axis 2 of current page
                        axis=ax2
                    else:
                        axis=ax1
                        ax2.remove()
                    for idet in range(self.norders):
                        # define the color
                        rr = (np.max(order_or_det) - order_or_det[idet]) \
                                / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
                        gg = 0.0
                        bb = (order_or_det[idet] - np.min(order_or_det)) \
                                / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
                        wave_gpm = self.sens['WAVE'][idet] > 1.0
                        axis.plot(self.sens[idet]['WAVE'][wave_gpm],
                                  self.sens[idet]['ZEROPOINT_FIT'][wave_gpm],
                                  color=(rr, gg, bb), linestyle='-', linewidth=2.5,
                                  label=thru_title[idet], zorder=5 * idet)
                    wave_zp_gpm = (self.wave_zp >= self.sens['WAVE_MIN'].min()) \
                                    & (self.wave_zp <= self.sens['WAVE_MAX'].max())
                    # If we are splicing, overplot the spliced zeropoint
                    if self.splice_multi_det:
                        axis.plot(self.wave_zp[wave_zp_gpm].flatten(),
                                  self.zeropoint[wave_zp_gpm].flatten(), color='black',
                                  linestyle='-', linewidth=2.5, label='Spliced Zeropoint',
                                  zorder=30, alpha=0.3)

                    axis.set_xlim((0.98 * self.sens['WAVE_MIN'].min(),
                                   1.02 * self.sens['WAVE_MAX'].max()))
                    axis.set_ylim((0.95*self.sens['ZEROPOINT_FIT'][self.sens['WAVE'] > 1.0].min(),
                                   1.05*self.sens['ZEROPOINT_FIT'][self.sens['WAVE'] > 1.0].max()))
                    axis.legend(fontsize=14)
                    axis.set_xlabel('Wavelength (Angstroms)')
                    axis.set_ylabel('Zeropoint (AB mag)')
                    axis.set_title('PypeIt Zeropoints for' + spec_str, fontsize=12)

                    pdf.savefig()
                    plt.close('all')

        # Plot throughput curve(s) for all orders/det
        fig = plt.figure(figsize=(12,8))
        axis = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        for idet in range(self.wave_zp.shape[1]):
            # define the color
            rr = (np.max(order_or_det) - order_or_det[idet]) \
                    / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
            gg = 0.0
            bb = (order_or_det[idet] - np.min(order_or_det)) \
                    / np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
            gpm = (self.throughput[:, idet] >= 0.0)
            axis.plot(self.wave_zp[gpm,idet], self.throughput[gpm,idet], color=(rr, gg, bb),
                      linestyle='-', linewidth=2.5, label=thru_title[idet], zorder=5*idet)

        axis.set_xlim((0.98*self.wave_zp[self.throughput >=0.0].min(),
                       1.02*self.wave_zp[self.throughput >=0.0].max()))
        axis.set_ylim((0.0,1.05*self.throughput[self.throughput >=0.0].max()))
        axis.legend()
        axis.set_xlabel('Wavelength (Angstroms)')
        axis.set_ylabel('Throughput')
        axis.set_title('PypeIt Throughput for' + spec_str)
        fig.savefig(self.thrufile)

# TODO Add a method which optionally merges sensfunc using the nsens > 1 logic

class OLDSensFunc:

    @classmethod
    def load(cls, sensfile):
        # Write to outfile
        msgs.info('Reading sensitivity function from file: {:}'.format(sensfile))
        hdulist = io.fits_open(sensfile)
        header = hdulist[0].header
        wave = hdulist['WAVE'].data
        zeropoint = hdulist['ZEROPOINT'].data
        meta_table = table.Table(hdulist['METADATA'].data)
        out_table  = table.Table(hdulist['OUT_TABLE'].data)
        return wave, zeropoint, meta_table, out_table, header

    def save(self):
        """
        Saves sensitivity to self.sensfile
        """

        # Write to outfile
        msgs.info('Writing sensitivity function results to file: {:}'.format(self.sensfile))

        # Standard init
        hdr = io.initialize_header()

        hdr['PYP_SPEC'] = (self.spectrograph.name, 'PypeIt: Spectrograph name')
        hdr['PYPELINE'] = self.spectrograph.pypeline
        #   - List the completed steps
        hdr['STEPS'] = (','.join(self.steps), 'Completed sensfunc steps')
        #   - Provide the file names
        hdr['SPC1DFIL'] = self.spec1dfile

        # Write the fits file
        data = [self.wave_zp, self.zeropoint, self.throughput]
        extnames = ['WAVE', 'ZEROPOINT', 'THROUGHPUT']
        # Write the fits file
        hdulist = fits.HDUList([fits.PrimaryHDU(header=hdr)] 
                               + [fits.ImageHDU(data=d, name=n) for d, n in zip(data, extnames)])
        hdu_meta = fits.table_to_hdu(self.meta_table)
        hdu_meta.name = 'METADATA'
        hdu_out = fits.table_to_hdu(self.out_table)
        hdu_out.name = 'OUT_TABLE'
        hdulist.append(hdu_meta)
        hdulist.append(hdu_out)
        hdulist.writeto(self.sensfile, overwrite=True, checksum=True)


class IRSensFunc(SensFunc):
    """
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

    algorithm = 'IR'
    """Algorithm used for the sensitivity calculation."""

    def compute_zeropoint(self):
        """
        Basic wrapper for :func:`pypeit.core.telluric.sensfunc_telluric`.
        """
        meta_table, out_table, self.telluric \
                = telluric.sensfunc_telluric(
                        self.wave, self.counts, self.counts_ivar, self.counts_mask,
                        self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'], self.std_dict,
                        self.par['IR']['telgridfile'], ech_orders=self.meta_spec['ECH_ORDERS'],
                        polyorder=self.par['polyorder'], mask_abs_lines=self.par['mask_abs_lines'],
                        delta_coeff_bounds=self.par['IR']['delta_coeff_bounds'],
                        minmax_coeff_bounds=self.par['IR']['minmax_coeff_bounds'],
                        sn_clip=self.par['IR']['sn_clip'], maxiter=self.par['IR']['maxiter'],
                        lower=self.par['IR']['lower'], upper=self.par['IR']['upper'], 
                        tol=self.par['IR']['tol'], popsize=self.par['IR']['popsize'],
                        recombination=self.par['IR']['recombination'],
                        polish=self.par['IR']['polish'], disp=self.par['IR']['disp'],
                        debug_init=self.debug, debug=self.debug)

        # TODO: This avoids messing with telluric.sensfunc_telluric much, but
        # this parsing should be consolidate with that method...

        # Parse the meta_table and out_table into the data model
        self.norders, self.nspec = out_table['TELLURIC'].shape
        ncoeff = out_table['SENS_COEFF'].shape[1]

        self.sens = self.__class__.empty_sensfunc_table(self.norders, self.nspec, ncoeff=ncoeff)

        self.stdname = meta_table['STD_NAME']
        self.stdfile = meta_table['STD_CALFILE']
        self.std_ra = meta_table['STD_RA']
        self.std_dec = meta_table['STD_DEC']
        self.airmass = meta_table['AIRMASS']
        self.exptime = self.meta_spec['EXPTIME']

        self.sens['WAVE'] = out_table['SENS_WAVE']
        self.sens['COUNTS_PER_ANG'] = out_table['SENS_COUNTS_PER_ANG']
        self.sens['ZEROPOINT'] = out_table['SENS_ZEROPOINT']
        self.sens['ZEROPOINT_GPM'] = out_table['SENS_ZEROPOINT_GPM']
        self.sens['ZEROPOINT_FIT'] = out_table['SENS_ZEROPOINT_FIT']
        self.sens['ZEROPOINT_FIT_GPM'] = out_table['SENS_ZEROPOINT_FIT_GPM']
        self.sens['COEFF'] = out_table['SENS_COEFF']
        self.sens['ECH_ORDER'] = self.meta_spec['ECH_ORDERS']
        self.sens['WAVE_MIN'] = out_table['WAVE_MIN']
        self.sens['WAVE_MAX'] = out_table['WAVE_MAX']

#        self.steps.append(inspect.stack()[0][3])

    def eval_zeropoint(self, wave, iorddet):
        """
        Evaluate the model zero-points

        Parameters
        ----------
        wave: ndarray shape (nspec)
           Wavelength array

        iorddet: int
           Order or detector

        Returns
        -------
        zeropoint: ndarray, shape (nspec,)
        """
        coeff = self.telluric.telluric['OBJ_THETA'][:self.TelObj.polyorder_vec[iorddet] + 2]
        return fitting.evaluate_fit(coeff, self.telluric.polyfunc, wave,
                                    minx=self.sens['WAVE_MIN'][iorddet],
                                    maxx=self.sens['WAVE_MAX'][iorddet])


class UVISSensFunc(SensFunc):

    algorithm = 'UVIS'
    """Algorithm used for the sensitivity calculation."""

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):
        super().__init__(spec1dfile, sensfile, par=par, debug=debug)

        # Add some cards to the meta spec. These should maybe just be added already in unpack object
        self.meta_spec['LATITUDE'] = self.spectrograph.telescope['latitude']
        self.meta_spec['LONGITUDE'] = self.spectrograph.telescope['longitude']

    def compute_zeropoint(self):
        """
        Calls routine to compute the sensitivity function.

        Returns
        -------
        meta_table: astropy table
               Table containing zeropoint meta data
        out_table: astropy table
               Table containing zeropoint information.
        """

        meta_table, out_table \
                = flux_calib.sensfunc(self.wave, self.counts, self.counts_ivar, self.counts_mask,
                                      self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'],
                                      self.std_dict, self.meta_spec['LONGITUDE'], 
                                      self.meta_spec['LATITUDE'], self.meta_spec['ECH_ORDERS'],
                                      polyorder=self.par['polyorder'],
                                      balm_mask_wid=self.par['UVIS']['balm_mask_wid'],
                                      nresln=self.par['UVIS']['nresln'],
                                      resolution=self.par['UVIS']['resolution'],
                                      trans_thresh=self.par['UVIS']['trans_thresh'],
                                      polycorrect=self.par['UVIS']['polycorrect'],
                                      polyfunc=self.par['UVIS']['polyfunc'], debug=self.debug)

        self.norders, self.nspec = out_table['SENS_ZEROPOINT'].shape

        self.sens = self.__class__.empty_sensfunc_table(self.norders, self.nspec)

        self.stdname = meta_table['STD_NAME'][0]
        self.stdfile = meta_table['CAL_FILE'][0]
        self.std_ra = meta_table['STD_RA'][0]
        self.std_dec = meta_table['STD_DEC'][0]
        self.airmass = meta_table['AIRMASS'][0]
        self.exptime = meta_table['EXPTIME'][0]

        self.sens['WAVE'] = out_table['SENS_WAVE']
        self.sens['ZEROPOINT'] = out_table['SENS_ZEROPOINT']
        self.sens['ZEROPOINT_GPM'] = out_table['SENS_ZEROPOINT_GPM']
        self.sens['ZEROPOINT_FIT'] = out_table['SENS_ZEROPOINT_FIT']
        self.sens['ZEROPOINT_FIT_GPM'] = out_table['SENS_ZEROPOINT_FIT_GPM']
        self.sens['ECH_ORDER'] = meta_table['ECH_ORDERS'][0]
        self.sens['WAVE_MIN'] = out_table['WAVE_MIN']
        self.sens['WAVE_MAX'] = out_table['WAVE_MAX']

#        self.steps.append(inspect.stack()[0][3])

    def eval_zeropoint(self, wave, iorddet):
        """

        Parameters
        ----------
        wave: ndarray shape (nspec)
           Wavelength array

        iorddet: int
           Order or detector

        Returns
        -------
        zeropoint: ndarray, shape (nspec,)

        """
        # This routine can extrapolate
        return scipy.interpolate.interp1d(self.sens['WAVE']wave[iorddet],
                                          self.sens['ZEROPOINT_FIT'][iorddet],
                                          bounds_error=False, fill_value='extrapolate')(wave)



