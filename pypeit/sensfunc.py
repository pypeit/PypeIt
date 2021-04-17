

import os
import inspect
import numpy as np
import scipy
from IPython import embed
import inspect

from matplotlib import pyplot as plt

from pypeit import msgs
from pypeit import specobjs
from pypeit.core import flux_calib
from pypeit.core import telluric, fitting
from pypeit.spectrographs.util import load_spectrograph
from astropy.io import fits
from astropy import table
from pypeit.core import coadd
from pypeit.core.wavecal import wvutils
from pypeit import utils
from pypeit import io
from pypeit.core import meta
from matplotlib.backends.backend_pdf import PdfPages


# TODO Add the data model up here as a standard thing using DataContainer.

#TODO Should this be a master frame? I think not.
#TODO Standard output location for sensfunc?

# TODO Add some QA plots, and plots to the screen if show is set.


class SensFunc(object):
    """
    Class to generate sensitivity function from a standard star spectrum.

    Args:
        spec1dfile (str):
            PypeIt spec1d file for the standard file.
        sensfile (str):
            File to write sensitivity function to.
        par (parset object):
            Parset containing parameters for sensitivity function computation.
        multi_spec_det (list):
            List of detectors for which to create a spliced sensitivity function.  If passed the sensitivity function
            code will merge together the sensitivity functions for the specified list of detectors. This option is
            required for instruments which have multiple detectors arranged in the spectral direction.
        debug:
            Run the sensitivity function codes in debug mode sending diagnostic information to the screen.
    """

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec1dfile, sensfile, par=None, debug=False):
        return next(c for c in cls.__subclasses__() if c.__name__ == par['algorithm'])(
            spec1dfile, sensfile, par=par, debug=debug)

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

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):
        """
        See docs above

        """
        # Arguments
        self.spec1dfile = spec1dfile
        self.sensfile = sensfile
        # Set spectrograph
        header = fits.getheader(self.spec1dfile)
        self.spectrograph = load_spectrograph(header['PYP_SPEC'])
        # TODO This line is necessary until we figure out a way to instantiate spectrograph objects with configuration
        # specific information from spec1d files.
        self.spectrograph.dispname = header['DISPNAME']
        self.par = self.spectrograph.default_pypeit_par()['sensfunc'] if par is None else par
        # QA and throughtput plot filenames
        self.qafile = sensfile.replace('.fits', '') + '_QA.pdf'
        self.thrufile = sensfile.replace('.fits', '') + '_throughput.pdf'
        self.debug = debug
        # Intermediate attributes
        self.wave = None
        self.counts = None
        self.conts_ivar = None
        self.counts_mask = None
        self.norderdet = None
        # Core attributes that will be output to file
        self.meta_table = None
        self.out_table = None
        self.wave_zp = None
        self.zeropoint = None
        self.throughput = None
        self.wave_zp_splice = None
        self.zeropoint_splice = None
        self.throughput_splice = None
        self.steps = []
        # Are we splicing together multiple detectors?
        self.splice_multi_det = True if self.par['multi_spec_det'] is not None else False

        # Read in the Standard star data
        sobjs_std = (specobjs.SpecObjs.from_fitsfile(self.spec1dfile)).get_std(multi_spec_det=self.par['multi_spec_det'])
        # Unpack standard
        wave, counts, counts_ivar, counts_mask, self.meta_spec, header = sobjs_std.unpack_object(ret_flam=False)
        # Perform any instrument tweaks
        wave_twk, counts_twk, counts_ivar_twk, counts_mask_twk = self.spectrograph.tweak_standard(
            wave, counts, counts_ivar, counts_mask, self.meta_spec)
        # Reshape to 2d arrays
        self.wave, self.counts, self.counts_ivar, self.counts_mask, nspec_in, self.norderdet = \
            utils.spec_atleast_2d(wave_twk, counts_twk, counts_ivar_twk, counts_mask_twk)

        # If the user provided RA and DEC use those instead of what is in meta
        star_ra = self.meta_spec['RA'] if self.par['star_ra'] is None else self.par['star_ra']
        star_dec = self.meta_spec['DEC'] if self.par['star_dec'] is None else self.par['star_dec']
        star_ra, star_dec = meta.convert_radec(star_ra, star_dec)  # Convert to decimal deg, as need be

        # Read in standard star dictionary
        self.std_dict = flux_calib.get_standard_spectrum(star_type=self.par['star_type'], star_mag=self.par['star_mag'],
                                                         ra=star_ra, dec=star_dec)

    @property
    def algorithm(self):
        return self.par['algorithm']


    def compute_zeropoint(self):
        """
        Dummy method overloaded by subclasses

        Returns:
            meta_table, out_table

        """
        pass
        return None, None

    def run(self):
        # Compute the sensitivity function
        self.meta_table, self.out_table = self.compute_zeropoint()
        # Extrapolate the zeropoint based on par['extrap_blu'], par['extrap_red']
        self.wave_zp, self.zeropoint = self.extrapolate(samp_fact = self.par['samp_fact'])
        if self.splice_multi_det:
            self.wave_zp_splice, self.zeropoint_splice = self.splice()

        # Compute the throughput
        self.throughput, self.throughput_splice = self.compute_throughput()

        # Write out QA and throughput plots
        self.write_QA()

        # TODO This is ugly code. Maybe it should go in save. Maybe we should just output all the zeropoint arrays
        # as 2d arrays rather than flattening them.
        # If the zeropoint has just one order, or detectors were spliced, flatten the output
        #if self.wave_zp.ndim == 2:
        #    # Always flatten for multi_spec_det
        #    if self.splice_multi_det:
        #        self.wave_zp_splice = self.wave_zp_splice.flatten()
        #        self.zeropoint_splice = self.zeropoint_splice.flatten()
        #        self.throughput_splice = self.throughput_splice.flatten()
        #    # Otherwise, flatten only if the number of det/orders = 1
        #    elif self.wave_zp.shape[1] == 1:
        #        self.wave_zp = self.wave_zp.flatten()
        #        self.zeropoint = self.zeropoint.flatten()
        #        self.throughput = self.throughput.flatten()

        return

    def eval_zeropoint(self, wave, iorddet):
        """
        Dummy method, overloaded by subclasses

        Returns:
        """
        pass

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
        if self.splice_multi_det:
            data = [self.wave_zp_splice, self.zeropoint_splice, self.throughput_splice]
        else:
            data = [self.wave_zp, self.zeropoint, self.throughput]

        extnames = ['WAVE', 'ZEROPOINT', 'THROUGHPUT']
        # Write the fits file
        hdulist = fits.HDUList([fits.PrimaryHDU(header=hdr)] + [fits.ImageHDU(data=d, name=n) for d, n in zip(data, extnames)])
        hdu_meta = fits.table_to_hdu(self.meta_table)
        hdu_meta.name = 'METADATA'
        hdu_out = fits.table_to_hdu(self.out_table)
        hdu_out.name = 'OUT_TABLE'
        hdulist.append(hdu_meta)
        hdulist.append(hdu_out)
        hdulist.writeto(self.sensfile, overwrite=True, checksum=True)



    def extrapolate(self, samp_fact=1.5):
        """
        Extrapolates the sensitivity function to cover an extra wavelength range set by the extrapl_blu extrap_red
        parameters. This is important for making sure that the sensfunc can be applied to data with slightly different
        wavelength coverage etc.

        Parameters
        ----------
        samp_fact: float
            Parameter governing the sampling of the wavelength grid used for the extrapolation.

        Returns
        -------
        wave_extrap: ndarray
            Extrapolated wavelength array
        zeropoint_extrap: ndarray
            Extrapolated sensfunc

        """

        # Create a new set of oversampled and padded wavelength grids for the extrapolation
        wave_extrap_min = self.out_table['WAVE_MIN'].data * (1.0 - self.par['extrap_blu'])
        wave_extrap_max = self.out_table['WAVE_MAX'].data * (1.0 + self.par['extrap_red'])
        nspec_extrap = 0
        # Find the maximum size of the wavewlength grids, since we want everything to have the same
        for idet in range(self.norderdet):
            wave = self.wave if self.wave.ndim == 1 else self.wave[:, idet]
            dwave_data, dloglam_data, resln_guess, pix_per_sigma = wvutils.get_sampling(wave)
            nspec_now = np.ceil(samp_fact * (wave_extrap_max[idet] - wave_extrap_min[idet]) / dwave_data).astype(int)
            nspec_extrap = np.max([nspec_now, nspec_extrap])
        # Create the wavelength grid
        wave_extrap = np.outer(np.arange(nspec_extrap), (wave_extrap_max - wave_extrap_min)/ (nspec_extrap - 1)) + \
                      np.outer(np.ones(nspec_extrap), wave_extrap_min)
        zeropoint_extrap = np.zeros_like(wave_extrap)
        # Evaluate extrapolated zerpoint for all orders detectors
        for iorddet in range(self.norderdet):
            zeropoint_extrap[:, iorddet] = self.eval_zeropoint(wave_extrap[:,iorddet], iorddet)

        self.steps.append(inspect.stack()[0][3])
        return wave_extrap, zeropoint_extrap

    def splice(self):
        """
        Routine to splice together sensitivity functions into one global sensitivity function for spectrographs
        with multiple detectors extending across the wavelength direction.

        Parameters
        ----------
        wave: ndarray, shape (nspec, norddet)

        Returns
        -------
        wave_splice: ndarray, shape (nspec_splice, 1)
        zeropoint_splice: ndarray, shape (nspec_splice, 1)


        """

        msgs.info('Merging sensfunc for {:d} detectors {:}'.format(self.norderdet, self.par['multi_spec_det']))
        wave_splice_min = self.wave_zp[self.wave_zp > 1.0].min()
        wave_splice_max = self.wave_zp[self.wave_zp > 1.0].max()
        wave_splice_1d, _, _ = coadd.get_wave_grid(self.wave_zp, wave_method='linear', wave_grid_min=wave_splice_min,
                                                wave_grid_max=wave_splice_max, spec_samp_fact=1.0)
        zeropoint_splice_1d = np.zeros_like(wave_splice_1d)
        for idet in range(self.norderdet):
            wave_min = self.out_table['WAVE_MIN'][idet]
            wave_max = self.out_table['WAVE_MAX'][idet]
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
            zeropoint_splice_1d[splice_wave_mask] = self.eval_zeropoint(wave_splice_1d[splice_wave_mask], idet)

        # Interpolate over gaps
        zeros = zeropoint_splice_1d == 0.
        if np.any(zeros):
            msgs.info("Interpolating over gaps (and extrapolating with fill_value=1, if need be)")
            interp_func = scipy.interpolate.interp1d(wave_splice_1d[np.invert(zeros)],
                                                 zeropoint_splice_1d[np.invert(zeros)],
                                                 kind='nearest', fill_value=0., bounds_error=False) #
            #kind='nearest', fill_value='extrapoloate', bounds_error=False) #  extrapolate fails for JXP, even on 1.4.1
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
        throughput (`np.ndarray`):
           Througphput, float shape (nspec, norders)

        """

        # Set the throughput to be -1 in places where it is not defined.
        throughput = np.full_like(self.zeropoint, -1.0)
        for idet in range(self.wave_zp.shape[1]):
            wave_gpm =  (self.wave_zp[:,idet] >= self.out_table[idet]['WAVE_MIN']) & \
                        (self.wave_zp[:,idet] <= self.out_table[idet]['WAVE_MAX']) & (self.wave_zp[:,idet] > 1.0)
            throughput[:,idet][wave_gpm] = flux_calib.zeropoint_to_throughput(
                self.wave_zp[:,idet][wave_gpm], self.zeropoint[:,idet][wave_gpm], self.spectrograph.telescope.eff_aperture())
        if self.splice_multi_det:
            wave_gpm = (self.wave_zp_splice >= self.out_table[:]['WAVE_MIN'].min()) & \
                        (self.wave_zp_splice <= self.out_table[:]['WAVE_MAX'].max()) & (self.wave_zp_splice > 1.0)
            throughput_splice = np.zeros_like(self.wave_zp_splice)
            throughput_splice[wave_gpm] = flux_calib.zeropoint_to_throughput(
                self.wave_zp_splice[wave_gpm], self.zeropoint_splice[wave_gpm], self.spectrograph.telescope.eff_aperture())
        else:
            throughput_splice = None

        return throughput, throughput_splice


    def write_QA(self):
        """
        Write out zeropoint QA files

        Returns
        -------

        """
        utils.pyplot_rcparams()

        # Plot QA for zeropoint
        if 'Echelle' in self.spectrograph.pypeline:
            order_or_det = self.spectrograph.orders[np.arange(self.norderdet)]
            order_or_det_str = 'order'
        else:
            order_or_det = np.arange(self.norderdet) + 1
            order_or_det_str = 'det'

        spec_str = ' {:s} {:s} {:s} '.format(self.spectrograph.name, self.spectrograph.pypeline, self.spectrograph.dispname)
        zp_title = ['PypeIt Zeropoint QA for' + spec_str + order_or_det_str +'={:d}'.format(order_or_det[idet]) for idet in range(self.norderdet)]
        thru_title = [order_or_det_str + '={:d}'.format(order_or_det[idet]) for idet in range(self.norderdet)]

        is_odd = self.norderdet % 2 != 0
        npages = int(np.ceil(self.norderdet/2)) if is_odd else self.norderdet//2 + 1
        # TODO PDF page logic is a bit complicated becauase we want to plot two plots per page, but the number of pages
        #  depends on the number of order/det. Consider just dumping out a set of plots or revamp once we have a dashboard
        with PdfPages(self.qafile) as pdf:
            for ipage in range(npages):
                figure, (ax1, ax2) = plt.subplots(2, figsize=(8.27, 11.69))
                if (2 * ipage) < self.norderdet:
                    flux_calib.zeropoint_qa_plot(
                        self.out_table[2*ipage]['SENS_WAVE'], self.out_table[2*ipage]['SENS_ZEROPOINT'],
                        self.out_table[2*ipage]['SENS_ZEROPOINT_GPM'], self.out_table[2*ipage]['SENS_ZEROPOINT_FIT'],
                        self.out_table[2*ipage]['SENS_ZEROPOINT_FIT_GPM'], title=zp_title[2*ipage], axis=ax1)
                if (2*ipage + 1) < self.norderdet:
                    flux_calib.zeropoint_qa_plot(
                        self.out_table[2*ipage+1]['SENS_WAVE'], self.out_table[2*ipage+1]['SENS_ZEROPOINT'],
                        self.out_table[2*ipage+1]['SENS_ZEROPOINT_GPM'], self.out_table[2*ipage+1]['SENS_ZEROPOINT_FIT'],
                        self.out_table[2*ipage+1]['SENS_ZEROPOINT_FIT_GPM'], title=zp_title[2*ipage+1], axis=ax2)
                if self.norderdet == 1:
                    # For single order/det just finish up after the first page
                    ax2.remove()
                    pdf.savefig()
                    plt.close('all')
                elif (self.norderdet > 1) & (ipage < npages-1):
                    # For multi order/det but not on the last page, finish up. No need to remove ax2 since there
                    # are always 2 plots per page except on the last page
                    pdf.savefig()
                    plt.close('all')
                else:
                    # For multi order/det but on the last page, add order/det summary plot to the final page
                    # Deal with even/odd page logic for axes
                    if is_odd:
                        # add order/det summary plot to axis 2 of current page
                        axis=ax2
                    else:
                        axis=ax1
                        ax2.remove()
                    for idet in range(self.norderdet):
                        # define the color
                        rr = (np.max(order_or_det) - order_or_det[idet]) / np.maximum(
                            np.max(order_or_det) - np.min(order_or_det), 1)
                        gg = 0.0
                        bb = (order_or_det[idet] - np.min(order_or_det)) / np.maximum(
                            np.max(order_or_det) - np.min(order_or_det), 1)
                        wave_gpm = self.out_table[idet]['SENS_WAVE'] > 1.0
                        axis.plot(self.out_table[idet]['SENS_WAVE'][wave_gpm],
                                 self.out_table[idet]['SENS_ZEROPOINT_FIT'][wave_gpm],
                                 color=(rr, gg, bb), linestyle='-', linewidth=2.5, label=thru_title[idet],
                                 zorder=5 * idet)
                    # If we are splicing, overplot the spliced zeropoint
                    if self.splice_multi_det:
                        wave_zp_gpm = (self.wave_zp_splice >= self.out_table['WAVE_MIN'].min()) & \
                                      (self.wave_zp_splice <= self.out_table['WAVE_MAX'].max()) & (self.wave_zp_splice > 1.0)
                        axis.plot(self.wave_zp_splice[wave_zp_gpm].flatten(), self.zeropoint_splice[wave_zp_gpm].flatten(),
                                 color='black', linestyle='-', linewidth=2.5, label='Spliced Zeropoint', zorder=30, alpha=0.3)

                    axis.set_xlim((0.98 * self.out_table['WAVE_MIN'].min(), 1.02 * self.out_table['WAVE_MAX'].max()))
                    axis.set_ylim((0.95*self.out_table['SENS_ZEROPOINT_FIT'][self.out_table['SENS_WAVE'] > 1.0].min(),
                                   1.05*self.out_table['SENS_ZEROPOINT_FIT'][self.out_table['SENS_WAVE'] > 1.0].max()))
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
            rr = (np.max(order_or_det) - order_or_det[idet])/np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
            gg = 0.0
            bb = (order_or_det[idet] - np.min(order_or_det))/np.maximum(np.max(order_or_det) - np.min(order_or_det), 1)
            gpm = (self.throughput[:, idet] >= 0.0)
            axis.plot(self.wave_zp[gpm,idet], self.throughput[gpm,idet], color=(rr, gg, bb), linestyle='-', linewidth=2.5, label=thru_title[idet], zorder=5*idet)
        if self.splice_multi_det:
            axis.plot(self.wave_zp_splice[wave_zp_gpm].flatten(), self.throughput_splice[wave_zp_gpm].flatten(),
                      color='black', linestyle='-', linewidth=2.5, label='Spliced Throughput', zorder=30, alpha=0.3)

        axis.set_xlim((0.98*self.wave_zp[self.throughput >=0.0].min(),1.02*self.wave_zp[self.throughput >=0.0].max()))
        axis.set_ylim((0.0,1.05*self.throughput[self.throughput >=0.0].max()))
        axis.legend()
        axis.set_xlabel('Wavelength (Angstroms)')
        axis.set_ylabel('Throughput')
        axis.set_title('PypeIt Throughput for' + spec_str)
        fig.savefig(self.thrufile)


# TODO Add a method which optionally merges sensfunc using the nsens > 1 logic



class IR(SensFunc):

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):
        super().__init__(spec1dfile, sensfile, par=par, debug=debug)

        self.TelObj = None

    def compute_zeropoint(self):
        """
        Calls routine to compute the sensitivity function.

        Returns
        -------
        meta_table: astropy table
               Table containing zeropoint meta data
        out_table: astropy table
               Table containing zerpoint information.
        """

        meta_table, out_table = telluric.sensfunc_telluric(
            self.wave, self.counts, self.counts_ivar, self.counts_mask, self.meta_spec['EXPTIME'],
            self.meta_spec['AIRMASS'], self.std_dict, self.par['IR']['telgridfile'], polyorder=self.par['polyorder'],
            ech_orders = self.meta_spec['ECH_ORDERS'],
            sn_clip=self.par['IR']['sn_clip'], mask_abs_lines=self.par['mask_abs_lines'],
            maxiter=self.par['IR']['maxiter'],lower=self.par['IR']['lower'],upper=self.par['IR']['upper'],
            delta_coeff_bounds=self.par['IR']['delta_coeff_bounds'],
            minmax_coeff_bounds=self.par['IR']['minmax_coeff_bounds'],
            tol=self.par['IR']['tol'], popsize=self.par['IR']['popsize'], recombination=self.par['IR']['recombination'],
            polish=self.par['IR']['polish'],
            disp=self.par['IR']['disp'], debug=self.debug, debug_init=self.debug)
        # Add the algorithm to the meta_table
        meta_table['ALGORITHM'] = self.par['algorithm']
        self.steps.append(inspect.stack()[0][3])

        return meta_table, out_table

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

        # Put this stuff in a function called eval_zeropoint for each algorithm
        wave_min = self.out_table[iorddet]['WAVE_MIN']
        wave_max = self.out_table[iorddet]['WAVE_MAX']
        polyorder_vec = self.meta_table['POLYORDER_VEC'][0]
        func = self.meta_table['FUNC'][0]
        coeff = self.out_table[iorddet]['OBJ_THETA'][0:polyorder_vec[iorddet] + 2]
        zeropoint = fitting.evaluate_fit(coeff, func, wave, minx=wave_min, maxx=wave_max)
        return zeropoint


class UVIS(SensFunc):
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

        meta_table, out_table = flux_calib.sensfunc(self.wave, self.counts, self.counts_ivar, self.counts_mask,
                                                    self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'], self.std_dict,
                                                    self.meta_spec['LONGITUDE'], self.meta_spec['LATITUDE'],
                                                    self.meta_spec['ECH_ORDERS'], polyorder=self.par['polyorder'],
                                                    balm_mask_wid=self.par['UVIS']['balm_mask_wid'],
                                                    nresln=self.par['UVIS']['nresln'],
                                                    resolution=self.par['UVIS']['resolution'],
                                                    trans_thresh=self.par['UVIS']['trans_thresh'],
                                                    polycorrect=self.par['UVIS']['polycorrect'],
                                                    polyfunc=self.par['UVIS']['polyfunc'],
                                                    debug=self.debug)
        # Add the algorithm to the meta_table
        meta_table['ALGORITHM'] = self.par['algorithm']

        self.steps.append(inspect.stack()[0][3])

        return meta_table, out_table


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
        zeropoint = scipy.interpolate.interp1d(self.out_table['SENS_WAVE'][iorddet,:], self.out_table['SENS_ZEROPOINT_FIT'][iorddet,:],
                                              bounds_error = False, fill_value='extrapolate')(wave)
        return zeropoint
