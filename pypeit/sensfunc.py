

import os
import inspect
import numpy as np
import scipy
from IPython import embed
import inspect

from pypeit import msgs
from pypeit import ginga
from pypeit import masterframe
from pypeit import specobjs
from pypeit.core import flux_calib
from pypeit.core import telluric
from pypeit.spectrographs.util import load_spectrograph
from astropy.io import fits
from astropy import table
from pypeit.core import coadd1d
from pypeit.core.wavecal import wvutils
from pypeit import utils
from pypeit.io import initialize_header


# TODO Add the data model up here as a standard thing, which is an astropy table.
# data model needs a tag on whether its merged or not. For merged specobjs, you only apply sensfunc directly coefficients
# are then nonsense. homogenize data model to be the same for both algorithms.

#TODO Should this be a master frame? I think not.
#TODO Standard output location for sensfunc?
#TODO What if the user wants to coadd exposures first for say a multi-detector instrument where detector runs in the wavelength
# dimension, then compute a single sensfunc, rather than one for each detector? In that case this routine would need to run
# on co-added data written in a different format? That is easy enough to do, but we need to require that that file
# have certain meta data in its header.


# TODO How do we deal with cases of multiple-detectors where detector runs in the spectral direction? This code currently
# knows nothing about detectors. This is a tricky case. The options are:
# 1) Flux calibrate each detector separately, and then apply detector specific sensfuncs. This works for longslit, but
# but will fail for multislit, since different slits have different wavelength coverage.
# 2) Co-add the standard star spectra in counts, and compute one global sensfunc for the data. This would work well for
# multislit, but it will fail in cases where the detectors have different througphut, since there will be jumps in the
# co-added spectrum across detector boundaries, making the sensunc discontinouus.

# TODO Define sensfunc data model here.



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
    def get_instance(cls, spec1dfile, sensfile, par, debug=False):
        return next(c for c in cls.__subclasses__() if c.__name__ == par['algorithm'])(spec1dfile, sensfile, par, debug=debug)

    @classmethod
    def load(cls, sensfile):
        # Write to outfile
        msgs.info('Reading sensitivity function from file: {:}'.format(sensfile))
        hdulist = fits.open(sensfile)
        header = hdulist[0].header
        wave = hdulist['WAVE'].data
        sensfunc = hdulist['SENSFUNC'].data
        meta_table = table.Table(hdulist['METADATA'].data)
        out_table  = table.Table(hdulist['OUT_TABLE'].data)
        return wave, sensfunc, meta_table, out_table, header

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):
        # Arguments
        self.spec1dfile = spec1dfile
        self.sensfile = sensfile
        self.par = par
        self.debug = debug
        # Core attributes that will be output to file
        self.meta_table = None
        self.out_table = None
        self.wave = None
        self.sensfunc = None
        self.steps = []
        # Are we splicing together multiple detectors?
        self.splice_multi_det = True if self.par['multi_spec_det'] is not None else False

        # Read in the Standard star data
        sobjs_std = (specobjs.SpecObjs.from_fitsfile(self.spec1dfile)).get_std(multi_spec_det=self.par['multi_spec_det'])
        # Put spectrograph info into meta
        self.wave, self.counts, self.counts_ivar, self.counts_mask, self.meta_spec, header = sobjs_std.unpack_object(ret_flam=False)
        self.norderdet = 1 if self.wave.ndim == 1 else self.wave.shape[1]
        # Set spectrograph
        self.spectrograph = load_spectrograph(self.meta_spec['PYP_SPEC'])

        # If the user provided RA and DEC use those instead of what is in meta
        star_ra = self.meta_spec['RA'] if self.par['star_ra'] is None else self.par['star_ra']
        star_dec = self.meta_spec['DEC'] if self.par['star_dec'] is None else self.par['star_dec']
        # Read in standard star dictionary
        self.std_dict = flux_calib.get_standard_spectrum(star_type=self.par['star_type'], star_mag=self.par['star_mag'],
                                                         ra=star_ra, dec=star_dec)

    @property
    def algorithm(self):
        return self.par['algorithm']


    def compute_sensfunc(self):
        """
        Dummy method overloaded by subclasses

        Returns:
            meta_table, out_table

        """
        pass
        return None, None

    def run(self):
        # Compute the sensitivity function
        self.meta_table, self.out_table = self.compute_sensfunc()
        # Extrapolate the sensfunc based on par['extrap_blu'], par['extrap_red']
        self.wave, self.sensfunc = self.extrapolate(samp_fact = self.par['samp_fact'])
        if self.splice_multi_det:
            self.wave, self.sensfunc = self.splice(self.wave)

        return

    def eval_sensfunc(self, wave, iorddet):
        """
        Dummy method, overloaded by subclasses

        Returns:
        """
        pass

    def save(self):

        # Write to outfile
        msgs.info('Writing sensitivity function results to file: {:}'.format(self.sensfile))

        # Standard init
        hdr = initialize_header()

        hdr['PYP_SPEC'] = (self.spectrograph.spectrograph, 'PypeIt: Spectrograph name')
        hdr['PYPELINE'] = self.spectrograph.pypeline
        #   - List the completed steps
        hdr['STEPS'] = (','.join(self.steps), 'Completed sensfunc steps')
        #   - Provide the file names
        hdr['SPC1DFIL'] = self.spec1dfile

        # Write the fits file
        data = [self.wave, self.sensfunc]
        extnames = ['WAVE', 'SENSFUNC']
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

        # Create a new set of oversampled and padded wavelength grids for the extrapolation
        wave_extrap_min = self.out_table['WAVE_MIN'].data * (1.0 - self.par['extrap_blu'])
        wave_extrap_max = self.out_table['WAVE_MAX'].data * (1.0 + self.par['extrap_red'])
        nspec_extrap = 0
        # Find the maximum size of the wavewlength grids, since we want everything to have the same
        for idet in range(self.norderdet):
            dwave_data, dloglam_data, resln_guess, pix_per_sigma = wvutils.get_sampling(self.out_table['WAVE'][idet])
            nspec_now = np.ceil(samp_fact * (wave_extrap_max[idet] - wave_extrap_min[idet]) / dwave_data).astype(int)
            nspec_extrap = np.max([nspec_now, nspec_extrap])
        # Create the wavelength grid
        wave_extrap = np.outer(np.arange(nspec_extrap), (wave_extrap_max - wave_extrap_min) / samp_fact / (nspec_extrap - 1)) + \
                      np.outer(np.ones(nspec_extrap), wave_extrap_min)
        sensfunc_extrap = np.zeros_like(wave_extrap)
        # Evaluate extrapolated sensfunc for all orders detectors
        for iorddet in range(self.norderdet):
            sensfunc_extrap[:, iorddet] = self.eval_sensfunc(wave_extrap[:,iorddet], iorddet)

        self.steps.append(inspect.stack()[0][3])

        return wave_extrap, sensfunc_extrap

    def splice(self, wave):

        msgs.info('Merging sensfunc for {:d} detectors {:}'.format(self.norderdet, self.par['multi_spec_det']))
        wave_splice_min = wave.min()
        wave_splice_max = wave.max()
        wave_splice, _, _ = coadd1d.get_wave_grid(wave, wave_method='linear', wave_grid_min=wave_splice_min,
                                                  wave_grid_max=wave_splice_max, samp_fact=1.0)
        sensfunc_splice = np.zeros_like(wave_splice)
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
            splice_wave_mask = (wave_splice >= wave_mask_min) & (wave_splice <= wave_mask_max)
            sensfunc_splice[splice_wave_mask] = self.eval_sensfunc(wave_splice[splice_wave_mask], idet)

        self.steps.append(inspect.stack()[0][3])

        return wave_splice, sensfunc_splice

    def show(self):
        pass

# TODO Add a method which optionally merges sensfunc using the nsens > 1 logic



class IR(SensFunc):

    def __init__(self, spec1dfile, sensfile, par, debug=False):
        super().__init__(spec1dfile, sensfile, par, debug=debug)

        self.TelObj = None

    def compute_sensfunc(self):

        meta_table, out_table = telluric.sensfunc_telluric(
            self.wave, self.counts, self.counts_ivar, self.counts_mask, self.meta_spec['EXPTIME'],
            self.meta_spec['AIRMASS'], self.std_dict, self.par['IR']['telgridfile'], polyorder=self.par['polyorder'],
            ech_orders = self.meta_spec['ECH_ORDERS'],
            sn_clip=self.par['IR']['sn_clip'], mask_abs_lines=self.par['mask_abs_lines'],
            # JFH Implement thease in parset?
            #delta_coeff_bounds=self.par['IR']['delta_coeff_bounds'],
            #minmax_coeff_bounds=self.par['IR']['min_max_coeff_bounds'],
            tol=self.par['IR']['tol'], popsize=self.par['IR']['popsize'], recombination=self.par['IR']['recombination'],
            polish=self.par['IR']['polish'],
            disp=self.par['IR']['disp'], debug=self.debug)
        # Add the algorithm to the meta_table
        meta_table['ALGORITHM'] = self.par['algorithm']
        self.steps.append(inspect.stack()[0][3])

        return meta_table, out_table

    def eval_sensfunc(self, wave, iorddet):
        # Put this stuff in a function called eval_sensfunc for each algorithm
        wave_min = self.out_table[iorddet]['WAVE_MIN']
        wave_max = self.out_table[iorddet]['WAVE_MAX']
        polyorder_vec = self.meta_table['POLYORDER_VEC'][0]
        func = self.meta_table['FUNC'][0]
        coeff = self.out_table[iorddet]['OBJ_THETA'][0:polyorder_vec[iorddet] + 2]
        sensfunc = np.exp(utils.func_val(coeff, wave, func, minx=wave_min, maxx=wave_max))
        return sensfunc


class UVIS(SensFunc):
    def __init__(self, spec1dfile, sensfile, par, debug=False):
        super().__init__(spec1dfile, sensfile, par, debug=debug)

        # Add some cards to the meta spec. These should maybe just be added already in unpack object
        self.meta_spec['LATITUDE'] = self.spectrograph.telescope['latitude']
        self.meta_spec['LONGITUDE'] = self.spectrograph.telescope['longitude']


    def compute_sensfunc(self):

        meta_table, out_table = flux_calib.sensfunc(self.wave, self.counts, self.counts_ivar, self.counts_mask,
                                                              self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'], self.std_dict,
                                                              self.meta_spec['LONGITUDE'], self.meta_spec['LATITUDE'],
                                                              self.meta_spec['ECH_ORDERS'], 
                                                              telluric=False, polyorder=self.par['polyorder'],
                                                              balm_mask_wid=self.par['UVIS']['balm_mask_wid'],
                                                              nresln=self.par['UVIS']['nresln'],
                                                              resolution=self.par['UVIS']['resolution'],
                                                              trans_thresh=self.par['UVIS']['trans_thresh'],
                                                              polycorrect=True, debug=self.debug)
        # Add the algorithm to the meta_table
        meta_table['ALGORITHM'] = self.par['algorithm']

        self.steps.append(inspect.stack()[0][3])

        return meta_table, out_table


    def eval_sensfunc(self, wave, iorddet):
        # This routine can extrapolate
        sensfunc = scipy.interpolate.interp1d(self.out_table['WAVE'][iorddet,:], self.out_table['SENSFUNC'][iorddet,:],
                                              bounds_error = False, fill_value='extrapolate')(wave)
        return sensfunc
