

import os
import inspect
import numpy as np

from IPython import embed

from pypeit import msgs
from pypeit import ginga
from pypeit import masterframe
from pypeit import specobjs
from pypeit.core import flux_calib
from pypeit.core import telluric
from pypeit.spectrographs.util import load_spectrograph
from astropy.io import fits
from astropy import table

from pypeit.par import pypeitpar
from pypeit.core import save
from pypeit.core import load
from pypeit.core import pixels
from pypeit.core import procimg

#TODO Should this be a master frame? I think not.
#TODO Standard output location for sensfunc?

class SensFunc(object):

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec1dfile, sensfile, par, debug=False):
        return next(c for c in cls.__subclasses__() if c.__name__ == par['algorithm'])(spec1dfile, sensfile, par, debug=debug)

    def __init__(self, spec1dfile, sensfile, par=None, debug=False):

        self.spec1dfile = spec1dfile
        self.sensfile = sensfile
        self.par = par
        self.debug = debug

        # Core attributes that will be output to file
        self.meta_table = None
        self.out_table = None

        # Read in the Standard star data
        sobjs_std = (specobjs.SpecObjs.from_fitsfile(self.spec1dfile)).get_std()
        # Put spectrograph info into meta
        self.wave, self.counts, self.counts_ivar, self.counts_mask, self.meta_spec, header = sobjs_std.unpack_object(ret_flam=False)
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


    def generate_sensfunc(self):
        pass


    def save(self):

        # Write to outfile
        msgs.info('Writing sensitivity function results to file: {:}'.format(self.sensfile))
        hdu_meta = fits.table_to_hdu(self.meta_table)
        hdu_meta.name = 'METADATA'
        hdu_out = fits.table_to_hdu(self.out_table)
        hdu_out.name = 'OUT_TABLE'
        hdulist = fits.HDUList()
        hdulist.append(hdu_meta)
        hdulist.append(hdu_out)
        hdulist.writeto(self.sensfile, overwrite=True)

    def load(self, sensfile):
        # Write to outfile
        msgs.info('Reading object and telluric models from file: {:}'.format(sensfile))
        meta_table = table.Table.read(sensfile, hdu=1)
        out_table = table.Table.read(sensfile, hdu=2)

        return meta_table, out_table

    def show(self):
        pass





class IR(SensFunc):

    def __init__(self, spec1dfile, sensfile, par, debug=False):
        super().__init__(spec1dfile, sensfile, par, debug=debug)

        self.TelObj = None

    def generate_sensfunc(self):
        self.TelObj = telluric.sensfunc_telluric(self.wave, self.counts, self.counts_ivar, self.counts_mask,
                                                 self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'], self.std_dict,
                                                 self.par['IR']['telgridfile'],
                                                 polyorder=self.par['polyorder'],
                                                 sn_clip=self.par['IR']['sn_clip'],
                                                 mask_abs_lines=self.par['mask_abs_lines'],
                                                 # JFH Implement thease in parset?
                                                 #delta_coeff_bounds=self.par['IR']['delta_coeff_bounds'],
                                                 #minmax_coeff_bounds=self.par['IR']['min_max_coeff_bounds'],
                                                 tol=self.par['IR']['tol'], popsize=self.par['IR']['popsize'],
                                                 recombination=self.par['IR']['recombination'], polish=self.par['IR']['polish'],
                                                 disp=self.par['IR']['disp'], debug=self.debug)

        self.meta_table, self.out_table = self.TelObj.meta_table, self.TelObj.meta_table
        # Add the algorithm to the meta_table
        self.meta_table['ALGORITHM'] = self.par['algorithm']
        return self.meta_table, self.out_table



class UVIS(SensFunc):
    def __init__(self, spec1dfile, sensfile, par, debug=False):
        super().__init__(spec1dfile, sensfile, par, debug=debug)

        # Add some cards to the meta spec. These should maybe just be added already in unpack object
        self.meta_spec['LATITUDE'] = self.spectrograph.telescope['latitude']
        self.meta_spec['LONGITUDE'] = self.spectrograph.telescope['longitude']


    def generate_sensfunc(self):
        self.meta_table, self.out_table = flux_calib.sensfunc(self.wave, self.counts, self.counts_ivar, self.counts_mask,
                                                              self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'], self.std_dict,
                                                              self.meta_spec['LONGITUDE'], self.meta_spec['LATITUDE'],
                                                              telluric=False, polyorder=self.par['polyorder'],
                                                              balm_mask_wid=self.par['UVIS']['balm_mask_wid'],
                                                              nresln=self.par['UVIS']['nresln'],
                                                              resolution=self.par['UVIS']['resolution'],
                                                              trans_thresh=self.par['UVIS']['trans_thresh'],
                                                              polycorrect=True, debug=self.debug)
        # Add the algorithm to the meta_table
        self.meta_table['ALGORITHM'] = self.par['algorithm']

        return self.meta_table, self.out_table


