

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
from pypeit.spectrograph.util import load_spectrograph
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
    def __init__(self, spec1dfile, sensfile, par=None, star_type=None, star_mag=None, star_ra=None, star_dec=None, debug=False):

        self.spec1dfile = spec1dfile
        self.sensfile = sensfile
        self.star_type = star_type
        self.star_mag = star_mag
        self.star_ra = star_ra
        self.star_dec = star_dec
        self.debug = debug

        # Read in the Standard star data
        sobjs_std = (specobjs.SpecObjs.from_fitsfile(self.spec1dfile)).get_std()
        # Put spectrograph info into meta
        self.wave, self.counts, self.counts_ivar, self.counts_mask, self.meta_spec, header = sobjs_std.unpack_object(ret_flam=False)
        # Set spectrograph
        self.spectrograph = load_spectrograph(self.meta_spec['PYP_SPEC'])

        # If the user provided RA and DEC use those instead of what is in meta
        star_ra = self.meta_spec['RA'] if self.star_ra is None else self.star_ra
        star_dec = self.meta_spec['DEC'] if self.star_dec is None else self.star_dec
        # Read in standard star dictionary
        self.std_dict = flux_calib.get_standard_spectrum(star_type=self.star_type, star_mag=self.star_mag, ra=star_ra, dec=star_dec)


    def generate_sensfunc(self):
        pass

    def save(self):
        pass

    def load(self, sensfile):
        pass

    def show(self):
        pass



class SensFuncIR(SensFunc):
    def __init__(self, spec1dfile, sensfile, par=None, star_type=None, star_mag=None, star_ra=None, star_dec=None, debug=False,
                 debug_init=False, telgridfile=None, only_orders=None):
        super().__init__(spec1dfile, sensfile, par=par, star_type=star_type, star_mag=star_mag, star_ra=star_ra, star_dec=star_dec, debug=debug)


        # Set parset
        self.par = pypeitpar.TelluricPar() if par is None else par
        self.only_orders = only_orders
        self.debug_init = debug_init
        if telgridfile is None:
            self.telgridfile = self.spectrograph.telluric_grid_file
        else:
            self.telgridfile=telgridfile
        self.TelObj = None

    def generate_sensfunc(self):

        self.TelObj = telluric.sensfunc_telluric(self.wave, self.counts, self.counts_ivar, self.counts_mask,
                                                 self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'],
                                                 self.std_dict, self.telgridfile,
                                                 polyorder=self.par['polyorder'],
                                                 mask_abs_lines=self.par['mask_abs_lines'],
                                                 delta_coeff_bounds=self.par['delta_coeff_bounds'],
                                                 minmax_coeff_bounds=self.par['min_max_coeff_bounds'],
                                                 sn_clip=self.par['sn_clip'], only_orders=self.only_orders,
                                                 tol=self.par['tol'], popsize=self.par['popsize'],
                                                 recombination=self.par['recombination'], polish=self.par['polish'],
                                                 disp=self.par['disp'], debug_init=self.debug_init, debug=self.debug)

        self.meta_table, self.out_table = self.TelObj.meta_table, self.TelObj.meta_table
        return self.meta_table, self.out_table

    def save(self):
        self.TelObj.save(self.sensfile)

    def load(self, sensfile):
        self.TelObj.load(sensfile)



class SensFuncUVIS(SensFunc):
    def __init__(self, spec1dfile, sensfile, par=None, star_type=None, star_mag=None, star_ra=None, star_dec=None, debug=False):
        super().__init__(spec1dfile, sensfile, par=par, star_type=star_type, star_mag=star_mag, star_ra=star_ra, star_dec=star_dec, debug=debug)

        # Set parset
        self.par = pypeitpar.SensfuncPar() if par is None else par

    def generate_sensfunc(self):
        self.meta_table, self.out_table = flux_calib.sensfunc(self.wave, self.counts, self.counts_ivar, self.counts_mask,
                                                              self.meta_spec['EXPTIME'], self.meta_spec['AIRMASS'], self.std_dict,
                                                              self.meta_spec['LONGITUDE'], self.meta_spec['LATITUDE'],
                                                              telluric=False, polyorder=self.par['polyorder'],
                                                              balm_mask_wid=self.par['balm_mask_wid'],
                                                              nresln=self.par['nresln'], resolution=self.par['resolution'],
                                                              trans_thresh=self.par['trans_thresh'],
                                                              polycorrect=True, debug=self.debug)

        return self.meta_table, self.out_table

    def save(self):

        # Write to outfile
        msgs.info('Writing object and telluric models to file: {:}'.format(self.sensfile))
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



