# Module for flux calibrating spectra
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits

from pypeit import msgs
from pypeit.core import flux_calib
from pypeit.core import load
from pypeit.core import save
from pypeit import sensfunc
from pypeit import specobjs
from astropy import table
from pypeit import debugger

from IPython import embed



class FluxCalibrate(object):

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec1dfiles, sensfiles, par, debug=False):
        meta_table = table.Table.read(sensfiles[0], hdu=1)
        algorithm = meta_table['ALGORITHM']
        return next(c for c in cls.__subclasses__() if c.__name__ == algorithm)(spec1dfiles, sensfiles, par, debug=debug)

    def __init__(self, spec1dfiles, sensfiles, par=None, debug=False):

        self.spec1dfiles = spec1dfiles
        self.sensfiles = sensfiles
        self.par = par
        self.debug = debug


        # Read in the data
        sobjs = (specobjs.SpecObjs.from_fitsfile(self.spec1dfile)).get_std()
        # Put spectrograph info into meta
        self.wave, self.counts, self.counts_ivar, self.counts_mask, self.meta_spec, header = sobjs.unpack_object(ret_flam=False)
        # Set spectrograph
        self.spectrograph = load_spectrograph(self.meta_spec['PYP_SPEC'])

        # If the user provided RA and DEC use those instead of what is in meta
        star_ra = self.meta_spec['RA'] if self.par['star_ra'] is None else self.par['star_ra']
        star_dec = self.meta_spec['DEC'] if self.par['star_dec'] is None else self.par['star_dec']
        # Read in standard star dictionary
        self.std_dict = flux_calib.get_standard_spectrum(star_type=self.par['star_type'], star_mag=self.par['star_mag'],
                                                         ra=star_ra, dec=star_dec)