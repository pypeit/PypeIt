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
    def get_instance(cls, spec1dfiles, sensfiles, spectrograph, par, debug=False):
        return next(c for c in cls.__subclasses__() if c.__name__ == spectrograph.pypeline)(
            spec1dfiles, sensfiles, spectrograph, par, debug=debug)

    def __init__(self, spec1dfiles, sensfiles, spectrograph, par, debug=False):

        self.spec1dfiles = spec1dfiles
        self.sensfiles = sensfiles
        self.spectrograph = spectrograph
        self.par = par
        self.debug = debug


        # Read in the data
        sobjs = (specobjs.SpecObjs.from_fitsfile(self.spec1dfile)).get_std()
        # Put spectrograph info into meta
        self.wave, self.counts, self.counts_ivar, self.counts_mask, self.meta_spec, header = sobjs.unpack_object(ret_flam=False)
