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

        sens_last = None
        for spec1, sens in zip(self.spec1dfiles,self.sensfiles):
            # Read in the data
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1)
            if sens != sens_last:
                wave, sensfunction, meta_table, out_table, header_sens = sensfunc.SensFunc.load(sens)
            self.flux_calib(sobjs, wave, sensfunction, meta_table)
            sobjs.write_to_fits(sobjs.header, spec1, overwrite=True)

    def flux_calib(self, sobjs, wave, sensfunction, meta_table):
        """
        Dummy method overloaded by subclass

        Args:
            sobjs:
            wave:
            sensfunction:
            meta_table:

        Returns:

        """
        pass

class MultiSlit(FluxCalibrate):
    """
    Child of FluxSpec for Multislit and Longslit reductions
    """

    def __init__(self, spec1dfiles, sensfiles, spectrograph, par, debug=False):
        super().__init__(spec1dfiles, sensfiles, spectrograph, par, debug=debug)


    def flux_calib(self, sobjs, wave, sensfunction, meta_table):
        """
        Apply sensitivity function to all the spectra in an sobjs object.

        Args:
            sobjs (object):
               SpecObjs object
            wave (ndarray):
               wavelength array for sensitivity function (nspec,)
            sensfunction (ndarray):
               sensitivity function
            meta_table (table):
               astropy table containing meta data for sensitivity function

        Returns:

        """

        # Run
        for sci_obj in sobjs:
            sci_obj.apply_flux_calib(wave, sensfunction,
                                     sobjs.header['EXPTIME'],
                                     extinct_correct=self.par['extinct_correct'],
                                     longitude=self.spectrograph.telescope['longitude'],
                                     latitude=self.spectrograph.telescope['latitude'],
                                     airmass=float(sobjs.header['AIRMASS']))




class Echelle(FluxCalibrate):
    """
    Child of FluxSpec for Echelle reductions
    """

    def __init__(self, spec1dfiles, sensfiles, spectrograph, par, debug=False):
        super().__init__(spec1dfiles, sensfiles, spectrograph, par, debug=debug)


    def flux_calib(self, sobjs, wave, sensfunction, meta_table):
        """
        Apply sensitivity function to all the spectra in an sobjs object.

        Args:
            sobjs (object):
               SpecObjs object
            wave (ndarray):
               wavelength array for sensitivity function (nspec,)
            sensfunction (ndarray):
               sensitivity function
            meta_table (table):
               astropy table containing meta data for sensitivity function

        Returns:

        """

        # Flux calibrate the orders that are mutually in the meta_table and in the sobjs. This allows flexibility
        # for applying to data for cases where not all orders are present in the data as in the sensfunc, etc.,
        # i.e. X-shooter with the K-band blocking filter.
        ech_orders = np.array(meta_table['ECH_ORDERS']).flatten()
        #norders = ech_orders.size
        for sci_obj in sobjs:
            # JFH Is there a more elegant pythonic way to do this without looping over both orders and sci_obj?
            indx = np.where(ech_orders == sci_obj.ECH_ORDER)[0]
            if indx.size==1:
                sci_obj.apply_flux_calib(wave[:, indx[0]],sensfunction[:,indx[0]],
                                         sobjs.header['EXPTIME'],
                                         extinct_correct=self.par['extinct_correct'],
                                         longitude=self.spectrograph.telescope['longitude'],
                                         latitude=self.spectrograph.telescope['latitude'],
                                         airmass=float(sobjs.header['AIRMASS']))
            elif indx.size == 0:
                msgs.info('Unable to flux calibrate order = {:} as it is not in your sensitivity function. '
                          'Something is probably wrong with your sensitivity function.'.format(sci_obj.ECH_ORDER))
            else:
                msgs.error('This should not happen')
