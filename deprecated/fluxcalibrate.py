# Module for flux calibrating spectra
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

from pypeit import msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import specobjs
from pypeit.history import History
from astropy import table
from IPython import embed



class FluxCalibrate:
    """
    Class for flux calibrating spectra.

    Args:
        spec1dfiles (list):
            List of PypeIt spec1d files that you want to flux calibrate
        sensfiles (list):
            List of sensitivity function files to use to flux calibrate the spec1d files. This list and the sensfiles
            list need to have the same length and be aligned
        par (pypeit.par.pypeitpar.FluxCalibrate, optional):
            Parset object containing parameters governing the flux calibration.
        outfiles (list, optional):
            Names of the output files.  If None, this is set to spec1dfiles and those are overwritten
    """
    def __init__(self, spec1dfiles, sensfiles, par=None, outfiles=None, debug=False):

        from pypeit import sensfunc

        self.spec1dfiles = spec1dfiles
        self.sensfiles = sensfiles

        # Output file names
        self.outfiles = spec1dfiles if outfiles is None else outfiles

        # Load the spectrograph
        header = fits.getheader(spec1dfiles[0])
        self.spectrograph = load_spectrograph(header['PYP_SPEC'])
        self.par = self.spectrograph.default_pypeit_par()['fluxcalib'] if par is None else par
        self.debug = debug

        sensf_last = None
        for spec1, sensf, outfile in zip(self.spec1dfiles, self.sensfiles, self.outfiles):
            # Read in the data
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1)
            history = History(sobjs.header)
            if sensf != sensf_last:
                sens = sensfunc.SensFunc.from_file(sensf)
                sensf_last = sensf
                history.append(f'PypeIt Flux calibration "{sensf}"')
            apply_flux_calib(self.par, self.spectrograph, sobjs, sens)
            sobjs.write_to_fits(sobjs.header, outfile, history=history, overwrite=True)


def apply_flux_calib(par, spectrograph, sobjs, sens):
    """
    Flux calibrate the provided object spectra (``sobjs``) using the
    provided sensitivity function (``sens``).


    Args:
        sobjs (:class:`~pypeit.specobjs.SpecObjs`):
            Object spectra
        sens (:class:`~pypeit.sensfunc.SensFunc`):
            Sensitivity function
    """

    _extinct_correct = (True if sens.algorithm == 'UVIS' else False) \
        if par['extinct_correct'] is None else par['extinct_correct']

    if spectrograph.pypeline == 'MultiSlit':
        for ii, sci_obj in enumerate(sobjs):
            if sens.wave.shape[1] == 1:
                sci_obj.apply_flux_calib(sens.wave[:, 0], sens.zeropoint[:, 0],
                                         sobjs.header['EXPTIME'],
                                         extinct_correct=_extinct_correct,
                                         longitude=spectrograph.telescope['longitude'],
                                         latitude=spectrograph.telescope['latitude'],
                                         extinctfilepar=par['extinct_file'],
                                         extrap_sens=par['extrap_sens'],
                                         airmass=float(sobjs.header['AIRMASS']))
            elif sens.wave.shape[1] > 1 and sens.splice_multi_det:
                # This deals with the multi detector case where the sensitivity function is spliced. Note that
                # the final sensitivity function written to disk is  the spliced one. This functionality is only
                # used internal to sensfunc.py for fluxing the standard for the QA plot.
                sci_obj.apply_flux_calib(sens.wave[:, ii], sens.zeropoint[:, ii],
                                         sobjs.header['EXPTIME'],
                                         extinct_correct=_extinct_correct,
                                         longitude=spectrograph.telescope['longitude'],
                                         latitude=spectrograph.telescope['latitude'],
                                         extinctfilepar=par['extinct_file'],
                                         extrap_sens=par['extrap_sens'],
                                         airmass=float(sobjs.header['AIRMASS']))
            else:
                msgs.error('This should not happen, there is a problem with your sensitivity function.')


    elif spectrograph.pypeline == 'Echelle':
        # Flux calibrate the orders that are mutually in the meta_table and in
        # the sobjs. This allows flexibility for applying to data for cases
        # where not all orders are present in the data as in the sensfunc, etc.,
        # i.e. X-shooter with the K-band blocking filter.
        ech_orders = np.array(sens.sens['ECH_ORDERS']).flatten()
        for sci_obj in sobjs:
            # JFH Is there a more elegant pythonic way to do this without looping over both orders and sci_obj?
            indx = np.where(ech_orders == sci_obj.ECH_ORDER)[0]
            if indx.size==1:
                sci_obj.apply_flux_calib(sens.wave[:, indx[0]], sens.zeropoint[:,indx[0]],
                                         sobjs.header['EXPTIME'],
                                         extinct_correct=_extinct_correct,
                                         extrap_sens = par['extrap_sens'],
                                         longitude=spectrograph.telescope['longitude'],
                                         latitude=spectrograph.telescope['latitude'],
                                         extinctfilepar=par['extinct_file'],
                                         airmass=float(sobjs.header['AIRMASS']))
            elif indx.size == 0:
                msgs.info('Unable to flux calibrate order = {:} as it is not in your sensitivity function. '
                          'Something is probably wrong with your sensitivity function.'.format(sci_obj.ECH_ORDER))
            else:
                msgs.error('This should not happen')

    else:
        msgs.error('Unrecognized pypeline: {0}'.format(spectrograph.pypeline))