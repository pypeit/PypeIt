# Module for flux calibrating spectra
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

from pypeit import msgs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import specobjs
from pypeit import sensfunc
from pypeit.history import History
from astropy import table
from IPython import embed


def flux_calibrate(spec1dfiles, sensfiles, par=None, outfiles=None, chk_version=True):
    """
    Function for flux calibrating spectra.

    Args:
        spec1dfiles (list):
            List of PypeIt spec1d files that you want to flux calibrate
        sensfiles (list):
            List of sensitivity function files to use to flux calibrate the spec1d files. This list and the sensfiles
            list need to have the same length and be aligned
        par (:class:`~pypeit.par.pypeitpar.FluxCalibratePar`, optional):
            Parset object containing parameters governing the flux calibration.
        outfiles (list, optional):
            Names of the output files.  If None, this is set to spec1dfiles and those are overwritten
        chk_version (bool, optional):
            Whether to check of the data model versions of spec1d and sens files. Defaults to True.
    """


    # Output file names
    outfiles = spec1dfiles if outfiles is None else outfiles

    # Load the spectrograph
    header = fits.getheader(spec1dfiles[0])
    spectrograph = load_spectrograph(header['PYP_SPEC'])
    par = spectrograph.default_pypeit_par()['fluxcalib'] if par is None else par

    sensf_last = None
    for spec1, sensf, outfile in zip(spec1dfiles, sensfiles, outfiles):
        # Read in the data
        sobjs = specobjs.SpecObjs.from_fitsfile(spec1, chk_version=chk_version)
        history = History(sobjs.header)
        if sensf != sensf_last:
            sens = sensfunc.SensFunc.from_file(sensf, chk_version=chk_version)
            sensf_last = sensf
            history.append(f'PypeIt Flux calibration "{sensf}"')
        sobjs.apply_flux_calib(par, spectrograph, sens)
        sobjs.write_to_fits(sobjs.header, outfile, history=history, overwrite=True)
