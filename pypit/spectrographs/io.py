""" Spectrograph specific methods for I/O (mainly I)
"""
import numpy as np

from astropy.io import fits

from pypit import msgs
from pypit.spectrographs import deimos
from pypit.spectrographs import lris
from pypit import arparse

from pypit import ardebug as debugger


def load_raw_frame(spectrograph, raw_file, det, dataext=None, disp_dir=0):
    """
    Load data frames, usually raw.

    Parameters
    ----------
    raw_file : str
       Full path to raw_file
    det : int
      Detector number requested, starts at 1
    dataext : int, optional
      Data extension for this detector in the HDU list
    disp_dir : int, optional
      if 1, Transpose the image to align spectral dimension with columns

    Returns
    -------
    frame : ndarray
      the raw_frame
    head : FITS header of the 0th HDU
    """
    msgs.info("Loading raw_file: {:s}".format(raw_file))
    if spectrograph in ['keck_lris_blue', 'keck_lris_red']:
        temp, head0, _ = lris.read_lris(raw_file, det=det)
    elif spectrograph in ['keck_deimos']:
        temp, head0, _ = deimos.read_deimos(raw_file, det=det)
    else:
        hdulist = fits.open(raw_file)
        temp = hdulist[dataext].data
        head0 = hdulist[0].header
    # Turn to float
    temp = temp.astype(np.float)
    if disp_dir == 1:
        temp = temp.T
    return temp, head0


def get_datasec(spectrograph, filename, det, settings_det):
    """  Determine the data and overscan sections of an image

    Currently only used for LRIS and DEIMOS (with their multiple detectors
    packed in funny ways).  Should consider another approach.

    Parameters
    ----------
    spectrograph : str
    filename : str, optional
    det : int
      Detector number, starts at 1
    settings_det : dict
      Used to grab the sections

    Returns
    -------
    datasec : list
    oscansec : list
    naxis0 : int
    naxis1 : int
    """
    # Get naxis0, naxis1, datasec, oscansec, ampsec for specific instruments
    datasec, oscansec, naxis0, naxis1 = [], [], 0, 0

    if spectrograph in ['keck_lris_blue', 'keck_lris_red']:
        msgs.info("Parsing datasec and oscansec from headers")
        temp, head0, secs = lris.read_lris(filename, det)
        for kk in range(settings_det['numamplifiers']):
            datasec.append(arparse.load_sections(secs[0][kk], fmt_iraf=False))
            oscansec.append(arparse.load_sections(secs[1][kk], fmt_iraf=False))
    elif spectrograph in ['keck_deimos']:
        msgs.info("Parsing datasec and oscansec from headers")
        temp, head0, secs = deimos.read_deimos(filename, det=det)
        datasec.append(arparse.load_sections(secs[0][0], fmt_iraf=False))
        oscansec.append(arparse.load_sections(secs[1][0], fmt_iraf=False))
    else:  # Other instruments are set in their settings file
        for i in range(settings_det['numamplifiers']):
            sdatasec = "datasec{0:02d}".format(i+1)
            datasec.append(settings_det[sdatasec])
            soscansec = "oscansec{0:02d}".format(i+1)
            oscansec.append(settings_det[soscansec])
        # Read the image for the shape (just in case)
        temp, _ = load_raw_frame(spectrograph, filename, det, dataext=settings_det['dataext01'],
                              disp_dir=settings_det['dispaxis'])

    # Need naxis0, naxis1 too
    naxis0 = temp.shape[0]
    naxis1 = temp.shape[1]

    # Return
    return datasec, oscansec, naxis0, naxis1


