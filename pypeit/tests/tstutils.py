"""
Odds and ends in support of tests
"""
import glob
import os
import pathlib
import pytest

from IPython import embed

import numpy as np
from astropy import time
from astropy.table import Table 
from pypeit import data

from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit.spectrographs.util import load_spectrograph
from pypeit.metadata import PypeItMetaData
from pypeit.inputfiles import PypeItFile 

# ----------------------------------------------------------------------
# pytest @decorators setting the tests to perform

# Tests require the PypeIt dev-suite
#dev_suite_required = pytest.mark.skipif(os.getenv('PYPEIT_DEV') is None
#                                        or not os.path.isdir(os.getenv('PYPEIT_DEV')),
#                                        reason='test requires dev suite')

# Tests require the Cooked data
cooked_required = pytest.mark.skipif(
    os.getenv('PYPEIT_DEV') is None
    or not (pathlib.Path(os.getenv('PYPEIT_DEV')) / 'Cooked').is_dir(),
    reason='no dev-suite cooked directory'
)

# Tests require the Telluric file (Mauna Kea)
par = Spectrograph.default_pypeit_par()
tell_test_grid = data.Paths.telgrid / 'TelFit_MaunaKea_3100_26100_R20000.fits'
telluric_required = pytest.mark.skipif(not tell_test_grid.is_file(),
                                       reason='no Mauna Kea telluric file')

# Tests require the bspline c extension
try:
    from pypeit.bspline import utilc
except:
    bspline_ext = False
else:
    bspline_ext = True
bspline_ext_required = pytest.mark.skipif(not bspline_ext, reason='Could not import C extension')
# ----------------------------------------------------------------------


def data_path(filename):
    data_dir = pathlib.Path(__file__).parent.absolute().joinpath('files')
    # TODO: This really should have the `.resolve()`, but it crashes the
    #       Windows/python3.9 CI test (only that one).  When PypeIt advances
    #       to python>=3.10, reinstate the last part of the following line:
    return str(data_dir.joinpath(filename))#.resolve())


def get_kastb_detector():
    """
    Pass back a shane_kast_blue detector when any old one will do

    Returns:
        :class:`pypeit.images.detector_container.DetectorContainer`:

    """
    return load_spectrograph('shane_kast_blue').get_detector_par(1)


def dummy_fitstbl(nfile=10, spectro_name='shane_kast_blue', directory='', notype=False):
    """
    Generate a dummy fitstbl for testing

    Parameters
    ----------
    nfile : int, optional
      Number of files to mimic
    spectro_name : str, optional
      Name of spectrograph to mimic
    notype : bool (optional)
      If True, do not add image type info to the fitstbl

    Returns
    -------
    fitstbl : PypeItMetaData

    """
    fitsdict = {}
    fitsdict['index'] = np.arange(nfile)
    fitsdict['directory'] = [directory]*nfile
    fitsdict['filename'] = ['b{:03d}.fits.gz'.format(i) for i in range(nfile)]
    # TODO: The below will fail at 60
    dates = ['2015-01-23T00:{:02d}:11.04'.format(i) for i in range(nfile)]
    ttime = time.Time(dates, format='isot')
    fitsdict['mjd'] = ttime.mjd
    fitsdict['target'] = ['Dummy']*nfile
    fitsdict['ra'] = ['00:00:00']*nfile
    fitsdict['dec'] = ['+00:00:00']*nfile
    fitsdict['exptime'] = [300.] * nfile
    fitsdict['dispname'] = ['600/4310'] * nfile
    fitsdict['dichroic'] = ['560'] * nfile
    fitsdict["binning"] = ['1,1']*nfile
    fitsdict["airmass"] = [1.0]*nfile

    if spectro_name == 'shane_kast_blue':
        fitsdict['numamplifiers'] = [1] * nfile
        # Lamps
        for i in range(1,17):
            fitsdict['lampstat{:02d}'.format(i)] = ['off'] * nfile
        fitsdict['exptime'][0] = 0        # Bias
        fitsdict['lampstat06'][1] = 'on'  # Arc
        fitsdict['exptime'][1] = 30       # Arc
        fitsdict['lampstat01'][2] = 'on'  # Trace, pixel, slit flat
        fitsdict['lampstat01'][3] = 'on'  # Trace, pixel, slit flat
        fitsdict['exptime'][2] = 30     # flat
        fitsdict['exptime'][3] = 30     # flat
        fitsdict['ra'][4] = '05:06:36.6'  # Standard
        fitsdict['dec'][4] = '52:52:01.0'
        fitsdict['airmass'][4] = 1.2
        fitsdict['ra'][5] = '07:06:23.45' # Random object
        fitsdict['dec'][5] = '+30:20:50.5'
        fitsdict['decker'] = ['0.5 arcsec'] * nfile

    # arrays
    for k in fitsdict.keys():
        fitsdict[k] = np.array(fitsdict[k])

    spectrograph = load_spectrograph(spectro_name)
    fitstbl = PypeItMetaData(spectrograph, spectrograph.default_pypeit_par(), data=fitsdict)
    fitstbl['instrume'] = spectro_name
    type_bits = np.zeros(len(fitstbl), dtype=fitstbl.type_bitmask.minimum_dtype())

    # Image typing
    if not notype:
        if spectro_name == 'shane_kast_blue':
            #fitstbl['sci_ID'] = 1  # This links all the files to the science object
            type_bits[0] = fitstbl.type_bitmask.turn_on(type_bits[0], flag='bias')
            type_bits[1] = fitstbl.type_bitmask.turn_on(type_bits[1], flag='arc')
            type_bits[1] = fitstbl.type_bitmask.turn_on(type_bits[1], flag='tilt')
            type_bits[2:4] = fitstbl.type_bitmask.turn_on(type_bits[2:4],
                                                          flag=['pixelflat', 'trace', 'illumflat'])
            type_bits[4] = fitstbl.type_bitmask.turn_on(type_bits[4], flag='standard')
            type_bits[5:] = fitstbl.type_bitmask.turn_on(type_bits[5:], flag='science')
            fitstbl.set_frame_types(type_bits)
            # Calibration groups
            cfgs = fitstbl.unique_configurations() #ignore_frames=['bias', 'dark'])
            fitstbl.set_configurations(configs=cfgs)
            fitstbl.set_calibration_groups() #global_frames=['bias', 'dark'])

    return fitstbl


def make_shane_kast_blue_pypeitfile():
    """ Generate a PypeItFile class """
    # Bits needed to generate a PypeIt file
    confdict = {'rdx': {'spectrograph': 'shane_kast_blue'}}

    data = Table()
    data['filename'] = [os.path.basename(item) for item in glob.glob(data_path('b2*fits.gz'))]
    data['frametype'] = ['science']*len(data)
    file_paths = [data_path('')]
    setup_dict = {'Setup A': ' '}

    # Return
    return PypeItFile(confdict, file_paths, data, setup_dict)
