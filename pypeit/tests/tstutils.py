"""
Odds and ends in support of tests
"""
import os
import pytest

from IPython import embed

import numpy as np
from astropy import time
from astropy.table import Table
import astropy.io.fits as fits

from pypeit import dataPaths
from pypeit.spectrographs.util import load_spectrograph
from pypeit.metadata import PypeItMetaData
from pypeit.inputfiles import PypeItFile 


# Tests require the bspline c extension
try:
    from pypeit.bspline import utilc
except:
    bspline_ext = False
else:
    bspline_ext = True
bspline_ext_required = pytest.mark.skipif(not bspline_ext, reason='Could not import C extension')
# ----------------------------------------------------------------------


# NOTE: Now that the test data files are kept remotely, we have to distinguish
# between paths to files that are *written* (data_output_path) as part of the
# tests and those that are *read* (data_input_path) as part of the tests.  I.e.,
# files to be written should not exist on GitHub, so we should not attempt to
# download them.
# def data_input_path(filename, to_pkg=None):
#    return str(dataPaths.tests.get_file_path(filename, to_pkg=to_pkg))


def data_output_path(filename):
    if len(filename) == 0:
        return str(dataPaths.tests.path)
    return str(dataPaths.tests.path / filename)


def get_kastb_detector():
    """
    Pass back a shane_kast_blue detector when any old one will do

    Returns:
        :class:`pypeit.images.detector_container.DetectorContainer`:

    """
    return load_spectrograph('shane_kast_blue').get_detector_par(1)


def install_shane_kast_blue_raw_data():
    # Download and move all the b*fits.gz files into the local package
    # installation
    files = [dataPaths.tests.get_file_path(f'b{i}.fits.gz', to_pkg='symlink') 
                for i in [1, 11, 12, 13, 21, 22, 23, 24, 27]]


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

    raw_files = [
        'b21.fits.gz',
        'b22.fits.gz',
        'b23.fits.gz',
        'b24.fits.gz',
        'b27.fits.gz'
    ]

    data = Table()
    data['filename'] = [os.path.basename(dataPaths.tests.get_file_path(f, to_pkg='symlink'))
                            for f in raw_files]
    data['frametype'] = ['science']*len(data)
    file_paths = [data_output_path('')]
    setup_dict = {'Setup A': ' '}

    # Return
    return PypeItFile(confdict, file_paths, data, setup_dict)


def make_fake_fits_files():
    """ Generate some raw files covering multiple setups
    """
    spectrograph = load_spectrograph("shane_kast_blue")
    filelist = []
    setups = ['grismA', 'grismB']  # GRISM_N
    nframes = dict({'bias':3, 'flat':2, 'arc':2, 'sci':1})
    # Make some bias frames (one setup only)
    for frmtyp in nframes.keys():
        for ss, setup in enumerate(setups):
            # Only have one set of bias frames, independent of setup
            if frmtyp == 'bias' and ss != 0:
                continue
            # Loop over frame types
            for ff in range(nframes[frmtyp]):
                frname = f"{frmtyp}_{ff+1}_{setup}.fits"
                filelist.append(frname)
                hdu = fits.PrimaryHDU(np.zeros((2,2)))  # Small fake image
                for key in spectrograph.meta.keys():
                    card = spectrograph.meta[key]['card']
                    if key == 'exptime':
                        if frmtyp == 'bias':
                            value = 0.0
                        elif frmtyp == 'sci':
                            value = 1800.0
                        else:
                            value = 60.0
                    elif key == 'mjd':
                        card, value = 'DATE', '2023-03-27T01:27:44.03'
                    elif key == 'binning':
                        continue
                    elif key == 'dispname':
                        value = setup
                    else:
                        value = 'None'
                    hdu.header[card] = value
                if frmtyp == 'flat':
                    hdu.header['LAMPSTA1'] = 'on'
                elif frmtyp == 'arc':
                    hdu.header['LAMPSTAC'] = 'on'
                # Save the fake fits file to disk
                hdu.writeto(frname, overwrite=True)

    # Return the filelist so that it can be later deleted
    return filelist
