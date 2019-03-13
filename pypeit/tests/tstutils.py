"""
Odds and ends in support of tests
"""
import os
import pytest
import numpy as np
import copy

from astropy import time

from pypeit import arcimage
from pypeit import traceslits
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit.masterframe import MasterFrame
from pypeit.core.wavecal import waveio
from pypeit.spectrographs.util import load_spectrograph
from pypeit.metadata import PypeItMetaData

# Create a decorator for tests that require the PypeIt dev suite
dev_suite_required = pytest.mark.skipif(os.getenv('PYPEIT_DEV') is None,
                                        reason='test requires dev suite')

cooked_required = pytest.mark.skipif(os.getenv('PYPEIT_DEV') is None or
                            not os.path.isdir(os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked')),
                            reason='no dev-suite cooked directory')

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

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
            type_bits[2:4] = fitstbl.type_bitmask.turn_on(type_bits[2:4], flag=['pixelflat', 'trace'])
            type_bits[4] = fitstbl.type_bitmask.turn_on(type_bits[4], flag='standard')
            type_bits[5:] = fitstbl.type_bitmask.turn_on(type_bits[5:], flag='science')
            fitstbl.set_frame_types(type_bits)
            # Calibration groups
            cfgs = fitstbl.unique_configurations(ignore_frames=['bias', 'dark'])
            fitstbl.set_configurations(cfgs)
            fitstbl.set_calibration_groups(global_frames=['bias', 'dark'])

    return fitstbl

# TODO: Need to split this into functions that do and do not require
# cooked.  We should remove the get_spectrograph option.
def load_kast_blue_masters(aimg=False, tslits=False, tilts=False, datasec=False, wvcalib=False):
    """
    Load up the set of shane_kast_blue master frames

    Args:
        get_spectrograph:
        aimg:
        tslits:
        tilts:
        datasec:
        wvcalib:

    Returns:

    """

    spectrograph = load_spectrograph('shane_kast_blue')
    spectrograph.naxis = (2112,350)     # Image shape with overscan

    master_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Shane_Kast_blue')
#    master_dir = root_path+'_'+spectrograph.spectrograph

    reuse_masters = True

    # Load up the Masters
    ret = []

#    if get_spectrograph:
#        ret.append(spectrograph)

    master_key = 'A_1_01'
    if aimg:
        AImg = arcimage.ArcImage(spectrograph, master_key=master_key, master_dir=master_dir,
                                 reuse_masters=reuse_masters)
        msarc = AImg.load()
        ret.append(msarc)

    if tslits:
        trace_file = os.path.join(master_dir, MasterFrame.construct_file_name('Trace', master_key))
        tslits_dict, mstrace = traceslits.TraceSlits.load_from_file(trace_file)
        ret.append(tslits_dict)
        ret.append(mstrace)

    if tilts:
        tilts_file = os.path.join(master_dir, MasterFrame.construct_file_name('Tilts', master_key))
        tilts_dict = wavetilts.WaveTilts.load_from_file(tilts_file)
        ret.append(tilts_dict)

    if datasec:
        datasec_img = spectrograph.get_datasec_img(data_path('b1.fits.gz'), 1)
        ret.append(datasec_img)

    if wvcalib:
        calib_file = os.path.join(master_dir,
                                  MasterFrame.construct_file_name('WaveCalib', master_key))
        wv_calib = waveio.load_wavelength_calibration(calib_file) 
        ret.append(wv_calib)

    # Return
    return ret

def instant_traceslits(mstrace_file, det=None):
    """
    Instantiate a TraceSlits object from the master file

    The loaded tslits_dict is set as the atribute

    Args:
        mstrace_file (str):
        det (int, optional):

    Returns:
        Spectrograph, TraceSlits:

    """
    # Load
    tslits_dict, mstrace = traceslits.TraceSlits.load_from_file(mstrace_file)
    # Instantiate
    spectrograph = load_spectrograph(tslits_dict['spectrograph'])
    par = spectrograph.default_pypeit_par()
    msbpm = spectrograph.bpm(shape=mstrace.shape, det=det)
    binning = tslits_dict['binspectral'], tslits_dict['binspatial']
    traceSlits = traceslits.TraceSlits(spectrograph, par['calibrations']['slits'],
                                       msbpm=msbpm, binning=binning)
    traceSlits.mstrace = copy.deepcopy(mstrace)
    traceSlits.tslits_dict = copy.deepcopy(tslits_dict)
    return spectrograph, traceSlits
