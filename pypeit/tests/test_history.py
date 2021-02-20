"""
Module to run tests on History class.
"""

import pytest
import os.path

from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
import numpy as np

from pypeit.io import initialize_header
from pypeit.tests.tstutils import cooked_required
from pypeit.history import History
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.framematch import FrameTypeBitMask

@cooked_required
def test_history_header_access():
    """Test building a History object from a fits header"""
    cooked_sci_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science')
    spec1d_file = os.path.join(cooked_sci_dir, 'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_2010Sep13T061231.334.fits')

    header = fits.getheader(spec1d_file)

    history = History(header)

    # TODO verify history once Cooked is generated with it

def test_write_and_append_history():
    """Test appending to a History object and writing the history to a fits header"""

    # Get a blank header
    header = initialize_header(primary=True)
    history = History(header)

    # Append a test history with and without a date
    test_history='Test history'
    history.append(test_history)
    history.append(test_history, add_date=False)

    # Convert to a fits header and verify the 'HISTORY' tag
    history.write_to_header(header)

    assert header['HISTORY'][0][17:] == test_history
    assert header['HISTORY'][1] == test_history

    # Convert date into an astropy time verifying it doesn't throw an exception
    t = Time(header['HISTORY'][0][0:16], format='isot')
    
def test_add_reduce():
    """Test adding history entries for reducing an exposure"""

    # Create a PypeItMetadata object to simulate
    # what would be created in PypeIt::reduce_all
    # This will contain two calibration groups, one
    # that combines two frames d has no calibration
    # frames, and another that combines and subtracts
    # frames and has calibration frames
    spectrograph = load_spectrograph('keck_deimos')
    spectrograph_def_par = spectrograph.default_pypeit_par()
    ftype_bitmask = FrameTypeBitMask()

    science_value = ftype_bitmask.turn_on(np.int_(0),'science')
    standard_value = ftype_bitmask.turn_on(np.int_(0),'standard')
    arc_value = ftype_bitmask.turn_on(np.int_(0),'arc')
    tilt_arc_science_value = ftype_bitmask.turn_on(np.int_(0), 'tilt') | arc_value | science_value

    t = Table(names=['filename',  'frametype',       'framebit',             'target', 'calib', 'comb_id', 'bkg_id'],
              rows=[['file1.fits','science',         science_value,          'targ1',  '1',     '-1',       '-1'],
                    ['file2.fits','standard',        standard_value,         'targ2',  '1',     '-1',       '-1'],
                    ['file3.fits','tilt/arc/science',tilt_arc_science_value, 'targ3',  '2',     '2',        '3'],
                    ['file4.fits','tilt/arc/science',tilt_arc_science_value, 'targ4',  '2',     '3',        '2'], 
                    ['file5.fits','tilt/arc/science',tilt_arc_science_value, 'targ3',  '2',     '2',        '3'],
                    ['file6.fits','tilt/arc/science',tilt_arc_science_value, 'targ4',  '2',     '3',        '2'],
                    ['file7.fits','arc',             arc_value,              'targ5',  '2',     '-1',       '-1']])

    metadata = PypeItMetaData(spectrograph, spectrograph_def_par, data=t)
    metadata.set_calibration_groups()

    #Create history entries for both calibration groups
    history = History()
    history.add_reduce(1, metadata, [0], [])
    history.add_reduce(2, metadata, [2, 4], [3, 5])

    # Verify we get the expected entries
    expected_history = ['PypeIt Reducing target targ1',
                        'Combining frames:', 
                        '"file1.fits"', 
                        'PypeIt Reducing target targ3',
                        'Combining frames:', 
                        '"file3.fits"', 
                        '"file5.fits"', 
                        'Subtracted background from frames:', 
                        '"file4.fits"', 
                        '"file6.fits"',
                        'Callibration frames:', 
                        'tilt/arc/science "file3.fits"', 
                        'tilt/arc/science "file4.fits"', 
                        'tilt/arc/science "file5.fits"', 
                        'tilt/arc/science "file6.fits"', 
                        'arc "file7.fits"']

    
    assert len(expected_history) == len(history.history)
    for i in range(len(expected_history)):
        if expected_history[i].startswith('PypeIt'):
            # The start of a reduce entry should have a 
            # date, construct a Time with it to
            # verify it doesn't raise an exception
            t = Time(history.history[i][0:16], format='isot')
            # Comprare the rest of the string
            assert history.history[i][16:] == ' ' + expected_history[i]
        else:
            assert history.history[i] == expected_history[i]

