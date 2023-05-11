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
from pypeit.history import History
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.framematch import FrameTypeBitMask

def verify_history(history, expected_history):

    if not len(history.history) == len(expected_history):
        return False

    for i in range(len(expected_history)):
        if expected_history[i].startswith(' '):
            # Verify first line has the date by building a Time
            t = Time(history.history[i][0:16], format='isot')
            if not history.history[i][16:] == expected_history[i]:
                return False
        else:
            if not history.history[i] == expected_history[i]:
                return False

        return True


def test_history_header_access():
    """Test building a History object from a fits header"""

    # Get a blank header
    header = initialize_header()

    history1 = History()

    expected_history = ['2021-02-23T21:17 PypeIt Reducing target CFHQS1',
                        'Combining frames:',
                        '"DE.20100913.22358.fits.gz"',
                        'Callibration frames:',
                        'arc,tilt "DE.20100913.56927.fits.gz"',
                        'pixelflat,illumflat,trace "DE.20100913.57161.fits.gz"',
                        'pixelflat,illumflat,trace "DE.20100913.57006.fits.gz"']

    # Write and then read history read from header
    for h in expected_history:
        history1.append(h, add_date=False)

    history1.write_to_header(header)

    history2 = History(header)


    assert verify_history(history2, expected_history)

def test_write_and_append_history():
    """Test appending to a History object and writing the history to a fits header"""

    # Get a blank header
    header = initialize_header()
    history = History(header)

    # Append a test history with and without a date
    test_history='Test history'
    history.append(test_history)
    history.append(test_history, add_date=False)

    # Convert to a fits header and verify the 'HISTORY' keyword
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
    # that combines two frames and has no calibration
    # frames, and another that combines and subtracts
    # frames and has calibration frames
    spectrograph = load_spectrograph('keck_deimos')
    spectrograph_def_par = spectrograph.default_pypeit_par()
    ftype_bitmask = FrameTypeBitMask()

    # TODO: Does this need to use np.int_?  Can we just use 'int'?
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
    expected_history = [' PypeIt Reducing target targ1',
                        'Combining frames:', 
                        '"file1.fits"', 
                        ' PypeIt Reducing target targ3',
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

    
    assert verify_history(history, expected_history)

def test_add_coadd1d(monkeypatch):
    """Test adding coadd1d entries into history"""

    # Create a mock getheader to return the SEMESTER/PROGID header keywords add_coadd1d will look for
    def mock_getheader(file, **kwargs):
        if file == "spec1d_file1.fits":
            return {"SEMESTER": "2021A",
                    "PROGID":   "test_prog" }
        else:
            return dict()

    with monkeypatch.context() as m:
        monkeypatch.setattr(fits, "getheader", mock_getheader)

        # First test when the spec1d_files and objids 
        # match length wise
        spec1d_files = ['spec1d_file1.fits', 
                        'spec1d_file1.fits', 
                        'spec1d_file2.fits',
                        'spec1d_file3.fits', 
                        'spec1d_file3.fits', 
                        'spec1d_file3.fits']

        objids = ['SPAT001-SLIT001-DET01',
                  'SPAT101-SLIT002-DET02',
                  'SPAT002-SLIT001-DET01',
                  'SPAT003-SLIT001-DET01',
                  'SPAT103-SLIT002-DET02',
                  'SPAT113-SLIT002-DET03']

        history = History()
        history.add_coadd1d(spec1d_files, objids)

        expected_history = [' PypeIt Coadded 6 objects from 3 spec1d files',
                            'From "spec1d_file1.fits"',
                            'Semester: 2021A Program ID: test_prog',
                            'SPAT001-SLIT001-DET01',
                            'SPAT101-SLIT002-DET02',
                            'From "spec1d_file2.fits"',
                            'SPAT002-SLIT001-DET01',                        
                            'From "spec1d_file3.fits"',
                            'SPAT003-SLIT001-DET01',
                            'SPAT103-SLIT002-DET02',
                            'SPAT113-SLIT002-DET03']

        assert verify_history(history, expected_history)

        # Now test when the # of object ids is shorter than the # of files
        # (This shouldn't happen but I want to make sure History behaves in a sane manner)
        short_objids = ['SPAT001-SLIT001-DET01',
                        'SPAT101-SLIT002-DET02',
                        'SPAT002-SLIT001-DET01',
                        'SPAT003-SLIT001-DET01']

        history = History()
        history.add_coadd1d(spec1d_files, short_objids)

        # THe expected result is for the extra files to be ignored
        expected_history = [' PypeIt Coadded 4 objects from 3 spec1d files',
                            'From "spec1d_file1.fits"',
                            'Semester: 2021A Program ID: test_prog',
                            'SPAT001-SLIT001-DET01',
                            'SPAT101-SLIT002-DET02',
                            'From "spec1d_file2.fits"',
                            'SPAT002-SLIT001-DET01',
                            'From "spec1d_file3.fits"',
                            'SPAT003-SLIT001-DET01']

        assert verify_history(history, expected_history)

        # Now test with fewer spec1d files than objects
        short_spec1d_files = ['spec1d_file1.fits', 
                            'spec1d_file1.fits', 
                            'spec1d_file2.fits',
                            'spec1d_file3.fits', 
                            'spec1d_file3.fits' ]

        history = History()
        history.add_coadd1d(short_spec1d_files, objids)

        # The expected result is to ignore the extra object ids
        expected_history = [' PypeIt Coadded 5 objects from 3 spec1d files',
                            'From "spec1d_file1.fits"',
                            'Semester: 2021A Program ID: test_prog',
                            'SPAT001-SLIT001-DET01',
                            'SPAT101-SLIT002-DET02',
                            'From "spec1d_file2.fits"',
                            'SPAT002-SLIT001-DET01',
                            'From "spec1d_file3.fits"',
                            'SPAT003-SLIT001-DET01',
                            'SPAT103-SLIT002-DET02']


        assert verify_history(history, expected_history)