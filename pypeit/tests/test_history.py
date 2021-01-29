"""
Module to run tests on History class.
"""

import pytest
import os.path

from astropy.io import fits
from astropy.time import Time

from pypeit.io import initialize_header
from pypeit.tests.tstutils import cooked_required
from pypeit.history import History

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
    