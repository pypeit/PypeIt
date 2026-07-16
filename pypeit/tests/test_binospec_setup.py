"""
Tests for Binospec setup grouping (configuration_keys) and frame typing
(check_frame_type), addressing GitHub issue #2137.
"""
import numpy as np
from astropy.table import Table

from pypeit.spectrographs.mmt_binospec import MMTBINOSPECSpectrograph


def test_configuration_keys_includes_decker():
    """Setups must be split by slit mask (decker), not just grating."""
    spec = MMTBINOSPECSpectrograph()
    keys = spec.configuration_keys()
    assert 'dispname' in keys, \
        f"dispname (grating) must be a configuration key; got {keys}"
    assert 'decker' in keys, \
        "decker (MASK) must be a configuration key so different slit masks " \
        "get separate setups"


def _frame_table(**columns):
    """Build a minimal metadata table for frame-typing tests.

    Every column is broadcast to the same length.
    """
    n = max(len(np.atleast_1d(v)) for v in columns.values())
    data = {k: np.broadcast_to(np.atleast_1d(v), n).copy()
            for k, v in columns.items()}
    return Table(data)


def test_incan_on_is_flat():
    """INCAN=on with screen deployed and arc lamp off is a flat/trace."""
    spec = MMTBINOSPECSpectrograph()
    fitstbl = _frame_table(exptime=10.0, lampstat01='off',
                           lampstat02='deployed', lampstat03='on')
    for ftype in ['pixelflat', 'trace', 'illumflat']:
        assert spec.check_frame_type(ftype, fitstbl)[0], \
            f"INCAN=on frame should be typed as {ftype}"


def test_incan_off_is_not_flat():
    """INCAN=off must NOT be typed as a flat/trace even if screen deployed."""
    spec = MMTBINOSPECSpectrograph()
    fitstbl = _frame_table(exptime=10.0, lampstat01='off',
                           lampstat02='deployed', lampstat03='off')
    for ftype in ['pixelflat', 'trace', 'illumflat']:
        assert not spec.check_frame_type(ftype, fitstbl)[0], \
            f"INCAN=off frame must not be typed as {ftype}"
