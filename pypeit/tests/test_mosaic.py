from pathlib import Path
from IPython import embed

import pytest

from astropy.io import fits

from pypeit.pypmsgs import PypeItDataModelError
from pypeit.tests.tstutils import data_output_path
from pypeit.images.mosaic import Mosaic
from pypeit.spectrographs.util import load_spectrograph


def test_io():
    # Create the mosaic
    spec = load_spectrograph('keck_deimos')
    mpar = spec.get_mosaic_par((1,5))

    # Write it
    ofile = data_output_path('tmp_mosaic.fits')
    mpar.to_file(ofile, overwrite=True)

    # Try to read it
    _mpar = Mosaic.from_file(ofile)

    # Change the version
    _ofile = data_output_path('tmp_mosaic_wrongver.fits')
    with fits.open(ofile) as hdu:
        hdu['MOSAIC'].header['DMODVER'] = '1.0.0'
        hdu.writeto(_ofile, overwrite=True)

    # Reading should fail because version is checked by default
    with pytest.raises(PypeItDataModelError):
        _mpar = Mosaic.from_file(_ofile)

    # Should not fail because skipping the version check    
    _mpar = Mosaic.from_file(_ofile, chk_version=False)

    # Remove files
    Path(ofile).unlink()
    Path(_ofile).unlink()
