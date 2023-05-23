"""
Module to test ArcImage.
"""
from pathlib import Path

from IPython import embed

import numpy as np

from pypeit.images import pypeitimage
from pypeit.images import buildimage
from pypeit.images import detector_container
from pypeit.tests.test_detector import def_det
from pypeit.tests.tstutils import data_path


def test_init():
    # Instantiate a simple pypeitImage
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.reinit_mask()
    pypeitImage.detector = detector_container.DetectorContainer(**def_det)
    # Now the arcimage
    arcImage = buildimage.ArcImage.from_pypeitimage(pypeitImage)


def test_io():
    # Instantiate a simple pypeitImage
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.reinit_mask()
    pypeitImage.detector = detector_container.DetectorContainer(**def_det)
    pypeitImage.PYP_SPEC = 'shane_kast_blue'
    # Now the arcimage
    arcImage = buildimage.ArcImage.from_pypeitimage(pypeitImage, calib_dir=data_path(''),
                                                    setup='A', calib_id=['1'],
                                                    detname='DET01')
    # Set paths and check name
    ofile = Path(arcImage.get_path()).resolve()
    assert str(ofile) == str(Path(data_path('Arc_A_1_DET01.fits')).resolve()), \
            'Calibration file name changed'
    # Write
    arcImage.to_file(overwrite=True)
    # Read
    _arcImage = buildimage.ArcImage.from_file(str(ofile))
    # Random set of checks to make sure the written and read versions of
    # arcImage are identical
    assert isinstance(_arcImage.detector, detector_container.DetectorContainer), \
            'detector has wrong type'
    assert arcImage.detector.version == detector_container.DetectorContainer.version, \
            'detector version changed'
    assert np.array_equal(arcImage.detector.gain, _arcImage.detector.gain), \
            'Detector properties changed'
    # Check that the image data itself did not change
    assert np.array_equal(_arcImage.image, arcImage.image), 'image data changed'
    assert np.array_equal(_arcImage.fullmask.mask, arcImage.fullmask.mask), 'mask changed'

    # Cleanup
    ofile.unlink()



