"""
Module to run tests on FlatField class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit.tests.tstutils import dev_suite_required, load_kast_blue_masters, cooked_required
from pypeit import flatfield
from pypeit import slittrace
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import pypeitimage
from pypeit import bspline

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


# TODO: Bring this test back in some way?
#def test_step_by_step():
#    if skip_test:
#        assert True
#        return
#    # Masters
#    spectrograph, TSlits, tilts, datasec_img \
#                = load_kast_blue_masters(get_spectrograph=True, tslits=True, tilts=True,
#                                         datasec=True)
#    # Instantiate
#    flatField = flatfield.FlatField(spectrograph, det=1, tilts=tilts,
#                                    tslits_dict=TSlits.tslits_dict.copy())
#    # Use mstrace
#    flatField.mspixelflat = TSlits.mstrace.copy()
#    # Normalize a slit
#    slit=0
#    flatField._prep_tck()
#    modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = flatField.slit_profile(slit)
#    assert np.isclose(iextrap_slit, 0.)
#    # Apply
#    word = np.where(flatField.tslits_dict['slitpix'] == slit + 1)
#    flatField.mspixelflatnrm = flatField.mspixelflat.copy()
#    flatField.mspixelflatnrm[word] /= nrmvals
#    assert np.isclose(np.median(flatField.mspixelflatnrm), 1.0267346)

def test_flatimages():
    tmp = np.ones((1000, 100)) * 10.
    x = np.random.rand(500)
    # Create bspline
    spat_bspline1 = bspline.bspline(x, bkspace=0.01*(np.max(x)-np.min(x)))
    spat_bspline2 = bspline.bspline(x, bkspace=0.01*(np.max(x)-np.min(x)))
    instant_dict = dict(pixelflat_raw=tmp,
                        pixelflat_norm=np.ones_like(tmp),
                        pixelflat_model=None,
                        pixelflat_spat_bsplines=np.asarray([spat_bspline1, spat_bspline2]),
                        pixelflat_spec_illum=None,
                        illumflat_raw=tmp,
                        illumflat_spat_bsplines=np.asarray([spat_bspline1, spat_bspline2]),
                        spat_id=np.asarray([100, 200]),
                        PYP_SPEC="specname")

    flatImages = flatfield.FlatImages(**instant_dict)
    assert flatImages.pixelflat_model is None
    assert flatImages.pixelflat_spec_illum is None
    assert flatImages.pixelflat_spat_bsplines is not None

    # I/O
    outfile = data_path('tst_flatimages.fits')
    flatImages.to_master_file(outfile)
    _flatImages = flatfield.FlatImages.from_file(outfile)

    # Test
    for key in instant_dict.keys():
        if key == 'pixelflat_spat_bsplines':
            np.array_equal(flatImages[key][0].breakpoints,
                           _flatImages[key][0].breakpoints)
            continue
        if key == 'illumflat_spat_bsplines':
            np.array_equal(flatImages[key][0].breakpoints,
                           _flatImages[key][0].breakpoints)
            continue
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(flatImages[key],_flatImages[key])
        else:
            assert flatImages[key] == _flatImages[key]

    os.remove(outfile)

    # Illumflat
#    left = np.full((1000,2), 90, dtype=float)
#    left[:,1] = 190.
#    right = np.full((1000,2), 110, dtype=float)
#    right[:,1] = 210
#    slits = slittrace.SlitTraceSet(left_init=left, right_init=right,
#                                   nspat=1000, PYP_SPEC='dummy')
#    illumflat = flatImages.generate_illumflat(slits)
#    pytest.set_trace()


#@cooked_required
#def test_run():
#    # Masters
#    spectrograph = load_spectrograph('shane_kast_blue')
#    edges, waveTilts = load_kast_blue_masters(edges=True, tilts=True)
#    # Instantiate
#    par = spectrograph.default_pypeit_par()
#    rawflatimg = pypeitimage.PypeItImage(edges.img.copy())
#    # TODO -- We would want to save the detector if we ever planned to re-run from EdgeTrace
#    hdul = fits.HDUList([])
#    rawflatimg.detector = spectrograph.get_detector_par(hdul, 1)
#    flatField = flatfield.FlatField(rawflatimg, spectrograph, par['calibrations']['flatfield'],
#                                    wavetilts=waveTilts, slits=edges.get_slits())
#
#    # Use the trace image
#    flatImages = flatField.run()
#    assert np.isclose(np.median(flatImages.pixelflat), 1.0)
