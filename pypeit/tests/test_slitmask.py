"""
Module to run tests on PypeItPar classes
"""
import os
import numpy

import pytest

from astropy.io import fits

from pypeit.images import buildimage
from pypeit import edgetrace, slittrace, specobjs
from pypeit.spectrographs.keck_deimos import KeckDEIMOSSpectrograph
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests.tstutils import dev_suite_required, cooked_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def deimos_flat_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '830G_M_8500', ifile)
            for ifile in ['DE.20100913.57161.fits.gz', 'DE.20100913.57006.fits.gz']]


@cooked_required
def test_assign_maskinfo():

    # Spectrograph
    keck_deimos = KeckDEIMOSSpectrograph()
    par = keck_deimos.default_pypeit_par()
    # working only on detector 3
    det = 3
    det_buffer = par['calibrations']['slitedges']['det_buffer']

    # Built trace image
    traceImage = buildimage.buildimage_fromlist(keck_deimos, det, par['calibrations']['traceframe'],
                                                deimos_flat_files())
    msbpm = keck_deimos.bpm(traceImage.files[0], det)

    # load specific config parameters
    par = keck_deimos.config_specific_par(traceImage.files[0])
    trace_par = par['calibrations']['slitedges']

    # Run edge trace
    edges = edgetrace.EdgeTraceSet(traceImage, keck_deimos, trace_par, bpm=msbpm, auto=True,
                                           debug=False, show_stages=False,qa_path=None)

    slits = edges.get_slits()

    slits_left, slits_right, _ = slits.select_edges()

    # Test that the maskfile is saved properly
    hdul = fits.open(slits.maskfile)
    det_par = keck_deimos.get_detector_par(hdul, det=det)

    specobjs_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits')
    # specobjs_file = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'keck_deimos',
    #                              '830G_M_8500', 'Science',
    #                              'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits')
    sobjs = specobjs.SpecObjs.from_fitsfile(specobjs_file)
    # Init at null
    for sobj in sobjs:
        sobj.MASKDEF_ID = None
        sobj.MASKDEF_OBJNAME = None
        sobj.RA = None
        sobj.DEC = None
        sobj.MASKDEF_EXTRACT = None

    # Run me
    slits.assign_maskinfo(sobjs, det_par['platescale'], slits_left, slits_right, det_buffer)

    # Test
    assert sobjs[sobjs.SLITID == 496].MASKDEF_OBJNAME == 'ero89', 'Wrong MASKDEF_OBJNAME'
    assert sobjs[sobjs.SLITID == 496].RA == 352.27471667, 'Wrong object RA'
    assert sobjs[sobjs.SLITID == 496].DEC == -3.09223056, 'Wrong object DEC'

    # Write sobjs
    sobjs.write_to_fits({}, data_path('tst_sobjs.fits'))
    os.remove(data_path('tst_sobjs.fits'))


@cooked_required
def test_add_missing_obj():

    # Spectrograph
    keck_deimos = KeckDEIMOSSpectrograph()
    par = keck_deimos.default_pypeit_par()
    # working only on detector 7
    det = 7
    det_buffer = par['calibrations']['slitedges']['det_buffer']

    # Built trace image
    traceImage = buildimage.buildimage_fromlist(keck_deimos, det, par['calibrations']['traceframe'],
                                                deimos_flat_files())
    msbpm = keck_deimos.bpm(traceImage.files[0], det)

    # load specific config parameters
    par = keck_deimos.config_specific_par(traceImage.files[0])
    trace_par = par['calibrations']['slitedges']

    # Run edge trace
    edges = edgetrace.EdgeTraceSet(traceImage, keck_deimos, trace_par, bpm=msbpm, auto=True,
                                           debug=False, show_stages=False,qa_path=None)

    slits = edges.get_slits()

    slits_left, slits_right, _ = slits.select_edges()

    # Test that the maskfile is saved properly
    hdul = fits.open(slits.maskfile)
    det_par = keck_deimos.get_detector_par(hdul, det=det)

    specobjs_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits')
    # specobjs_file = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'keck_deimos',
    #                              '830G_M_8500', 'Science',
    #                              'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_2010Sep13T061231.334.fits')
    sobjs = specobjs.SpecObjs.from_fitsfile(specobjs_file)
    cut_sobjs = sobjs[sobjs.DET == det]
    # Init at null
    idx_remove = []
    for i, sobj in enumerate(cut_sobjs):
        if sobj.MASKDEF_EXTRACT:
            idx_remove.append(i)
        else:
            sobj.MASKDEF_ID = None
            sobj.MASKDEF_OBJNAME = None
            sobj.RA = None
            sobj.DEC = None
            sobj.MASKDEF_EXTRACT = None
    cut_sobjs.remove_sobj(idx_remove)

    # Run assign mask info
    expected_objpos = slits.assign_maskinfo(cut_sobjs, det_par['platescale'], slits_left, slits_right, det_buffer)

    # Run add missing objects
    median_off = 0.
    fwhm = par['reduce']['findobj']['find_fwhm']
    cut_sobj = slits.mask_add_missing_obj(cut_sobjs, expected_objpos, fwhm, median_off, slits_left, slits_right)

    # Test that undetected object are found at the correct location (the correct location is
    # verified by visual inspection)
    assert round(cut_sobj[cut_sobj.MASKDEF_OBJNAME == 'ero884'].SPAT_PIXPOS[0]) == 2022, \
        'Wrong object (ero884) location on the slit'
    assert round(cut_sobj[cut_sobj.MASKDEF_OBJNAME == 'ero191'].SPAT_PIXPOS[0]) == 1125, \
        'Wrong object (ero191) location on the slit'

@dev_suite_required
def test_deimosslitmask():
    f = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_deimos', '830G_M_8500',
                     'DE.20100913.22358.fits.gz')
    spec = KeckDEIMOSSpectrograph()
    spec.get_slitmask(f)
    assert spec.slitmask.nslits == 106, 'Incorrect number of slits read!'

