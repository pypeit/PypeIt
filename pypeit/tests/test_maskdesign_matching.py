import os

from pypeit.images import buildimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet


# This test check that `maskdesign_matching` method is properly assigning DEIMOS slit-mask design IDs
# to traced slits

# Load flats files
def deimos_flat_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '1200G_Cooper', ifile)
                for ifile in ['d0115_0023.fits.gz', 'd0115_0024.fits.gz', 'd0115_0025.fits.gz']]

@dev_suite_required
def test_maskdef_id():
    # Load instrument
    keck_deimos = load_spectrograph('keck_deimos')
    par = keck_deimos.default_pypeit_par()
    par['calibrations']['traceframe']['process']['use_biasimage'] = False

    # working only on detector 1
    det=1

    # Built trace image
    traceImage = buildimage.buildimage_fromlist(keck_deimos, det, par['calibrations']['traceframe'],
                                                deimos_flat_files())
    msbpm = keck_deimos.bpm(traceImage.files[0], det)

    # load specific config parameters
    par = keck_deimos.config_specific_par(traceImage.files[0])
    trace_par = par['calibrations']['slitedges']


    # Run edge trace
    edges = EdgeTraceSet(traceImage, keck_deimos, trace_par, bpm=msbpm, auto=True,
                                           debug=False, show_stages=False,qa_path=None)

    slits = edges.get_slits()
    # Check that the `maskdef_id` assigned to the first and last slits is correct
    assert slits.maskdef_id[0] == 1098247, 'maskdef_id not found or wrong'
    assert slits.maskdef_id[-1] == 1098226, 'maskdef_id not found or wrong'
    # `maskdef_id` is the slit_id ("dSlitId") from the DEIMOS slit-mask design. If the matching is done
    # correctly these are the slit_id values that the first and last slits should have. These values are
    # coming from a reduction run with DEEP2 IDL-based pipeline.

