import os

import pytest
import numpy as np

from astropy.io import fits
from astropy.table import Table

from pypeit.images import buildimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs.util import load_spectrograph
from pypeit.edgetrace import EdgeTraceSet
from pypeit.slittrace import SlitTraceSet
from pypeit import masterframe


# This test check that `maskdesign_matching` method is properly assigning DEIMOS slit-mask design IDs
# to traced slits

# Load flats files
@dev_suite_required
def deimos_flat_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '1200G_Cooper', ifile)
                for ifile in ['d0115_0023.fits.gz', 'd0115_0024.fits.gz', 'd0115_0025.fits.gz']]

def test_maskdef_id():
    # Load instrument
    keck_deimos = load_spectrograph('keck_deimos')
    par = keck_deimos.default_pypeit_par()
    par['calibrations']['traceframe']['process']['use_biasimage'] = False
    trace_par = par['calibrations']['slitedges']

    # working only on detector 1
    det=1
    master_key = 'A_1_{}'.format(str(det).zfill(2))
    redux_path = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'keck_deimos', '1200G_Cooper')
    master_dir = os.path.join(redux_path, 'Masters')

    # Built trace image
    traceImage = buildimage.buildimage_fromlist(keck_deimos, det, par['calibrations']['traceframe'],
                                                deimos_flat_files())

    msbpm = keck_deimos.bpm(traceImage.files[0], det)

    # Run edge trace
    edges = EdgeTraceSet(traceImage, keck_deimos, trace_par, bpm=msbpm, auto=True, maskdesign=True,
                                           debug=False, show_stages=False,qa_path=None)

    # Save masterEdge and masterSlit to file
    edge_masterframe_name = masterframe.construct_file_name(EdgeTraceSet, master_key, master_dir=master_dir)
    edges.to_master_file(edge_masterframe_name)

    slit_masterframe_name = masterframe.construct_file_name(SlitTraceSet,master_key, master_dir=master_dir)
    edges.get_slits().to_master_file(slit_masterframe_name)

    # Check that the `maskdef_id` assigned to the first and last lists is correct
    tab=Table(fits.getdata(slit_masterframe_name))
    assert tab['maskdef_id'][0] == 1098247, 'maskdef_id not found or wrong'
    assert tab['maskdef_id'][-1] == 1098226, 'maskdef_id not found or wrong'