# Module to run tests on TraceSlits class
#   Requires files in Development suite and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
import glob
import numpy as np

from pypeit import traceslits
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs import util
from pypeit.core import trace_slits

def chk_for_files(root):
    files = glob.glob(root+'*')
    return len(files) != 0

def instantiate(mstrace, tslits_dict, det=None):
    # Instantiate
    spectrograph = util.load_spectrograph(tslits_dict['spectrograph'])
    par = spectrograph.default_pypeit_par()
    msbpm = spectrograph.bpm(shape=mstrace.shape, det=det)
    binning = tslits_dict['binspectral'], tslits_dict['binspatial']
    traceSlits = traceslits.TraceSlits(mstrace, spectrograph, par['calibrations']['slits'],
                                       msbpm=msbpm, binning=binning)
    return spectrograph, traceSlits

@dev_suite_required
def test_addrm_slit():
    """ This tests the add and remove methods for user-supplied slit fussing. """

    # Check for files
    mstrace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace',
                                'MasterTrace_KeckLRISr_400_8500_det1.fits')
    assert chk_for_files(mstrace_file)
    tslits_dict, mstrace = traceslits.load_tslits(mstrace_file)
    norig = tslits_dict['nslits']
    # Instantiate
    _, traceSlits = instantiate(mstrace, tslits_dict)
    traceSlits.lcen = tslits_dict['slit_left']
    traceSlits.rcen = tslits_dict['slit_righ']

    #  Add dummy slit
    #  y_spec, left edge, right edge on image
    add_user_slits = [[1024, 140, 200]]
    traceSlits.lcen, traceSlits.rcen = trace_slits.add_user_edges(
        traceSlits.lcen, traceSlits.rcen, add_user_slits)
    # Test
    assert traceSlits.nslit == (norig+1)
    xcen0 = np.median((traceSlits.lcen[:,0] + traceSlits.rcen[:,0])/2.)
    assert np.abs(xcen0-170) < 3

    # Remove it
    rm_user_slits = [[1024, 170]]
    traceSlits.lcen, traceSlits.rcen = trace_slits.rm_user_edges(
        traceSlits.lcen, traceSlits.rcen, rm_user_slits)
    assert traceSlits.nslit == norig

@dev_suite_required
def test_chk_kast_slits():
    """ This tests finding the longslit for Kast blue """
    # Red, blue
    for root in ['MasterTrace_ShaneKastred_600_7500_d55.fits', 'MasterTrace_ShaneKastblue_600_4310_d55.fits']:
        # Load
        mstrace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
        assert chk_for_files(mstrace_file)
        tslits_dict, mstrace = traceslits.load_tslits(mstrace_file)
        norig = tslits_dict['nslits']
        # Instantiate
        _, traceSlits = instantiate(mstrace, tslits_dict)
        # Run me
        traceSlits.run(trim_slits=False, write_qa=False)  # Don't need plate_scale for longslit
        # Test
        assert traceSlits.nslit == norig

@dev_suite_required
def test_chk_lris_blue_slits():
    """ This tests slit finding for LRISb """
    for nslit, binning, det, root in zip([1, 14, 15],
                                         [(2,2), (2,2), (2,2)],
                                         [2,1,2],
                                         ['MasterTrace_KeckLRISb_long_600_4000_det2.fits',
                                          'MasterTrace_KeckLRISb_multi_600_4000_det1.fits',
                                          'MasterTrace_KeckLRISb_multi_600_4000_det2.fits',
                                          ]):
        mstrace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
        assert chk_for_files(mstrace_file)
        # Instantiate
        tslits_dict, mstrace = traceslits.load_tslits(mstrace_file)
        ipos = mstrace_file.find('det')  # Update to using the trace dict eventually
        spectrograph, traceSlits = instantiate(mstrace, tslits_dict, det=int(mstrace_file[ipos+3]))
        norig = tslits_dict['nslits']
        #assert norig == nslit
        # Run me
        plate_scale = binning[0]*spectrograph.detector[det-1]['platescale']
        traceSlits.run(show=False, plate_scale=plate_scale, write_qa=False)
        # Test
        assert traceSlits.nslit == norig

@dev_suite_required
def test_chk_lris_red_slits():
    """ This tests slit finding for LRISr """
    for nslit, binning, det, root in zip([1, 13, 14],
                                         [(2,2), (2,2), (2,2)],
                                         [2,1,2],
                           ['MasterTrace_KeckLRISr_long_600_7500_d560.fits',
                            'MasterTrace_KeckLRISr_400_8500_det1.fits',
                            'MasterTrace_KeckLRISr_400_8500_det2.fits',
                               ]):
        mstrace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
        assert chk_for_files(mstrace_file)

        # Instantiate
        tslits_dict, mstrace = traceslits.load_tslits(mstrace_file)
        spectrograph, traceSlits = instantiate(mstrace, tslits_dict)
        norig = tslits_dict['nslits']
        assert norig == nslit
        # Run me
        plate_scale = binning[0]*spectrograph.detector[det-1]['platescale']
        traceSlits.run(show=False, plate_scale=plate_scale, write_qa=False)
        # Test
        assert traceSlits.nslit == norig

@dev_suite_required
def test_chk_deimos_slits():
    """ This tests slit finding for DEIMOS """
    for nslit, binning, det, root in zip([27],
                                         [(1,1)],
                                         [3],
                                         ['MasterTrace_KeckDEIMOS_830G_8600_det3.fits',
                                          ]):
        mstrace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
        assert chk_for_files(mstrace_file)
        # Instantiate
        tslits_dict, mstrace = traceslits.load_tslits(mstrace_file)
        ipos = mstrace_file.find('det')  # Update to using the trace dict eventually
        spectrograph, traceSlits = instantiate(mstrace, tslits_dict, det=int(mstrace_file[ipos+3]))
        norig = tslits_dict['nslits']
        assert norig == nslit
        # Run me
        plate_scale = binning[0]*spectrograph.detector[det-1]['platescale']
        traceSlits.run(show=False, plate_scale=plate_scale, write_qa=False)
        # Test
        assert traceSlits.nslit == norig

