"""
Module to run tests on TraceSlits class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.tests.tstutils import cooked_required
from pypeit.spectrographs import util

from pypeit import edgetrace

@cooked_required
def test_addrm_slit():
    """ This tests the add and remove methods for user-supplied slit fussing. """

    # Check for files
    trace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace',
                              'MasterEdges_KeckLRISr_400_8500_det1.fits.gz')
    assert os.path.isfile(trace_file), 'Trace file does not exist!'
    # Instantiate
    edges = edgetrace.EdgeTraceSet.from_file(trace_file)
    assert edges.is_synced, 'Error in constructed Cooked trace file: not synced.'
    nslits = edges.ntrace//2

    #  Add dummy slit
    #  y_spec, left edge, right edge on image
    add_user_slits = [[1024, 140, 200]]
    edges.add_user_traces(add_user_slits)
    assert edges.is_synced, 'Added slit caused object to de-synchronize.'
    assert edges.ntrace//2 == nslits + 1, 'Number of slits incorrect.'
    xcen0 = np.median(np.mean(edges.edge_fit[:,:2], axis=0))
    assert np.absolute(xcen0-170) < 3, 'Unexpected slit center coordinate'

    # Remove it
    rm_user_slits = [[1024, 170]]
    edges.rm_user_traces(rm_user_slits)
    assert edges.ntrace//2 == nslits, 'Did not remove trace.'


# TODO: Can we (and is it useful to) get these tests back?

'''
@cooked_required
def test_chk_kast_slits():
    """ This tests finding the longslit for Kast blue """
    # (valuable) and that it didn't change between when we run the dev
    # suite and when we run the unit tests (not so valuable).  We should be
    # doing the more difficult thing of checking that the tracing is
    # providing consistent results...
    # Red, blue
    for root in ['MasterEdges_ShaneKastred_600_7500_d55_ret.fits.gz',
                 'MasterEdges_ShaneKastblue_600_4310_d55.fits.gz']:
        # Load
        trace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
        assert os.path.isfile(trace_file), 'Trace file does not exist!'
        edges = edgetrace.EdgeTraceSet.from_file(trace_file)
        assert edges.is_synced, 'Error in constructed Cooked trace file: not synced.'
        nslits = edges.ntrace//2
        # Rerun (but don't write QA plots)
        edges.qa_path = None
        edges.auto_trace(bpm=edges.bpm, det=edges.det, binning=edges.binning)
        assert edges.ntrace//2 == nslits, 'Did not regain the same slits!'

@cooked_required
def test_chk_lris_blue_slits():
    """ This tests slit finding for LRISb """
    for root in ['MasterEdges_KeckLRISb_long_600_4000_det2.fits.gz',
                 'MasterEdges_KeckLRISb_multi_600_4000_det2.fits.gz']:
        trace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
        assert os.path.isfile(trace_file), 'Trace file does not exist!'
        edges = edgetrace.EdgeTraceSet.from_file(trace_file)
        assert edges.is_synced, 'Error in constructed Cooked trace file: not synced.'
        nslits = edges.ntrace//2
        # Rerun (but don't write QA plots)
        edges.qa_path = None
        edges.auto_trace(bpm=edges.bpm, det=edges.det, binning=edges.binning)
        assert edges.ntrace//2 == nslits, 'Did not regain the same slits!'

@cooked_required
def test_chk_lris_red_slits():
    """ This tests slit finding for LRISr """
    for root in ['MasterEdges_KeckLRISr_long_600_7500_d560.fits.gz',
                 'MasterEdges_KeckLRISr_400_8500_det1.fits.gz',
                 'MasterEdges_KeckLRISr_400_8500_det2.fits.gz']:
        trace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
        assert os.path.isfile(trace_file), 'Trace file does not exist!'
        edges = edgetrace.EdgeTraceSet.from_file(trace_file)
        assert edges.is_synced, 'Error in constructed Cooked trace file: not synced.'
        nslits = edges.ntrace//2
        # Rerun (but don't write QA plots)
        edges.qa_path = None
        edges.auto_trace(bpm=edges.bpm, det=edges.det, binning=edges.binning)
        assert edges.ntrace//2 == nslits, 'Did not regain the same slits!'

@cooked_required
def test_chk_deimos_slits():
    """ This tests slit finding for DEIMOS """
    root = 'MasterEdges_KeckDEIMOS_830G_8500_det3.fits.gz'
    trace_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Trace', root)
    assert os.path.isfile(trace_file), 'Trace file does not exist!'
    edges = edgetrace.EdgeTraceSet.from_file(trace_file)
    assert edges.is_synced, 'Error in constructed Cooked trace file: not synced.'
    nslits = edges.ntrace//2
    # Rerun (but don't write QA plots)
    edges.qa_path = None
    edges.auto_trace(bpm=edges.bpm, det=edges.det, binning=edges.binning)
    assert edges.ntrace//2 == nslits, 'Did not regain the same slits!'

'''

