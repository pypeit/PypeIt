"""
Tests for the spec2d display preparation functions.
"""
from astropy.io import fits
import numpy as np

from pypeit.display import spec2d_data


class FakeSpecObj:
    def __init__(self, name, det, trace, maskdef_objname=None, maskdef_extract=None,
                 hand_extract_flag=None):
        self.NAME = name
        self.DET = det
        self.MASKDEF_OBJNAME = maskdef_objname
        self.MASKDEF_EXTRACT = maskdef_extract
        self.hand_extract_flag = hand_extract_flag
        self._trace = trace

    def __getitem__(self, key):
        assert key == 'TRACE_SPAT', 'Only the trace is read from the object'
        return self._trace


class FakeSpecObjs:
    def __init__(self, objs):
        self.objs = objs
        self.DET = np.array([o.DET for o in objs])

    def __getitem__(self, i):
        return self.objs[i]


def test_trace_overlay_args():
    nspec = 50
    sobjs = FakeSpecObjs([
        FakeSpecObj('SPAT0100-SLIT0001-DET01', 'DET01', np.full(nspec, 100.)),
        FakeSpecObj('SPAT0300-SLIT0001-DET02', 'DET02', np.full(nspec, 300.)),
        FakeSpecObj('SPAT0200-SLIT0001-DET01', 'DET01', np.full(nspec, 200.),
                    maskdef_objname='gal1', maskdef_extract=True),
    ])
    args = spec2d_data.trace_overlay_args(sobjs, 'DET01')
    assert args['trace'].shape == (nspec, 2), 'Only traces on the detector should be kept'
    assert np.all(args['trace'][:,1] == 200.), 'Traces should be in columns'
    assert args['trc_name'][1].endswith('OBJNAME:gal1'), \
        'The slit-mask object name should be added to the trace name'
    assert args['maskdef_extr'].tolist() == [False, True], \
        'Slit-mask extractions should be flagged'
    assert spec2d_data.trace_overlay_args(sobjs, 'DET03') is None, \
        'A detector without objects should have no traces'
    assert spec2d_data.trace_overlay_args(None, 'DET01') is None, \
        'Without a spec1d file there should be no traces'


def test_spec2d_detname(tmp_path):
    ofile = tmp_path / 'spec2d.fits'
    hdr = fits.Header()
    hdr['HIERARCH ALLSPEC2D_DETS'] = 'MSC02,DET03'
    fits.PrimaryHDU(header=hdr).writeto(ofile)
    assert spec2d_data.spec2d_detname(ofile) == 'MSC02', \
        'Without --det, the first detector in the file should be used'
    assert spec2d_data.spec2d_detname(ofile, det='3') == 'DET03', \
        'A detector number should be converted to a detector name'
    assert spec2d_data.spec2d_detname(ofile, det='MSC01') == 'MSC01', \
        'A detector name should be used as given'
