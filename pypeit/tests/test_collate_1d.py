"""
Module to run tests on collate_1d code.
"""

import pytest
import os, os.path

import numpy as np
import astropy.units as u
from astropy.io import fits
from pypeit import specobjs
from pypeit.spec2dobj import AllSpec2DObj
from pypeit.core.collate import collate_spectra_by_source, SourceObject
from pypeit.scripts.collate_1d import find_spec2d_from_spec1d,find_slits_to_exclude, exclude_source_objects
from pypeit.scripts.collate_1d import flux, coadd, build_coadd_file_name, get_report_metadata, refframe_correction
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit.pypmsgs import PypeItError
from pypeit.images.detector_container import DetectorContainer
from pypeit import fluxcalibrate
from pypeit import coadd1d
from pypeit.core import wave


class MockSpecObj:
    def __init__(self, MASKDEF_OBJNAME, MASKDEF_ID, DET, RA, DEC, SPAT_PIXPOS, NAME, WAVE_RMS, OPT_FLAM=None, OPT_COUNTS=None, BOX_COUNTS=None, VEL_CORR=None):
        self.MASKDEF_OBJNAME = MASKDEF_OBJNAME
        self.MASKDEF_ID = MASKDEF_ID
        self.DET = DetectorContainer.get_name(DET)
        self.RA = RA
        self.SPAT_PIXPOS = SPAT_PIXPOS
        self.DEC = DEC
        self.NAME = NAME
        self.OPT_FLAM = OPT_FLAM
        self.OPT_COUNTS = OPT_COUNTS
        self.BOX_COUNTS = BOX_COUNTS
        self.WAVE_RMS = WAVE_RMS
        self.VEL_CORR = VEL_CORR

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)

    def apply_helio(self, *args, **kwargs):
        self.VEL_CORR = 1.0


def mock_header(file):
    if os.path.basename(file) == 'spec1d_file1':
        return {'MJD': '58878.0',
                'PYP_SPEC': 'keck_deimos',
                'DISPNAME': '830G',
                'DECKER': 'Z6CL01B',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0',
                'RA':       '201.1517', 
                'DEC':      '+27.3246',
                'FILENAME': 'DE.20100913.22358',
                'SEMESTER': '2019B',
                'PROGID':   'TEST1'}
    elif os.path.basename(file) == 'spec1d_file3':
        # Invalid header w/o MJD
        return {'PYP_SPEC': 'keck_deimos',
                'DISPNAME': '830G',
                'DECKER': 'Z6CL01B',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0',
                'FILENAME': 'DE.20100913.22358',
                'SEMESTER': '2019B',
                'PROGID':   'TEST1'}
    else:
        # Return a different decker to make sure it's properly ignored
        return {'MJD': '58879.0',
                'PYP_SPEC': 'keck_deimos',
                'DISPNAME': '830G',
                'DECKER': 'foo',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0',
                'RA':      '201.0052', 
                'DEC':     '+27.2418',
                'FILENAME': 'DE.20100914.12358'}

class MockSpecObjs:
    def __init__(self, file):
        self.header = mock_header(file)

        # specobjs for testing group_spectra_by_source
        # object1 entries test a match within a single file
        # object2 entry tests a spectra that doesn't match anything
        # object3 entries test a source that matches four spectra across two files
        # object3 also has a bad entry with no counts
        # The SERENDIP in file1 tests two serendipitous detections that match within a .0003d tolerance
        # object4 in file 2 along with the SERENDIP test an object that matches two spectra but narrowly misses
        # a SERENDIP detection with a tollerance of 0.0003d, but will match at 0.0004d.
        # object4 also has boxcar counts and no opt_counts

        if file == "spec1d_file1":
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=1, RA=201.1517, DEC=27.3246, SPAT_PIXPOS=1234.0, NAME='SPAT1234_SLIT1234_DET01', WAVE_RMS=0.01, OPT_COUNTS=np.zeros(100), OPT_FLAM=np.zeros(100), BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=1, RA=201.1522, DEC=27.3250, SPAT_PIXPOS=1334.0, NAME='SPAT1334_SLIT1234_DET01', WAVE_RMS=0.02, VEL_CORR=2.0, OPT_COUNTS=np.zeros(100), OPT_FLAM=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object2',  MASKDEF_ID='3002', DET=2, RA=201.0051, DEC=27.2228, SPAT_PIXPOS=5334.0, NAME='SPAT5334_SLIT4934_DET02', WAVE_RMS=0.01, OPT_COUNTS=np.zeros(100), OPT_FLAM=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3233.0, NAME='SPAT3233_SLIT3235_DET03', WAVE_RMS=0.01, OPT_COUNTS=np.zeros(100), OPT_FLAM=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3232.0, NAME='SPAT3232_SLIT3235_DET03', WAVE_RMS=0.03),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3236.0, NAME='SPAT3236_SLIT3245_DET05', WAVE_RMS=0.01, OPT_COUNTS=np.zeros(100), OPT_FLAM=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=7, RA=201.1517, DEC=27.3246, SPAT_PIXPOS=1233.0, NAME='SPAT1233_SLIT1235_DET07', WAVE_RMS=0.11, OPT_COUNTS=np.zeros(100), OPT_FLAM=np.zeros(100), BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=7, RA=201.1520, DEC=27.3249, SPAT_PIXPOS=1336.0, NAME='SPAT1336_SLIT1235_DET07', WAVE_RMS=0.01, OPT_COUNTS=np.zeros(100), OPT_FLAM=np.zeros(100))]
        elif file == "spec1d_file4":
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=None, DEC=None, SPAT_PIXPOS=3234.0, NAME='SPAT3234_SLIT3236_DET03', WAVE_RMS=0.01, VEL_CORR=2.0, OPT_FLAM=np.zeros(100), OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=3, RA=None, DEC=None, SPAT_PIXPOS=6250.0, NAME='SPAT6250_SLIT6235_DET03', WAVE_RMS=0.02, BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=5, RA=None, DEC=None, SPAT_PIXPOS=6256.0, NAME='SPAT6256_SLIT6245_DET05', WAVE_RMS=0.01, BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='4004', DET=5, RA=None, DEC=None, SPAT_PIXPOS=6934.0, NAME='SPAT6934_SLIT6245_DET05', WAVE_RMS=0.20, BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=None, DEC=None, SPAT_PIXPOS=3237.0, NAME='SPAT3237_SLIT3246_DET05', WAVE_RMS=0.01, OPT_COUNTS=np.zeros(100))]
        else:
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3234.0, NAME='SPAT3234_SLIT3236_DET03', WAVE_RMS=0.01, OPT_FLAM=np.zeros(100), OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=3, RA=201.0052, DEC=27.2418, SPAT_PIXPOS=6250.0, NAME='SPAT6250_SLIT6235_DET03', WAVE_RMS=0.02, BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=5, RA=201.0052, DEC=27.2418, SPAT_PIXPOS=6256.0, NAME='SPAT6256_SLIT6245_DET05', WAVE_RMS=0.01, BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='4004', DET=5, RA=201.0056, DEC=27.2419, SPAT_PIXPOS=6934.0, NAME='SPAT6934_SLIT6245_DET05', WAVE_RMS=0.20, BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3237.0, NAME='SPAT3237_SLIT3246_DET05', WAVE_RMS=0.01, OPT_COUNTS=np.zeros(100))]

    def __getitem__(self, idx):
        return self.specobjs[idx]

    def write_to_fits(self, *args, **kwargs):
        pass

def mock_specobjs(file):
    return MockSpecObjs(file)

def test_group_spectra_by_radec(monkeypatch):
    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", mock_specobjs)

    file_list = ['spec1d_file1', 'spec1d_file2']
    uncollated_list = SourceObject.build_source_objects(file_list, 'ra/dec')
    source_list = collate_spectra_by_source(uncollated_list, 0.0003, u.deg)

    assert len(source_list) == 6
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file1','spec1d_file1','spec1d_file1','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT3233_SLIT3235_DET03', 'SPAT3232_SLIT3235_DET03', 'SPAT3236_SLIT3245_DET05', 'SPAT3234_SLIT3236_DET03', 'SPAT3237_SLIT3246_DET05']

    assert source_list[4].spec1d_file_list == ['spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6256_SLIT6245_DET05']

    assert source_list[5].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[5].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']

    source_list = collate_spectra_by_source(uncollated_list, 1.44)

    assert len(source_list) == 5
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file1','spec1d_file1','spec1d_file1','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT3233_SLIT3235_DET03', 'SPAT3232_SLIT3235_DET03', 'SPAT3236_SLIT3245_DET05', 'SPAT3234_SLIT3236_DET03', 'SPAT3237_SLIT3246_DET05']

    assert source_list[4].spec1d_file_list == ['spec1d_file2','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6256_SLIT6245_DET05','SPAT6934_SLIT6245_DET05']

def test_group_spectra_by_pixel(monkeypatch):
    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", mock_specobjs)

    file_list = ['spec1d_file1', 'spec1d_file4']
    spectrograph = load_spectrograph('keck_deimos')
    # Test matching by pixel and that unit argument is ignored
    uncollated_list = SourceObject.build_source_objects(file_list, 'pixel')
    source_list = collate_spectra_by_source(uncollated_list, 5.0, u.deg)

    assert len(source_list) == 7
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file1','spec1d_file1','spec1d_file1','spec1d_file4','spec1d_file4']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT3233_SLIT3235_DET03', 'SPAT3232_SLIT3235_DET03', 'SPAT3236_SLIT3245_DET05', 'SPAT3234_SLIT3236_DET03', 'SPAT3237_SLIT3246_DET05']

    assert source_list[4].spec1d_file_list == ['spec1d_file4']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03']
    
    assert source_list[5].spec1d_file_list == ['spec1d_file4']
    assert [x.NAME for x in source_list[5].spec_obj_list] == ['SPAT6256_SLIT6245_DET05']

    assert source_list[6].spec1d_file_list == ['spec1d_file4']
    assert [x.NAME for x in source_list[6].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']

    source_list = collate_spectra_by_source(uncollated_list, 10.0)

    assert len(source_list) == 6
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file1','spec1d_file1','spec1d_file1','spec1d_file4','spec1d_file4']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT3233_SLIT3235_DET03', 'SPAT3232_SLIT3235_DET03', 'SPAT3236_SLIT3245_DET05', 'SPAT3234_SLIT3236_DET03', 'SPAT3237_SLIT3246_DET05']

    assert source_list[4].spec1d_file_list == ['spec1d_file4','spec1d_file4']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6256_SLIT6245_DET05']
    
    assert source_list[5].spec1d_file_list == ['spec1d_file4']
    assert [x.NAME for x in source_list[5].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']


def test_config_key_match():

    file_list = ['spec1d_file1', 'spec1d_file2']

    spectrograph = load_spectrograph('keck_deimos')

    header1 = {'dispname': '830G',
               'MJD': '58878.0',
               'PYP_SPEC': 'keck_deimos',
               'decker': 'Z6CL01B',
               'binning': '1,1'}
    header2 = {'dispname': '830G',
               'MJD': '58878.0',
               'PYP_SPEC': 'keck_deimos',
               'decker': 'foo',
               'binning': '1,1'}
    header3 = {'dispname': '830L',
               'MJD': '58878.0',
               'PYP_SPEC': 'keck_deimos',
               'decker': 'Z6CL01B',
               'binning': '1,1'}
    header4 = {'dispname': '830G',
               'MJD': '58878.0',
               'PYP_SPEC': 'keck_deimos',
               'decker': 'Z6CL01B',
               'binning': '1,1',
               'amp':  'foo'}
    header5 = {'dispname': '830G',
               'MJD': '58878.0',
               'decker': 'foo',
               'binning': '1,1'}
    header6 = {'dispname': '830G',
               'MJD': '58878.0',
               'PYP_SPEC': 'keck_nires',
               'decker': 'foo',
               'binning': '1,1'}

    sobjs = mock_specobjs('spec1d_file1')
    sobj = sobjs.specobjs[0]

    so1 = SourceObject(sobj, header1, 'spec1d_file1', spectrograph, 'ra/dec')
    so4 = SourceObject(sobj, header4, 'spec1d_file1', spectrograph, 'ra/dec')

    # Test that config keys match regardless of decker and even if neither
    # header has amp
    assert so1._config_key_match(header2) is True

    # Test that config keys don't match if dispname doesn't match
    assert so1._config_key_match(header3) is False

    # Test that config keys don't match if amp doesn't match
    # Also tests that the calling method match calls conifg_key_match correclty
    assert so1.match(sobj, header4, 3.0) is False

    # Test that config_key matching doesn't depend on order
    assert so4._config_key_match(header1) is False

    # Test that config keys don't match if there's no PYP_SPEC,
    # or the wrong PY_SPEC
    assert so1._config_key_match(header5) is False
    assert so1._config_key_match(header6) is False


def test_build_coadd_file_name():
    mock_sobjs = mock_specobjs('spec1d_file1')
    spectrograph = load_spectrograph('keck_deimos')

    source = SourceObject(mock_sobjs.specobjs[0], mock_sobjs.header, 'spec1d_file1',
                          spectrograph, 'ra/dec')

    assert build_coadd_file_name(source) == 'J132436.41+271928.56_DEIMOS_20200130_20200130.fits'

    source2 = SourceObject(mock_sobjs.specobjs[0], mock_sobjs.header, 'spec1d_file1',
                           spectrograph, 'pixel')
    assert build_coadd_file_name(source2) == 'SPAT1234_DEIMOS_20200130_20200130.fits'

def test_find_slits_to_exclude(monkeypatch):

    # Return mock data structure that contains what
    # find_slits_to_exclude needs from a AllSpec2DObj
    def mock_from_fits(*args, **kwargs):

        class MockSpec2dObj:
            def __init__(self, detector):
                if detector == 1:
                    self.slit_info = [  [    30,      0, 654112],
                                        [   813,      2, 654144],
                                        [   882,      0, 654125],
                                        [   952,      2, 654145],
                                        [  1024,      0, 654057]]

                else:
                    self.slit_info = [  [    22,      0, 654101],
                                        [   397,      2, 654143],
                                        [   463,      0, 654073],
                                        [   531,      6, 654142],
                                        [   601,      0, 654071],
                                        [  1101,      8, 654104],
                                        [  1840,     16, 654046],
                                        [  1932,      0, 654128]]
            # So that sobj2d['slits'].slit_info works
            def __getitem__(self, key):
                return self


        class MockAllSpec2dObj:
            def __init__(self):
                self.detectors = [1, 4]
        
            def __getitem__(self, key):
                return MockSpec2dObj(key)

        return MockAllSpec2dObj()

    par = pypeitpar.PypeItPar()
    par['collate1d'] = pypeitpar.Collate1DPar()
    par['collate1d']['exclude_slit_trace_bm'] = ['BOXSLIT', 'USERIGNORE']

    monkeypatch.setattr(AllSpec2DObj, 'from_fits', mock_from_fits)

    exclude_map = find_slits_to_exclude(['spec2d_file_name'], par)
    assert len(exclude_map) == 4
    assert exclude_map[654142] == {'BOXSLIT','USERIGNORE'}
    assert exclude_map[654143] == {'BOXSLIT'}
    assert exclude_map[654144] == {'BOXSLIT'}
    assert exclude_map[654145] == {'BOXSLIT'}

    par['collate1d']['exclude_slit_trace_bm'] = 'USERIGNORE'
    exclude_map = find_slits_to_exclude(['spec2d_file_name'], par)
    assert len(exclude_map) == 1
    assert exclude_map[654142] == {'USERIGNORE'}

def test_exclude_source_objects(monkeypatch):
    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", mock_specobjs)

    file_list = ['spec1d_file1', 'spec1d_file2']
    uncollated_list = SourceObject.build_source_objects(file_list, 'ra/dec')
    par = pypeitpar.PypeItPar()
    par['collate1d']['exclude_serendip'] = True
    par['collate1d']['wv_rms_thresh'] = 0.1
    filtered_list, excluded_msgs = exclude_source_objects(uncollated_list, {'3003': 'Test Exclude`'}, par)
    assert [so.spec_obj_list[0].NAME for so in filtered_list] == ['SPAT1234_SLIT1234_DET01', 'SPAT5334_SLIT4934_DET02']
    assert [so.spec1d_file_list[0] for so in filtered_list] == ['spec1d_file1', 'spec1d_file1']

    par['collate1d']['exclude_serendip'] = False
    par['coadd1d']['ex_value'] = 'BOX'
    par['collate1d']['wv_rms_thresh'] = None

    filtered_list, excluded_msgs = exclude_source_objects(uncollated_list, dict(), par)
    assert [so.spec_obj_list[0].NAME for so in filtered_list] == ['SPAT1234_SLIT1234_DET01', 'SPAT1233_SLIT1235_DET07', 'SPAT6250_SLIT6235_DET03', 'SPAT6256_SLIT6245_DET05', 'SPAT6934_SLIT6245_DET05']
    assert [so.spec1d_file_list[0] for so in filtered_list] == ['spec1d_file1', 'spec1d_file1', 'spec1d_file2', 'spec1d_file2', 'spec1d_file2' ]


def test_find_spec2d_from_spec1d(tmp_path):

    # Write dummy files to avoid dependency on Cooked
    spec1d_names = [str(tmp_path / 'spec1d_test-name.1.fits'), str(tmp_path / 'spec1d_test-name.2.fits')]
    spec2d_names = [str(tmp_path / 'spec2d_test-name.1.fits'), str(tmp_path / 'spec2d_test-name.2.fits')]

    for file in spec1d_names + spec2d_names:
        with open(file, "w") as f:
            print("Dummy fits file", file=f)
    
    result = find_spec2d_from_spec1d(spec1d_names)

    assert len(result) == 2
    assert spec2d_names[0] in result
    assert spec2d_names[1] in result

    # Test when one of the spec1d files does not have a corresponding spec2d
    spec1d_names.append(str(tmp_path/'spec1d_test-name.3.fits'))
    with pytest.raises(PypeItError):
        result = find_spec2d_from_spec1d(spec1d_names)

def test_get_report_metadata(monkeypatch):

    spectrograph = load_spectrograph('keck_deimos')
    filenames = ['spec1d_file1', 'spec1d_file2']
    specobjs_file1 = mock_specobjs(filenames[0])
    specobjs_file2 = mock_specobjs(filenames[1])
    header_file1 = mock_header(filenames[0])
    header_file2 = mock_header(filenames[1])
    file1_objects = [3,  # 'SPAT3233_SLIT3235_DET03'
                     5,] # 'SPAT3236_SLIT3245_DET05'
                     
    file2_objects = [0,  # 'SPAT3234_SLIT3236_DET03'
                     4,] # 'SPAT3237_SLIT3246_DET05'
    
    source_object = SourceObject(specobjs_file1.specobjs[file1_objects[0]], 
                                 header_file1, 
                                 filenames[0], 
                                 spectrograph, 
                                 'ra/dec')

    for object in file1_objects[1:]:
        source_object.spec_obj_list.append(specobjs_file1.specobjs[object])
        source_object.spec1d_file_list.append(filenames[0])
        source_object.spec1d_header_list.append(header_file1)

    for object in file2_objects:
        source_object.spec_obj_list.append(specobjs_file2.specobjs[object])
        source_object.spec1d_file_list.append(filenames[1])
        source_object.spec1d_header_list.append(header_file2)

    monkeypatch.setattr(fits, "getheader", mock_header)

    (metadata_rows, files_to_copy) = get_report_metadata(['DISPNAME','MJD', 'GUIDFHWM'],
                                                         ['MASKDEF_OBJNAME', 'NAME'],
                                                         source_object)
    dest_file = 'J132500.41+271959.88_DEIMOS_20200130_20200131.fits'
    assert len(metadata_rows) == 4
    assert metadata_rows[0] == [dest_file, 'object3', 'SPAT3233_SLIT3235_DET03', 'spec1d_file1', '830G', '58878.0', None]
    assert metadata_rows[1] == [dest_file, 'object3', 'SPAT3236_SLIT3245_DET05', 'spec1d_file1', '830G', '58878.0', None]
    assert metadata_rows[2] == [dest_file, 'object3', 'SPAT3234_SLIT3236_DET03', 'spec1d_file2', '830G', '58879.0', None]
    assert metadata_rows[3] == [dest_file, 'object3', 'SPAT3237_SLIT3246_DET05', 'spec1d_file2', '830G', '58879.0', None]
    assert files_to_copy is None
    
    assert (None, None) ==  get_report_metadata(['DISPNAME', 'MJD', 'GUIDFHWM'],
                                                ['MASKDEF_OBJNAME', 'NAME'],
                                                "afilename")

def test_flux(monkeypatch):
    
    def mock_get_header(file):
        if file == "fail_grating.fits":
            return {"DISPNAME": "Unknown"}
        else: 
            return {"DISPNAME": "600ZD" }

    def mock_get_flux_calib_instance(spec1d_files, sens_files, par):
        if spec1d_files[0] == "fail_flux.fits":
            raise PypeItError("Test failure")
        else:
            # The collate_1d caller doesn't use the output, it just
            # depends on the side effect of fluxing
            return None 

    # Test success
    with monkeypatch.context() as m:
        monkeypatch.setattr(fits, "getheader", mock_get_header)
        monkeypatch.setattr(fluxcalibrate.FluxCalibrate, "get_instance", mock_get_flux_calib_instance)

        spectrograph = load_spectrograph('keck_deimos')
        # This could break if we start supporting it
        unsupported_spectrograph = load_spectrograph('shane_kast_red')

        par = pypeitpar.PypeItPar()
        par['fluxcalib'] = pypeitpar.FluxCalibratePar()

        # Test success
        failed_messages = []
        flux(par, spectrograph, ["no_fail.fits"], failed_messages)

        assert len(failed_messages) == 0

        # Test failure due to no archived sensfunc
        flux(par, spectrograph, ["fail_grating.fits"], failed_messages)
        assert failed_messages[0] == "Could not find archived sensfunc to flux fail_grating.fits, skipping it."

        failed_messages = []

        # Test failure in fluxing
        flux(par, spectrograph, ["fail_flux.fits"], failed_messages)
        assert failed_messages[0] == "Failed to flux calibrate fail_flux.fits, skipping it."
        
        # Test failure due to unsupported spectrograph
        with pytest.raises(PypeItError):
            flux(par, unsupported_spectrograph, ["600ZD_file.fits"], failed_messages)            

def test_coadd(monkeypatch):

    class mock_coadd:

        def run(self):
            return

        def save(self, file):
            return

    def mock_get_instance(*args, **kwargs):
        return mock_coadd()

    with monkeypatch.context() as m:
        monkeypatch.setattr(coadd1d.CoAdd1D, "get_instance", mock_get_instance)
        par = pypeitpar.PypeItPar()
        par['collate1d'] = pypeitpar.Collate1DPar()
        par['coadd1d'] = pypeitpar.Coadd1DPar()
        spectrograph = load_spectrograph('keck_deimos')

        filenames = ['spec1d_file1', 'spec1d_file2']
        specobjs_file1 = mock_specobjs(filenames[0])
        specobjs_file2 = mock_specobjs(filenames[1])
        header_file1 = mock_header(filenames[0])
        header_file2 = mock_header(filenames[1])

        # Both source object 1's SpecObj objects will have OPT_FLAM
        file1_objects = [3,  # 'SPAT3233_SLIT3235_DET03'
                         5,] # 'SPAT3236_SLIT3245_DET05'
                            
        # Only one of source object 2's SpecObj objects will have OPT_FLAM
        file2_objects = [0,  # 'SPAT3234_SLIT3236_DET03'
                         4,] # 'SPAT3237_SLIT3246_DET05'

        source_object1 = SourceObject(specobjs_file1.specobjs[file1_objects[0]], 
                                      header_file1, 
                                      filenames[0], 
                                      spectrograph, 
                                      'ra/dec')

        for object in file1_objects[1:]:
            source_object1.spec_obj_list.append(specobjs_file1.specobjs[object])
            source_object1.spec1d_file_list.append(filenames[0])
            source_object1.spec1d_header_list.append(header_file1)

        source_object2 = SourceObject(specobjs_file2.specobjs[file2_objects[0]], 
                                      header_file2, 
                                      filenames[1], 
                                      spectrograph, 
                                      'ra/dec')

        for object in file2_objects:
            source_object2.spec_obj_list.append(specobjs_file2.specobjs[object])
            source_object2.spec1d_file_list.append(filenames[1])
            source_object2.spec1d_header_list.append(header_file2)

        # Test with out using fluxed data
        par['collate1d']['ignore_flux'] = True
        par['coadd1d']['flux_value'] = True
        test_file = "test_coadd_name"
        coadd(par, test_file, source_object1)

        assert par['coadd1d']['flux_value'] == False

        # Test using fluxed data
        par['collate1d']['ignore_flux'] = False
        par['coadd1d']['flux_value'] = False
        coadd(par, test_file, source_object1)

        assert par['coadd1d']['flux_value'] == True

        # Test not using fluxed data because not all SpecObj objects had flux data
        par['collate1d']['ignore_flux'] = False
        par['coadd1d']['flux_value'] = True
        coadd(par, test_file, source_object2)

        assert par['coadd1d']['flux_value'] == False

def test_refframe_correction(monkeypatch):
    def mock_geomotion_correct(*args, **kwargs):
        return 1.0, 1.0

    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", mock_specobjs)
    monkeypatch.setattr(wave, "geomotion_correct", mock_geomotion_correct)

    par = pypeitpar.PypeItPar()
    par['collate1d'] = pypeitpar.Collate1DPar()
    par['collate1d']['refframe'] = 'heliocentric'
    spectrograph = load_spectrograph('keck_deimos')

    # Test that should fail due to no RA/DEC nor mjd in header
    spec1d_failure_msgs = []
    spec1d_files = ["spec1d_file3"]
    refframe_correction(par, spectrograph, spec1d_files, spec1d_failure_msgs)
    assert len(spec1d_failure_msgs) == 1
    assert spec1d_failure_msgs[0].startswith('Failed to perform heliocentric correction on spec1d_file3')

    # Test where onf of the SpecObjs already has a VEL_CORR that should not be overwritten
    spec1d_failure_msgs = []
    spec1d_files = ["spec1d_file4"]

    # Test where one VEL_CORR is already set, and the SpecObj objects have no RA/DEC so the header RA/DEC must be used instead
    sobjs = MockSpecObjs("spec1d_file4")
    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", lambda x: sobjs)

    refframe_correction(par, spectrograph, spec1d_files, spec1d_failure_msgs)
    assert len(spec1d_failure_msgs) == 1
    assert spec1d_failure_msgs[0].startswith('Not performing heliocentric correction for spec1d_file4 object SPAT3234_SLIT3236_DET03 because it has already been corrected')
    assert sobjs[0].VEL_CORR == 2.0 # Original value, should not have been overwritten
    assert sobjs[1].VEL_CORR == 1.0 # New value, from apply_helio


