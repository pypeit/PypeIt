"""
Module to run tests on collate_1d code.
"""

import pytest
from collections import namedtuple
import os, os.path
import filecmp

from astropy.coordinates import Angle, SkyCoord
import numpy as np
import astropy.units as u
from astropy.io import fits
from pypeit import specobjs
from pypeit.spec2dobj import AllSpec2DObj
from pypeit.core.collate import collate_spectra_by_source, SourceObject
from pypeit.archive import ArchiveDir
from pypeit.scripts.collate_1d import find_spec2d_from_spec1d,find_slits_to_exclude, exclude_source_objects, extract_id, get_metadata_by_id, get_object_based_metadata
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit.pypmsgs import PypeItError

from pypeit.tests.tstutils import data_path, cooked_required

class mock_coadd:
    pass

def mock_get_instance():
    pass


class MockSpecObj:
    def __init__(self, MASKDEF_OBJNAME, MASKDEF_ID, DET, RA, DEC, SPAT_PIXPOS, NAME, OPT_COUNTS=None, BOX_COUNTS=None):
        self.MASKDEF_OBJNAME = MASKDEF_OBJNAME
        self.MASKDEF_ID = MASKDEF_ID
        self.DET = DET
        self.RA = RA
        self.SPAT_PIXPOS = SPAT_PIXPOS
        self.DEC = DEC
        self.NAME = NAME
        self.OPT_COUNTS = OPT_COUNTS
        self.BOX_COUNTS = BOX_COUNTS

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)


def mock_header(file):
    if file == 'spec1d_file1':
        return {'MJD': '58878.0',
                'PYP_SPEC': 'keck_deimos',
                'DISPNAME': '830G',
                'DECKER': 'Z6CL01B',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0',
                'FILENAME': 'DE.20100913.22358'}
    else:
        # Return a different decker to make sure it's properly ignored
        return {'MJD': '58879.0',
                'PYP_SPEC': 'keck_deimos',
                'DISPNAME': '830G',
                'DECKER': 'foo',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0',
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
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=1, RA=201.1517, DEC=27.3246, SPAT_PIXPOS=1234.0, NAME='SPAT1234_SLIT1234_DET01', OPT_COUNTS=np.zeros(100), BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=1, RA=201.1522, DEC=27.3250, SPAT_PIXPOS=1334.0, NAME='SPAT1334_SLIT1234_DET01', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object2',  MASKDEF_ID='3002', DET=2, RA=201.0051, DEC=27.2228, SPAT_PIXPOS=5334.0, NAME='SPAT5334_SLIT4934_DET02', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3233.0, NAME='SPAT3233_SLIT3235_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3232.0, NAME='SPAT3232_SLIT3235_DET03'),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3236.0, NAME='SPAT3236_SLIT3245_DET05', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=7, RA=201.1517, DEC=27.3246, SPAT_PIXPOS=1233.0, NAME='SPAT1233_SLIT1235_DET07', OPT_COUNTS=np.zeros(100), BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=7, RA=201.1520, DEC=27.3249, SPAT_PIXPOS=1336.0, NAME='SPAT1336_SLIT1235_DET07', OPT_COUNTS=np.zeros(100))]
        else:
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3234.0, NAME='SPAT3234_SLIT3236_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=3, RA=201.0052, DEC=27.2418, SPAT_PIXPOS=6250.0, NAME='SPAT6250_SLIT6235_DET03', BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=5, RA=201.0052, DEC=27.2418, SPAT_PIXPOS=6256.0, NAME='SPAT6256_SLIT6245_DET05', BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='4004', DET=5, RA=201.0056, DEC=27.2419, SPAT_PIXPOS=6934.0, NAME='SPAT6934_SLIT6245_DET05', BOX_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3237.0, NAME='SPAT3237_SLIT3246_DET05', OPT_COUNTS=np.zeros(100))]

    def __getitem__(self, idx):
        return self.specobjs[idx]

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

    file_list = ['spec1d_file1', 'spec1d_file2']
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

    assert source_list[3].spec1d_file_list == ['spec1d_file1','spec1d_file1','spec1d_file1','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT3233_SLIT3235_DET03', 'SPAT3232_SLIT3235_DET03', 'SPAT3236_SLIT3245_DET05', 'SPAT3234_SLIT3236_DET03', 'SPAT3237_SLIT3246_DET05']

    assert source_list[4].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03']
    
    assert source_list[5].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[5].spec_obj_list] == ['SPAT6256_SLIT6245_DET05']

    assert source_list[6].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[6].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']

    source_list = collate_spectra_by_source(uncollated_list, 10.0)

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

    assert source.coaddfile == 'J132436.41+271928.56_DEIMOS_20200130.fits'

    source2 = SourceObject(mock_sobjs.specobjs[0], mock_sobjs.header, 'spec1d_file1',
                           spectrograph, 'pixel')
    assert source2.coaddfile == 'SPAT1234_DEIMOS_20200130.fits'

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
 
    filtered_list = exclude_source_objects(uncollated_list, {'3003': 'Test Exclude`'}, par)
    assert [so.spec_obj_list[0].NAME for so in filtered_list] == ['SPAT1234_SLIT1234_DET01', 'SPAT5334_SLIT4934_DET02', 'SPAT1233_SLIT1235_DET07']
    assert [so.spec1d_file_list[0] for so in filtered_list] == ['spec1d_file1', 'spec1d_file1', 'spec1d_file1']

    par['collate1d']['exclude_serendip'] = False
    par['coadd1d']['ex_value'] = 'BOX'

    filtered_list = exclude_source_objects(uncollated_list, dict(), par)
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

def test_extract_id():

    mock_header_with_koaid = {'KOAID': 'DE.20121017.20019.fits'}
    mock_header_with_koaid_in_filename = {'FILENAME': 'DE.20121017.20019blah.fits'}
    mock_header_with_short_non_koa_filename = {'FILENAME':'short.fits'}
    mock_header_with_long_non_koa_filename  =  {'FILENAME': 'file_name_longer_than_17.fits'}

    assert extract_id(mock_header_with_koaid) == 'DE.20121017.20019.fits'
    assert extract_id(mock_header_with_koaid_in_filename) == 'DE.20121017.20019.fits'
    assert extract_id(mock_header_with_short_non_koa_filename) == 'short.fits'
    assert extract_id(mock_header_with_long_non_koa_filename) == 'file_name_longer_than_17.fits'

def test_get_metadata_by_id(monkeypatch):

    monkeypatch.setattr(fits, "getheader", mock_header)
    (metadata_rows, file_info, filename) = get_metadata_by_id(['DISPNAME', 'TESTKEY'], 'spec1d_file1')
    header = mock_header('spec1d_file1')
    assert len(metadata_rows) == 1
    assert metadata_rows[0] == [header['FILENAME'] + ".fits", 'spec1d_file1', '830G', None]
    assert file_info == 'spec1d_file1'
    assert filename == 'spec1d_file1'    

    spectrograph = load_spectrograph('keck_deimos')
    source_object = SourceObject(mock_specobjs(file_info).specobjs[0], header, file_info, spectrograph, 'ra/dec')

    assert (None, None, None) == get_metadata_by_id(['DISPNAME', 'TESTKEY'], source_object)

def test_get_object_based_metadata(monkeypatch):

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
    source_object.coaddfile = "/user/test/coaddfile.fits"
    for object in file1_objects[1:]:
        source_object.spec_obj_list.append(specobjs_file1.specobjs[object])
        source_object.spec1d_file_list.append(filenames[0])
        source_object.spec1d_header_list.append(header_file1)

    for object in file2_objects:
        source_object.spec_obj_list.append(specobjs_file2.specobjs[object])
        source_object.spec1d_file_list.append(filenames[1])
        source_object.spec1d_header_list.append(header_file2)

    (metadata_rows, file_info, filename) = get_object_based_metadata(['DISPNAME','MJD', 'GUIDFHWM'],
                                                                     ['MASKDEF_OBJNAME', 'NAME'],
                                                                     source_object)

    assert len(metadata_rows) == 4
    assert metadata_rows[0] == ['coaddfile.fits', 'object3', 'SPAT3233_SLIT3235_DET03', 'DE.20100913.22358.fits', '830G', '58878.0', None]
    assert metadata_rows[1] == ['coaddfile.fits', 'object3', 'SPAT3236_SLIT3245_DET05', 'DE.20100913.22358.fits', '830G', '58878.0', None]
    assert metadata_rows[2] == ['coaddfile.fits', 'object3', 'SPAT3234_SLIT3236_DET03', 'DE.20100914.12358.fits', '830G', '58879.0', None]
    assert metadata_rows[3] == ['coaddfile.fits', 'object3', 'SPAT3237_SLIT3246_DET05', 'DE.20100914.12358.fits', '830G', '58879.0', None]
    assert file_info == source_object.coaddfile
    assert filename == source_object.coaddfile
    
    assert (None, None, None) ==  get_object_based_metadata(['DISPNAME', 'MJD', 'GUIDFHWM'],
                                                            ['MASKDEF_OBJNAME', 'NAME'],
                                                            "afilename")
