"""
Module to run tests on collate_1d script.
"""

import pytest
from collections import namedtuple
import os, os.path
import filecmp

from astropy.coordinates import Angle, SkyCoord
import numpy as np
import astropy.units as u

from pypeit import specobjs
from pypeit.spec2dobj import AllSpec2DObj
from pypeit.scripts.collate_1d import group_spectra_by_source, SourceObject
from pypeit.scripts.collate_1d import ArchiveDir, find_slits_to_exclude, find_spec2d_from_spec1d
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit.pypmsgs import PypeItError

from pypeit.tests.tstutils import data_path, cooked_required

class mock_coadd:
    pass

def mock_get_instance():
    pass


class MockSpecObj:
    def __init__(self, MASKDEF_OBJNAME, MASKDEF_ID, DET, RA, DEC, SPAT_PIXPOS, NAME, OPT_COUNTS):
        self.MASKDEF_OBJNAME = MASKDEF_OBJNAME
        self.MASKDEF_ID = MASKDEF_ID
        self.DET = DET
        self.RA = RA
        self.SPAT_PIXPOS = SPAT_PIXPOS
        self.DEC = DEC
        self.NAME = NAME
        self.OPT_COUNTS = OPT_COUNTS

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)


def dummy_header(file):
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
        return {'MJD': '58878.0',
                'PYP_SPEC': 'keck_deimos',
                'DISPNAME': '830G',
                'DECKER': 'foo',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0',
                'FILENAME': 'DE.20100913.22358'}

class MockSpecObjs:
    def __init__(self, file):
        self.header = dummy_header(file)

        # specobjs for testing group_spectra_by_source
        # object1 entries test a match within a single file
        # object2 entry tests a spectra that doesn't match anything
        # object3 entries test a source that matches four spectra across two files
        # The SERENDIP in file1 tests two serendipitous detections that match within a .0003d tolerance
        # object4 in file 2 along with the SERENDIP test an object that matches two spectra but narrowly misses
        # a SERENDIP detection with a tollerance of 0.0003d, but will match at 0.0004d.

        if file == "spec1d_file1":
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=1, RA=201.1517, DEC=27.3246, SPAT_PIXPOS=1234.0, NAME='SPAT1234_SLIT1234_DET01', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=1, RA=201.1522, DEC=27.3250, SPAT_PIXPOS=1334.0, NAME='SPAT1334_SLIT1234_DET01', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object2',  MASKDEF_ID='3002', DET=2, RA=201.0051, DEC=27.2228, SPAT_PIXPOS=5334.0, NAME='SPAT5334_SLIT4934_DET02', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3233.0, NAME='SPAT3233_SLIT3235_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3236.0, NAME='SPAT3236_SLIT3245_DET05', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=7, RA=201.1517, DEC=27.3246, SPAT_PIXPOS=1233.0, NAME='SPAT1233_SLIT1235_DET07', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=7, RA=201.1520, DEC=27.3249, SPAT_PIXPOS=1336.0, NAME='SPAT1336_SLIT1235_DET07', OPT_COUNTS=np.zeros(100))]
        else:
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3234.0, NAME='SPAT3234_SLIT3236_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=3, RA=201.0052, DEC=27.2418, SPAT_PIXPOS=6250.0, NAME='SPAT6250_SLIT6235_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=5, RA=201.0052, DEC=27.2418, SPAT_PIXPOS=6256.0, NAME='SPAT6256_SLIT6245_DET05', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='4004', DET=5, RA=201.0056, DEC=27.2419, SPAT_PIXPOS=6934.0, NAME='SPAT6934_SLIT6245_DET05', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, SPAT_PIXPOS=3237.0, NAME='SPAT3237_SLIT3246_DET05', OPT_COUNTS=np.zeros(100))]


def mock_specobjs(file):
    return MockSpecObjs(file)

def test_group_spectra_by_radec(monkeypatch):
    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", mock_specobjs)

    file_list = ['spec1d_file1', 'spec1d_file2']

    source_list = group_spectra_by_source(file_list, dict(), 'ra/dec', 0.0003, u.deg)

    assert len(source_list) == 6
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file1','spec1d_file1','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT3233_SLIT3235_DET03','SPAT3236_SLIT3245_DET05','SPAT3234_SLIT3236_DET03','SPAT3237_SLIT3246_DET05']

    assert source_list[4].spec1d_file_list == ['spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6256_SLIT6245_DET05']

    assert source_list[5].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[5].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']

    exclude_map = {'3003': 'TEST_FLAG'}
    source_list = group_spectra_by_source(file_list, exclude_map, 'ra/dec', 1.44)

    assert len(source_list) == 4
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file2','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6256_SLIT6245_DET05','SPAT6934_SLIT6245_DET05']

def test_group_spectra_by_pixel(monkeypatch):
    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", mock_specobjs)

    file_list = ['spec1d_file1', 'spec1d_file2']
    spectrograph = load_spectrograph('keck_deimos')
    # Test matching by pixel and that unit argument is ignored
    source_list = group_spectra_by_source(file_list, dict(), 'pixel', 5.0, u.deg)

    assert len(source_list) == 7
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file1','spec1d_file1','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT3233_SLIT3235_DET03','SPAT3236_SLIT3245_DET05','SPAT3234_SLIT3236_DET03','SPAT3237_SLIT3246_DET05']

    assert source_list[4].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03']
    
    assert source_list[5].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[5].spec_obj_list] == ['SPAT6256_SLIT6245_DET05']

    assert source_list[6].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[6].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']

    exclude_map = {'3003': 'TEST_FLAG'}
    source_list = group_spectra_by_source(file_list, exclude_map, 'pixel', 10.0)

    assert len(source_list) == 5
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6256_SLIT6245_DET05']
    
    assert source_list[4].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']


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
    assert so1._config_key_match(header4) is False

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
                                        [   531,      2, 654142],
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
    par['collate1d']['slit_exclude_flags'] = 'BOXSLIT'

    monkeypatch.setattr(AllSpec2DObj, 'from_fits', mock_from_fits)

    exclude_map = find_slits_to_exclude(['spec2d_file_name'], par)
    assert len(exclude_map) == 4
    assert exclude_map[654142] == {'BOXSLIT'}
    assert exclude_map[654143] == {'BOXSLIT'}
    assert exclude_map[654144] == {'BOXSLIT'}
    assert exclude_map[654145] == {'BOXSLIT'}

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


@cooked_required
def test_archive(tmp_path):
    cooked_sci_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science')

    spec1d_name = 'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_2010Sep13T061231.334.fits'
    spec2d_name = 'spec2d_KB.20191219.56886-BB1245p4238_KCWI_2019Dec19T154806.538.fits'
    spec1d_files = [os.path.join(cooked_sci_dir, spec1d_name)]
    spec2d_files = [os.path.join(cooked_sci_dir, spec2d_name)]

    # Don't actually need a real coadd file, just any file to see if it gets copied to the archive directory    
    coadd_name = 'spec1d_b24-Feige66_KASTb_2015May20T041246.960.fits'
    coadd_file = os.path.join(cooked_sci_dir, coadd_name)

    # Build source for add_coadd_sources
    mock_sobjs = mock_specobjs('spec1d_file1')
    spectrograph = load_spectrograph('keck_deimos')
    source_list = [SourceObject(mock_sobjs.specobjs[0], mock_sobjs.header, 'spec1d_file1',
                                spectrograph,'ra/dec')]

    source_list[0].coaddfile = coadd_file
    archive_dir = str(tmp_path)
    archive = ArchiveDir(archive_dir)
    archive.add_files(spec1d_files)
    archive.add_files(spec2d_files)
    archive.add_coadd_sources(source_list)
    archive.save()
    assert(os.path.exists(os.path.join(archive_dir, spec1d_name)))
    assert(os.path.exists(os.path.join(archive_dir, spec2d_name)))
    assert(os.path.exists(os.path.join(archive_dir, coadd_name)))

    # Use filecmp to compare .dat files
    good_dat_dir = data_path('ipac')
    good_by_id_dat =  os.path.join(good_dat_dir, 'by_id_meta.dat')
    good_by_obj_dat = os.path.join(good_dat_dir, 'by_object_meta.dat')
    test_by_id_dat =  os.path.join(archive_dir, 'by_id_meta.dat')
    test_by_obj_dat = os.path.join(archive_dir, 'by_object_meta.dat')
    assert(filecmp.cmp(good_by_id_dat, test_by_id_dat, shallow=False))
    assert(filecmp.cmp(good_by_obj_dat, test_by_obj_dat, shallow=False))
 
    # Now test opening and adding to an existing archive
    archive2 = ArchiveDir(archive_dir)
    spec2d_name = 'spec2d_KB.20191219.57662-BB1245p4238_KCWI_2019Dec19T160102.755.fits'
    spec2d_list = [os.path.join(cooked_sci_dir, spec2d_name)]

    coadd_name = 'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits'
    mock_sobjs = mock_specobjs('spec1d_file2')
    spectrograph = load_spectrograph('keck_deimos')
    source_list = [SourceObject(mock_sobjs.specobjs[0], mock_sobjs.header, 'spec1d_file1',
                                spectrograph,'ra/dec')]
    source_list[0].coaddfile = os.path.join(cooked_sci_dir, coadd_name)

    archive2.add_files(spec2d_list)
    archive2.add_coadd_sources(source_list)
    archive2.save()
    assert(os.path.exists(os.path.join(archive_dir, spec2d_name)))
    assert(os.path.exists(os.path.join(archive_dir, coadd_name)))

    good_by_id_dat =  os.path.join(good_dat_dir, 'by_id_meta2.dat')
    good_by_obj_dat = os.path.join(good_dat_dir, 'by_object_meta2.dat')
    assert(filecmp.cmp(good_by_id_dat, test_by_id_dat, shallow=False))
    assert(filecmp.cmp(good_by_obj_dat, test_by_obj_dat, shallow=False))
