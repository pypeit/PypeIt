"""
Module to run tests on collate_1d script.
"""

import pytest
from collections import namedtuple
import os, os.path
import filecmp

from astropy.coordinates import Angle, SkyCoord
import numpy as np

from pypeit import specobjs
from pypeit.scripts.collate_1d import group_spectra_by_source, SourceObject, build_coadd_file_name
from pypeit.scripts.collate_1d import  config_key_match, ADAPArchive, find_slits_to_exclude, find_spec2d_from_spec1d
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit.pypmsgs import PypeItError

from pypeit.tests.tstutils import cooked_required, data_path

class mock_coadd:
    pass

def mock_get_instance():
    pass


class MockSpecObj:
    def __init__(self, MASKDEF_OBJNAME, MASKDEF_ID, DET, RA, DEC, NAME, OPT_COUNTS):
        self.MASKDEF_OBJNAME = MASKDEF_OBJNAME
        self.MASKDEF_ID = MASKDEF_ID
        self.DET = DET
        self.RA = RA
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
                'DISPNAME': '830G',
                'DECKER': 'Z6CL01B',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0'}
    else:
        # Return a different decker to make sure it's properly ignored
        return {'MJD': '58878.0',
                'DISPNAME': '830G',
                'DECKER': 'foo',
                'BINNING': '1,1',
                'AIRMASS': '1.0',
                'EXPTIME': '1200.0'}

class MockSpecObjs:
    def __init__(self, file):
        self.header = dummy_header('spec1d_file1')

        # specobjs for testing group_spectra_by_source
        # object1 entries test a match within a single file
        # object2 entry tests a spectra that doesn't match anything
        # object3 entries test a source that matches four spectra across two files
        # The SERENDIP in file1 tests two serendipitous detections that match within a .0003d threshold
        # object4 in file 2 along with the SERENDIP test an object that matches two spectra but narrowly misses
        # a SERENDIP detection with a threshold of 0.0003d, but will match at 0.0004d.

        if file == "spec1d_file1":
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=1, RA=201.1517, DEC=27.3246, NAME='SPAT1234_SLIT1234_DET01', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=1, RA=201.1522, DEC=27.3250, NAME='SPAT1334_SLIT1234_DET01', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object2',  MASKDEF_ID='3002', DET=2, RA=201.0051, DEC=27.2228, NAME='SPAT5334_SLIT4934_DET02', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, NAME='SPAT3233_SLIT3235_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, NAME='SPAT3236_SLIT3245_DET05', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object1',  MASKDEF_ID='1001', DET=7, RA=201.1517, DEC=27.3246, NAME='SPAT1233_SLIT1235_DET07', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='1001', DET=7, RA=201.1520, DEC=27.3249, NAME='SPAT1336_SLIT1235_DET07', OPT_COUNTS=np.zeros(100))]
        else:
            self.specobjs = [MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=3, RA=201.2517, DEC=27.3333, NAME='SPAT3234_SLIT3236_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=3, RA=201.0052, DEC=27.2418, NAME='SPAT6250_SLIT6235_DET03', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object4',  MASKDEF_ID='4004', DET=5, RA=201.0052, DEC=27.2418, NAME='SPAT6255_SLIT6245_DET05', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='SERENDIP', MASKDEF_ID='4004', DET=5, RA=201.0056, DEC=27.2419, NAME='SPAT6934_SLIT6245_DET05', OPT_COUNTS=np.zeros(100)),
                             MockSpecObj(MASKDEF_OBJNAME='object3',  MASKDEF_ID='3003', DET=5, RA=201.2517, DEC=27.3333, NAME='SPAT3237_SLIT3246_DET05', OPT_COUNTS=np.zeros(100))]


def mock_specobjs(file):
    return MockSpecObjs(file)

def test_group_spectra_by_source(monkeypatch):
    monkeypatch.setattr(specobjs.SpecObjs, "from_fitsfile", mock_specobjs)

    file_list = ['spec1d_file1', 'spec1d_file2']
    spectrograph = load_spectrograph('keck_deimos')

    source_list = group_spectra_by_source(file_list, spectrograph, dict(), Angle('.0003d'))

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
    assert [x.NAME for x in source_list[4].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6255_SLIT6245_DET05']

    assert source_list[5].spec1d_file_list == ['spec1d_file2']
    assert [x.NAME for x in source_list[5].spec_obj_list] == ['SPAT6934_SLIT6245_DET05']

    exclude_map = {'3003': 'TEST_FLAG'}
    source_list = group_spectra_by_source(file_list, spectrograph, exclude_map, Angle('.0004d'))

    assert len(source_list) == 4
    assert source_list[0].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[0].spec_obj_list] == ['SPAT1234_SLIT1234_DET01','SPAT1233_SLIT1235_DET07']

    assert source_list[1].spec1d_file_list == ['spec1d_file1','spec1d_file1']
    assert [x.NAME for x in source_list[1].spec_obj_list] == ['SPAT1334_SLIT1234_DET01','SPAT1336_SLIT1235_DET07']

    assert source_list[2].spec1d_file_list == ['spec1d_file1']
    assert [x.NAME for x in source_list[2].spec_obj_list] == ['SPAT5334_SLIT4934_DET02']

    assert source_list[3].spec1d_file_list == ['spec1d_file2','spec1d_file2','spec1d_file2']
    assert [x.NAME for x in source_list[3].spec_obj_list] == ['SPAT6250_SLIT6235_DET03','SPAT6255_SLIT6245_DET05','SPAT6934_SLIT6245_DET05']



def test_config_key_match():
    spectrograph = load_spectrograph('keck_deimos')

    header1 = {'dispname': '830G',
               'decker': 'Z6CL01B',
               'binning': '1,1'}
    header2 = {'dispname': '830G',
               'decker': 'foo',
               'binning': '1,1'}
    header3 = {'dispname': '830L',
               'decker': 'Z6CL01B',
               'binning': '1,1'}
    header4 = {'dispname': '830G',
               'decker': 'Z6CL01B',
               'binning': '1,1',
               'amp':  'foo'}

    assert config_key_match(spectrograph, header1, header2) is True
    assert config_key_match(spectrograph, header1, header3) is False
    assert config_key_match(spectrograph, header1, header4) is False
    assert config_key_match(spectrograph, header4, header1) is False

    """
    for source in source_list:
        print(f'Source: {source.coord}')
        for i in range(len(source.spec1d_file_list)):
            print(f'    {source.spec1d_file_list[i]} {source.spec_obj_list[i].NAME} {source.spec_obj_list[i].MASKDEF_OBJNAME}')
    """

def test_build_coadd_file_name():
    mock_sobjs = mock_specobjs('spec1d_file1')

    source = SourceObject(SkyCoord(ra=mock_sobjs.specobjs[0].RA,
                                   dec=mock_sobjs.specobjs[0].DEC, unit='deg'),
                          mock_sobjs.specobjs[0], mock_sobjs.header, 'spec1d_file1')
    assert build_coadd_file_name(source) == 'J132436.41+271928.56_DEIMOS_20200130.fits'

@cooked_required
def test_find_slits_to_exclude():
    cooked_sci_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science')
    spec2d_file_name =  os.path.join(cooked_sci_dir, 'spec2d_DE.20100913.22358-CFHQS1_DEIMOS_2010Sep13T061231.334.fits')

    par = pypeitpar.PypeItPar()
    par['collate1d'] = pypeitpar.Collate1DPar()
    par['collate1d']['slit_exclude_flags'] = 'BOXSLIT'
    exclude_map = find_slits_to_exclude([spec2d_file_name], par)
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
def test_adap_archive(tmp_path):
    cooked_sci_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science')

    spec1d_name = 'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560'
    spec2d_name = 'spec2d_b27-J1217p3905_KASTb_2015May20T045733.560'
    spec1d_files = [os.path.join(cooked_sci_dir, spec1d_name + '.fits')]
    spec2d_files = [os.path.join(cooked_sci_dir, spec2d_name + '.fits')]

    # Don't actually need a real coadd file, just any file to see if it gets copied to the archive directory    
    coadd_name = 'spec1d_b24-Feige66_KASTb_2015May20T041246.960'
    coadd_file = os.path.join(cooked_sci_dir, coadd_name + '.fits')

    # Build source for add_coadd_sources
    mock_sobjs = mock_specobjs('spec1d_file1')

    source_list = [SourceObject(SkyCoord(ra=mock_sobjs.specobjs[0].RA,
                                         dec=mock_sobjs.specobjs[0].DEC, unit='deg'),
                                mock_sobjs.specobjs[0], mock_sobjs.header, 'spec1d_file1')]
    source_list[0].coaddfile = coadd_file
    archive_dir = str(tmp_path)
    archive = ADAPArchive(archive_dir)
    archive.add_spec1d_files(spec1d_files)
    archive.add_spec2d_files(spec2d_files)
    archive.add_coadd_sources(source_list)
    spec_archive_dir = os.path.join(archive_dir, '2015', '05', '20')
    coadd_archive_dir = os.path.join(archive_dir, '2020', '01', '30')
    assert(os.path.exists(os.path.join(spec_archive_dir, spec1d_name + '.fits')))
    assert(os.path.exists(os.path.join(spec_archive_dir, spec2d_name + '.fits')))
    assert(os.path.exists(os.path.join(coadd_archive_dir, coadd_name + '.fits')))

    # Use filecmp to compare .dat files
    good_dat_dir = data_path('adap')
    spec1d_dat = os.path.join(good_dat_dir, spec1d_name + '.dat')
    spec2d_dat = os.path.join(good_dat_dir, spec2d_name + '.dat')
    coadd_dat = os.path.join(good_dat_dir, coadd_name + '.dat')
    assert(filecmp.cmp(spec1d_dat,os.path.join(spec_archive_dir, spec1d_name + '.dat'), shallow=False))
    assert(filecmp.cmp(spec2d_dat,os.path.join(spec_archive_dir, spec2d_name + '.dat'), shallow=False))
    assert(filecmp.cmp(coadd_dat,os.path.join(coadd_archive_dir, coadd_name + '.dat'), shallow=False))
