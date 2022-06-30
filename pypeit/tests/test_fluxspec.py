"""
Module to run tests on SensFunc and FluxCalibrate classes
Requires files in PypeIt/pypeit/data
"""
import os

import pytest

from IPython import embed

import numpy as np
import configobj

from astropy.table import Table
from astropy.io import fits

from pypeit import fluxcalibrate
from pypeit import sensfunc
from pypeit.par import pypeitpar
from pypeit.tests.tstutils import data_path
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs import keck_deimos
from pypeit import specobjs, specobj
from pypeit.tests import tstutils
from pypeit import inputfiles 

from pypeit import fluxcalibrate
from pypeit.pypmsgs import PypeItError
from pypeit import scripts

def test_input_flux_file():
    """Tests for generating and reading fluxing input files
    """
    # Generate an input file
    flux_input_file = tstutils.data_path('test.flux')
    if os.path.isfile(flux_input_file):
        os.remove(flux_input_file)

    cfg_lines = ['[fluxcalib]']
    cfg_lines += ['  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm\n']
    cfg_lines += ['# Please add your SENSFUNC file name below before running pypeit_flux_calib']

    # These files need to be in tests/files/
    data = Table()
    data['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data['sensfile'] = 'sens_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    # 
    paths = [tstutils.data_path('')]

    fluxFile = inputfiles.FluxFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data)
    # Write
    fluxFile.write(flux_input_file)

    # Read
    fluxFile2 = inputfiles.FluxFile.from_file(flux_input_file)
    assert np.all(fluxFile2.data['filename'] == data['filename'])

    # Test path
    assert fluxFile2.file_paths[0] == paths[0]
    assert fluxFile2.filenames[0] == os.path.join(paths[0], data['filename'][0])

    # #################
    # Tickle the other ways to do sensfiles
    data3 = Table()
    data3['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data3['sensfile'] = ['sens_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits',
                         '']

    fluxFile3 = inputfiles.FluxFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data3)
    assert fluxFile3.sensfiles[1] == os.path.join(paths[0], data['sensfile'][0])
    
    data4 = Table()
    data4['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data4['sensfile'] = ''

    fluxFile4 = inputfiles.FluxFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data4)
    assert len(fluxFile4.sensfiles) == 0

    # Clean up
    os.remove(flux_input_file)

def test_flux_calib(tmp_path, monkeypatch):
    """ Some of these items are also tested in test_fluxspec.py"""

    # Change to the tmp_path so the fluxing.par file is written there
    os.chdir(tmp_path)

    # Test the flux_calib script (but not fluxing itself)
    def mock_get_header(*args, **kwargs):
        return {"DISPNAME": "600ZD",
                "PYP_SPEC": "keck_deimos" }

    def mock_get_flux_calib_instance(*args, **kwargs):
        # The flux_calib caller doesn't use the output, it just
        # depends on the side effect of fluxing
        return None 


    with monkeypatch.context() as m:
        monkeypatch.setattr(fits, "getheader", mock_get_header)
        monkeypatch.setattr(fluxcalibrate.FluxCalibrate, "get_instance", mock_get_flux_calib_instance)

        # Test with a flux file missing "flux end"

        config_file_missing_end = str(tmp_path / "test_flux_calib_missing_end.flux")

        with open(config_file_missing_end, "w") as f:
            print("flux read", file=f)
            print("filename | sensfile", file=f)
            print("spec1d_file1.fits | sens_file1.fits", file=f)
            print("spec1d_file2.fits | sens_file2.fits", file=f)


        with pytest.raises(PypeItError, match="Missing 'flux end'"):
            parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_missing_end])
            scripts.flux_calib.FluxCalib.main(parsed_args)

        # Test with a flux file missing the flux block entirely
        config_file_missing_flux = str(tmp_path / "test_flux_calib_missing_flux.flux")
        with open(config_file_missing_flux, "w") as f:
            print("filename | sensfile", file=f)
            print("spec1d_file1.fits | sens_file1.fits", file=f)
            print("spec1d_file2.fits | sens_file2.fits", file=f)
        
        with pytest.raises(PypeItError, match="You have not specified the data block!"):
            parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_missing_flux])
            scripts.flux_calib.FluxCalib.main(parsed_args)

        # Test with no sensfunc, but it's an error because an archive sensfunc
        # was not requested
        config_file_no_sens = str(tmp_path / "test_flux_calib_no_sens.flux")
        with open(config_file_no_sens, "w") as f:
            print("flux read", file=f)
            print(f"path {data_path('')}", file=f)
            print("filename", file=f)
            print("spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits", file=f)
            print("spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits", file=f)
            print("flux end", file=f)

        with pytest.raises(PypeItError, match = 'Invalid format for .flux'):
            parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_no_sens])
            scripts.flux_calib.FluxCalib.main(parsed_args)
        

# TODO: Include tests for coadd2d, sensfunc


def test_extinction_correction_uvis():
    extinction_correction_tester('UVIS')

def test_extinction_correction_ir():
    extinction_correction_tester('IR')

def extinction_correction_tester(algorithm):
    spec1d_file = data_path('spec1d_test.fits')
    sens_file = data_path('sens_test.fits')

    if os.path.isfile(spec1d_file):
        os.remove(spec1d_file)
    if os.path.isfile(sens_file):
        os.remove(sens_file)

    # make a bogus spectrum that has N_lam = 1
    wave = np.linspace(4000, 6000)
    counts = np.ones_like(wave)
    ivar = np.ones_like(wave)
    sobj = specobj.SpecObj.from_arrays('MultiSlit', wave, counts, ivar)
    sobjs = specobjs.SpecObjs([sobj])

    # choice of PYP_SPEC, DISPNAME and EXPTIME are unimportant here
    # AIRMASS must be > 1
    sobjs.write_to_fits({
        'PYP_SPEC': 'p200_dbsp_blue',
        'DISPNAME': '600/4000',
        'EXPTIME': 1.0,
        'AIRMASS': 1.1}, spec1d_file)

    par = pypeitpar.PypeItPar()

    # set the senfunc algorithm
    par['sensfunc']['algorithm'] = algorithm
    # needed to initiate SensFunc (dummy standard star Feige34)
    par['sensfunc']['star_ra'] = 159.9042
    par['sensfunc']['star_dec'] = 43.1025

    sensobj = sensfunc.SensFunc.get_instance(spec1d_file, sens_file, par['sensfunc'])

    sensobj.wave = np.linspace(3000, 6000, 300).reshape((300, 1))
    sensobj.sens = sensobj.empty_sensfunc_table(*sensobj.wave.T.shape)
    # make the zeropoint such that the sensfunc is flat
    sensobj.zeropoint = 30 - np.log10(sensobj.wave ** 2) / 0.4

    sensobj.to_file(sens_file)

    # now flux our N_lam = 1 specobj
    par['fluxcalib']['extinct_correct'] = None
    fluxCalibrate = fluxcalibrate.MultiSlitFC([spec1d_file], [sens_file], par=par['fluxcalib'])
    # without extinction correction, we should get constant F_lam
    # with extinction correction, the spectrum will be blue

    # make sure that the appropriate default behavior occurred
    sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
    print(sobjs[0].keys())
    if algorithm == 'UVIS':
        assert sobjs[0]['OPT_FLAM'][0] > sobjs[0]['OPT_FLAM'][-1], \
            "UVIS sensfunc was not extinction corrected by default, but should have been"
    elif algorithm == 'IR':
        assert np.isclose(sobjs[0]['OPT_FLAM'][0], sobjs[0]['OPT_FLAM'][-1]), \
            "IR sensfunc was extinction corrected by default, but shouldn't have been"

    # clean up
    os.remove(spec1d_file)
    os.remove(sens_file)


