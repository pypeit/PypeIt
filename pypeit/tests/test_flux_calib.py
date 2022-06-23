"""
Test pypeit_flux_calib file parsing and command line argument parsing logic.
"""
import pytest

import os

from astropy.io import fits

from pypeit import fluxcalibrate
from pypeit.pypmsgs import PypeItError
from pypeit import scripts


def test_flux_calib(tmp_path, monkeypatch):

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
            print(" spec1d_file1.fits sens_file1.fits", file=f)
            print(" spec1d_file2.fits sens_file2.fits", file=f)


        with pytest.raises(PypeItError, match="Missing 'flux end'"):
            parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_missing_end])
            scripts.flux_calib.FluxCalib.main(parsed_args)

        # Test with a flux file missing the flux block entirely
        config_file_missing_flux = str(tmp_path / "test_flux_calib_missing_flux.flux")
        with open(config_file_missing_flux, "w") as f:
            print(" spec1d_file1.fits sens_file1.fits", file=f)
            print(" spec1d_file2.fits sens_file2.fits", file=f)
        
        with pytest.raises(PypeItError, match="Missing flux block in"):
            parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_missing_flux])
            scripts.flux_calib.FluxCalib.main(parsed_args)

        # Test 1 sens file with multiple spec1ds
        config_file_one_to_many = str(tmp_path / "test_flux_calib_1_to_many.flux")
        with open(config_file_one_to_many, "w") as f:
            print("flux read", file=f)
            print(" spec1d_file1.fits sens_file1.fits", file=f)
            print(" spec1d_file2.fits", file=f)
            print(" spec1d_file3.fits", file=f)
            print("flux end", file=f)

        parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_one_to_many])
        assert scripts.flux_calib.FluxCalib.main(parsed_args) == 0

        # Test 1 sens file per spec1d
        config_file_one_to_one = str(tmp_path / "test_flux_calib_one_to_one.flux")
        with open(config_file_one_to_one, "w") as f:
            print("flux read", file=f)
            print(" spec1d_file1.fits sens_file1.fits", file=f)
            print(" spec1d_file2.fits sens_file2.fits", file=f)
            print(" spec1d_file3.fits sens_file1.fits", file=f)
            print("flux end", file=f)

        parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_one_to_one])
        assert scripts.flux_calib.FluxCalib.main(parsed_args) == 0
        
        # Test with no sensfunc, but using an archived sensfunc
        config_file_use_arxiv = str(tmp_path / "test_flux_calib_use_arxiv.flux")
        with open(config_file_use_arxiv, "w") as f:
            print("[fluxcalib]", file=f)
            print(" use_archived_sens = True", file=f)
            print("flux read", file=f)
            print(" spec1d_file1.fits", file=f)
            print(" spec1d_file2.fits", file=f)
            print(" spec1d_file3.fits", file=f)
            print("flux end", file=f)

        parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_use_arxiv])
        assert scripts.flux_calib.FluxCalib.main(parsed_args) == 0
        
        
        # Test with no sensfunc, but it's an error because an archive sensfunc
        # was not requested
        config_file_no_sens = str(tmp_path / "test_flux_calib_no_sens.flux")
        with open(config_file_no_sens, "w") as f:
            print("flux read", file=f)
            print(" spec1d_file1.fits", file=f)
            print(" spec1d_file2.fits", file=f)
            print(" spec1d_file3.fits", file=f)
            print("flux end", file=f)

        with pytest.raises(PypeItError, match = 'Invalid format for .flux'):
            parsed_args = scripts.flux_calib.FluxCalib.parse_args([config_file_no_sens])
            scripts.flux_calib.FluxCalib.main(parsed_args)
        
