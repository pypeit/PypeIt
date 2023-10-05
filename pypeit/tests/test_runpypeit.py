"""
Module to do a full run of PypeIt on good ole Kastb
"""
from pathlib import Path
import shutil
from IPython.terminal.embed import embed
import numpy as np

import matplotlib
matplotlib.use('agg')  

from pypeit.scripts.setup import Setup
from pypeit.scripts.run_pypeit import RunPypeIt
from pypeit.scripts.sensfunc import SensFunc
from pypeit.scripts.flux_calib import FluxCalib
from pypeit.tests.tstutils import data_path
from pypeit import specobjs, sensfunc
from pypeit.par import pypeitpar 


def test_run_pypeit():

    # Just get a few files
    testrawdir = Path(data_path('')).resolve()
    outdir = testrawdir / 'REDUX_OUT_TEST'

    # For previously failed tests
    if outdir.exists():
        shutil.rmtree(outdir)

    # Run the setup
    sargs = Setup.parse_args(['-r', str(testrawdir / 'b'), '-s', 
                              'shane_kast_blue', '-c all', '-o', 
                              '--output_path', str(outdir)])
    Setup.main(sargs)

    # Change to the configuration directory and set the pypeit file
    configdir = outdir / 'shane_kast_blue_A'
    pyp_file = configdir / 'shane_kast_blue_A.pypeit'
    assert pyp_file.exists(), 'PypeIt file not written.'

    # Try to run with -m and -o
    pargs = RunPypeIt.parse_args([str(pyp_file), '-o', '-m', '-r', str(configdir)])
    RunPypeIt.main(pargs)

    # TODO: Add some code that will try to open the QA HTML and check that it
    # has the correct PNGs in it.

    # #########################################################
    # Test!!
    # Files exist
    spec1d_file = configdir / 'Science' / 'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    assert spec1d_file.exists(), 'spec1d file missing'

    # .par file was written and loads
    par_file = str(list(configdir.glob('shane_kast_blue*.par'))[0])
    par = pypeitpar.PypeItPar.from_cfg_file(par_file)
    assert isinstance(par, pypeitpar.PypeItPar)
                               
    # Now re-use those calibration files
    pargs = RunPypeIt.parse_args([str(pyp_file), '-o', '-r', str(configdir)])
    RunPypeIt.main(pargs)

    # Generate a sensitivity function from the standard star spec1d
    std_file = configdir / 'Science' / 'spec1d_b24-Feige66_KASTb_20150520T041246.960.fits'
    assert std_file.exists(), 'std spec1d file missing'

    sens_file = outdir / "sens_shane_kast_blue.fits"

    pargs = SensFunc.parse_args(["--algorithm", "UVIS", "-o", str(sens_file), str(std_file)])
    SensFunc.main(pargs)

    # Assert sensfunc file exists and contains data
    assert sens_file.exists(), "sensfunc file missing"
    sf = sensfunc.SensFunc.from_file(str(sens_file))
    assert len(sf.wave) > 0
    assert len(sf.zeropoint) > 0

    # Flux calibrate
    with open(outdir / "fluxfile", "w") as f:
        print("flux read", file=f)
        print("filename | sensfile", file=f)
        print(f"{spec1d_file} | {sens_file}", file=f)
        print("flux end", file=f)
    pargs = FluxCalib.parse_args([str(outdir / "fluxfile")])
    FluxCalib.main(pargs)

    # Verify spec 1d
    # spec1d
    specObjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
    
    # Check RMS
    assert specObjs[0].WAVE_RMS < 0.08  # RMS must be less than 0.08 pixels

    # Flexure
    assert abs(-0.03 - specObjs[0].FLEX_SHIFT_TOTAL) < 0.1  # difference must be less than 0.1 pixels

    # Helio
    assert abs(specObjs[0].VEL_CORR - 0.9999251762866389) < 1.0E-10

    # Flux
    mask = specObjs[0].OPT_MASK
    assert len(np.nonzero(mask)) != 0

    assert len(specObjs[0].OPT_FLAM[mask]) != 0
    assert len(specObjs[0].OPT_FLAM_IVAR[mask]) != 0
    assert len(specObjs[0].OPT_FLAM_SIG[mask]) != 0

    mask = specObjs[0].BOX_MASK
    assert len(np.nonzero(mask)) != 0

    assert len(specObjs[0].BOX_FLAM[mask]) != 0
    assert len(specObjs[0].BOX_FLAM_IVAR[mask]) != 0
    assert len(specObjs[0].BOX_FLAM_SIG[mask]) != 0

    # Clean-up
    shutil.rmtree(outdir)


