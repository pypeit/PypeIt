"""
Module to do a full run of PypeIt on good ole Kastb
"""
import os
from pathlib import Path
import shutil
import sys
from IPython.terminal.embed import embed
import numpy as np
import pytest

import matplotlib
matplotlib.use('agg')  

from pypeit import dataPaths
from pypeit.scripts.setup import Setup
from pypeit.scripts.run_pypeit import RunPypeIt
from pypeit.scripts.sensfunc import SensFunc
from pypeit.scripts.flux_calib import FluxCalib
from pypeit import specobjs, sensfunc
from pypeit import log
from pypeit.par import pypeitpar 
from pypeit.tests import tstutils

@pytest.mark.remote_data
def test_run_pypeit():

    # Just get a few files
    testrawdir = dataPaths.tests.path
    tstutils.install_shane_kast_blue_raw_data()
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


@pytest.mark.remote_data
def test_run_pypeit_ncpu():
    """
    Reduce shane_kast_blue twice, serially and with --ncpu 2, into separate
    output directories and require the products to match.

    shane_kast_blue has a single detector, so :func:`~pypeit.parallel.nworkers`
    collapses to one worker; :func:`~pypeit.parallel.map_over_detectors` still
    enters the process pool because ``ncpu>1``, which is the point.  Starting
    the worker, the BLAS thread pinning, the buffered worker logging and the
    result marshalling all get exercised by this reduction.

    On Linux the products must be bit-identical.  Elsewhere they are only
    required to agree to within a tolerance: on macOS, threadpoolctl cannot
    control the Accelerate BLAS, so the serial run is not pinned, and long-slit
    reductions there have been seen to differ even between two serial runs.
    """
    if os.environ.get('PYPEIT_SKIP_SLOW_TESTS'):
        pytest.skip('slow test')

    testrawdir = dataPaths.tests.path
    tstutils.install_shane_kast_blue_raw_data()
    outdir = testrawdir / 'REDUX_OUT_NCPU_TEST'

    # For previously failed tests
    if outdir.exists():
        shutil.rmtree(outdir)

    sargs = Setup.parse_args(['-r', str(testrawdir / 'b'), '-s', 'shane_kast_blue',
                              '-c all', '-o', '--output_path', str(outdir)])
    Setup.main(sargs)

    configdir = outdir / 'shane_kast_blue_A'
    pyp_file = configdir / 'shane_kast_blue_A.pypeit'
    assert pyp_file.exists(), 'PypeIt file not written.'

    # Reduce into two separate directories, serially and in parallel.  Both read
    # the same .pypeit file; only the redux path and --ncpu differ.
    serial_dir = outdir / 'serial'
    par_dir = outdir / 'parallel'
    for redux_dir, ncpu in ((serial_dir, '1'), (par_dir, '2')):
        redux_dir.mkdir(parents=True, exist_ok=True)
        RunPypeIt.main(RunPypeIt.parse_args([str(pyp_file), '-o', '-r', str(redux_dir),
                                             '--ncpu', ncpu]))

    basename = 'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits'
    serial_spec1d = serial_dir / 'Science' / basename
    par_spec1d = par_dir / 'Science' / basename
    assert serial_spec1d.exists(), 'serial spec1d file missing'
    assert par_spec1d.exists(), 'parallel spec1d file missing'

    s = specobjs.SpecObjs.from_fitsfile(str(serial_spec1d))
    p = specobjs.SpecObjs.from_fitsfile(str(par_spec1d))
    assert s.nobj == p.nobj, 'serial and parallel runs found a different number of objects'

    exact = sys.platform.startswith('linux')
    for key in ('OPT_COUNTS', 'OPT_COUNTS_IVAR', 'OPT_WAVE', 'BOX_COUNTS', 'TRACE_SPAT'):
        sval, pval = getattr(s[0], key), getattr(p[0], key)
        msg = f'{key} differs between the serial and the --ncpu 2 reduction'
        if exact:
            # Exact comparison, not allclose: any difference means the BLAS
            # thread pinning is not doing its job, since the sigma-clipped
            # iterative fits amplify a last-ulp change into a visibly different
            # sky model.
            np.testing.assert_array_equal(sval, pval, err_msg=msg)
        else:
            # Scientific equivalence: differences relative to the array's
            # maximum, judged at the 99.9th percentile so that a handful of
            # pixels at a mask boundary do not dominate.
            scale = max(float(np.max(np.absolute(sval))), np.finfo(float).tiny)
            rel = np.absolute(np.asarray(sval, dtype=float) - pval) / scale
            assert np.quantile(rel, 0.999) < 1e-3, msg
    if exact:
        assert s[0].WAVE_RMS == p[0].WAVE_RMS, \
                'wavelength solution RMS differs between the serial and the --ncpu 2 reduction'
    else:
        assert np.isclose(s[0].WAVE_RMS, p[0].WAVE_RMS, rtol=1e-3), \
                'wavelength solution RMS differs between the serial and the --ncpu 2 reduction'

    # Clean-up.  The reduction leaves its log file open, and Windows refuses to
    # remove a file that is still held, so release it before deleting the tree.
    log.close_file()
    shutil.rmtree(outdir)
