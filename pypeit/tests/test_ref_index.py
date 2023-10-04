"""
Module to run tests on pyidl functions
"""
import numpy as np

from pypeit.core import ref_index
import pytest


def test_nist_ciddor_1():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix t at 20, p at 101325, rh at 50
    """
    wave = [321.456, 500, 600.1234, 633.0, 700, 1000.987, 1500.8, 1700.0]
    nist_n = [1.000283543, 1.000273781, 1.000271818, 1.000271373,
              1.000270657, 1.000269038, 1.00026819, 1.000268041]
    nist_w = [321.364879, 499.863147, 599.96032, 632.828268,
              699.810591, 1000.717769, 1500.397608, 1699.544453]

    xv = ref_index.rh2mole_fraction(50, 101325, 20)

    n = [ref_index.ciddor_ri(i, 20, 101325, xv) for i in wave]
    wave_n = [ref_index.vac2air(i, t=20, p=101325, rh=50.0) for i in wave]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6

    n = [ref_index.ciddor(i, 20, 101325, 50.0) for i in wave]
    wave_n = [ref_index.vac2air(i, t=20, p=101325, rh=50.0) for i in wave]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def test_nist_ciddor_2():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix wave at 633.0 p at 101325 rh at 50
    """
    t = [-20.0, 0.0, 20, 26.7982, 40.123, 60.45]
    nist_w = [632.800737, 632.815441, 632.828268, 632.832303, 632.839872,
              632.850953]
    nist_n = [1.00031489, 1.000291647, 1.000271373, 1.000264994, 1.000253031,
              1.000235516]

    xv = [ref_index.rh2mole_fraction(50, 101325, i) for i in t]
    n = [ref_index.ciddor_ri(633.0, i, 101325, j) for i, j in zip(t, xv)]

    wave_n = [ref_index.vac2air(633.0, i, 101325, 50) for i in t]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def test_nist_ciddor_3():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix wave at 633.0, t at 20, rh at 50. vary p
    """
    p = [1000 * i for i in [10, 50.123, 100.1234, 140.0]]

    nist_n = [1.000026385, 1.000133999, 1.000268148, 1.000375169]
    nist_w = [632.983299, 632.91519, 632.830308, 632.762607]

    xv = [ref_index.rh2mole_fraction(50, i, 20) for i in p]
    n = [ref_index.ciddor_ri(633.0, 20, i, j) for i, j in zip(p, xv)]

    wave_n = [ref_index.vac2air(633.0, 20, i, 50) for i in p]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def test_nist_ciddor_4():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix wave at 633.0, t at 20, p at 101325, vary rh.
    """
    rh = [0.0, 20.123, 40, 50.9876, 70, 90.7432, 100.0]
    nist_n = [1.0002718, 1.000271627, 1.000271458, 1.000271364,
              1.000271203, 1.000271027, 1.000270949]
    nist_w = [632.827997, 632.828106, 632.828214, 632.828273,
              632.828375, 632.828486, 632.828535]

    xv = [ref_index.rh2mole_fraction(i, 101325, 20) for i in rh]
    n = [ref_index.ciddor_ri(633.0, 20, 101325, j) for j in xv]

    wave_n = [ref_index.vac2air(633.0, 20, 101325, i) for i in rh]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def test_air2vac():
    """Test reversibility with vac2air."""
    wave = np.array([321.456, 500, 600.1234, 633.0, 700, 1000.987, 1500.8, 1700.0])
    wave_o = ref_index.air2vac(ref_index.vac2air(wave))
    assert np.allclose(wave, wave_o)


def test_idlastro():
    # Using IDLASTRO downloaded on 2011/10/07. The vac2air.pro uses a
    # formulation of the Ciddor equation. Previous versions used a
    # different equation.

    # The REVISION HISTORY from the vac2air.pro file is:
    # ; REVISION HISTORY
    # ;       Written W. Landsman                November 1991
    # ;       Use Ciddor (1996) formula for better accuracy in the infrared
    # ;           Added optional output vector, W Landsman Mar 2011
    # ;       Iterate for better precision W.L./D. Schlegel  Mar 2011

    # The REVISION HISTORY from air2vac.pro file is:
    # ; REVISION HISTORY
    # ;	Written, D. Lindler 1982
    # ;	Documentation W. Landsman  Feb. 1989
    # ;       Use Ciddor (1996) formula for better accuracy in the infrared
    # ;           Added optional output vector, W Landsman Mar 2011

    # Velocity errors in m/s for different wave length errors, at
    # different wave lengths.
    # >>> 1e-5/330.0 * 299792458
    # 9.0846199393939404
    # >>> 1e-5/200.0 * 299792458
    # 14.989622900000001
    # >>> 1e-5/1000.0 * 299792458
    # 2.9979245800000003

    # nm
    wave = np.array([200.0, 300.0, 500.0, 800.0, 1200.0, 1600.0, 1700.0])

    # angstrom
    wave_idl_vactoair = np.array([1999.3526550081103323, 2999.1255923046301177,
                                  4998.6055889614663101, 7997.8003315140686027,
                                  11996.7167708424640296, 15995.6298776736693981,
                                  16995.3579139663052047])
    wave_vac2air = ref_index.vac2air(wave, t=15, rh=0)

    # values in wave_idl_vactoair was fed to airtovac idl procedure.
    wave_idl_airtovac = np.array([1999.3526550081103323,
                                  3000.0000371189012185,
                                  5000.0000183785432455,
                                  8000.0000108292333607,
                                  12000.0000070745754783,
                                  16000.0000052688483265,
                                  17000.0000049538284657])
    # Have to convert angstrom to nm.
    wave_air2vac = ref_index.air2vac(wave_idl_vactoair / 10.0, t=15, rh=0)

    assert np.allclose(wave_vac2air, wave_idl_vactoair/10.0)

    # IDL code ignores values under 2000 angstrom.
    assert np.allclose(wave_air2vac[1:], wave_idl_airtovac[1:]/10.0)
