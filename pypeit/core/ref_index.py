"""
Module containing the core methods for calculating the refractive index
of the atmosphere based on current conditions.

Note that this code is directly copied (on 25 September 2023) from the following
repository: https://github.com/phn/ref_index

CREDIT :: All equations used in this module come from the documentation for the
NIST online refractive index calculator, written by Jack A. Stone and Jay H. Zimmerman,
and is available here:
https://emtoolbox.nist.gov/Wavelength/Documentation.asp

Credit to the original source
################
Refractive index of air.

NIST provides an online calculator for calculating refractive index of
air, for light of a certain wave length, under varying atmospheric
conditions. This module implements the equations provided in the
documentation for the online calculator.

In addition to calculating the refractive index, this module also has
functions for converting wave length of light in vacuum to that in air,
and vice-versa.

The documentation for the online calculator is provided at
http://emtoolbox.nist.gov/Wavelength/Documentation.asp, and includes a
link to the online calculator.

The following comments are based on the discussions presented in the
NIST documentation. It is intended as a brief overview. See
http://emtoolbox.nist.gov/Wavelength/Documentation.asp, for detailed
discussions.

Refractive index of air can be caclulated using two different
algorithms: one due to Edlén (updated by Birch and Down), and one due
to Ciddor. The latter has been adopted by the International Association
of Geodesy (IAG) as the reference equation for calculating refractive
index of air. Functions for calculating refractive index using either
of these are defined in this module.

The vacuum to air and air to vacuum wave length conversion functions in
this module use the Ciddor equation, in the form presented in the NIST
documentation.

Uncertainities in refractive index, and hence in wave length
conversions, due to uncertanities in measured values of temperature,
pressure, and humidity exceeds that due to the intrinsic uncertainity
in the equations used.

An uncertainty of 1e-6 in refractive index can result from a
combination of:

  + an error of 1°C (1.8 °F) in air temperature

  + an error of 0.4kPa (3mm of Hg) in air pressure

  + an error of 50% in relative humidity at sufficiently high air
    temperatures (near 35°C)

Valid range for input parameters for the refractive index calculations
are presented below. The online calculator issues a warning if input
parameters are outside a smaller interval within the maximum
range. Functions in this module do not raise a warning by default. But
they accept a keyword ``warn``, which when set to ``True`` will result
in warnings, when the input parameters are outside the accepted range.

  + Wavelength [300nm - 1700nm]

    Warning is issued if value is outside [350nm - 1600nm].

  + Pressure [10kPa - 140kPa]

    Warning is issued if value is outside [60kPa - 120kPa].

  + Temperature [-40∘C - 100∘C].

    Warning is issued if value is outside [0∘C - 40∘C].

  + Humidity [0 - 100]

    Can be given as relative humidity, dew point, frost point or
    partial pressure of water vapour. A warning is given if the mole
    fraction of water vapour exceeds 20% or, equivalently, relative
    humidity exceeds 85%. A warning is issued if relative humidity is
    less than 1%.

  + CO2 concentration [0µmole/mole - 2000µmole/mole]

    The common value to use is 450. Outdoor values are rarely below 300
    and indoor can be as high as 600. A difference of 150 will lead to
    a difference of only ~ 2e-8 in index of refraction.

    A warning is issued if a value other than 450 is used.


In astronomy, the convention is to use the refraction correction for
wave length greater than 200nm, eventhough the equations are not
strictly valid at wave lengths shorter than 300nm. For example, the
popular IDLASTRO IDL code vactoair.pro and airtovac.pro will accept any
wave length greater than 2000Å.

To accomodate this type of usage, instead of limiting the possible
input wave lengths, functions in this module will accept any wave
length value. It is up to the user to decide if a particular wave
length is to be used as an input to the equations.

Comparison with the IDLASTRO vactoair.pro and airtovac.pro algorithms
show that the equivalent functions in this module, vac2air and air2vac,
give results that agree to within 1e-4nm, over a range of wavelengths
from 200nm to 1700nm. This uncertainty translates to a velocity
difference of 150m/s to 17m/s, over the wave length range 1700nm to
200nm.

The IDLASTRO code uses a fixed value of temperature and humidity which
is not documented in the code. The above comparison was carried out at
a temperature of 15∘C and a relative humidity of 0.

The IDL code used for testing was downloaded on 2011/10/07. The
revision history indicates that the IDL code in vactoair.pro and
airtovac.pro were last modified in March 2011.

The PypeIt developers have made some minor adjustments to the code.

Original author details:
:author: Prasanth Nair
:contact: prasanthhn@gmail.com
:license: BSD (http://www.opensource.org/licenses/bsd-license.php)
###################

.. include:: ../include/links.rst

"""

from __future__ import division
from __future__ import print_function
import numpy as np
from pypeit import msgs


def f2k(f):
    """Converts Fahrenheit to Kelvin."""
    return (f - 32.0) * (100.0 / 180.0) + 273.15


def k2f(k):
    """Converts Kelvin to Fahrenheit."""
    return (k - 273.15) * (180.0 / 100.0) + 32.0


def c2k(c):
    """Converts Celsius to Kelvin."""
    return c + 273.15


def k2c(k):
    """Converts Kelvin to Celsius."""
    return k - 273.15


def c2f(c):
    """Converts Celsius to Fahrenheit."""
    return c * (180.0 / 100.0) - 32.0


def f2c(f):
    """Converts Fahrenheit to Celsius."""
    return (f - 32.0) * (100.0 / 180.0)


def svp_water(t):
    """Saturation vapour pressure over water at given temperature.

    Parameters
    ----------
    t : float
        Air temperature in degree Celsius.

    Returns
    -------
    p_sv : float
        Saturation vapour pressure over water, at the given
        temperature, in Pascal.

    Notes
    -----
    From section A-I of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    K1 = 1.16705214528e+03
    K2 = -7.24213167032e+05
    K3 = -1.70738469401e+01
    K4 = 1.20208247025e+04
    K5 = -3.23255503223e+06
    K6 = 1.49151086135e+01
    K7 = -4.82326573616e+03
    K8 = 4.05113405421e+05
    K9 = -2.38555575678e-01
    K10 = 6.50175348448e+02

    T = t + 273.15
    omega = T + K9 / (T - K10)
    A = omega ** 2 + K1 * omega + K2
    B = K3 * omega ** 2 + K4 * omega + K5
    C = K6 * omega ** 2 + K7 * omega + K8
    X = -B + np.sqrt(B ** 2 - 4 * A * C)

    p_sv = 1.0e6 * ((2.0 * C / X) ** 4)

    return p_sv


def svp_ice(t):
    """Saturation vapour pressure over ice at given temperature.


    Parameters
    ----------
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    p_sv : float
        Saturation vapour pressure over ice, at the given
        temperature, in Pascal.

    Notes
    -----
    From section A-I of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    A1 = -13.928169
    A2 = 34.7078238

    t += 273.15
    theta = t / 273.16
    Y = A1 * (1 - theta ** -1.5) + A2 * (1 - theta ** -1.25)

    p_sv = 611.657 * np.exp(Y)

    return p_sv


def dew_point_wvpp(td):
    """Water vapour saturation pressure, given dew point temperature."""
    return svp_water(td)


def frost_point_wvpp(tf):
    """Water vapour saturation pressure, given frost point temperature."""
    return svp_ice(tf)


def rh2wvpp(rh, t):
    """Convert relative humidity to water vapour partial pressure.

    Parameters
    ----------
    rh : float
        Relative humidity as a number between 0 and 100.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    p_sv : float
        Water vapour partial pressure, in Pascal.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    # t > 0 according to documentation.
    if t >= 0:
        p_sv = svp_water(t)
    elif t < 0:
        p_sv = svp_ice(t)

    return (rh / 100.0) * p_sv


def f_factor(p, t):
    """Enhancement factor for calculating mole fraction.

    Parameters
    ----------
    p : float
        Pressure in Pascal.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    f : float
        Enhancement factor needed in calculation of mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    alpha = 1.00062
    beta = 3.14e-8
    gamma = 5.60e-7

    return alpha + beta * p + gamma * (t ** 2)


def dew_point_mole_fraction(p, t):
    """Water vapour mole fraction for given dew point temperature.
    Parameters
    ----------
    p : float
        Pressure in Pascal.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * dew_point_wvpp(t) / p


def frost_point_mole_fraction(p, t):
    """Water vapour mole fraction for given frost point temperature.
    Parameters
    ----------
    p : float
        Pressure in Pascal.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * frost_point_wvpp(t) / p


def rh2mole_fraction(rh, p, t):
    """Water vapour mole fraction from relative humidity.

    Parameters
    ----------
    rh : float
        Relative humidity as a number between 0 and 100.
    p : float
        Pressure in Pascal.
    t : float
        Temperature in Kelvin.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * rh2wvpp(rh, t) / p


def pp2mole_fraction(pv, p, t):
    """Water vapour mole fraction from partial pressure.

    Parameters
    ----------
    rh : float
        Relative humidity as a number between 0 and 100.
    p : float
        Pressure in Pascal.
    t : float
        Temperature in Kelvin.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * pv / p


def _check_range(**kwargs):
    """Return True if value is inside accepted range."""
    if not (350 <= kwargs.get('wave', 633) <= 1600):
        msgs.warn("Wave length outside [350nm, 1600nm].")
    if not (60000 <= kwargs.get('p', 101325) <= 120000):
        msgs.warn("Pressure outside [60000Pa - 120000Pa].")
    if not (0 <= kwargs.get('t', 20) <= 40):
        msgs.warn("Temperature outside [0C - 40C].")
    if not (1 < kwargs.get('rh', 50) <= 85):
        msgs.warn("Relative humidity outside (1 - 85].")
    if not (kwargs.get('xv', 0.4) >= 0.2):
        msgs.warn("Mole fraction less than 0.2.")
    if kwargs.get('co2', 450) != 450:
        msgs.warn("CO2 concentration is not 450.")


def ciddor_ri(wave, t, p, xv, co2=450, warn=False):
    """Refractive index of air according to the Ciddor equation.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    xv : float
        Water vapour mole fraction, as a number between 0 and
        1. Default is set to 0.
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 should be enough for most purposes. Valid range is from
        0 - 2000 µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-III of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    See
    """
    if warn:
        _check_range(wave, t, p, xv)

    w0 = 295.235
    w1 = 2.6422
    w2 = -0.03238
    w3 = 0.004028
    k0 = 238.0185
    k1 = 5792105
    k2 = 57.362
    k3 = 167917
    a0 = 1.58123e-6
    a1 = -2.9331e-8
    a2 = 1.1043e-10
    b0 = 5.707e-6
    b1 = -2.051e-8
    c0 = 1.9898e-4
    c1 = -2.376e-6
    d = 1.83e-11
    e = -0.765e-8
    pr1 = 101325
    tr1 = 288.15
    Za = 0.9995922115
    rhovs = 0.00985938
    R = 8.314472
    Mv = 0.018015

    wave = wave * 1.0e-3
    S = 1.0 / wave ** 2

    ras = 1e-8 * ((k1 / (k0 - S)) + (k3 / (k2 - S)))
    rvs = 1.022e-8 * (w0 + w1 * S + w2 * S ** 2 + w3 * S ** 3)

    Ma = 0.0289635 + 1.2011e-8 * (co2 - 400.0)

    raxs = ras * (1 + 5.34e-7 * (co2 - 450.0))

    T = t + 273.15

    Zm = a0 + a1 * t + a2 * t ** 2 + (b0 + b1 * t) * xv + \
        (c0 + c1 * t) * xv ** 2
    Zm *= -(p / T)
    Zm += (p / T ) ** 2 * (d + e * xv ** 2)
    Zm += 1

    rhoaxs = pr1 * Ma / (Za * R * tr1)

    rhov = xv * p * Mv / (Zm * R * T)

    rhoa = (1 - xv) * p * Ma / (Zm * R * T)

    n = 1.0 + (rhoa / rhoaxs) * raxs + (rhov / rhovs) * rvs

    return n


def ciddor(wave, t, p, rh, co2=450, warn=False):
    """Refractive index of air according to the Ciddor equation.

    Accepts relative humidity instead of mole fraction, as done in
    ``ciddor_ri()``.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    rh : float
        Relative humidity [0 - 100].
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 should be enough for most purposes. Valid range is from
        0 - 2000 µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-III of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    if warn:
        _check_range(wave, t, p, rh)
        # turn off warning, so that ciddor_ri doesn't issue duplicate
        # warning.
        warn = False

    xv = rh2mole_fraction(rh=rh, p=p, t=t)
    return ciddor_ri(wave=wave, t=t, p=p, xv=xv, co2=co2, warn=warn)


def edlen_ri(wave, t, p, pv, warn=False):
    """Refractive index of air according to the Edlén equation.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    pv : float
        Water vapour partial pressure, in Pascal.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-IV of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    if warn:
        _check_range(wave, t, p)

    A = 8342.54
    B = 2406147
    C = 15998
    D = 96095.43
    E = 0.601
    F = 0.00972
    G = 0.003661

    wave = wave * 1.0e-3
    S = 1.0 / wave ** 2

    ns = 1 + 1e-8 * (A + B / (130.0 - S) + C / (38.9 - S))

    X = (1 + 1e-8 * (E - F * t) * p) / (1 + G * t)

    ntp = 1 + p * (ns - 1) * X / D

    n = ntp - 1e-10 * ((292.75 / (t + 273.15)) * \
                           (3.7345 - 0.0401 * S)) * pv

    return n


def edlen(wave, t, p, rh, warn=False):
    """Refractive index of air according to the Edlén equation.

    Accepts relative humidity instead of water vapour partial pressure,
    as in ``edlen_ri()``.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    rh : float
        Relative humidity in [0 - 100].
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-IV of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    if warn:
        _check_range(wave, t, p)
        # turn off warning so that edlen_ri() doesn't raise duplicate
        # warning.
        warn = False

    pv = rh2wvpp(rh=rh, t=t)
    return edlen_ri(wave=wave, t=t, p=p, pv=pv, warn=warn)


def vac2air(wave, t=15.0, p=101325, rh=0.0, co2=450, warn=False):
    """Wavelength of light in air, using Ciddor refractive index.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in nano-meters. Valid range is 300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid range is -40 - 100 degree
        Celsius. Default is 15 degree Celsius (288.15 Kelvin).
    p : float
        Pressure in Pascal. Valid range is 10kPa - 140kPa. Default is
        101325 Pa (1 atmosphere).
    rh : float
        Relative humidity as a number between 0 and 100. Default is 0.
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 is sufficient for most purposes. Valid range is 0 - 2000
        µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Returns
    -------
    w : float
        Wavelength in air, in nm.

    """
    if warn:
        _check_range(wave, t, p, rh, co2)

    n = ciddor(wave, t, p, rh, co2)
    return wave / n


def air2vac(wave, t=15.0, p=101325, rh=0.0, co2=450, warn=False):
    """Wavelength of light in vacuum, using Ciddor refractive index.

    The refractive index calculation needs wavelength in vacuum. In
    this function, the wavelength in air is used. The errors are on the
    order of 1e-5 nm.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in nano-meters. Valid range is 300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid range is -40 - 100 degree
        Celsius. Default is 15 degree Celsius (288.15 Kelvin).
    p : float
        Pressure in Pascal. Valid range is 10kPa - 140kPa. Default is
        101325 Pa (1 atmosphere).
    rh : float
        Relative humidity as a number between 0 and 100. Default is 0.
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 is sufficient for most purposes. Valid range is 0 - 2000
        µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Returns
    -------
    w : float
        Wavelength in vacuum, in nm.

    """
    if warn:
        _check_range(wave=wave, t=t, p=p, rh=rh, co2=co2)

    n = ciddor(wave, t, p, rh, co2)
    return wave * n
