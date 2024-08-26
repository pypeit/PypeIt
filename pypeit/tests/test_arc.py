"""
Module to run tests on ararclines
"""
import pytest
import numpy as np

from pypeit.core import arc
from pypeit import io


def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    arx_sky = io.load_sky_spectrum('paranal_sky.fits')
    arx_amp_true, arx_amp, arx_cent, arx_wid, arx_centerr, arx_w, arx_yprep, _ \
            = arc.detect_lines(arx_sky.flux.value)
    assert (len(arx_w) > 3275)


def test_mask_around_peaks():
    # Generate a fake spectrum
    np.random.seed(0)
    npks = 5   # Number of peaks to generate
    npix = 20  # Number of pixels in each segment
    sigma = 2.0  # Fluctuations
    linewid = 2.0  # Width of the peaks
    noiselev = [np.random.normal(0.0, sigma, npix) for i in range(npks+1)]
    ampl = np.random.randint(100*sigma, 200*sigma, npks)
    spec = np.array([])
    outbpm = np.array([], dtype=bool)
    for i in range(npks):
        spec = np.append(spec, noiselev[i])
        outbpm = np.append(outbpm, np.zeros(npix, dtype=bool))
        gauss = ampl[i]*np.exp(-0.5*(np.arange(npix)-npix/2)**2/linewid**2)
        spec = np.append(spec, np.random.normal(gauss, sigma))
        if i == 0 or i == 1:
            tmpbpm = np.array(2*[False] + 16*[True] + 2*[False])
        elif i == 2:
            tmpbpm = np.array(1*[False] + 18*[True] + 1*[False])
        elif i == 3:
            tmpbpm = np.array(3*[False] + 15*[True] + 2*[False])
        elif i == 4:
            tmpbpm = np.array(20*[True])
        else:
            assert (False, "Input test data have changed")
        outbpm = np.append(outbpm, tmpbpm.copy())
    outbpm = np.append(outbpm, np.zeros(npix, dtype=bool))
    spec = np.append(spec, noiselev[npks])
    inbpm = np.zeros(len(spec), dtype=bool)
    # This is the expected bpm
    bpm = arc.mask_around_peaks(spec, inbpm=inbpm)
    # Check that all peaks were identified and masked
    ind = arc.detect_peaks(bpm)
    assert (len(ind) == npks)
    # Check the BPM matches exactly
    assert (np.array_equal(bpm, outbpm))


# Many more functions in pypeit.core.arc that need tests!

