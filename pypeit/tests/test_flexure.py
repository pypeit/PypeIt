"""
Module to run tests on simple fitting routines for arrays
"""
import os

import pytest

import numpy as np

from linetools.spectra.io import readspec

import pypeit
from pypeit.core import flexure, arc

from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests.tstutils import cooked_required

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

from IPython import embed

#dev_path = os.getenv('PYPEIT_DEV')
#lris_path =  os.path.join(dev_path, 'REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560')
#trc_file = os.path.join(lris_path, 'Masters', 'MasterTrace_A_1_01.fits')
#tilts_file = os.path.join(lris_path, 'Masters', 'MasterTilts_A_1_01.fits')
#spec2dfile = os.path.join(lris_path, 'Science', 'spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits')

#def test_flex_lris():
#    dpath = '/scratch/REDUX/Keck/LRIS/2020jan28/keck_lris_red_C'
#    spec2dfile = os.path.join(dpath, 'Science', 'spec2d_LR.20200128.57449-frb1907B_LRISr_2020Jan28T155729.779.fits')
#    tilts_file = os.path.join(dpath, 'Masters', 'MasterTilts_C_1_01.fits')
#    slits_file = os.path.join(dpath, 'Masters', 'MasterSlits_C_1_01.fits.gz')
#    # Load
#    hdul = fits.open(spec2dfile)
#    sciimg = hdul[1].data
#    slits = slittrace.SlitTraceSet.from_file(slits_file)
#    tilts = wavetilts.WaveTilts.from_file(tilts_file)
#
#    flexure.spat_flexure_shift(sciimg, slits)


def test_flex_shift():
    # Dummy slf
    # Read spectra
    obj_spec = readspec(data_path('obj_lrisb_600_sky.fits'))
    arx_file = pypeit.__path__[0]+'/data/sky_spec/sky_LRISb_600.fits'
    arx_spec = readspec(arx_file)
    arx_lines = arc.detect_lines(arx_spec.flux.value)

    # Call
    flex_dict = flexure.spec_flex_shift(obj_spec, arx_spec, arx_lines, mxshft=60)

#    # Apply
#    from scipy import interpolate
#    print(flex_dict['shift'])
#    npix = len(obj_spec.wavelength)
#    x = np.linspace(0., 1., npix)
#    f = interpolate.interp1d(x, obj_spec.wavelength.value, bounds_error=False,
#                             fill_value="extrapolate")
#    new_wave = f(x+flex_dict['shift']/(npix-1))
#
#    from matplotlib import pyplot
#    pyplot.plot(arx_spec.wavelength, arx_spec.flux)
#    pyplot.plot(obj_spec.wavelength, obj_spec.flux)
#    pyplot.plot(new_wave, obj_spec.flux)
#    pyplot.show()
    assert np.abs(flex_dict['shift'] - 43.7) < 0.1

@cooked_required
def test_flex_multi():
   
    spec1d_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science',
                            'spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits')

    msFlex = flexure.MultiSlitFlexure(s1dfile=spec1d_file) 
    # Parameters
    keck_deimos = load_spectrograph('keck_deimos')
    par = keck_deimos.default_pypeit_par()
    msFlex.init(keck_deimos, par['flexure'])
    # Init                    
    outfile = data_path('tst_multi_flex.fits')
    # INITIAL SKY LINE STUFF
    msFlex.measure_sky_lines()
    # FIT SURFACES
    msFlex.fit_mask_surfaces()
    # Apply
    msFlex.update_fit()
    # QA
    #mask = header['TARGET'].strip()
    #fnames = header['FILENAME'].split('.')
    #root = mask+'_'+fnames[2]
    #mdFlex.qa_plots('./', root)

    # Write
    msFlex.to_file(outfile, overwrite=True)

    # Read
    msFlex2 = flexure.MultiSlitFlexure.from_file(outfile)
    msFlex2.to_file(outfile, overwrite=True)

    # Check
    assert np.all(np.isclose(msFlex2.fit_b, msFlex.fit_b))

    # Clean up
    if os.path.isfile(outfile):
        os.remove(outfile)