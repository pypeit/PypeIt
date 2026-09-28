"""Regression tests for propagation of gain-corrected FITS units."""
import numpy as np
import pytest
from astropy.io import fits

from pypeit import dataPaths, exposure, outputfiles, spec2dobj, specobj, specobjs
from pypeit.images import buildimage, pypeitimage, rawimage
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests import tstutils


@pytest.fixture(params=[None, 'ADU'])
def kast_raw_file(request, tmp_path):
    """Copy a real Kast exposure with missing or ADU-valued input units."""
    path = tmp_path / 'b24.fits'
    with fits.open(dataPaths.tests.get_file_path('b24.fits.gz')) as hdus:
        for hdu in hdus:
            hdu.header.pop('BUNIT', None)
            if request.param is not None:
                hdu.header['BUNIT'] = request.param
        hdus.writeto(path)
    return path


@pytest.mark.remote_data
def test_gain_headers(kast_raw_file):
    spectrograph = load_spectrograph('shane_kast_blue')
    original_file = kast_raw_file.read_bytes()
    raw = rawimage.RawImage(str(kast_raw_file), spectrograph, 1)
    try:
        original_image = raw.rawimage.copy()
        original_headers = [hdu.header.copy() for hdu in raw.hdu]
        assert raw.headarr == original_headers, 'Initial headers should match the raw input.'
        raw.apply_gain()
        assert all(header['BUNIT'] == 'electron' for header in raw.headarr), \
            'Gain correction should set electron units in every processed header.'
        for amp, gain in enumerate(raw.detector[0].gain, start=1):
            pixels = (raw.rawdatasec_img == amp) | (raw.oscansec_img == amp)
            assert np.any(pixels), 'Each Kast amplifier should have pixels in the raw image.'
            assert np.allclose(raw.image[pixels], original_image[pixels] * gain), \
                'Science and overscan pixels should be multiplied by their amplifier gain.'
        corrected_image = raw.image.copy()
        raw.apply_gain()
        assert np.array_equal(raw.image, corrected_image), \
            'Calling apply_gain twice should not apply the correction a second time.'
        assert np.array_equal(raw.rawimage, original_image), \
            'Gain correction should leave the original raw pixels unchanged.'
        assert [hdu.header for hdu in raw.hdu] == original_headers, \
            'Gain correction should modify copied headers, not the input HDU headers.'
    finally:
        raw.hdu.close()
    assert kast_raw_file.read_bytes() == original_file, \
        'Gain correction should not modify the raw FITS file on disk.'


@pytest.mark.parametrize('units, bunit', [('e-', 'electron'), ('ADU', 'ADU')])
@pytest.mark.parametrize('image_class', [pypeitimage.PypeItImage, buildimage.BiasImage])
def test_image_units(tmp_path, units, bunit, image_class):
    img = image_class(np.ones((4, 4)), units=units)
    img.process_steps = ['apply_gain'] if units == 'e-' else []
    img.rawheadlist = [fits.Header({'BUNIT': bunit})]
    difference = img.sub(image_class(np.zeros((4, 4)), units=units))
    assert difference.process_steps == img.process_steps, \
        'Subtraction should preserve the science image processing history.'
    assert difference.rawheadlist[0]['BUNIT'] == bunit, \
        'Subtraction should preserve the science image header units.'
    assert difference.to_hdu(add_primary=True)[0].header['BUNIT'] == bunit, \
        'The subtracted image primary header should retain its units.'
    path = tmp_path / 'image.fits'
    img.to_file(path)
    with fits.open(path) as hdus:
        assert hdus[0].header['BUNIT'] == bunit, \
            'The primary header should record image units.'
        assert hdus[1].header['BUNIT'] == bunit, \
            'The image extension should record image units.'
    restored = image_class.from_file(path)
    assert restored.units == units, \
        'Image units should survive a FITS round trip.'
    assert restored.to_hdu(add_primary=True)[0].header['BUNIT'] == bunit, \
        'Rewriting a restored image should preserve BUNIT.'


@pytest.mark.remote_data
@pytest.mark.parametrize('apply_gain', [False, True])
def test_spectral_units(tmp_path, kast_raw_file, apply_gain):
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()
    par['rdx']['redux_path'] = str(tmp_path)
    fitstbl = PypeItMetaData(spectrograph, par, files=[str(kast_raw_file)], strict=True)
    spec = spec2dobj.Spec2DObj(sciimg=np.ones((4, 4)),
                              detector=tstutils.get_kastb_detector(),
                              ivarraw=None, skymodel=None, bkg_redux_skymodel=None,
                              objmodel=None, ivarmodel=None, scaleimg=None, waveimg=None,
                              bpmmask=None, sci_spat_flexure=None, sci_spec_flexure=None,
                              vel_type=None, vel_corr=None, slits=None, wavesol=None,
                              tilts=None, maskdef_designtab=None)
    spec.process_steps = ['apply_gain'] if apply_gain else []
    spectra = spec2dobj.AllSpec2DObj()
    spectra[spec.detname] = spec
    raw_header = fits.getheader(kast_raw_file)
    raw_units = raw_header.get('BUNIT')
    expected = 'electron' if apply_gain else raw_units
    primary = spectra.build_primary_hdr(raw_header, spectrograph)
    assert primary.get('BUNIT') == expected, \
        'The spec2d primary header should reflect whether gain was applied.'
    subheader = spectrograph.subheader_for_spec(fitstbl[0], primary)
    assert subheader.get('BUNIT') == expected, \
        'Spectral subheaders should preserve BUNIT.'

    # Exercise the real metadata, naming, and writing paths in save_exposure.
    obj = specobj.SpecObj('MultiSlit', spec.detname, SLITID=0)
    obj.BOX_WAVE = np.arange(4, dtype=float)
    obj.BOX_COUNTS = np.ones(4)
    obj.BOX_R_ASEC = 1.
    obj.SPAT_PIXPOS = 2.
    obj.SPAT_PIXPOS_ID = 2
    obj.SPAT_FRACPOS = 0.5
    obj.SPAT_FWHM = 1.
    obj.S2N = 1.
    obj.DETECTOR = spec.detector
    extracted = specobjs.SpecObjs([obj])
    exposure.save_exposure(spectrograph, fitstbl, par, 0, spectra, extracted,
                           str(tmp_path))
    outfile1d = outputfiles.spec_output_file(fitstbl, par, 0)
    outfile2d = outputfiles.spec_output_file(fitstbl, par, 0, twod=True)
    assert fits.getheader(outfile1d).get('BUNIT') == expected, \
        'The saved spec1d primary header should reflect the gain-corrected units.'
    assert outputfiles.spec_output_file(fitstbl, par, 0, ext='.txt').is_file(), \
        'Saving an extracted spectrum should also write its summary file.'
    with fits.open(outfile2d) as hdus:
        assert hdus[0].header.get('BUNIT') == expected, \
            'The saved spec2d primary header should reflect the gain-corrected units.'
        if apply_gain:
            assert hdus[1].header['BUNIT'] == 'electron', \
                'The gain-corrected science image extension should have electron units.'
    assert fits.getheader(kast_raw_file).get('BUNIT') == raw_units, \
        'Saving spectra should leave the original raw file units unchanged.'
