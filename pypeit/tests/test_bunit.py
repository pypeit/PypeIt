"""Regression tests for propagation of gain-corrected FITS units."""
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from astropy.io import fits

from pypeit import exposure, spec2dobj, specobj, specobjs
from pypeit.images import buildimage, pypeitimage, rawimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests import tstutils


@pytest.mark.parametrize('nimg', [1, 2])
def test_gain_headers(nimg):
    detector = tstutils.get_kastb_detector()
    detector.gain = np.array([2.])
    shape = (4, 4) if nimg == 1 else (nimg, 4, 4)
    image = np.ones(shape)
    hdus = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU()])
    hdus[1].header['BUNIT'] = 'ADU'
    spectrograph = SimpleNamespace(
        get_rawimage=lambda filename, det: (
            detector if nimg == 1 else SimpleNamespace(detectors=[detector]*nimg),
            image, hdus, 1., np.ones(shape, dtype=int), np.zeros(shape, dtype=int)),
        get_headarr=lambda hdu: [h.header for h in hdu])
    raw = rawimage.RawImage('unused.fits', spectrograph, 1 if nimg == 1 else (1, 2))
    assert 'BUNIT' not in raw.headarr[0]
    assert raw.headarr[1]['BUNIT'] == 'ADU'
    raw.apply_gain()
    assert np.all(raw.image == 2.)
    assert all(header['BUNIT'] == 'electron' for header in raw.headarr)
    raw.apply_gain()
    assert np.all(raw.image == 2.)
    # The original input headers and pixels are not modified.
    assert 'BUNIT' not in hdus[0].header
    assert hdus[1].header['BUNIT'] == 'ADU'
    assert np.all(raw.rawimage == 1.)


@pytest.mark.parametrize('units, bunit', [('e-', 'electron'), ('ADU', 'ADU')])
@pytest.mark.parametrize('image_class', [pypeitimage.PypeItImage, buildimage.BiasImage])
def test_image_units(tmp_path, units, bunit, image_class):
    img = image_class(np.ones((4, 4)), units=units)
    img.process_steps = ['apply_gain'] if units == 'e-' else []
    img.rawheadlist = [fits.Header({'BUNIT': bunit})]
    difference = img.sub(image_class(np.zeros((4, 4)), units=units))
    assert difference.process_steps == img.process_steps
    assert difference.rawheadlist[0]['BUNIT'] == bunit
    assert difference.to_hdu(add_primary=True)[0].header['BUNIT'] == bunit
    path = tmp_path / 'image.fits'
    img.to_file(path)
    with fits.open(path) as hdus:
        assert hdus[0].header['BUNIT'] == bunit
        assert hdus[1].header['BUNIT'] == bunit
    restored = image_class.from_file(path)
    assert restored.units == units
    assert restored.to_hdu(add_primary=True)[0].header['BUNIT'] == bunit


@pytest.mark.parametrize('apply_gain', [False, True])
@pytest.mark.parametrize('raw_units', [None, 'ADU'])
def test_spectral_units(tmp_path, monkeypatch, apply_gain, raw_units):
    spectrograph = load_spectrograph('shane_kast_blue')
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
    raw_header = fits.Header()
    if raw_units is not None:
        raw_header['BUNIT'] = raw_units
    expected = 'electron' if apply_gain else raw_units
    primary = spectra.build_primary_hdr(raw_header, spectrograph)
    assert primary.get('BUNIT') == expected
    subheader = spectrograph.subheader_for_spec(
        {'filename': 'raw.fits'}, primary, allow_missing=True)
    assert subheader.get('BUNIT') == expected

    # Exercise save_exposure, which rereads the untouched input file header.
    raw_path = tmp_path / 'raw.fits'
    fits.PrimaryHDU(header=raw_header).writeto(raw_path)
    class Fitstbl:
        def __getitem__(self, frame):
            return {'filename': 'raw.fits'}

        def frame_paths(self, frame):
            return str(raw_path)

    original_subheader = spectrograph.subheader_for_spec
    monkeypatch.setattr(spectrograph, 'subheader_for_spec',
                        lambda row, hdr: original_subheader(row, hdr, allow_missing=True))
    monkeypatch.setattr(exposure.outputfiles, 'science_path', lambda par: tmp_path)
    monkeypatch.setattr(exposure.outputfiles, 'spec_output_file',
                        lambda *args, twod=False, **kwargs:
                        tmp_path / ('spec2d.fits' if twod else 'spec1d.fits'))
    obj = specobj.SpecObj('MultiSlit', spec.detname, SLITID=0)
    obj.BOX_WAVE = np.arange(4, dtype=float)
    obj.BOX_COUNTS = np.ones(4)
    obj.DETECTOR = spec.detector
    extracted = specobjs.SpecObjs([obj])
    monkeypatch.setattr(specobjs.SpecObjs, 'write_info', Mock())
    par = spectrograph.default_pypeit_par()
    exposure.save_exposure(spectrograph, Fitstbl(), par, 0, spectra, extracted,
                           str(tmp_path))
    assert fits.getheader(tmp_path / 'spec1d.fits').get('BUNIT') == expected
    with fits.open(tmp_path / 'spec2d.fits') as hdus:
        assert hdus[0].header.get('BUNIT') == expected
        if apply_gain:
            assert hdus[1].header['BUNIT'] == 'electron'
    assert fits.getheader(raw_path).get('BUNIT') == raw_units
