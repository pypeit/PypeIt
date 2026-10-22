"""
Module to run tests on FlatField class
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import flatfield
from pypeit.core.bspline import Knots
from pypeit.containers.bspline import BSplineContainer
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests.tstutils import data_output_path
# from pypeitdev.flats.dev_flats import pixelflat


def test_flatimages():
    tmp = np.ones((1000, 100)) * 10.
    x = np.random.rand(500)
    # Create bspline
    bsp = 0.01 * (np.max(x) - np.min(x))
    spat_bspline1 = BSplineContainer(x=x, knots=Knots(spacing=bsp), nord=4)
    spat_bspline2 = BSplineContainer(x=x, knots=Knots(spacing=bsp), nord=4)
    instant_dict = dict(pixelflat_raw=tmp,
                        pixelflat_norm=np.ones_like(tmp),
                        pixelflat_model=None,
                        pixelflat_spat_bsplines=np.asarray([spat_bspline1, spat_bspline2]),
                        pixelflat_spec_illum=None,
                        pixelflat_spat_flexure=None,
                        illumflat_raw=tmp,
                        illumflat_spat_bsplines=np.asarray([spat_bspline1, spat_bspline2]),
                        illumflat_spat_flexure=None,
                        spat_id=np.asarray([100, 200]),
                        PYP_SPEC="specname")

    flatImages = flatfield.FlatImages(**instant_dict)
    assert flatImages.pixelflat_model is None
    assert flatImages.pixelflat_spec_illum is None
    assert flatImages.pixelflat_spat_bsplines is not None
    flatImages.set_paths(data_output_path(''), 'A', '1', 'DET01')

    # I/O
    ofile = Path(flatImages.get_path()).absolute()
    flatImages.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    _flatImages = flatfield.FlatImages.from_file(ofile)

    # Test
    for key in instant_dict.keys():
        if key == 'pixelflat_spat_bsplines':
            assert np.array_equal(flatImages[key][0].breakpoints,
                                  _flatImages[key][0].breakpoints), 'pixelflat breakpoints changed'
            continue
        if key == 'illumflat_spat_bsplines':
            assert np.array_equal(flatImages[key][0].breakpoints,
                                  _flatImages[key][0].breakpoints), 'illumflat breakpoints changed'
            continue
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(flatImages[key],_flatImages[key])
        else:
            assert flatImages[key] == _flatImages[key]

    ofile.unlink()


def test_fit_det_response():
    spec = load_spectrograph('keck_kcwi')
    # Generate a good pixel mask
    frsize = 4100
    gpm = np.ones((frsize,frsize), dtype=bool)
    # Generate a fake image
    sinemodel = lambda xx, yy, amp, scl, phase, wavelength, angle: 1 + (amp + xx * scl) * np.sin(
                2 * np.pi * (xx * np.cos(angle*np.pi / 180.0) + yy * np.sin(angle*np.pi / 180.0)) / wavelength + phase)
    x = np.arange(frsize)
    y = np.arange(frsize)
    xx, yy = np.meshgrid(x, y, indexing='ij')
    amp, scale, wavelength, phase, angle = 0.02, 0.0, 1.41*frsize/31.5, 0.0, -45.34
    img = sinemodel(xx, yy, amp, scale, phase, wavelength, angle)
    model = spec.fit_2d_det_response(img, gpm)
    assert np.allclose(img, model, atol=0.001), 'structure fitting failed.'


def test_load_pixflat_preserves_calib_key(tmp_path):
    """
    ``load_pixflat`` merges a user/archived pixel flat into the existing
    ``FlatImages``.  The merge rebuilds the object from datamodel fields
    only, so the calibration identifiers must be carried over explicitly;
    otherwise ``calib_key`` is None and ``get_path()`` raises (this broke
    the flats state recording for keck_lris_blue, which uses an archived
    pixel flat).
    """
    spec = load_spectrograph('keck_lris_blue')
    cal_dir = str(tmp_path)
    # Minimal archived pixel-flat file: one extension per detector, named
    # '<detname>-PIXELFLAT_NORM', as load_pixflat expects.
    pixflat_file = str(tmp_path / 'mypixflat.fits')
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(np.full((10, 10), 0.5, dtype=float),
                      name='DET01-PIXELFLAT_NORM')]).writeto(pixflat_file,
                                                             overwrite=True)
    # An existing FlatImages with calibration identifiers set
    flats = flatfield.FlatImages(pixelflat_norm=np.ones((10, 10)))
    flats.calib_key = 'B_1_DET01'
    flats.calib_dir = cal_dir
    merged = flatfield.load_pixflat(pixflat_file, spec, 1, flats,
                                    calib_dir=cal_dir, chk_version=False)
    # The merged result keeps the user pixel flat...
    assert np.allclose(merged.pixelflat_norm, 0.5)
    # ...and, crucially, the calibration identifiers so get_path() works.
    assert merged.calib_key == 'B_1_DET01'
    assert merged.calib_dir == cal_dir
    assert Path(merged.get_path()).name == 'Flat_B_1_DET01.fits'




def test_flatfield_spec_illum_flexure():
    # The relative spectral illumination is defined on the slits of the pixel
    # flat, and must be remapped to the spatial flexure of the frame
    from pypeit.core import flexure
    from pypeit.images.rawimage import RawImage
    from pypeit.par.pypeitpar import ProcessImagesPar
    from pypeit.slittrace import SlitTraceSet
    nspec, nspat = 40, 80
    left = np.full((nspec, 2), 8.) + np.array([0., 35.])[None, :]
    right = left + 28.
    slits = SlitTraceSet(left_init=left, right_init=right, pypeline='SlicerIFU', nspat=nspat,
                         PYP_SPEC='dummy')
    flat_flexure = np.full((2, 2), 0.4)
    sci_flexure = np.array([[2.2, 2.6], [1.8, 2.0]])
    # A different spectral illumination for each slit, varying across the slit
    spec_illum = np.ones((nspec, nspat))
    slitidx_img = slits.slit_img(pad=0, initial=True, spat_flexure=flat_flexure, use_spatial=False)
    for slit_idx in range(slits.nslits):
        spec, spat = np.where(slitidx_img == slit_idx)
        spec_illum[spec, spat] = 1.2 + 0.3 * slit_idx + 0.002 * spat
    flatimages = flatfield.FlatImages(pixelflat_norm=np.ones((nspec, nspat)),
                                      pixelflat_spec_illum=spec_illum,
                                      pixelflat_spat_flexure=flat_flexure,
                                      PYP_SPEC='dummy', spat_id=slits.spat_id)
    # Build a minimal processed image (without reading a raw file)
    raw = object.__new__(RawImage)
    raw.par = ProcessImagesPar(use_overscan=False, use_illumflat=False, use_specillum=True)
    raw.steps = dict(flatfield=False)
    raw.nimg = 1
    raw.image = np.ones((1, nspec, nspat))
    raw.spat_flexure_shift = sci_flexure
    raw.flatfield(flatimages, slits=slits)
    expected = flexure.spat_flexure_remap(spec_illum, slits, spat_flexure=sci_flexure,
                                          src_spat_flexure=flat_flexure, initial=True)
    assert not np.array_equal(expected, spec_illum), 'Test requires a change of flexure'
    assert np.allclose(raw.image[0] * expected, 1.), \
        'The spectral illumination should be applied at the spatial flexure of the frame'


def test_flatfield_fit_padded_slits(tmp_path):
    # The spectral model of each slit must be evaluated over the padded slit,
    # so that the illumination profile close to (and beyond) the slit edges,
    # and therefore the tweaked slit edges, are correct
    from scipy.special import erf, erfinv
    from pypeit.images.detector_container import DetectorContainer
    from pypeit.images.pypeitimage import PypeItImage
    from pypeit.par.pypeitpar import FlatFieldPar
    from pypeit.slittrace import SlitTraceSet
    from pypeit.tests import tstutils
    from pypeit.wavetilts import WaveTilts
    nspec, nspat = 200, 100
    spec = np.arange(nspec)
    left = np.column_stack([12.3 + 0.004 * spec, 56.6 + 0.004 * spec])
    right = left + 31.
    slits = SlitTraceSet(left_init=left, right_init=right, pypeline='MultiSlit', nspat=nspat,
                         PYP_SPEC='dummy')
    slits.set_paths(tmp_path, 'A', '1', 'DET01')
    slits.to_file()
    # The tilts are the spectral coordinate
    coeffs = np.zeros((2, 2, 2))
    coeffs[:2, 0, :] = 0.5
    wvtilts = WaveTilts(coeffs=coeffs, nslit=slits.nslits, spat_id=slits.spat_id,
                        spat_order=np.array([1, 1]), spec_order=np.array([1, 1]),
                        func2d='legendre2d', slits_filename=Path(slits.get_path()).name)
    wvtilts.set_paths(tmp_path, 'A', '1', 'DET01')
    # The flat is not illuminated in the first 40 spectral rows, so the
    # spectral model is not defined there
    blaze = 1e4 * np.clip((spec - 40) / 50., 0, None) * (1.5 - spec / nspec)
    sig = 0.8
    xx = np.arange(nspat)[None, :]
    illum = sum(0.5 * (erf((xx - left[:, i, None]) / (np.sqrt(2) * sig))
                       - erf((xx - right[:, i, None]) / (np.sqrt(2) * sig)))
                for i in range(slits.nslits))
    model = blaze[:, None] * illum + 5.
    rng = np.random.default_rng(1)
    rawflat = PypeItImage(rng.normal(model, np.sqrt(model)), ivar=1 / model)
    rawflat.detector = DetectorContainer(**tstutils.default_detector())
    flatpar = FlatFieldPar()
    flatField = flatfield.FlatField(rawflat, load_spectrograph('shane_kast_blue'), flatpar, slits,
                                    wavetilts=wvtilts)
    flatField.fit(doqa=False)
    # The illumination profile beyond the slit edges (in the padding) must be
    # close to 0
    width = np.median(right - left)
    outside = np.array([-3., width + 3.]) / width
    for bspl in flatField.list_of_spat_bsplines:
        assert np.all(np.absolute(bspl.value(outside)[0]) < 0.1), \
            'Illumination profile beyond the slit edges is wrong'
    # The tweaked edges must be where the illumination profile reaches the
    # threshold
    inset = sig * np.sqrt(2) * erfinv(2 * flatpar['tweak_slits_thresh'] - 1)
    assert np.allclose(np.median(slits.left_tweak - left, axis=0), inset, atol=0.35), \
        'Incorrect tweaked left edges'
    assert np.allclose(np.median(slits.right_tweak - right, axis=0), -inset, atol=0.35), \
        'Incorrect tweaked right edges'
