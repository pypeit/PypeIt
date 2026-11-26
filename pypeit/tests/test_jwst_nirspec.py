"""
Test the JWST/NIRSpec-specific processing with synthetic data: the measurement
of the readnoise missing from the JWST ``VAR_RNOISE`` arrays, the conversion
between image counts and electrons, and the conversion to flux per pixel.
"""
import types

import numpy as np
import pytest
from astropy.io import fits

from pypeit import PypeItError
from pypeit.calibrations import NIRSpecSlitCalibrations
from pypeit.images import rawimage
from pypeit.spectrographs import jwst_nirspec
from pypeit.spectrographs.util import load_spectrograph


def _write_rate(path, rng, structure, readnoise, excess, detector='NRS1', exposure=1):
    # A synthetic _assign_wcs.fits file: no signal, fixed structure, and noise
    # with the VAR_RNOISE readnoise plus an excess (DN/s)
    shape = (256, 256)
    var_rnoise = np.full(shape, readnoise**2, dtype=np.float32)
    sci = structure + rng.normal(scale=np.sqrt(readnoise**2 + excess**2), size=shape)
    hdr = fits.Header()
    hdr.update(INSTRUME='NIRSPEC', PROGRAM='01967', OBSERVTN='012', VISIT='001', VISITGRP='03',
               SEQ_ID='1', ACT_ID='02', DETECTOR=detector, EXP_TYPE='NRS_FIXEDSLIT',
               GRATING='G395M', FILTER='F290LP', READPATT='NRSIRS2RAPID', NGROUPS=45, NINTS=1,
               SUBARRAY='FULL', FXD_SLIT='S200A2', EXPOSURE=str(exposure))
    fits.HDUList([fits.PrimaryHDU(header=hdr),
                  fits.ImageHDU(sci.astype(np.float32), name='SCI'),
                  fits.ImageHDU(np.zeros(shape, dtype=np.uint32), name='DQ'),
                  fits.ImageHDU(var_rnoise, name='VAR_RNOISE')]).writeto(path)


class _MetaData:
    # Stand-in for PypeItMetaData with only what cache_metadata uses
    def __init__(self, files):
        self.files = [str(f) for f in files]

    def __getitem__(self, key):
        return np.array([f.split('/')[-1] for f in self.files])

    def frame_paths(self, indx):
        return [f for f, i in zip(self.files, indx) if i]


def test_measure_rn_excess(tmp_path):
    # The excess readnoise measured on exposure differences ignores the fixed
    # structure of the pixels without signal
    rng = np.random.default_rng(seed=11)
    readnoise, excess = 0.0013, 0.002
    # Faint fixed structure, below the threshold used to select the pixels
    # without signal, with a pixel-to-pixel scatter comparable to the noise
    structure = rng.normal(scale=0.002, size=(256, 256))
    files = [tmp_path / f'jw01967012001_03102_0000{i}_nrs1_assign_wcs.fits' for i in range(1, 4)]
    for i, f in enumerate(files):
        _write_rate(f, rng, structure, readnoise, excess, exposure=i+1)
    groups = jwst_nirspec.group_exposures(files)
    assert len(groups) == 1, 'Exposures with the same configuration should be grouped together.'
    rn_excess, npairs = jwst_nirspec.measure_rn_excess(groups[0], min_pixels=1000)
    assert npairs == 2, 'The consecutive exposures of the group should be paired.'
    assert np.isclose(rn_excess, excess, rtol=0.05), \
        'The excess readnoise should be measured without the fixed structure.'


def test_get_rn_excess(tmp_path):
    # The spectrograph measures the excess readnoise of the exposures listed in
    # the metadata table, grouping them by configuration; an exposure without
    # another exposure of the same configuration is not measured
    rng = np.random.default_rng(seed=12)
    structure = np.zeros((256, 256))
    nrs1 = [tmp_path / f'jw01967012001_03102_0000{i}_nrs1_assign_wcs.fits' for i in range(1, 3)]
    for i, f in enumerate(nrs1):
        _write_rate(f, rng, structure, 0.0013, 0.002, exposure=i+1)
    nrs2 = tmp_path / 'jw01967012001_03102_00001_nrs2_assign_wcs.fits'
    _write_rate(nrs2, rng, structure, 0.0013, 0.002, detector='NRS2')
    # A file that is not a full-frame science file is not used
    other = tmp_path / 'jw01967012001_03102_00001_nrs1_cal.fits'
    other.touch()

    spec = load_spectrograph('jwst_nirspec')
    assert spec.get_rn_excess(nrs1[0]) is None, \
        'Without the metadata, the excess readnoise should not be measured.'
    spec.cache_metadata(_MetaData(nrs1 + [nrs2, other]))
    groups = jwst_nirspec.group_exposures(spec._rn_excess_files)
    assert sorted(len(g) for g in groups) == [1, 2], \
        'Exposures of different detectors should be in different groups.'
    assert np.isclose(spec.get_rn_excess(nrs1[0]), 0.002, rtol=0.05), \
        'The excess readnoise of a pair of exposures should be measured.'
    assert spec.get_rn_excess(nrs1[1]) == spec.get_rn_excess(nrs1[0]), \
        'All the exposures of a group should have the same excess readnoise.'
    assert spec.get_rn_excess(nrs2) is None, \
        'An exposure without another exposure of the same configuration should not be measured.'


def test_counts_per_electron():
    # The conversion from electrons to image counts (the count rate times
    # EFFEXPTM) uses the gain of each detector and the ramp time, (NGROUPS-1)
    # * TGROUP * NINTS, and is 1 in the gap of a mosaic
    img = rawimage.NIRSpecRawImage.__new__(rawimage.NIRSpecRawImage)
    img.headarr = [fits.Header({'NGROUPS': 11, 'TGROUP': 100., 'NINTS': 2})]
    img.detector = [types.SimpleNamespace(det=1, gain=np.array([1.0])),
                    types.SimpleNamespace(det=2, gain=np.array([1.25]))]
    img.exptime = 11 * 100. * 2
    img.nimg = 1
    img.image = np.zeros((1, 10, 20))
    img.det_img = np.zeros((1, 10, 20), dtype=int)
    img.det_img[..., :8] = 1
    img.det_img[..., 12:] = 2
    assert np.isclose(img.electrons_per_rate(1), 1.25 * 10 * 100. * 2), \
        'The electrons per unit rate should be the gain times the ramp time.'
    scale = img.counts_per_electron()
    assert np.allclose(scale[..., :8], 2200 / 2000), 'Wrong conversion for the first detector.'
    assert np.allclose(scale[..., 12:], 2200 / 2500), 'Wrong conversion for the second detector.'
    assert np.allclose(scale[..., 8:12], 1.), 'The conversion in the mosaic gap should be 1.'


def test_flux_per_pixel_scale():
    # The slits of the JWST _cal.fits files are converted to MJy per pixel:
    # surface brightnesses are multiplied by the pixel solid angle
    hdr = fits.Header({'BUNIT': 'MJy', 'PIXAR_SR': 4.8e-13})
    assert NIRSpecSlitCalibrations.flux_per_pixel_scale(hdr) == 1., \
        'A slit already in MJy should not be rescaled.'
    hdr['BUNIT'] = 'MJy/sr'
    assert NIRSpecSlitCalibrations.flux_per_pixel_scale(hdr) == 4.8e-13, \
        'A slit in MJy/sr should be multiplied by the pixel solid angle.'
    hdr['BUNIT'] = 'DN/s'
    with pytest.raises(PypeItError):
        NIRSpecSlitCalibrations.flux_per_pixel_scale(hdr)
