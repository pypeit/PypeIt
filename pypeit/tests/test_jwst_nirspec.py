"""
Test the measurement of the readnoise missing from the JWST ``VAR_RNOISE``
arrays, used by the JWST/NIRSpec noise model, with synthetic files.
"""
import numpy as np
from astropy.io import fits

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
