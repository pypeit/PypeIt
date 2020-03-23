"""
Module to test TracePCA object.
"""
import os
import numpy as np
import pytest

from astropy.io import fits

from pypeit.tracepca import TracePCA

@pytest.fixture
def vec_coo():
    nvec = 50
    return np.linspace(0,1,nvec)

@pytest.fixture
def bogus_vectors(vec_coo):
    # Generate some bogus vectors
    nspec = 1000
    spec_coo = np.linspace(0,1,nspec)
    coeff0 = 0.1*np.square(vec_coo) + vec_coo + 10
    coeff1 = 2*vec_coo - 3
    base_vector = 3*spec_coo + np.power(spec_coo,3)
    return coeff0[None,:] + coeff1[None,:]*base_vector[:,None]


def test_build(vec_coo, bogus_vectors):
    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])

    assert pca.npca == 2, 'Incorrect number of components'
    assert pca.nspec == bogus_vectors.shape[0], 'Incorrect number of pixels'
    assert np.array_equal(vec_coo, pca.trace_coo), 'Coordinates do not match'
    # TODO: More checks?


def test_prediction(vec_coo, bogus_vectors):
    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])

    pred = pca.predict(0.5)
    assert pred.size == bogus_vectors.shape[0], 'Bad prediction'
    # TODO: More checks?


def test_write(vec_coo, bogus_vectors):

    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])

    ofile = 'junkpca.fits'
    if os.path.isfile(ofile):
        os.remove(ofile)
    fits.HDUList([fits.PrimaryHDU(), pca.to_hdu()]).writeto(ofile)
    os.remove(ofile)


def test_read(vec_coo, bogus_vectors):

    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])

    ofile = 'junkpca.fits'
    if os.path.isfile(ofile):
        os.remove(ofile)
    fits.HDUList([fits.PrimaryHDU(), pca.to_hdu()]).writeto(ofile)
    readpca = TracePCA.from_file(ofile)

    assert np.array_equal(pca.trace_coo, readpca.trace_coo), 'Bad read'
    assert np.array_equal(pca.pca_mean, readpca.pca_mean), 'Bad read'
    assert np.array_equal(pca.pca_coeffs, readpca.pca_coeffs), 'Bad read'
    assert np.array_equal(pca.pca_components, readpca.pca_components), 'Bad read'
    assert np.array_equal(pca.pca_bpm, readpca.pca_bpm), 'Bad read'
    assert pca.npca == readpca.npca, 'Bad read'
    assert pca.nspec == readpca.nspec, 'Bad read'
    for i in range(pca.npca):
        assert np.array_equal(pca.fit_coeff[i], readpca.fit_coeff[i]), 'Bad read'

    os.remove(ofile)

