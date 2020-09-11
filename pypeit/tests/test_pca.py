"""
Module to test TracePCA object.
"""
import os
from IPython import embed
import numpy as np
import pytest

from astropy.io import fits

from pypeit.tracepca import TracePCA
from pypeit.core.fitting import PypeItFit
from pypeit.tests.tstutils import data_path

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
    assert pca.pca_coeffs_model.size == pca.npca, 'Incorrect number of models'
    assert isinstance(pca.pca_coeffs_model[0], PypeItFit)
    # TODO: More checks?


def test_prediction(vec_coo, bogus_vectors):
    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])

    pred = pca.predict(0.5)
    assert pred.size == bogus_vectors.shape[0], 'Bad prediction'
    assert pred[0] > 10 and pred[-1] < 3, 'Prediction changed'
    # TODO: More checks?


def test_rms():
    """Test on some real data."""
    center = np.load(data_path('example_trace_deimos_1200G_M_7750.npz'))['center']

    pca = TracePCA(trace_cen=center, npca=2)
    pca.build_interpolator(np.array([3,1]), function='legendre') #, debug=True)
    pca_center = pca.predict(center[pca.reference_row,:])

    rms = np.std(center-pca_center, axis=0)
    assert np.sum(rms > 0.2) == 1, 'Change in the accuracy of the PCA.'


def test_to_hdu(vec_coo, bogus_vectors):
    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])
    hdu = pca.to_hdu(add_primary=True)
    hdu_names = [h.name for h in hdu][1:]
    assert hdu_names == ['PCA', 'PCA_MODEL_1', 'PCA_MODEL_2'], 'Bad HDU extensions'

def test_write(vec_coo, bogus_vectors):
    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])
    ofile = 'junkpca.fits'
    if os.path.isfile(ofile):
        os.remove(ofile)
    pca.to_file(ofile)
    os.remove(ofile)

def test_read(vec_coo, bogus_vectors):

    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])
    ofile = 'junkpca.fits'
    if os.path.isfile(ofile):
        os.remove(ofile)
    pca.to_file(ofile)
    readpca = TracePCA.from_file(ofile)

    assert np.array_equal(pca.trace_coo, readpca.trace_coo), 'Bad read'
    assert np.array_equal(pca.pca_mean, readpca.pca_mean), 'Bad read'
    assert np.array_equal(pca.pca_coeffs, readpca.pca_coeffs), 'Bad read'
    assert np.array_equal(pca.pca_components, readpca.pca_components), 'Bad read'
    assert pca.npca == readpca.npca, 'Bad read'
    assert pca.nspec == readpca.nspec, 'Bad read'
    for i in range(pca.npca):
        assert np.array_equal(pca.pca_coeffs_model[i].fitc, readpca.pca_coeffs_model[i].fitc), \
                    'Bad read'
    os.remove(ofile)

def test_write_two(vec_coo, bogus_vectors):
    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])
    pca2 = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    ofile = 'junkpca.fits'
    if os.path.isfile(ofile):
        os.remove(ofile)
    fits.HDUList([fits.PrimaryHDU()]
                  + pca.to_hdu(hdu_prefix='LEFT_')
                  + pca2.to_hdu(hdu_prefix='RIGHT_')).writeto(ofile)
    os.remove(ofile)

def test_read_two(vec_coo, bogus_vectors):
    pca = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    pca.build_interpolator([1,2])
    pca2 = TracePCA(trace_cen=bogus_vectors, npca=2, coo=vec_coo)
    ofile = 'junkpca.fits'
    if os.path.isfile(ofile):
        os.remove(ofile)
    hdu = fits.HDUList([fits.PrimaryHDU()]
                        + pca.to_hdu(hdu_prefix='LEFT_')
                        + pca2.to_hdu(hdu_prefix='RIGHT_'))
    _pca = TracePCA.from_hdu(hdu, hdu_prefix='LEFT_')
    _pca2 = TracePCA.from_hdu(hdu, hdu_prefix='RIGHT_')

    assert _pca.pca_coeffs_model is not None, 'Should find LEFT_PCA model'
    assert _pca2.pca_coeffs_model is None, 'Should not find RIGHT_PCA model'
    for i in range(pca.npca):
        assert np.array_equal(pca.pca_coeffs_model[i].fitc, _pca.pca_coeffs_model[i].fitc), \
                'Bad read'
