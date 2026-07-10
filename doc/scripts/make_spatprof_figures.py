"""
Generate example spatial-profile QA figures for doc/spatprof.rst.

Produces two PNGs in doc/figures/ using the same synthetic data construction
as :func:`~pypeit.tests.test_spatialprofile.test_fourier_flux_and_curved_trace`:
a Fourier-generated spectrum with rapid flux variations and a quadratic curved
trace.

Usage
-----
Run from the repo root after installing pypeit in development mode::

    python doc/scripts/make_spatprof_figures.py
"""
from importlib import resources
import os

import matplotlib
matplotlib.use('Agg')
import numpy as np

from pypeit.core import spatialprofile
from pypeit.tests.test_spatialprofile import make_fourier_spectrum, make_profile_inputs

def _make_inputs(sn_ratio, seed=42):
    """Build synthetic inputs matching test_fourier_flux_and_curved_trace."""
    nspec, nspat = 200, 100
    fwhm = 4.0

    flux_level = make_fourier_spectrum(
        nspec, nmodes=40, mean=1000.0, std=50.0, min_period=20., power=2.
    )

    t = np.linspace(-1.0, 1.0, nspec)
    trace_offset = 10.0 * (2.0 * t ** 2 - 1.0)

    image, ivar, waveimg, thismask, spat_img, _, wave, flux, \
        fluxivar, inmask, true_trace = make_profile_inputs(
            nspec=nspec, nspat=nspat, fwhm=fwhm, sn_ratio=sn_ratio,
            flux_level=flux_level, trace_offset=trace_offset, seed=seed,
        )
    return (image, ivar, waveimg, thismask, spat_img, true_trace,
            wave, flux, fluxivar, inmask, fwhm)


def make_bspline_figure(out_path):
    """B-spline path: S/N=20, Fourier flux, quadratic curved trace."""
    image, ivar, waveimg, thismask, spat_img, trace_in, \
        wave, flux, fluxivar, inmask, fwhm = _make_inputs(sn_ratio=20.0)
    spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask,
        thisfwhm=fwhm, sn_gauss=4.0,
        obj_string='Synthetic point source, S/N=20',
        generate_qa=str(out_path),
    )


def make_gaussian_figure(out_path):
    """Gaussian fallback path: S/N=2, Fourier flux, quadratic curved trace."""
    image, ivar, waveimg, thismask, spat_img, trace_in, \
        wave, flux, fluxivar, inmask, fwhm = _make_inputs(sn_ratio=2.0)
    spatialprofile.fit_profile(
        image=image, ivar=ivar, waveimg=waveimg, thismask=thismask,
        spat_img=spat_img, trace_in=trace_in, wave=wave, flux=flux,
        fluxivar=fluxivar, inmask=inmask,
        thisfwhm=fwhm, sn_gauss=4.0,
        obj_string='Synthetic point source, S/N=2 (Gaussian fallback)',
        generate_qa=str(out_path),
    )


def main():
    root = resources.files('pypeit').parent
    out_dir = root / 'doc' / 'figures'


    os.makedirs(out_dir, exist_ok=True)

    out_bsp = out_dir / 'spatprof_example_bspline.png'
    make_bspline_figure(out_bsp)
    print(f'Wrote {out_bsp}')

    out_gauss = out_dir / 'spatprof_example_gaussian.png'
    make_gaussian_figure(out_gauss)
    print(f'Wrote {out_gauss}')


if __name__ == '__main__':
    main()
