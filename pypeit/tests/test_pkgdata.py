"""
Module to test the various routines in `pypeit.data.utils`
"""

import requests

from linetools.spectra import xspectrum1d

from pypeit import data
from pypeit import __version__
from pypeit.core.wavecal import waveio


def test_cloud_url():

    # The telgrid files live on a cloud server.  Test for file existance (or URL change)
    telgrid_file = "TelFit_MaunaKea_3100_26100_R20000.fits"
    _, telgrid_src = data.utils._build_remote_url(telgrid_file, "telluric/atm_grids",
                                               remote_host="s3_cloud")

    # Get Url; status code == 200 is success
    get = requests.head(telgrid_src[0])
    assert (get.status_code == requests.codes.ok), \
           f"Got status {get.status_code} (!= 200) for URL {telgrid_src[0]}"


def test_fetch_github_files():

    # These are commonly used files, do all three in one test; the test just ensures
    #   the routines don't crash
    # First test a `reid_arxiv` file
    data.fetch_remote_file("keck_deimos_600ZD.fits", "arc_lines/reid_arxiv",
                           force_update=True)
    # Next, test a `skisim` file
    data.fetch_remote_file("mktrans_zm_10_10.dat", "skisim",
                           force_update=True)
    # Finally, test a `sensfunc` file
    data.fetch_remote_file("keck_deimos_600ZD_sensfunc.fits", "sensfuncs",
                           force_update=True)


def test_filepath_routines():

    # Test each of the get_*_filepath() routines to ensure they return a valid file path

    # reid_arxiv (returns tuple):
    filepath, _ = data.get_reid_arxiv_filepath("keck_deimos_600ZD.fits")
    assert filepath.is_file()

    # others (return just the filepath):
    assert data.get_skisim_filepath("mktrans_zm_10_10.dat").is_file()
    assert data.get_sensfunc_filepath("keck_deimos_600ZD_sensfunc.fits").is_file()
    assert data.get_linelist_filepath("ArI_lines.dat").is_file()


def test_load_sky_spectrum():

    # Load in the most common sky spectrum, check that the return is valid
    skyspec = data.load_sky_spectrum("paranal_sky.fits")
    assert isinstance(skyspec, xspectrum1d.XSpectrum1D)


def test_search_cache():

    # Make sure a junk search returns an empty list (and not None or something else)
    assert data.search_cache('junkymcjunkface.txt') == []

    # Place a file in the cache, and retrieve it
    data.write_file_to_cache(
        data.Paths.linelist / 'ArI_lines.dat',
        'totally_special_argon_lines.dat',
        'arc_lines/reid_arxiv'
    )
    assert data.search_cache('totally_special')[0].is_file()


def test_waveio_load_reid_arxiv():

    # Test the extension logic, given the download/cache system
    waveio.load_reid_arxiv("vlt_xshooter_vis1x1.fits")
    waveio.load_reid_arxiv("vlt_xshooter_vis1x1.json")
