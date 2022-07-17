"""
Module to run fetching remote package data from GitHub
"""

import requests

from pypeit import data
from pypeit.core.wavecal import waveio
from pypeit import __version__ 


def test_cloud_url():

    # The telgrid files live on a cloud server.  Test for file existance (or URL change)
    telgrid_file = "TelFit_MaunaKea_3100_26100_R20000.fits"
    _, telgrid_src = data.utils._build_remote_url(telgrid_file, "telluric/atm_grids",
                                               remote_host="s3_cloud")

    # Get Url; status code == 200 is success
    get = requests.head(telgrid_src[0])
    try:
        assert get.status_code == requests.codes.ok
    except AssertionError as err:
        raise Exception(f"Got status {get.status_code} (!= 200) for URL {telgrid_src[0]}") from err


def test_fetch_github_files():

    # These are commonly used files, do all three in one test
    # First test a `reid_arxiv` file
    data.fetch_remote_file("keck_deimos_600ZD.fits", "arc_lines/reid_arxiv",
                           force_update=True)

    # Next, try a `skisim` file
    data.fetch_remote_file("mktrans_zm_10_10.dat", "skisim",
                           force_update=True)
 
    # Finally, try a `sensfunc` file
    data.fetch_remote_file("keck_deimos_600ZD_sensfunc.fits", "sensfuncs",
                           force_update=True)


def test_load_sky_spectrum():

    # Load in the most common sky spectrum to ensure this function works 
    data.load_sky_spectrum("paranal_sky.fits")

def test_waveio_load_reid_arxiv():

    # Test the extension logic, given the download/cache system
    waveio.load_reid_arxiv("vlt_xshooter_vis1x1.fits")
    waveio.load_reid_arxiv("vlt_xshooter_vis1x1.json")
