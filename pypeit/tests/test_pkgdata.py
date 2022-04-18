"""
Module to run fetching remote package data from GitHub
"""

import pkg_resources
import requests

from pypeit import data
from pypeit import __version__ 


def test_cloud_url():

    # The telgrid files live on a cloud server.  Test for file existance (or URL change)
    telgrid_file = "TelFit_MaunaKea_3100_26100_R20000.fits"
    telgrid_url = data.utils._build_remote_url(telgrid_file, "telluric/atm_grids",
                                               remote_host="s3_cloud")

    # Get Url; status code == 200 is success
    get = requests.head(telgrid_url)
    try:
        assert get.status_code == requests.codes.ok
    except AssertionError as err:
        raise Exception(f"Got status {get.status_code} (!= 200) for URL {telgrid_url}") from err


def test_fetch_github_files():

    # Define the version as the most recent stable release for GitHub fetch
    pv = pkg_resources.parse_version(__version__)
    if pv.is_devrelease:
        prev_release = f"{pv.major}.{pv.minor}.{pv.micro - 1}"
    else:
        prev_release = __version__

    # These are commonly used files, do all three in one test
    # First test a `reid_arxiv` file
    data.fetch_remote_file("keck_deimos_600ZD.fits", "arc_lines/reid_arxiv",
                           test_version=prev_release, force_update=True)

    # Next, try a `skisim` file
    data.fetch_remote_file("mktrans_zm_10_10.dat", "skisim",
                           test_version=prev_release, force_update=True)
 
    # Finally, try a `sensfunc` file
    data.fetch_remote_file("keck_deimos_600ZD_sensfunc.fits", "sensfuncs",
                           test_version=prev_release, force_update=True)


def test_load_sky_spectrum():

    # Load in the most common sky spectrum to ensure this function works 
    data.load_sky_spectrum("paranal_sky.fits")
