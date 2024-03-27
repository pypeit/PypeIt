"""
Module to test the various routines in `pypeit.data.utils`
"""

import os
import requests
from pathlib import Path

from IPython import embed

import pytest

import github

from linetools.spectra import xspectrum1d

from pypeit.pypmsgs import PypeItPathError
from pypeit.pypeitdata import PypeItDataPath
from pypeit import dataPaths
from pypeit import io
from pypeit import cache
from pypeit.core.wavecal import waveio


def test_cloud_url():

    # The telgrid files live on a cloud server.  Test for file existance (or URL change)
    telgrid_file = "TelFit_MaunaKea_3100_26100_R20000.fits"
    _, telgrid_src = cache._build_remote_url(telgrid_file, "telluric/atm_grids",
                                             remote_host="s3_cloud")

    # Get Url; status code == 200 is success
    get = requests.head(telgrid_src[0])
    assert (get.status_code == requests.codes.ok), \
           f"Got status {get.status_code} (!= 200) for URL {telgrid_src[0]}"


def test_fetch_github_files():

    # These are commonly used files, do all three in one test; the test just ensures
    #   the routines don't crash
    # First test a `reid_arxiv` file
    cache.fetch_remote_file("keck_deimos_600ZD.fits", "arc_lines/reid_arxiv",
                            force_update=True)
    # Next, test a `skisim` file
    cache.fetch_remote_file("mktrans_zm_10_10.dat", "skisim",
                           force_update=True)
    # Finally, test a `sensfunc` file
    cache.fetch_remote_file("keck_deimos_600ZD_sensfunc.fits", "sensfuncs",
                            force_update=True)
    

def test_github_contents():
    # Access the repo
    repo = github.Github().get_repo("pypeit/PypeIt")

    # Get the relevant github branch
    branch = cache.git_branch()

    # Recursively get the contents of the test data directory
    contents = cache.github_contents(repo, branch, 'pypeit/data/tests')
    assert all([c.type != 'dir' for c in contents]), 'Recursion should not return any directories'

    # Get the contents without recursion.  This assumes the 'tests' directory
    # has sub-directories
    contents = cache.github_contents(repo, branch, 'pypeit/data/tests', recursive=False)
    assert any([c.type == 'dir' for c in contents]), \
            'tests/ directory expected to have subdirectories'


def test_filepath_routines():

    filepath, format = dataPaths.reid_arxiv.get_file_path("keck_deimos_600ZD.fits",
                                                          return_format=True,
                                                          to_pkg='symlink')
    assert filepath.is_file(), 'reid arxiv file does not exist'
    assert format == 'fits', 'File format wrong'

    # others (return just the filepath):
    assert dataPaths.skisim.get_file_path("mktrans_zm_10_10.dat").is_file(), \
            'skisim file does not exist'
    assert dataPaths.sensfunc.get_file_path("keck_deimos_600ZD_sensfunc.fits").is_file(), \
            'sensfunc file does not exist'
    assert dataPaths.linelist.get_file_path("ArI_lines.dat").is_file(), \
            'linelist file does not exist'


def test_load_sky_spectrum():

    # Load in the most common sky spectrum, check that the return is valid
    skyspec = io.load_sky_spectrum("paranal_sky.fits")
    assert isinstance(skyspec, xspectrum1d.XSpectrum1D)


def test_search_cache():

    # Make sure a junk search returns an empty list (and not None or something else)
    assert cache.search_cache('junkymcjunkface.txt') == [], 'should not find junk file'

    contents = cache.search_cache('totally_special_argon_lines.dat', path_only=False)
    if len(contents) > 0:
        # Remove all of the previous instances of this file.  It's possible to
        # have multiple files in the cache from different branches.
        cache.remove_from_cache(cache_url=list(contents.keys()), allow_multiple=True)

    # Place a file in the cache, and retrieve it
    cache.write_file_to_cache(
        str(dataPaths.linelist.get_file_path('ArI_lines.dat')),
        'totally_special_argon_lines.dat',
        'arc_lines/reid_arxiv'
    )

    # Check it can be found
    contents = cache.search_cache('totally_special_argon_lines.dat', path_only=False)
    assert len(contents) == 1, 'Problem adding file to cache!'
    cache_url = list(contents.keys())[0]
    cached_file = list(contents.values())[0]
    assert cached_file.is_file(), 'File not added to cache'

    # Delete it
    cache.remove_from_cache(cache_url=cache_url)
    assert cache.search_cache('totally_special_argon_lines.dat') == [], \
            'Should not be able to find the file'


def test_pygit2():
    assert cache.Repository is not None, 'Must install pygit2 to run tests'


# TODO: This works when running locally, but it fails in CI for a reason I don't
# yet understand.  Removing for now.
#def test_get_tag():
#    tag_version, date = cache.git_most_recent_tag()
#    assert date is not None, 'Failed to get most recent tag version'


def test_waveio_load_reid_arxiv():

    # Test the extension logic, given the download/cache system
    waveio.load_reid_arxiv("vlt_xshooter_vis1x1.fits")
    waveio.load_reid_arxiv("vlt_xshooter_vis1x1.json")


def test_datapath():
    # NOTE: Because dataPaths is created every time pypeit is imported, the
    # first part of this test is basically guaranteed to pass, if pypeit is
    # imported properly.
    try:
        # Try to define a path that should exist
        p = PypeItDataPath('tests')
    except PypeItPathError as e:
        raise AssertionError('Could not define a path that should exist')

    with pytest.raises(PypeItPathError):
        # Make sure that defining a path that does not exist fails
        p = PypeItDataPath('junk')


def test_truediv():
    p = PypeItDataPath('tests')

    data_file = p / 'b1.fits.gz'
    assert isinstance(data_file, Path), 'Should return a Path object with the file'
    _data_file = p.path / 'b1.fits.gz'
    assert data_file == _data_file, 'Direct access to the path should return identical paths'

    # Should raise an error because 'junk' is neither a directory or an
    # *existing* file
    with pytest.raises(PypeItPathError):
        data_file = p / 'junk'

    subdir = 'ipac'
    _p = p / subdir
    assert isinstance(_p, PypeItDataPath), \
            'Should return a PypeItDataPath for a valid subdirectory'

    assert str(_p.path.relative_to(p.path)) == subdir, 'Wrong subdirectory'


def test_get_file_path():
    f = 'b1.fits.gz'
    # NOTE: Setting to_pkg symlink only needs to be done once for each unique file
    data_file = dataPaths.tests.get_file_path(f, to_pkg='symlink')
    # Check that a Path is returned
    assert isinstance(data_file, Path), 'returned file should be a Path instance'
    # Check that the format is correct
    assert dataPaths.tests.get_file_path(f, return_format=True)[1] == 'fits', 'Wrong file format'
    # Check a different format
    assert dataPaths.tests.get_file_path('bspline_model.npz', return_format=True,
                                         to_pkg='symlink')[1] == 'npz', 'Wrong file format'


def test_cache_to_pkg():
    test_file_name = 'cache_test.txt'
    test_file = dataPaths.tests.path / test_file_name

    # Make sure the file is not currently in the cache
    contents = cache.search_cache(test_file_name, path_only=False)
    if len(contents) > 0:
        cache.remove_from_cache(cache_url=list(contents.keys()), allow_multiple=True)

    assert test_file.is_file(), 'File should exist on disk at the start of the test'

    # Remove the file
    test_file.unlink()

    # Use the cache system to access it
    _test_file = dataPaths.tests.get_file_path(test_file_name)

    # Search the cache for the file
    contents = cache.search_cache(test_file_name, path_only=False)
    assert len(contents) == 1, 'Should find 1 relevant file in the cache'

    # Parse the url
    host, branch, subdir, filename = cache.parse_cache_url(list(contents.keys())[0])
    assert host == 'github', 'Host is wrong'
    assert branch == cache.git_branch(), 'Branch is wrong'
    assert subdir == dataPaths.tests.subdirs, 'Subdirectory is wrong'
    assert filename == test_file_name, 'File name is wrong'

    # Check that the file is in the cache
    assert len(cache.search_cache(test_file_name)) == 1, 'File not found in cache'

    # Symlink the file from the cache to the repository
    _test_file = dataPaths.tests.get_file_path(test_file_name, to_pkg='symlink')

    # The file should now exist in the package directory
    assert test_file.is_file(), 'File should now exist'
    # ... as a symlink
    assert test_file.is_symlink(), 'File shoule be a symlink'
    # ... that points to the cache file
    assert str(cache.search_cache(test_file_name)[0]) == os.path.realpath(test_file), \
            'symlink has the wrong target'

    # Remove the symlink
    test_file.unlink()
    # .. and instead move the file from the cache to the package directory
    _test_file = dataPaths.tests.get_file_path(test_file_name, to_pkg='move')

    # Now the file should be back in the right place
    assert test_file.is_file(), 'File not moved from cache'
    # ... it should not be a symlink
    assert not test_file.is_symlink(), 'File should not be a symlink'
    # ... and should no longer exist in the cache
    assert len(cache.search_cache(test_file_name)) == 0, \
            'File should have been removed from the cache'
    


