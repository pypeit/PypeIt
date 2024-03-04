
from IPython import embed

from pypeit import dataPaths
from pypeit import scripts
from pypeit import cache


def test_install_telluric():
    scripts.install_telluric.InstallTelluric.main(
        scripts.install_telluric.InstallTelluric.parse_args(
            ['TellPCA_3000_26000_R25000.fits', '--force_update']
        )
    )

def test_install_extinctfile():

    test_file_name = 'cache_test.txt'
    test_file = dataPaths.tests.path / test_file_name

    # Make sure the file is not currently in the cache
    contents = cache.search_cache(test_file_name, path_only=False)
    if len(contents) > 0:
        cache.remove_from_cache(cache_url=list(contents.keys()), allow_multiple=True)

    assert test_file.is_file(), 'File should exist on disk at the start of the test'

    # Use the install script to install the file as an extinction file
    scripts.install_extinctfile.InstallExtinctfile.main(
        scripts.install_extinctfile.InstallExtinctfile.parse_args([str(test_file)])
    )

    # Search the cache for the file
    contents = cache.search_cache(test_file_name, path_only=False)
    assert len(contents) == 1, 'Should find 1 relevant file in the cache'
    cache_url = list(contents.keys())[0]

    # Parse the url
    host, branch, subdir, filename = cache.parse_cache_url(cache_url)
    assert host == 'github', 'Host is wrong'
    assert branch == cache.git_branch(), 'Branch is wrong'
    assert subdir == dataPaths.extinction.subdirs, 'Subdirectory is wrong'
    assert filename == test_file_name, 'File name is wrong'

    # Delete it from the cache
    cache.remove_from_cache(cache_url=cache_url)


def test_install_linelist():

    test_file_name = 'cache_test.txt'
    test_file = dataPaths.tests.path / test_file_name

    # Make sure the file is not currently in the cache
    contents = cache.search_cache(test_file_name, path_only=False)
    if len(contents) > 0:
        cache.remove_from_cache(cache_url=list(contents.keys()), allow_multiple=True)

    assert test_file.is_file(), 'File should exist on disk at the start of the test'

    # Use the install script to install the file as a line list
    scripts.install_linelist.InstallLinelist.main(
        scripts.install_linelist.InstallLinelist.parse_args([str(test_file)])
    )

    # Search the cache for the file
    contents = cache.search_cache(test_file_name, path_only=False)
    assert len(contents) == 1, 'Should find 1 relevant file in the cache'
    cache_url = list(contents.keys())[0]

    # Parse the url
    host, branch, subdir, filename = cache.parse_cache_url(cache_url)
    assert host == 'github', 'Host is wrong'
    assert branch == cache.git_branch(), 'Branch is wrong'
    assert subdir == dataPaths.linelist.subdirs, 'Subdirectory is wrong'
    assert filename == test_file_name, 'File name is wrong'

    # Delete it from the cache
    cache.remove_from_cache(cache_url=cache_url)

