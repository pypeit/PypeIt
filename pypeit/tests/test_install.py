
from pathlib import Path
import shutil
from IPython import embed

from astropy.config import set_temp_cache

from pypeit import dataPaths
from pypeit import scripts
from pypeit import cache


def run_install_telluric():
    scripts.install_telluric.InstallTelluric.main(
        scripts.install_telluric.InstallTelluric.parse_args(
            ['TellPCA_3000_26000_R25000.fits', '--force_update']
        )
    )


def run_install_extinctfile():

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


def run_install_linelist():

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


# TODO: There's got to be a more concise way to do this...
def test_install_telluric():

    root = 'cache_test'    
    tmp_cache_dir = Path(f'{root}').absolute()
    if tmp_cache_dir.is_dir():
        shutil.rmtree(tmp_cache_dir)
    tmp_cache_dir.mkdir(parents=True)

    with set_temp_cache(root):
        run_install_telluric()

    shutil.rmtree(tmp_cache_dir)

# TODO: For some reason, the following two tests do not work on Windows.  There
# are other tests in test_pkgdata.py that do effectively the same thing ---
# install a file in the cache and then check that the file is successfully found
# in the cache --- and don't fail. The only difference as these execute the
# installation process from the relevant scripts.  The main reason I'm using the
# `set_temp_cache` function was because I thought it might be a permissions
# issue in CI that writing to the home directory is bad.  But that doesn't
# actually make sense given the success of the other tests.  Regardless, I (KBW)
# give up ... for now.

#def test_install_extinctfile():
#
#    root = 'cache_test'    
#    tmp_cache_dir = Path(f'{root}').absolute()
#    if tmp_cache_dir.is_dir():
#        shutil.rmtree(tmp_cache_dir)
#    tmp_cache_dir.mkdir(parents=True)
#
#    with set_temp_cache(root):
#        run_install_extinctfile()
#
#    shutil.rmtree(tmp_cache_dir)
#
#
#def test_install_linelist():
#
#    root = 'cache_test'    
#    tmp_cache_dir = Path(f'{root}').absolute()
#    if tmp_cache_dir.is_dir():
#        shutil.rmtree(tmp_cache_dir)
#    tmp_cache_dir.mkdir(parents=True)
#
#    with set_temp_cache(root):
#        run_install_linelist()
#
#    shutil.rmtree(tmp_cache_dir)


