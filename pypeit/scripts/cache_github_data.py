"""
Script to install telluric model grids into the user's pypeit installation.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase
from pypeit import data
from pypeit import msgs
from pypeit.spectrographs import available_spectrographs
from pypeit import __version__ as pypeit_version

class CacheGithubData(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Script to download/cache PypeIt github data',
                                    width=width)
        parser.add_argument('spectrograph', type=str, nargs='+',
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))
        parser.add_argument('--force_update', action='store_true',
                            help='Force download of GitHub file')
        return parser

    @staticmethod
    def main(args):
        import github

        # First, get the list of reid_arxiv files in GitHub for the present PypeIt version
        # Look in the current `develop` branch if the code is not a tagged release
        tag = "develop" if ".dev" in pypeit_version else pypeit_version
        
        # Use PyGithub to get the download URLs from the repo, with error checking
        repo = github.Github().get_repo("pypeit/PypeIt")
        try:
            arxiv_listing = repo.get_contents(f"pypeit/data/arc_lines/reid_arxiv", tag)
            sensfunc_listing = repo.get_contents(f"pypeit/data/sensfuncs", tag)
            skisim_listing = repo.get_contents(f"pypeit/data/skisim", tag)
        except github.GithubException as err:
            raise ValueError(f"Directory not found in the '{tag}' GitHub tree") from err

        # Loop through the files passed; check against dir_listing and download
        for instrument in args.instrument:

            # Make the list of files to get
            dload_arxiv = [listing for listing in arxiv_listing if instrument in listing.name]
            dload_sensfunc = [listing for listing in sensfunc_listing 
                              if instrument in listing.name]

            # If not matches, warn and continue
            if not dload_arxiv:
                msgs.warn(f"No reid_arxiv files for spectrograph `{instrument}` "
                          f"found on GitHub in tree '{tag}'.")
                continue

            # Loop through found files, using AstroPy's download_file() to cache them
            for file in dload_arxiv:
                data.fetch_remote_file(file.name, "arc_lines/reid_arxiv",
                                       force_update=args.force_update,
                                       full_url=file.download_url)
            for file in dload_sensfunc:
                data.fetch_remote_file(file.name, "sensfuncs",
                                        force_update=args.force_update,
                                        full_url=file.download_url)

        # Download the currently used skisim files
        skisim_files = ['rousselot2000.dat','atm_transmission_secz1.5_1.6mm.dat',
                        'HITRAN.txt','mktrans_zm_10_10.dat']
        dload_skisim = [listing for listing in skisim_listing for skisim in skisim_files 
                        if skisim in listing.name]
        for file in dload_skisim:
            data.fetch_remote_file(file.name, "skisim",
                                   force_update=args.force_update,
                                   full_url=file.download_url)
