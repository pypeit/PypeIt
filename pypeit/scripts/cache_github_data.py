"""
Script to install files from the github repo into the user's cache.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit import dataPaths
from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs

class CacheGithubData(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        valid_paths = list(dataPaths.github_paths().keys())

        parser = super().get_parser(description='Script to download/cache PypeIt github data',
                                    width=width)
        parser.add_argument('spectrograph', type=str, nargs='+',
                            help='A valid spectrograph identifier: '
                                 f'{", ".join(available_spectrographs)}')
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--exclude', type=str, default=['tests'], nargs='+',
                           help='A subset of the directories to *exclude* from the list of '
                                'files to download. Options are: '
                                f'{", ".join(valid_paths)}.  This option is '
                                'mutually exclusive with --include.')
        group.add_argument('--include', type=str, default=None, nargs='+',
                           help='The directories to *include* in the list of '
                                'files to download.  Use "--include all" to include all '
                                'directories. Options are: '
                                f'{", ".join(["all"] + valid_paths)}.  This '
                                'option is mutually exclusive with --exclude.')
        parser.add_argument('--spec_dependent_only', default=False, action='store_true',
                            help='Only include files that are specific to the provided list of '
                                 'spectrographs.  By default, the script also includes any files '
                                 'in the selected directories that are *not* specific to a given '
                                 'spectrograph (e.g., atmospheric extinction curves).')
        parser.add_argument('--force_update', action='store_true',
                            help='Force re-download of existing files')
        return parser

    @staticmethod
    def main(args):
        import os
        import pathlib
        from IPython import embed

        import numpy as np

        import github

        from pypeit import msgs
        from pypeit import cache
        from pypeit.pypeitdata import PypeItDataPath

        # First check the input spectrograph list
        if any([inst not in available_spectrographs for inst in args.spectrograph]):
            raise ValueError('Provided invalid spectrograph name. Options are: '
                             f'{", ".join(available_spectrographs)}')

        github_paths = dataPaths.github_paths()
        valid_paths = list(github_paths.keys())
        if any([path not in valid_paths for path in args.exclude]):
            raise ValueError('Provided path to exclude not valid.  Options are: '
                             f'{", ".join(valid_paths)}')
        include_paths = ['all'] + valid_paths
        if args.include is not None and any([path not in include_paths for path in args.include]):
            raise ValueError('Provided path to include not valid.  Options are: '
                             f'{", ".join(include_paths)}')
        selected_paths = np.setdiff1d(valid_paths, args.exclude).tolist() \
                            if args.include is None else args.include
        if 'all' in selected_paths:
            selected_paths = valid_paths
        
        # Get the relevant github branch
        branch = cache.git_branch()

        # Access the repo; use a token if one is available
        if os.getenv('GITHUB_TOKEN') is None:
            msgs.warn('GITHUB_TOKEN environmental variable is not defined, meaning script will '
                      'not authenticate a GitHub user via an OAuth token.  Beware of rate limits!')
            auth = None
        else:
            auth = github.Auth.Token(os.getenv('GITHUB_TOKEN'))
        repo = github.Github(auth=auth).get_repo("pypeit/PypeIt")

        # Cycle through all the data paths that use github as their host and
        # collect the contents of the directory
        msgs.info('Searching github repository ... ')
        contents = {}
        for datadir, meta in github_paths.items():
            if datadir not in selected_paths:
                continue
            # Recursively get the directory contents
            contents[meta['path']] \
                = cache.github_contents(repo, branch, f'pypeit/data/{meta["path"]}')
        msgs.info('Searching github repository ... done.')
        
        # Determine which files should be in the cache (or in the repo)
        # TODO: This is currently broken because not all spectrograph-dependent
        # files are named after their *exact* pypeit spectrograph names.  E.g.,
        # the reid_arxiv file for gemini_flamingos2 is
        # ``Flamingos2_HK_HK.fits``, not ``gemini_flamingos2_HK_HK.fits``.  It's
        # a major effort to rename all those files.  The other (better?) option
        # would be to put all these reid_arxiv files into subdirectories named
        # after each spectrograph.
        msgs.info('Parsing which files to include in cache.')
        to_download = {}
        for path, files in contents.items():
            nfiles = len(files)
            to_download[path] = np.zeros(nfiles, dtype=bool)
            for i, f in enumerate(files):
                # NOTE: These check the spectrograph name against the *full
                # path* (e.g.,
                # pypeit/data/arc_lines/reid_arxiv/Flamingos2_HK_HK.fits), not
                # just the file name (e.g., Flamingos2_HK_HK.fits)

                spec_dependent = any([spec in f.path for spec in available_spectrographs])
                to_download[path][i] \
                    = (spec_dependent and any([spec in f.path for spec in args.spectrograph])) \
                        or (not spec_dependent and not args.spec_dependent_only)

        # Use the `get_file_path`` function to find each file locally or pull it
        # into the cache
        msgs.info(f'Number of files to check against package installation/cache:')
        for path in contents.keys():
            msgs.info(f'    {path}: {np.sum(to_download[path])}')
            files = np.array(contents[path])[to_download[path]]
            if len(files) == 0:
                continue
            data_path = PypeItDataPath(path, remote_host="github")
            # NOTE: I'm using POSIX path here because I'm unsure what will
            # happen on Windows if the file is in a subdirectory.
            root = pathlib.PurePosixPath(f'pypeit/data/{path}')
            for f in files:
                # We need the POSIX path relative to the source path.
                rel_path = str(pathlib.PurePosixPath(f.path).relative_to(root))
                data_path.get_file_path(rel_path, force_update=args.force_update, quiet=True)
            
