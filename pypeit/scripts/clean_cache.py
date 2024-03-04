"""
Script to clean cache of specified files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase

class CleanCache(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='View/Remove fils in the PypeIt data cache',
                                    width=width)
        parser.add_argument('-p', '--pattern', type=str, nargs='+',
                            help='Remove any files matching the provided pattern.  If combined '
                                 'with --version, this selects only files downloaded from the '
                                 'identified GitHub versoin.  If the version is not specified, '
                                 'any file matching the provided pattern(s) are removed.')
        parser.add_argument('-v', '--version', type=str, nargs='+',
                            help='Remove files associated one or more provided tags, branches, '
                                 'or commit references on GitHub.  These must be an exact match '
                                 'to the relevant GitHub reference.  If combined with --pattern, '
                                 'this selects the GitHub reference for the files found.  If no '
                                 'files are specified, all files associated with the given '
                                 'reference are removed.  Note this is only relevant for the '
                                 'files on GitHub, not s3.  For files on s3, do not specify the '
                                 'version.')
        parser.add_argument('--remove_all', default=False, action='store_true',
                            help='BEWARE: Removes all data from the pypeit cache.  Use of this '
                                 'option ignores the --pattern and --version options.')
        parser.add_argument('-l', '--list', default=False, action='store_true',
                            help='Only list the contents of the cache.')

        return parser

    @staticmethod
    def main(args):
        from IPython import embed
        import numpy as np
        import astropy.utils.data

        from pypeit import msgs
        from pypeit import cache

        if args.list:
            # Print the full contents
            contents = cache.search_cache(None, path_only=False)
            print(f' {"HOST":>10} {"BRANCH":>20} {"SUBDIR":>20} {"FILE":<30}')
            for url in contents.keys():
                head, branch, subdir, f = cache.parse_cache_url(url)
                print(f' {head:>10} {"..." if branch is None else branch:>20}'
                      f' {subdir:>20} {f:<30}')
            return

        if args.pattern is None and args.version is None and not args.remove_all:
            msgs.error('Arguments provided not sufficient to find files for deletion.')

        if args.remove_all:
            # Removes the entire cache
            astropy.utils.data.clear_download_cache(pkgname='pypeit')
            return
        
        # Get *all* of the contents of the cache
        if args.pattern is None:
            contents = cache.search_cache(None, path_only=False)
        else:
            contents = {}
            for p in args.pattern:
                contents.update(cache.search_cache(pattern=p, path_only=False))

        # TODO: For symlinked files, is there a way to follow the symlinks?  Or
        # should we search for broken symlinks in the package directory
        # structure after the cache contents are removed?

        # For now, we only need the urls.
        contents = list(contents.keys())

        # If versions are set, down select to files on github *and* in the selected versions
        if args.version is not None:
            versions = np.array([cache.parse_cache_url(c)[1] for c in contents])
            contents = np.array(contents)[np.isin(versions, args.version)].tolist()

        if len(contents) == 0:
            msgs.warn('No files to remove.')
            return

        # Report
        msgs.info('Removing the following files from the cache:')
        for c in contents:
            msgs.info(f'    {c}')

        # Remove the selected contents.  cache_url argument must be a list
        cache.remove_from_cache(cache_url=contents, allow_multiple=True)


