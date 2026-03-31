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
        parser.add_argument('--all', default=False, action='store_true',
                            help='By default, the presence of any of the listed patterns yields '
                                 'a match.  This flag requires all patterns to be present for a '
                                 'match.')
        parser.add_argument('--clear', default=False, action='store_true',
                            help='BEWARE: Removes all data from the pypeit cache.  Use of this '
                                 'option ignores the --pattern options.')
        parser.add_argument('-l', '--list', default=False, action='store_true',
                            help='Only list the contents of the cache.')
        return parser

    @staticmethod
    def main(args):
        from IPython import embed
        import astropy.utils.data

        from pypeit import msgs
        from pypeit import cache

        if args.list:
            # Print the full contents
            contents = cache.search_cache(None, path_only=False)
            if len(contents) == 0:
                msgs.info('Cache is empty!')
                return
            cache.list_cache_contents(contents)
            return

        if args.pattern is None and not args.clear:
            msgs.error('Arguments provided not sufficient to find files for deletion.')

        if args.clear:
            # Removes the entire cache
            msgs.info('Clearing the cache!')
            astropy.utils.data.clear_download_cache(pkgname='pypeit')
            return
        
        if args.pattern is None:
            # Get *all* of the contents of the cache
            contents = cache.search_cache(None, path_only=False)
        else:
            # Match cache contents to multiple patterns
            for i, p in enumerate(args.pattern):
                new_contents = cache.search_cache(pattern=p, path_only=False)
                if i == 0:
                    contents = new_contents
                    continue
                if args.all:
                    contents = {k:v for k, v in contents.items() if k in new_contents}
                else:
                    contents.update(new_contents)

        # TODO: For symlinked files, is there a way to follow the symlinks?  Or
        # should we search for broken symlinks in the package directory
        # structure after the cache contents are removed?

        # For now, we only need the urls.
        contents = list(contents.keys())
        if len(contents) == 0:
            msgs.warn('No files to remove.')
            return

        # Report
        msgs.info('Removing the following files from the cache:')
        for c in contents:
            msgs.info(f'    {c}')
        # TODO: Require confirmation?

        # Remove the selected contents.  cache_url argument must be a list
        cache.remove_from_cache(cache_url=contents, allow_multiple=True)


