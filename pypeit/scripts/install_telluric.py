"""
Script to install telluric model grids into the user's pypeit installation.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase

class InstallTelluric(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Script to download/install PypeIt telluric files',
                                    width=width)
        parser.add_argument('files', type=str, nargs='+',
                            help='Exact paths to TelFits files to be downloaded '
                                 'from the Cloud and installed in the PypeIt cache')
        parser.add_argument('--force_update', action='store_true',
                            help='Force download of latest version of the telluric grid')
        parser.add_argument('--local_file', action='store_true',
                            help='This is a local file to be installed in the cache')
        return parser

    @staticmethod
    def main(args):
        import os
        from pypeit import cache

        # Loop through the files passed
        for file in args.files:

            if args.local_file:
                # Copy the previously downloaded or power-user-created file to the cache
                cache.write_file_to_cache(file, os.path.basename(file),
                                         'telluric/atm_grids', remote_host="s3_cloud")

            else:
                # Download the file into the cache if not already there (unless force_update)
                cache.fetch_remote_file(file, 'telluric/atm_grids', remote_host="s3_cloud",
                                        install_script=True, force_update=args.force_update)
