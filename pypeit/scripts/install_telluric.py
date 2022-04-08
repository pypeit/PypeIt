"""
Script to install telluric model grids into the user's pypeit installation.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase
from pypeit import data

class InstallTelluric(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Script to download/install PypeIt telluric files',
                                    width=width)
        parser.add_argument('files', type=str, nargs='+',
                            help='Filename(s) of the TelFits files to be downloaded '
                                 'from the Cloud and installed in the PypeIt cache')
        parser.add_argument('--force_update', action='store_true',
                            help='Force download of latest version of the telluric grid')
        return parser

    @staticmethod
    def main(args):

        # Loop through the files passed
        for file in args.files:

            # Download the file into the cache if not already there (unless force_update)
            data.fetch_remote_file(file, 'telgrid', remote_host='s3_cloud',
                                   install_script=True, force_update=args.force_update)
