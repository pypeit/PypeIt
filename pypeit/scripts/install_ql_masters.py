"""
Script to install quick-look master files into the user's pypeit installation.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class InstallQLMasters(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import os

        parser = super().get_parser(description='Script to install PypeIt QL Master files',
                                    width=width)
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--zip', type=str, default=None,
                           help='Zip file of the full QL_MASTERS directory downloaded from the '
                                'PypeIt Google Drive')
        group.add_argument('--ql_path', type=str, default=None,
                           help='An existing directory to symlink as the QL_MASTERS directory.')
        parser.add_argument('--odir', type=str, default=os.getcwd(),
                            help='The directory in which to extract the zip file.  Ignored if a '
                                 'direct path is provided using --ql_path.')
        parser.add_argument('--rmzip', default=False, action='store_true',
                            help='Remove the downloaded zip file')
        return parser

    @staticmethod
    def main(args):
        import os
        import zipfile
        from pkg_resources import resource_filename

        from pypeit.io import create_symlink

        # Check that either the zip file or the directory is provided
        if args.zip is None and args.ql_path is None:
            raise ValueError('Must provide either the zip file or the path to an existing '
                             'QL_MASTERS directory.')

        if args.zip is None:
            # Directory should already exist.  Check that
            ql_dir = os.path.abspath(args.ql_path)
            if not os.path.isdir(ql_dir):
                raise NotADirectoryError(f'{args.ql_path} does not exist!')
        else:
            # Check that the zip file exists
            if not os.path.isfile(args.zip):
                raise FileExistsError(f'{args.zip} does not exist!')

            # Check the output path
            _odir = os.path.abspath(args.odir)
            if not os.path.isdir(_odir):
                os.makedirs(_odir)

            # Unzip the file
            with zipfile.ZipFile(args.zip, 'r') as zip_ref:
                zip_ref.extractall()

            # Check that it created the expected directory
            ql_dir = os.path.join(_odir, 'QL_MASTERS')
            if not os.path.isdir(ql_dir):
                raise NotADirectoryError('Zip file should have created a QL_MASTERS directory in '
                                        f'{_odir}, but that directory does not exist.  Check '
                                        'and/or re-download the zip file.')

        # Get the pypeit data directory and make sure it exists
        data_dir = resource_filename('pypeit', 'data')
        if not os.path.isdir(data_dir):
            raise NotADirectoryError(f'Unable to find {data_dir}.  Check your installation.')

        # Create a symlink to the QL_MASTERS directory in the pypeit/data
        # directory.
        create_symlink(ql_dir, data_dir, overwrite=True)

        # Remove the zip file if the user requested
        if args.rmzip:
            print(f'Removing: {args.zip}')
            os.remove(args.zip)


