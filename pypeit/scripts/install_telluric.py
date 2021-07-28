"""
Script to install telluric model grids into the user's pypeit installation.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class InstallTelluric(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import os

        parser = super().get_parser(description='Script to install PypeIt telluric files',
                                    width=width)
        parser.add_argument('--path', type=str, default=os.getcwd(),
                            help='Path to directory with TelFit files downloaded from the PypeIt'
                                 ' Google Drive Telluric/ folder')
        return parser

    @staticmethod
    def main(args):
        import os
        import glob
        from pkg_resources import resource_filename

        from pypeit.io import create_symlink

        # Check the input path
        if not os.path.isdir(args.path):
            raise NotADirectoryError(f'{args.path} is not a directory!')
        _path = os.path.abspath(args.path)

        # Find all the telluric grid files
        tell_files = glob.glob(os.path.join(_path, 'TelFit*'))

        # Check that files were found
        if len(tell_files) == 0:
            raise ValueError(f'No telluric grid files found in {args.path}')

        # Get the directory for the grids and make sure it exists
        grid_dir = resource_filename('pypeit', 'data/telluric/atm_grids/')
        if not os.path.isdir(grid_dir):
            raise NotADirectoryError(f'Unable to find {grid_dir}.  Check your installation.')

        # Create a symlink for each file
        for f in tell_files:
            create_symlink(f, grid_dir, overwrite=True)


