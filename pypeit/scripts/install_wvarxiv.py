"""
Script to install user wavelength templates (arxivs) into the PypeIt cache.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase
from pypeit import cache

class InstallWvArxiv(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Script to install user-created wavelength templates',
                                    width=width)
        parser.add_argument('files', type=str, nargs='+',
                            help='Filename(s) of the template files to be installed '
                                 'in the PypeIt cache')
        return parser

    @staticmethod
    def main(args):
        import numpy as np
        from pypeit import msgs

        # Grab all the files
        files = np.concatenate([sorted(scriptbase.ScriptBase.expandpath(f)) for f in args.files])
        # Remove any that are *not* files (i.e., directories or symlinks)
        files = [f for f in files if f.is_file()]

        # Loop through the files passed
        for f in files:
            if not f.is_file():
                msgs.warn(f'{f} is not a file.')
                continue
            # Copy the user-created file to the cache
            cache.write_file_to_cache(str(f), f.name, 'arc_lines/reid_arxiv')
