"""
This script displays the Trace image and the traces in an RC Ginga window.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import numpy as np
from pypeit.scripts import scriptbase
from IPython import embed


class ChkTilts(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display Tiltimg image and 2D fitted tilts in Ginga viewer '
                                                'or Matplotlib window. Tiltimg file must be in the same '
                                                'directory as Tilts.',
                                    width=width)

        parser.add_argument('file', type=str, default=None,
                            help='PypeIt Tilts file [e.g. Tilt_A_1_01.fits]')
        parser.add_argument('--mpl', default=False, action='store_true',
                            help='Use a matplotlib window instead of ginga to show the tilts. Faster plotting.')
        parser.add_argument('--show_traces', default=False, action='store_true',
                            help='Show the traced tilts. This slows down the plotting (mostly in Ginga). If not set, '
                                 'only the fitted, masked and rejected in the fit tilts are shown.')
        parser.add_argument('--try_old', default=True, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(pargs):
        from pypeit import wavetilts

        # Load
        tilts = wavetilts.WaveTilts.from_file(pargs.file, chk_version=(not pargs.try_old))
        tilts.show(in_ginga=np.logical_not(pargs.mpl), show_traces=pargs.show_traces)



