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
        parser = super().get_parser(description='Display Tilts image and traced tilts in Ginga viewer',
                                    width=width)

        parser.add_argument('file', type=str, default=None,
                            help='PypeIt Tilts file [e.g. Tilt_A_1_01.fits]')
        parser.add_argument('--tilts_image', type=str, default=None,
                            help='PypeIt Tiltimg file [e.g. Tiltimg_A_1_01.fits]. If not provided,'
                                 'PypeIt will attempt to load the file from the same directory as the Tilts file.')
        parser.add_argument('--mpl', default=False, action='store_true',
                            help='Use a matplotlib window instead of ginga to show the tilts')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(pargs):
        from pypeit import wavetilts

        # Load
        tilts = wavetilts.WaveTilts.from_file(pargs.file, chk_version=(not pargs.try_old))
        tilts.show(tilt_img=pargs.tilts_image, in_ginga=np.logical_not(pargs.mpl))
        # embed()



