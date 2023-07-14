"""
Wrapper to matplotlib to show an arc spectrum

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import numpy as np
from matplotlib import pyplot as plt

from pypeit.scripts import scriptbase

from pypeit import wavecalib
from pypeit import slittrace

from IPython import embed

class ShowWvCalib(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show the result of wavelength calibration',
                                    width=width)
        parser.add_argument("file", type=str, help="WaveCalib file")
        parser.add_argument("slit_order", type=int, help="Slit or Order number")
        parser.add_argument("--slit_file", type=str, help="Slit file")
        parser.add_argument("--is_order", default=False, action="store_true",
                            help="Input slit/order is an order")
        return parser

    @staticmethod
    def main(pargs, unit_test=False):
        """ Shows the spectrum
        """

        from matplotlib import pyplot as plt

        # Load
        wvcalib = wavecalib.WaveCalib.from_file(pargs.file)
        if pargs.slit_file is not None:
            slits = slittrace.SlitTraceSet.from_file(pargs.slit_file)

        # Parse
        if pargs.is_order:
            idx = np.where(slits.ech_order == pargs.slit_order)[0][0]
        else:
            idx = np.where(wvcalib.spat_ids == pargs.slit_order)[0][0]

        # Grab it
        spec = wvcalib.arc_spectra[:,idx]
        nspec = len(spec)

        # Generate wavelengths
        wave = wvcalib.wv_fits[idx].wave_soln

        # Plot
        plt.clf()
        ax = plt.gca()
        ax.plot(wave, spec)
        plt.show()


