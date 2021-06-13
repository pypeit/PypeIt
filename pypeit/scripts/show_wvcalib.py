"""
Wrapper to matplotlib to show an arc spectrum

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


# TODO: Should this script be deprecated?
class ShowWvCalib(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show the result of wavelength calibration',
                                    width=width)
        parser.add_argument("file", type=str, help="WaveCalib JSON file")
        parser.add_argument("slit", type=str)
        return parser

    @staticmethod
    def main(pargs, unit_test=False):
        """ Shows the spectrum
        """

        from matplotlib import pyplot as plt
        from linetools import utils as ltu

        wvcalib = ltu.loadjson(pargs.file)

        # Grab it
        spec = wvcalib[pargs.slit]['spec']

        plt.clf()
        ax = plt.gca()
        ax.plot(spec)
        plt.show()


