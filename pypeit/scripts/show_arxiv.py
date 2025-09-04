"""
Wrapper to matplotlib to show an archived arc spectrum

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import argparse

from pypeit.scripts import scriptbase
from pypeit import data


class ShowArxiv(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show an archived arc spectrum located in '
                                                'pypeit/data/arc_liens/reid_arxiv',
                                    width=width)
        parser.add_argument("file", type=str, help="Arxiv filename, e.g. gemini_gmos_r831_ham.fits")
        parser.add_argument('--det', default=1, type=int, help='Detector number')
        parser.add_argument('--test', default=False, action='store_true',
                            help=argparse.SUPPRESS)
        return parser

    @staticmethod
    def main(args):
        """ Shows the spectrum
        """
        from matplotlib import pyplot as plt
        from pypeit.core.wavecal import waveio

        # NOTE: Path is checked within load_template()
        wave, flux = waveio.load_template(args.file, args.det)[:2]

        # We're done if this is a test
        if args.test:
            return

        # Show
        plt.clf()
        ax = plt.gca()
        ax.plot(wave, flux)
        plt.show()

