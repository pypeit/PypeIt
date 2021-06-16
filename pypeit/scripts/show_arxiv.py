"""
Wrapper to matplotlib to show an archived arc spectrum

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ShowArxiv(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show an archived arc spectrum located in '
                                                'pypeit/data/arc_liens/reid_arxiv',
                                    width=width)
        parser.add_argument("file", type=str, help="Arxiv filename, e.g. gemini_gmos_r831_ham.fits")
        parser.add_argument('--det', default=1, type=int, help='Detector number')
        return parser

    @staticmethod
    def main(args):
        """ Shows the spectrum
        """
        import os
        from pkg_resources import resource_filename

        from matplotlib import pyplot as plt
        from pypeit.core.wavecal import waveio

        # Path
        if os.path.basename(args.file) == args.file:
            args.file = os.path.join(resource_filename('pypeit', 'data'),
                                     'arc_lines', 'reid_arxiv', args.file)

        wave, flux, binspec = waveio.load_template(args.file, args.det)

        plt.clf()
        ax = plt.gca()
        ax.plot(wave, flux)
        plt.show()


