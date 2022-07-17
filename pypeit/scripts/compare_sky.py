"""
Plots an extracted sky spectrum with an archived one.  Probably most useful for
exploring sky spectra in the blue

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class CompareSky(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Compare the extracted sky spectrum against an '
                                                'archived sky model maintained by PypeIt.',
                                    width=width)
        parser.add_argument('file', type=str, help='spec1d Spectral file')
        parser.add_argument('skyfile', type=str,
                            help='Archived PypeIt sky file (e.g. paranal_sky.fits)')
        parser.add_argument('--exten', type=int, help='FITS extension')
        parser.add_argument('--optimal', default=False, action='store_true',
                            help='Show Optimal? Default is boxcar')
        parser.add_argument('--scale_user', default=1., type=float,
                            help='Scale user spectrum by a factor')
        parser.add_argument('--test', default=False, action='store_true',
                            help='Load files but do not show plot')
        return parser

    # Script to run XSpec from the command line or ipython
    @staticmethod
    def main(args):

        import os

        from matplotlib import pyplot as plt

        from linetools.spectra.io import readspec
        from pypeit import data


        # Extension
        exten = args.exten if args.exten is not None else 1

        # Read spec keywords
        ikwargs = {}
        if args.optimal:
            ikwargs['wave_tag'] = 'OPT_WAVE'
            ikwargs['flux_tag'] = 'OPT_COUNTS_SKY'
            ikwargs['ivar_tag'] = 'OPT_COUNTS_IVAR'
            ikwargs['sig_tag'] = 'OPT_COUNTS_SIG'
        else:
            ikwargs['wave_tag'] = 'BOX_WAVE'
            ikwargs['flux_tag'] = 'BOX_COUNTS_SKY'
            ikwargs['ivar_tag'] = 'BOX_COUNTS_IVAR'
            ikwargs['sig_tag'] = 'BOX_COUNTS_SIG'

        # Load user file
        user_sky = readspec(args.file, exten=exten, **ikwargs)
        # Load sky spec
        arx_sky = data.load_sky_spectrum(args.skyfile)

        # Plot
        plt.clf()
        plt.plot(user_sky.wavelength, user_sky.flux*args.scale_user, 'k-', label='user')
        plt.plot(arx_sky.wavelength, arx_sky.flux, 'b-', label='archive')
        legend = plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                            handletextpad=0.3, fontsize='small', numpoints=1)
        if not args.test:
            plt.show()


