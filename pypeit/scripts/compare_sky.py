"""
Plots an extracted sky spectrum with an archived one.  Probably most useful for
exploring sky spectra in the blue

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import argparse
from pypeit.scripts import scriptbase

from IPython import embed

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
                            help=argparse.SUPPRESS)
        return parser

    # Script to run XSpec from the command line or ipython
    @classmethod
    def main(cls, args):

        import matplotlib.pyplot as plt

        from pypeit import specobjs
        from pypeit.core import skyspec

        # Initialize the log
        cls.init_log(args)

        # Extension
        exten = args.exten if args.exten is not None else 1

        # Load user file
        user_sobjs = specobjs.SpecObjs.from_fitsfile(args.file)
        user_sobj = user_sobjs[exten-1]

        # Load sky spec
        arx_sky = skyspec.load_sky_spectrum(args.skyfile)

        # Plot
        plt.clf()
        if args.optimal:
            usr_wave = user_sobj.OPT_WAVE 
            usr_flux = user_sobj.OPT_COUNTS_SKY
        else:
            usr_wave = user_sobj.BOX_WAVE 
            usr_flux = user_sobj.BOX_COUNTS_SKY
        plt.plot(usr_wave, usr_flux*args.scale_user, 'k-', label='user')
        plt.plot(arx_sky.wave, arx_sky.flux, 'b-', label='archive')
        legend = plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                            handletextpad=0.3, fontsize='small', numpoints=1)
        if not args.test:
            plt.show()


