"""
This script enables the user to convert a MasterWaveCalib wavelength solution fits file
into a PypeIt arxiv solution that can be used with the full_template method.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
from pypeit import msgs
from pypeit import par
from pypeit import inputfiles
from pypeit import utils
from pypeit.scripts import scriptbase


class ArxivSolution(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Read in a MasterWaveCalib solution and convert it into the '
                                                'format required for the PypeIt full template archive', width=width)
        parser.add_argument('file', type = str, default=None, help='MasterWaveCalib file')
        parser.add_argument('binning', type=int, help="Spectral binning")
        parser.add_argument('-s', '--slit', default=0, type=int, help='Slit number to use')
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename make_arxiv_solution_YYYYMMDD-HHMM.log')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):
        import os
        from pypeit.wavecalib import WaveCalib
        from pypeit.core.wavecal import wvutils

        chk_version = not args.try_old

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('arxiv_solution', args.verbosity)

        # Check that a file has been provided
        if args.file is None:
            msgs.error('You must input a MasterWaveCalib file')
        elif not os.path.exists(args.file):
            msgs.error("The following MasterWaveCalib file does not exist:" + msgs.newline() + args.file)

        # Load the wavelength calibration file
        wv_calib = WaveCalib.from_file(args.file, chk_version=chk_version)
        # Check if a wavelength solution exists
        if wv_calib['wv_fits'][args.slit]['wave_soln'] is None:
            gd_slits = []
            for slit in range(len(wv_calib['wv_fits'])):
                if wv_calib['wv_fits'][slit]['wave_soln'] is not None:
                    gd_slits.append(f"{slit}")
            # Prepare the message
            thismsg = f"A wavelength solution does not exist for slit {args.slit}. "
            if len(gd_slits) == 0:
                thismsg += "There are no good slits - the WaveCalib file is bad."
            else:
                thismsg += "Try one of the following slits, instead: " + msgs.newline() + ", ".join(gd_slits)
            msgs.error(thismsg)
        wave = wv_calib['wv_fits'][args.slit]['wave_soln'].flatten()
        spec = wv_calib['wv_fits'][args.slit]['spec'].flatten()
        outname = args.file.replace(".fits", "_arXiv.fits")
        wvutils.write_template(wave, spec, args.binning, './', outname, cache=True)
