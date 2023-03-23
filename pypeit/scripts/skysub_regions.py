"""
This script enables the user to view a 2D FITS file
and define the sky background regions interactively.
Run above the Science/ folder.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class SkySubRegions(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display a Raw science image and interactively '
                                                'define the sky regions using a GUI. Run in the '
                                                'same folder as your .pypeit file',
                                    width=width)
        parser.add_argument('file', type = str, default=None, help='PypeIt file')
        parser.add_argument('--det', default=1, type=int, help="Detector")
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')
        parser.add_argument('-i', '--initial', default=False, action='store_true',
                            help='Use initial slit edges?')
        parser.add_argument('-f', '--flexure', default=False, action='store_true',
                            help='Use flexure corrected slit edges?')
        parser.add_argument('-s', '--standard', default=False, action='store_true',
                            help='List standard stars as well?')
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename skysub_regions_YYYYMMDD-HHMM.log')
        return parser

    @staticmethod
    def main(args):

        import os

        from pypeit import msgs
        from pypeit.core.gui.skysub_regions import SkySubGUI
        from pypeit.core import flexure
        from pypeit.scripts import utils
        from pypeit.images import buildimage

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('skysub_regions', args.verbosity)

        # Generate a utilities class
        info = utils.Utilities(pypeit_file=args.file, det=args.det)

        # Interactively select a science frame
        sciIdx = info.select_science_frame(standard=args.standard)

        # Load the spectrograph and parset
        info.load_par(iFile=sciIdx)
        info.load_metadata()

        # Load the image data
        frame = info.load_frame()

        # Load the slits information
        slits = info.get_slits()
        spat_flexure = flexure.spat_flexure_shift(frame, slits) if args.flexure else None

        # Derive an appropriate output filename
        file_base = info.get_basename()
        prefix = os.path.splitext(file_base)
        outname = os.path.splitext(prefix[0])[0] if prefix[1] == ".gz" else prefix[0]

        info.load_calib_dir()
        calib_key = info.get_calib_key(iFile=iFile)
        regfile = buildimage.SkyRegions.construct_file_name(calib_key, calib_dir=info.calib_dir,
                                                            basename=outname)

        # Finally, initialise the GUI
        skyreg = SkySubGUI.initialize(args.det, frame, slits, info.spectrograph.pypeline,
                                      info.spectrograph.name, outname=regfile,
                                      overwrite=args.overwrite, runtime=False, printout=True,
                                      initial=args.initial, flexure=spat_flexure)

        # Get the results
        skyreg.get_result()

        # Reset the defaults
        skyreg.finalize()
