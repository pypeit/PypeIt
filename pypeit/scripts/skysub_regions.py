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
        parser.add_argument('file', type=str, default=None, help='spec2d file')
        parser.add_argument('--det', default='1', type=str, help="Detector")
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
        from IPython import embed
        from pypeit import spec2dobj
        import os
        import astropy.io.fits as fits
        from pypeit import msgs
        from pypeit import io
        from pypeit.core.gui.skysub_regions import SkySubGUI
        from pypeit.images import buildimage
        from pypeit.images.detector_container import DetectorContainer
        from pypeit.edgetrace import EdgeTraceSet

        # Parse the detector name
        try:
            det = int(args.det)
        except:
            detname = args.det
        else:
            detname = DetectorContainer.get_name(det)

        # Load it up
        spec2DObj = spec2dobj.Spec2DObj.from_file(args.file, detname, chk_version=True)
        frame = spec2DObj.sciimg
        hdr = fits.open(args.file)[0].header
        fname = hdr['FILENAME']
        calib_dir = hdr['CALIBDIR']
        pypeline = hdr['PYPELINE']
        specname = hdr['PYP_SPEC']

        # Use the edges calibration frame to set the calibration key
        key = EdgeTraceSet.calib_type.upper()
        if key not in spec2DObj.calibs:
            # TODO: Until I can figure out a better approach...
            msgs.error(f'EdgeTrace calibration frame not recorded in {args.file}!')
        calib_key, _ = EdgeTraceSet.parse_key_dir(spec2DObj.calibs[key], from_filename=True)

        # Use the appropriate class to get the "detector" number
        det = spec2DObj.detector.parse_name(detname)

        # Setup for PypeIt imports
        msgs.reset(verbosity=args.verbosity)

        # Grab the slit edges
        slits = spec2DObj.slits

        # Get the spatial flexure
        spat_flexure = None
        if args.flexure:
            spat_flexure = spec2DObj.sci_spat_flexure

        # Derive an appropriate output filename
        file_base = os.path.basename(fname)
        regfile = buildimage.SkyRegions.construct_file_name(calib_key, calib_dir=calib_dir,
                                                            basename=io.remove_suffix(file_base))

        # Finally, initialise the GUI
        skyreg = SkySubGUI.initialize(det, frame, slits, pypeline, specname, outname=regfile,
                                      overwrite=args.overwrite, runtime=False, printout=True,
                                      initial=args.initial, flexure=spat_flexure)

        # Get the results
        skyreg.get_result()

        # Reset the defaults
        skyreg.finalize()
