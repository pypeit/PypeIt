"""
This script displays the flat images in an RC Ginga window.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase

class ChkScattLight(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display the scattered light image in a Ginga viewer',
                                    width=width)
        parser.add_argument('file', type=str,
                            help='PypeIt Scattered Light file [e.g. ScatteredLight_A_0_DET01.fits.gz]')
        parser.add_argument('slits', type=str,
                            help='Slits calibration file [e.g. Slits_A_0_DET01.fits.gz]')
        parser.add_argument('--spec2d', type=str, default=None,
                            help='PypeIt science spec2d file')
        parser.add_argument('--det', default='1', type=str,
                            help='Detector name or number.  If a number, the name is constructed '
                                 'assuming the reduction is for a single detector.  If a string, '
                                 'it must match the name of the detector object (e.g., DET01 for '
                                 'a detector, MSC01 for a mosaic).')
        parser.add_argument('--mask', type=bool, default=False,
                            help='If True, the detector pixels that are considered on the slit will be '
                                 'masked to highlight the scattered light regions')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        from pypeit import scattlight, spec2dobj, slittrace
        from pypeit import msgs
        from pypeit.pypmsgs import PypeItError, PypeItDataModelError
        from pypeit.images.detector_container import DetectorContainer
        from pypeit import io

        chk_version = not args.try_old

        # Parse the detector name
        try:
            det = int(args.det)
        except:
            detname = args.det
        else:
            detname = DetectorContainer.get_name(det)

        # Load scattered light calibration frame
        ScattLightImage = scattlight.ScatteredLight.from_file(args.file, chk_version=chk_version)

        # Load slits information
        slits = slittrace.SlitTraceSet.from_file(args.slits, chk_version=chk_version)

        # Load the alternate file if requested
        display_frame = None  # The default is to display the frame used to calculate the scattered light model
        if args.spec2d is not None:
            msgs.error("displaying the spec2d scattered light is not currently supported")
            try:
                # TODO :: the spec2d file may have already had the scattered light removed, so this is not correct. This script only works when the scattered light is turned off for the spec2d file
                spec2D = spec2dobj.Spec2DObj.from_file(args.spec2d, detname,
                                                       chk_version=chk_version)
            except PypeItDataModelError:
                msgs.warn(f"Error loading spec2d file {args.spec2d} - attempting to load science image from fits")
                spec2D = None

            # Now set the frame to be displayed
            display_frame = io.fits_open(args.spec2d)[detname+'-SCIIMG'].data if spec2D is None else spec2D.sciimg

        # Show
        ScattLightImage.show(image=display_frame, slits=slits, mask=args.mask)
