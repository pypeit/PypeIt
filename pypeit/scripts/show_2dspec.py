"""
This script displays the images in a spec2d file, with the slit edges and
object traces, in a ginga viewer.

The data are read and displayed inside the ginga viewer by the
``PypeItRemote`` plugin; see :mod:`pypeit.display.remote` and
:func:`pypeit.display.spec2d_data.load_spec2d`.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pypeit.scripts import scriptbase


class Show2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display sky subtracted, spec2d image in a '
                                                'ginga viewer.',
                                    width=width, default_log_file=True)

        parser.add_argument('file', type=str, default=None, help='Path to a PypeIt spec2d file')
        parser.add_argument('--list', default=False, action='store_true',
                            help='List the extensions only?')
        parser.add_argument('--det', default=None, type=str,
                            help='Detector name or number.  If a number, the name is constructed '
                                 'assuming the reduction is for a single detector.  If a string, '
                                 'it must match the name of the detector object (e.g., DET01 for '
                                 'a detector, MSC01 for a mosaic). If not set, the first available '
                                 'detector in the spec2d file will be shown')
        parser.add_argument('--spat_id', type=int, default=None,
                            help='Restrict plotting to this slit (PypeIt ID notation)')
        parser.add_argument('--maskID', type=int, default=None,
                            help='Restrict plotting to this maskID')
        parser.add_argument('--showmask', nargs='*', default=None,
                            help='Include a channel showing the mask.  If no arguments are '
                                 'provided, the mask bit values are provided directly.  You can '
                                 'also specify one or more mask flags used to construct an '
                                 'image identifying which pixels are flagged by any of these '
                                 'issues.  E.g., to show pixels flagged by the instrument '
                                 'specific bad-pixel mask or cosmic arrays, use --showmask BPM CR '
                                 '.  See https://pypeit.readthedocs.io/en/release/out_masks.html '
                                 'for the list of flags.')
        parser.add_argument('--removetrace', default=False, action='store_true',
                            help='Do not overplot traces in the skysub, sky_resid, and resid '
                                 'channels')
        parser.add_argument('--ignore_extract_mask', default=False, action='store_true',
                            help='Ignore the extraction mask')
        parser.add_argument('--show_wave_rms', default=False, action='store_true',
                            help='Show the RMS (in pixels) of the wavelength solution for the '
                                 'slit under the cursor in the ginga readout. Use of this option '
                                 'on many-slit frames may take a long time')
#        parser.add_argument('--sensfunc', type=str, default=None,
#                            help='Pass in a sensfunc to display the sky-subtracted image with a '
#                                 'flux calibration')
        parser.add_argument('--channels', type=str,
                            help='Only show a subset of the channels (0-indexed), e.g. 1,3')
        parser.add_argument('--prefix', type=str, default='',
                            help='Channel name prefix [lets you display more than one set]')
        parser.add_argument('--no_clear', dest='clear', default=True, 
                            action='store_false',
                            help='Do *not* clear all existing tabs')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):

        from pypeit import io
        from pypeit.display import remote

        # List only?
        if args.list:
            io.fits_open(args.file).info()
            return

        # Initialize the log
        cls.init_log(args)

        # Read and show the spec2d file in ginga
        remote.send_command(
            'show_2dspec', remote.args_to_json(args, path_args=['file']), timeout=args.timeout
        )
