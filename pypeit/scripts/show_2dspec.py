"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pathlib import Path

from IPython import embed

from pypeit import log
from pypeit import PypeItError
from pypeit import io
from pypeit import utils
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
                                 'a detector, MSC01 for a mosaic). If not set, the first available detector'
                                 'in the spec2d file will be shown')
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
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        parser.add_argument('--ignore_extract_mask', default=False, action='store_true',
                            help='Ignore the extraction mask')
#        parser.add_argument('--sensfunc', type=str, default=None,
#                            help='Pass in a sensfunc to display the sky-subtracted image with a '
#                                 'flux calibration')
        parser.add_argument('--channels', type=str,
                            help='Only show a subset of the channels (0-indexed), e.g. 1,3')
        parser.add_argument('--prefix', type=str, default='',
                            help='Channel name prefix [lets you display more than one set]')
        parser.add_argument('--no_clear', dest='clear', default=True,
                            action='store_false',
                            help='Do *not* clear all existing tabs.  If a PypeIt viewer is '
                                 'already open (e.g., from a previous call to this script), the '
                                 'file is shown in that viewer, in addition to what it already '
                                 'shows, and this script returns immediately.  Use --prefix to '
                                 'avoid replacing channels with the same name.')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @classmethod
    def main(cls, args):

        # List only?
        if args.list:
            io.fits_open(args.file).info()
            return

        # Initialize the log
        cls.init_log(args)

        if not Path(args.file).is_file():
            raise PypeItError(f'{args.file} does not exist!')

        # NOTE: The import of the viewer is deferred so that the --list option
        # does not need to import ginga.
        from pypeit.display import viewer

        # Find the set of channels to show
        show_channels = None if args.channels is None \
                            else [int(item) for item in args.channels.split(',')]

        kwargs = dict(file=str(Path(args.file).absolute()), det=args.det, channels=show_channels,
                      prefix=args.prefix, showmask=args.showmask, removetrace=args.removetrace,
                      ignore_extract_mask=args.ignore_extract_mask, spat_id=args.spat_id,
                      maskID=args.maskID, try_old=args.try_old)

        # Show the file in an already open viewer?
        if not args.clear:
            if viewer.send_command('spec2d', clear=False, **kwargs):
                log.info(f'Sent {args.file} to the open PypeIt viewer.')
                return
            log.info('No open PypeIt viewer found; opening a new one.')

        # Open the viewer and show the file
        v = viewer.get_viewer()
        plugin = v.get_global_plugin('Spec2dView')
        # Allow other calls to this script (and pypeit_view_fits) to add files
        # to this viewer
        v.listen()
        plugin.load_spec2d(clear=True, **kwargs)

        if args.embed:
            viewer.embed(header=utils.embed_header())

        # Wait for the user to close the window
        v.show(message='close the window to exit')
