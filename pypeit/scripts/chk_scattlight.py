"""
This script displays the scattered-light data, model and residuals in a ginga
viewer.

The data are read and displayed inside the viewer by the ``PypeItRemote``
plugin; see :mod:`pypeit.display.remote`.

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
        scriptbase.add_ginga_timeout_argument(parser)
        return parser

    @classmethod
    def main(cls, args):

        from pypeit.display import remote

        # Initialize the log
        cls.init_log(args)

        # Read and show the scattered light in ginga
        remote.send_command(
            'chk_scattlight', remote.args_to_json(args, path_args=['file', 'slits', 'spec2d']),
            timeout=args.timeout
        )
