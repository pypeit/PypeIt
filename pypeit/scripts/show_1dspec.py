"""
Wrapper for 1D spectrum viewer in ginga.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

from pypeit.scripts import scriptbase


class Show1DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show a 1D spectrum', width=width)
        parser.add_argument(
            'file', type=str,
            help='PypeIt spec1d file, or a file with a single spectrum (OneSpec), e.g. from '
                 'pypeit_coadd_1dspec.  For OneSpec files, the extension, object, and '
                 'extraction options are ignored.'
        )
        parser.add_argument(
            '--list', default=False, action='store_true',
            help='Instead of plotting any spectra, simply list the extensions with spectra'
        )
        grp = parser.add_mutually_exclusive_group()
        grp.add_argument(
            '--exten', type=int, default=1, help='Number of the extension to plot'
        )
        grp.add_argument(
            '--obj', type=str, help='Extension (object) name to plot, e.g. SPAT0424-SLIT0000-DET01'
        )
        parser.add_argument(
            '--extract', type=str, default='OPT', choices=['BOX', 'OPT'],
            help='Method used to extract the spectrum'
        )
        parser.add_argument(
            '--flux', default=False, action='store_true',
            help='Show the flux-calibrated spectrum (if available)'
        )
        parser.add_argument(
            '-m', '--unmasked', dest='masked', default=True, action='store_false',
            help='Only show the unmasked data.'
        )
        parser.add_argument(
            '--chname', type=str, default=None,
            help='Name of the ginga channel (tab).  By default, the name of the file is used, '
                 'and showing the same file again reuses its channel.'
        )
        parser.add_argument(
            '--no_clear', dest='clear', default=True, action='store_false',
            help='If a PypeIt viewer is already open (e.g., from a previous call to this script '
                 'or pypeit_show_2dspec), show the spectrum in a new channel of that viewer and '
                 'return immediately.'
        )
#        parser.add_argument('--jdaviz', default=False, action='store_true',
#                            help='Open the spectrum in jdaviz (requires specutils and jdaviz '
#                                 'to be installed)')
        return parser

    @classmethod
    def main(cls, args):
        """ Runs the XSpecGui on an input file
        """
        from pathlib import Path

        import numpy as np

        from pypeit import log
        from pypeit import PypeItError
        from pypeit import specobjs
        from pypeit import onespec

        # Initialize the log
        cls.init_log(args)

        # A single spectrum (e.g., from pypeit_coadd_1dspec)?
        if onespec.is_onespec_file(args.file):
            spec = onespec.OneSpec.from_file(args.file, chk_version=False)
            if args.list:
                print(f'{args.file} contains a single spectrum (OneSpec): '
                      f'extraction={spec.ext_mode}, fluxed={spec.fluxed}, npix={spec.npix}')
                return
            if args.obj is not None or args.exten != 1:
                log.warning('The file contains a single spectrum; ignoring --exten/--obj.')
            cls.show(args, dict(filename=str(Path(args.file).absolute()), ext=0,
                                masked=args.masked, chname=args.chname))
            return

        sobjs = specobjs.SpecObjs.from_fitsfile(args.file, chk_version=False)

        # List only?
        if args.list:
            print("Showing object names for input file...")
            for ii in range(len(sobjs)):
                line = "EXT{:07d} = {}".format(ii + 1, sobjs[ii].NAME)
                if sobjs[ii].RA is not None:
                    line += " {:0.5f} {:0.5f} {:s}".format(
                        sobjs[ii].RA,
                        sobjs[ii].DEC,
                        sobjs[ii].MASKDEF_OBJNAME)
                if sobjs[ii].MASKDEF_EXTRACT is not None and sobjs[ii].MASKDEF_EXTRACT is True:
                    line += " maskdef_extract"
                if sobjs[ii].hand_extract_flag is True:
                    line += " manual_extract"
                #
                print(line)
            return

        # TODO: Keep this for now, assuming jdaviz ever allows users to
        # instantiate from within a python script.
#        if args.jdaviz:
#            from pypeit.specutils import Spectrum1D, SpectrumList
#            if Spectrum1D is None:
#                raise PypeItError('specutils package must be installed.')
#            try:
#                from jdaviz import Specviz
#            except ModuleNotFoundError:
#                raise PypeItError('jdaviz package must be installed.')
#
#            # First try reading it as a list
#            try:
#                spec = SpectrumList.read(args.file, extract=args.extract, fluxed=args.flux)
#            except:
#                pass
#            else:
#                specviz = Specviz()
#                specviz.load_spectrum(spec)
#                specviz.show()
#                return
#
#            # Maybe it's a OneSpec?
#            try:
#                # TODO: add "grid" as a command-line argument
#                spec = Spectrum1D.read(args.file)
#            except:
#                pass
#            else:
#                specviz = Specviz()
#                specviz.load_spectrum(spec)
#                specviz.show()
#                return
#
#            # If we get here, the file couldn't be parsed
#            raise PypeItError(f'Could not parse input file: {args.file}')


        if args.obj is not None:
            indx = np.where(sobjs.NAME == args.obj)[0]
            if len(indx) == 0:
                raise PypeItError(f"Bad input object name: {args.obj}")
            exten = int(indx[0])
        else:
            exten = args.exten-1 # 1-index in FITS file
        if exten < 0 or exten >= len(sobjs):
            raise PypeItError(f"Bad input extension: {args.exten}")

        # Check Extraction
        if args.extract == 'OPT':
            if sobjs[exten]['OPT_WAVE'] is None: #not in sobjs[exten]._data.keys():
                    raise PypeItError("Spectrum not extracted with OPT.  Try --extract BOX")

        cls.show(args, dict(filename=str(Path(args.file).absolute()), ext=exten,
                            masked=args.masked, extraction=args.extract, fluxed=args.flux,
                            chname=args.chname))

    @classmethod
    def show(cls, args, kwargs):
        """
        Show the spectrum in ginga and wait for the user to close the window.

        If requested (``--no_clear``), the spectrum is instead sent to an open
        PypeIt viewer, if there is one, and the method returns immediately.

        Parameters
        ----------
        args : `argparse.Namespace`_
            The parsed command-line arguments.
        kwargs : :obj:`dict`
            Keyword arguments passed to
            :func:`~pypeit.display.display.show_1dspec`.  Must be JSON
            serializable.
        """
        from pypeit import log
        from pypeit.display import viewer
        from pypeit.display.display import show_1dspec

        # Show the spectrum in an already open viewer?
        if not args.clear:
            if viewer.send_command('spec1d', **kwargs):
                log.info(f'Sent {args.file} to the open PypeIt viewer.')
                return
            log.info('No open PypeIt viewer found; opening a new one.')

        v = viewer.get_viewer()
        # Allow other calls to this script (and pypeit_show_2dspec,
        # pypeit_view_fits) to add spectra to this viewer
        v.listen()
        show_1dspec(**kwargs)

        # Wait for the user to close the window
        v.show(message='close the window to exit')
