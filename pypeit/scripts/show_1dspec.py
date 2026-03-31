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
        parser.add_argument("file", type=str, help="Spectral file")
        parser.add_argument("--list", default=False, help="List the extensions only?",
                            action="store_true")
        parser.add_argument("--exten", type=int, default=1, help="FITS extension")
        parser.add_argument("--obj", type=str,
                            help="Object name in lieu of extension, e.g. SPAT0424-SLIT0000-DET01")
        parser.add_argument("--extract", type=str, default='OPT',
                            help="Extraction method. Default is OPT. ['BOX', 'OPT']")
        parser.add_argument("--flux", default=False, action="store_true",
                            help="Show fluxed spectrum?")
        parser.add_argument('-m', '--unmasked', dest='masked', default=True, action='store_false',
                            help='Only show unmasked data.')
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

        from pypeit import PypeItError
        from pypeit import specobjs
        from pypeit.display.display import show_1dspec

        # Initialize the log
        cls.init_log(args)

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
            exten = np.where(sobjs.NAME == args.obj)[0][0]
            if exten < 0:
                raise PypeItError(f"Bad input object name: {args.obj}")
        else:
            exten = args.exten-1 # 1-index in FITS file

        # Check Extraction
        if args.extract == 'OPT':
            if sobjs[exten]['OPT_WAVE'] is None: #not in sobjs[exten]._data.keys():
                    raise PypeItError("Spectrum not extracted with OPT.  Try --extract BOX")

        # Pre-pend cwd to filename (in case RC Ginga was launched already)
        full_file = Path(args.file).absolute()

        # Ginga
        show_1dspec(
            str(full_file), ext=exten, masked=args.masked, extraction=args.extract,
            fluxed=args.flux
        )
