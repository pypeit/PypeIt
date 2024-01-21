"""
This script allows the user to read a spec3D FITS file (DataCube)
from IFU instruments, and extract a 1D spectrum of the brightest
object. This script is primarily used to extract a spectrum of a
point source from a DataCube, and save it as a spec1d file. A
common usage is to extract a spectrum of a standard star from a
DataCube, and use it to flux calibrate the science DataCubes.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
from pypeit import msgs
from pypeit import par
from pypeit import inputfiles
from pypeit import utils
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts import scriptbase
from pypeit.coadd3d import DataCube


class ExtractDataCube(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Read in a datacube, extract a spectrum of a point source,'
                                                'and save it as a spec1d file.', width=width)
        parser.add_argument('file', type = str, default=None, help='spec3d.fits DataCube file')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename extract_datacube_YYYYMMDD-HHMM.log')
        return parser

    @staticmethod
    def main(args):
        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('extract_datacube', args.verbosity)

        # Check that a file has been provided
        if args.file is None:
            msgs.error('You must input a spec3d (i.e. PypeIt DataCube) fits file')

        # Read in the relevant information from the .extract file
        ext3dfile = inputfiles.ExtractFile.from_file(args.file)
        spectrograph = load_spectrograph(ext3dfile.config['rdx']['spectrograph'])

        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                              merge_with=(ext3dfile.cfg_lines,))

        # Load the DataCube
        tstart = time.time()
        extcube = DataCube.from_file(args.file)

        # Extract the spectrum
        extcube.extract_spec(parset['reduce']['extraction'], overwrite=args.overwrite)

        # Report the extraction time
        msgs.info(utils.get_time_string(time.time()-tstart))
