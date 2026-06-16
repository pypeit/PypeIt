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
from pypeit.scripts import scriptbase


class ExtractDataCube(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Read in a datacube, extract a spectrum of a point source, and save it as '
                        'a spec1d file.',
            width=width, default_log_file=True
        )
        parser.add_argument('file', type = str, default=None, help='spec3d.fits DataCube file')
        parser.add_argument("-e", "--ext_file", type=str,
                            help='Configuration file with extraction parameters')
        parser.add_argument("-s", "--save", type=str,
                            help='Output spec1d filename')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')
        parser.add_argument('-b', '--boxcar_radius', type=float, default=None,
                            help='Radius of the circular boxcar (in arcseconds) to use for the extraction.')
        return parser

    @classmethod
    def main(cls, args):
        import time

        from pypeit import log
        from pypeit import PypeItError
        from pypeit import par
        from pypeit import inputfiles
        from pypeit import utils
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.coadd3d import DataCube

        # Initialize the log
        cls.init_log(args)

        # Check that a file has been provided
        if args.file is None:
            raise PypeItError('You must input a spec3d (i.e. PypeIt DataCube) fits file')
        extcube = DataCube.from_file(args.file)
        spectrograph = load_spectrograph(extcube.PYP_SPEC)

        if args.ext_file is None:
            parset = spectrograph.default_pypeit_par()
        else:
            # Read in the relevant information from the .extract file
            ext3dfile = inputfiles.ExtractFile.from_file(args.ext_file)

            # Parameters
            spectrograph_def_par = spectrograph.default_pypeit_par()
            parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                  merge_with=(ext3dfile.cfg_lines,))

        # Set the boxcar radius
        boxcar_radius = args.boxcar_radius

        # Set the output name
        outname = None if args.save is None else args.save

        # Load the DataCube
        tstart = time.time()

        # Extract the spectrum
        extcube.extract_spec(parset['reduce'], outname=outname, boxcar_radius=boxcar_radius, overwrite=args.overwrite)

        # Report the extraction time
        log.info(utils.get_time_string(time.time()-tstart))
