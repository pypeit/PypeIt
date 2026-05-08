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
from pathlib import Path
from IPython import embed

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
                            help='Basename for output files, i.e. outputs will be written to'
                            'spec1d_basename.fits and spec2d_basename.fits')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')
        parser.add_argument('-b', '--boxcar_radius', type=float, default=None,
                            help='Radius of the circular boxcar (in arcseconds) to use for the extraction.')
        parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
        return parser

    @classmethod
    def main(cls, args):
        import time

        from pypeit import log
        from pypeit import PypeItError
        from pypeit.par import pypeitpar
        from pypeit import inputfiles
        from pypeit import utils
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.coadd3d import DataCube, CoAdd3D

        # Initialize the log
        cls.init_log(args)

        # Check that a file has been provided
        if args.file is None:
            raise PypeItError('You must input a spec3d (i.e. PypeIt DataCube) fits file')
        extcube = DataCube.from_file(args.file)
        spectrograph = load_spectrograph(extcube.PYP_SPEC)

        if args.ext_file is None:
            par = spectrograph.default_pypeit_par()
        else:
            # Read in the relevant information from the .extract file
            ext3dfile = inputfiles.ExtractFile.from_file(args.ext_file)
            # Parameters
            spectrograph_def_par = spectrograph.default_pypeit_par()
            par = pypeitpar.PypeItPar.from_cfg_lines(
                cfg_lines=spectrograph_def_par.to_config(), merge_with=(ext3dfile.cfg_lines,)
            )

        # Set the boxcar radius
        boxcar_radius = args.boxcar_radius

        # Set the output name. If one was provided by the user 
        if args.save is not None:
            par['reduce']['cube_extraction']['output_filename'] = args.save
        if args.boxcar_radius is not None:
            par['reduce']['cube_extraction']['boxcar_radius'] = args.boxcar_radius
        
        # Load the DataCube
        tstart = time.time()
        
        # Get the paths
        coadd_scidir, qa_path = map(lambda x : Path(x).absolute(),
                CoAdd3D.output_paths(args.file, par, coadd_dir=par['rdx']['redux_path']))

        # Extract the spectrum
        extcube.extract_spec(
            par['reduce']['cube_extraction'], output_dir=str(coadd_scidir), overwrite=args.overwrite, 
            debug=args.debug)

        # Report the extraction time
        log.info(utils.get_time_string(time.time()-tstart))
