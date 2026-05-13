"""
This script enables the user to convert spec2D FITS files
from SlicerIFU instruments into a 3D cube with a defined WCS.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pypeit.scripts import scriptbase
from IPython import embed

class CoAddDataCube(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Read in an array of spec2D files and convert them into a datacube',
            width=width, default_log_file=True
        )
        parser.add_argument('file', type=str, default=None, help='filename.coadd3d file')
        parser.add_argument('--det', default=1, type=int, help="Detector")
        parser.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Overwrite any existing files/directories'
        )
        parser.add_argument(
            "--debug", default=False, action="store_true", help="Run in debugging mode"
        )
        return parser

    @classmethod
    def main(cls, args):
        import time

        from pathlib import Path
        from pypeit import log
        from pypeit import PypeItError
        from pypeit import par
        from pypeit import inputfiles
        from pypeit import utils
        from pypeit.coadd3d import CoAdd3D
        from pypeit.spectrographs.util import load_spectrograph

        # Initialize the log
        cls.init_log(args)

        # Check that a file has been provided
        if args.file is None:
            raise PypeItError('You must input a coadd3d file')

        # Read in the relevant information from the .coadd3d file
        coadd3dfile = inputfiles.Coadd3DFile.from_file(args.file)
        spectrograph = load_spectrograph(
            coadd3dfile.config['rdx']['spectrograph'], pypeit_fits=True
        )

        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        parset = par.PypeItPar.from_cfg_lines(
            cfg_lines=spectrograph_def_par.to_config(), merge_with=(coadd3dfile.cfg_lines,)
        )

        # If detector was passed as an argument override whatever was in the coadd3d file
        if args.det is not None:
            log.info("Restricting to detector={}".format(args.det))
            parset['rdx']['detnum'] = int(args.det)

        # Extract the options
        ra_offsets = coadd3dfile.options['ra_offset']
        dec_offsets = coadd3dfile.options['dec_offset']
        skysub_frame = coadd3dfile.options['skysub_frame']
        scale_corr = coadd3dfile.options['scale_corr']
        sensfile = coadd3dfile.options['sensfile']
        grating_corr = coadd3dfile.options['grating_corr']
        
#        # Get the paths
#        coadd_scidir, qa_path = map(lambda x : Path(x).absolute(),
#                CoAdd3D.output_paths(coadd3dfile.filenames, parset, coadd_dir=parset['rdx']['redux_path']))
#
#        # Write the par to disk
#        par_outfile = coadd_scidir.parent / f"{parset['reduce']['cube']['output_filename']}_datacube.par"
#        log.info(f'Writing full parameter set to {par_outfile}.')
#        parset.to_config(par_outfile, exclude_defaults=True, include_descr=False)

        # Instantiate CoAdd3d
        tstart = time.time()
        coadd = CoAdd3D.get_instance(
            coadd3dfile.filenames, parset, output_dir=parset['rdx']['redux_path'], #str(coadd_scidir),
            skysub_frame=skysub_frame, sensfile=sensfile, scale_corr=scale_corr,
            grating_corr=grating_corr, ra_offsets=ra_offsets, dec_offsets=dec_offsets,
            spectrograph=spectrograph, det=args.det, overwrite=args.overwrite, debug=args.debug
        )

        # Coadd the files
        coadd.run()

        # Report the execution time
        log.info(utils.get_time_string(time.time()-tstart))
