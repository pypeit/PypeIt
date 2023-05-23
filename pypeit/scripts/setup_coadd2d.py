"""
Script for preparing a 2d-coadd configuration file.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class SetupCoAdd2D(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
                description='Prepare a configuration file for performing 2D coadds', width=width)

        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument('-f', '--pypeit_file', type=str, default=None,
                                 help='PypeIt reduction file')
        input_group.add_argument('-d', '--science_dir', type=str, default=None,
                                 help='Directory with spec2d files to stack')

        parser.add_argument('--keep_par', dest='clean_par', default=True, action='store_false',
                            help='Propagate all parameters from the pypeit file to the coadd2d '
                                 'file(s).  If not set, only the required parameters and their '
                                 'default values are included in the output file(s).')
        parser.add_argument('--obj', type=str, nargs='+',
                            help='Limit the coadd2d files created to observations of the '
                                 'specified target.  If not provided, a coadd2D file is written '
                                 'for each target found in the science directory.  The target '
                                 'names are included in the PypeIt spec2d file names.'
                                 'For example, the target for spec2d file '
                                 '"spec2d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits" '
                                 'is "pisco".')
        parser.add_argument('--det', type=str, nargs='+',
                            help='A space-separated set of detectors or detector mosaics to '
                                 'coadd.  By default, *all* detectors or default mosaics for '
                                 'this instrument will be coadded.  Detectors in a mosaic must '
                                 'be a mosaic "allowed" by PypeIt and should be provided as '
                                 'comma-separated integers (with no spaces).  For example, to '
                                 'separately coadd detectors 1 and 5 for Keck/DEIMOS, you would '
                                 'use --det 1 5; to coadd mosaics made up of detectors 1,5 and '
                                 '3,7, you would use --det 1,5 3,7')
        parser.add_argument('--only_slits', type=str, nargs='+',
                            help='A space-separated set of slits to coadd.  If not provided, all '
                                 'slits are coadded.')
        parser.add_argument('--offsets', type=str, default=None,
                            help='Spatial offsets to apply to each image; see the '
                                 '[coadd2d][offsets] parameter.  Options are restricted here to '
                                 'either maskdef_offsets or auto.  If not specified, the '
                                 '(spectrograph-specific) default is used.  Other options exist '
                                 'but must be entered by directly editing the coadd2d file.')
        parser.add_argument('--weights', type=str, default=None,
                            help='Weights used to coadd images; see the [coadd2d][weights] '
                                 'parameter.  Options are restricted here to '
                                 'either uniform or auto.  If not specified, the '
                                 '(spectrograph-specific) default is used.  Other options exist '
                                 'but must be entered by directly editing the coadd2d file.')

        return parser

    @staticmethod
    def main(args):

        from pathlib import Path

        from IPython import embed

        import numpy as np

        from astropy.table import Table

        from pypeit import msgs
        from pypeit import io
        from pypeit import utils
        from pypeit import inputfiles
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.coadd2d import CoAdd2D

        if args.pypeit_file is None:
            pypeitFile = None
            par = None
            spec_name = None
            sci_dir = Path(args.science_dir).resolve()
        else:
            # Read the pypeit file
            pypeitFile = inputfiles.PypeItFile.from_file(args.pypeit_file)
            # Get the spectrograph instance and the parameters used
            spec, par, _ = pypeitFile.get_pypeitpar()
            spec_name = spec.name
            # Get the Science directory used.
            # NOTE: When [rdx][redux_path] isn't defined, the get_pypeitpar()
            # function will use the current working directory.  If the pypeit_file
            # being source is in a different directory, this means that the Science
            # directory won't be found.  So first try to set the science directory
            # based on the parameter value, then try to base it on the parent
            # directory of the provided pypeit file.  The latter is critical to the
            # vet_test in the dev-suite.
            sci_dir = Path(par['rdx']['redux_path']).resolve() / par['rdx']['scidir']
            if not sci_dir.exists():
                sci_dir = Path(args.pypeit_file).resolve().parent / par['rdx']['scidir']

        if not sci_dir.exists():
            msgs.error(f'Science directory not found: {sci_dir}\n')

        # Find all the spec2d files:
        spec2d_files = sorted(sci_dir.glob('spec2d*'))
        if len(spec2d_files) == 0:
            msgs.error(f'No spec2d files found in {sci_dir}.')

        if spec_name is None:
            with io.fits_open(spec2d_files[0]) as hdu:
                spec_name = hdu[0].header['PYP_SPEC']

        # Get the set of objects
        # TODO: Direct parsing of the filenames will be wrong if any of the
        # reduced files have dashes in them or if the objects have underscores.
        objects = np.unique(pypeitFile.data['target'].data 
                                if pypeitFile is not None and 'target' in pypeitFile.data.keys()
                                else [f.name.split('-')[1].split('_')[0] for f in spec2d_files])
        if args.obj is not None:
            # Limit to the selected objects
            _objects = [o for o in objects if o == args.obj]
            # Check some were found
            if len(_objects) == 0:
                msgs.error('Unable to find relevant objects.  Unique objects are '
                           f'{objects.tolist()}; you requested {args.obj}.')
            objects = _objects

        # Match spec2d files to objects
        object_spec2d_files = {}
        for obj in objects:
            object_spec2d_files[obj] = [f for f in spec2d_files if obj in f.name]
            if len(object_spec2d_files[obj]) == 0:
                msgs.warn(f'No spec2d files found for target={obj}.')
                del object_spec2d_files[obj]

        # Check spec2d files exist for the selected objects
        if len(object_spec2d_files.keys()) == 0:
            msgs.error('Unable to match any spec2d files to objects.')

        # Add the paths to make sure they match the pypeit file.
        # NOTE: cfg does *not* need to include the spectrograph parameter in
        # rdx.  This is added by `CoAdd2D.default_par`.
        cfg = {} if args.clean_par or pypeitFile is None else dict(pypeitFile.config)
        if par is None:
            # This is only used to get the default offsets and weights (see below)
            par = load_spectrograph(spec_name).default_pypeit_par()
            utils.add_sub_dict(cfg, 'rdx')
            cfg['rdx']['scidir'] = sci_dir.name
        else:
            utils.add_sub_dict(cfg, 'rdx')
            cfg['rdx']['redux_path'] = par['rdx']['redux_path']
            cfg['rdx']['scidir'] = par['rdx']['scidir']
            cfg['rdx']['qadir'] = par['rdx']['qadir']
            utils.add_sub_dict(cfg, 'calibrations')
            cfg['calibrations']['calib_dir'] = par['calibrations']['calib_dir']
        utils.add_sub_dict(cfg, 'coadd2d')
        cfg['coadd2d']['offsets'] = par['coadd2d']['offsets'] \
                if args.offsets is None else args.offsets
        cfg['coadd2d']['weights'] = par['coadd2d']['weights'] \
                if args.weights is None else args.weights

        # Build the default parameters
        cfg = CoAdd2D.default_par(spec_name, inp_cfg=cfg, det=args.det, slits=args.only_slits)

        # Create a coadd2D file for each object
        # NOTE: Below expect all spec2d files have the same path
        for obj, files in object_spec2d_files.items():
            tbl = Table()
            tbl['filename'] = [f.name for f in files]
            ofile = args.pypeit_file.replace('.pypeit', f'_{obj}.coadd2d')
            inputfiles.Coadd2DFile(config=cfg, file_paths=[str(sci_dir)],
                                   data_table=tbl).write(ofile)

