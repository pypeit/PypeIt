"""
Script to preprocess JWST data through the JWST pipeline (Stage 1 and Stage 2)
prior to ingestion by PypeIt.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pypeit.scripts import scriptbase


class PreprocessJWST(scriptbase.ScriptBase):

    @staticmethod
    def default_det1_params():
        """
        Return the default parameter dictionary for the JWST Detector1Pipeline
        (Stage 1).

        Returns
        -------
        :obj:`dict`:
            Default parameters for the Detector1Pipeline.
        """
        return {
            'clean_flicker_noise': {
                'skip': False,
                'fit_method': 'fft',
                'n_sigma': 2,
                'mask_science_regions': False,
                'background_method': None,
            },
        }

    @staticmethod
    def default_spec2_params(source_type='EXTENDED'):
        """
        Return the default parameter dictionary for the JWST Spec2Pipeline
        (Stage 2).

        Parameters
        ----------
        source_type : :obj:`str`, optional
            Source type for the srctype step.  Must be ``'POINT'`` or
            ``'EXTENDED'``.

        Returns
        -------
        :obj:`dict`:
            Default parameters for the Spec2Pipeline.
        """
        return {
            'assign_wcs': {'save_results': True},
            'extract_2d': {'save_results': True},
            'bkg_subtract': {'skip': True},
            'imprint_subtract': {'save_results': True},
            'master_background_mos': {'skip': True},
            'srctype': {'source_type': source_type},
            'resample_spec': {'skip': True},
            'extract_1d': {'skip': True},
            'flat_field': {'save_interpolated_flat': True},
            # Forces the barshadow step to always run. Default is to apply
            # barshadow only for extended sources, which means the slit
            # won't be flat.
            'barshadow': {'source_type': 'EXTENDED'},
            'nsclean': {'skip': True, 'save_results': False},
        }

    @staticmethod
    def _filter_jwst_files(file_list, program=None, observation=None, visit=None,
                           visit_group=None, parallel_seq=None, activity=None,
                           exposure=None):
        """
        Filter a list of JWST files by fields in the JWST file naming
        convention:
        ``jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>_<detector>_<prodType>.fits``.

        Parameters
        ----------
        file_list : :obj:`list`
            List of `Path` objects for the JWST files.
        program : :obj:`list`, optional
            List of program ID numbers (e.g., ``['01234', '01235']``).
            Matched against the 5-digit program field in the filename.
        observation : :obj:`list`, optional
            List of observation numbers (e.g., ``['003', '004']``).
            Matched against the 3-digit observation field in the filename.
        visit : :obj:`list`, optional
            List of visit numbers (e.g., ``['001', '002']``).  Matched
            against the 3-digit visit field in the filename.
        visit_group : :obj:`list`, optional
            List of visit group numbers (e.g., ``['01', '02']``).  Matched
            against the 2-digit visit group field in the filename.
        parallel_seq : :obj:`list`, optional
            List of parallel sequence IDs (e.g., ``['1', '2']``).  1 means
            prime, 2-5 are parallel sequences.  Matched against the 1-digit
            sequence field in the filename.
        activity : :obj:`list`, optional
            List of activity numbers in base-36 (e.g., ``['01', 'a1']``).
            Matched case-insensitively against the 2-character activity field
            in the filename.
        exposure : :obj:`list`, optional
            List of exposure numbers (e.g., ``['00001', '00002']``).
            Matched against the 5-digit exposure field in the filename.

        Returns
        -------
        :obj:`list`:
            Filtered list of `Path` objects.
        """
        import re
        # Match the JWST naming pattern:
        #   jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>_<detector>_<prodType>.fits
        pattern = re.compile(
            r'^jw(\d{5})(\d{3})(\d{3})_(\d{2})(\d)([0-9a-z]{2})_(\d{5})_',
            re.IGNORECASE)
        # Zero-pad / normalise the filter values
        pid_set = {p.zfill(5) for p in program} if program is not None else None
        obs_set = {o.zfill(3) for o in observation} if observation is not None else None
        vis_set = {v.zfill(3) for v in visit} if visit is not None else None
        vg_set = {g.zfill(2) for g in visit_group} if visit_group is not None else None
        pseq_set = set(parallel_seq) if parallel_seq is not None else None
        act_set = {a.lower().zfill(2) for a in activity} if activity is not None else None
        exp_set = {e.zfill(5) for e in exposure} if exposure is not None else None
        filtered = []
        for f in file_list:
            match = pattern.match(f.name)
            if match is None:
                # Skip files that don't match the JWST naming convention
                # (they can't be filtered)
                continue
            prog, obs, vis, vg, pseq, act, exp = match.groups()
            if pid_set is not None and prog not in pid_set:
                continue
            if obs_set is not None and obs not in obs_set:
                continue
            if vis_set is not None and vis not in vis_set:
                continue
            if vg_set is not None and vg not in vg_set:
                continue
            if pseq_set is not None and pseq not in pseq_set:
                continue
            if act_set is not None and act.lower() not in act_set:
                continue
            if exp_set is not None and exp not in exp_set:
                continue
            filtered.append(f)
        return filtered

    @staticmethod
    def _deep_merge(base, override):
        """
        Recursively merge ``override`` into ``base``.  Values in ``override``
        take precedence.

        Parameters
        ----------
        base : :obj:`dict`
            Base dictionary with default values.
        override : :obj:`dict`
            Dictionary with user-provided overrides.

        Returns
        -------
        :obj:`dict`:
            Merged dictionary.
        """
        merged = base.copy()
        for key, value in override.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = PreprocessJWST._deep_merge(merged[key], value)
            else:
                merged[key] = value
        return merged

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Preprocess JWST data through the JWST Detector1Pipeline (Stage 1) and '
                        'Spec2Pipeline (Stage 2) prior to ingestion by PypeIt.',
            width=width, formatter=scriptbase.SmartFormatter, default_log_file=True
        )
        parser.add_argument('raw_dir', type=str,
                            help='Directory containing the raw uncalibrated JWST files '
                                 '(*_uncal.fits)')
        parser.add_argument('-o', '--output_dir', type=str, default='jwst_output',
                            help='Directory for the JWST pipeline output files.  '
                                 'Default is jwst_output in the current working directory.')
        parser.add_argument('--crds_path', type=str, default=None,
                            help='Path to the CRDS cache directory. If not provided, the '
                                 'CRDS_PATH environment variable must be set.')
        parser.add_argument('--pid', type=str, nargs='+', default=None,
                            help='Select only files matching these JWST program ID number(s) '
                                 '(e.g., 01234 01235).  Based on the JWST file naming convention: '
                                 'jw<ppppp><ooo><vvv>_...  If not provided, files are not '
                                 'filtered by program ID.')
        parser.add_argument('--obs', type=str, nargs='+', default=None,
                            help='Select only files matching these JWST observation number(s) '
                                 '(e.g., 003 004).  Based on the JWST file naming convention: '
                                 'jw<ppppp><ooo><vvv>_...  If not provided, files are not '
                                 'filtered by observation number.')
        parser.add_argument('--vis', type=str, nargs='+', default=None,
                            help='Select only files matching these JWST visit number(s) '
                                 '(e.g., 001 002).  Based on the JWST file naming convention: '
                                 'jw<ppppp><ooo><vvv>_...  If not provided, files are not '
                                 'filtered by visit number.')
        parser.add_argument('--vg', type=str, nargs='+', default=None,
                            help='Select only files matching these JWST visit group number(s) '
                                 '(e.g., 01 02).  Based on the JWST file naming convention: '
                                 'jw<ppppp><ooo><vvv>_<gg><s><aa>_...  If not provided, files '
                                 'are not filtered by visit group.')
        parser.add_argument('--s', type=str, nargs='+', default=None,
                            help='Select only files matching these JWST parallel sequence ID(s) '
                                 '(e.g., 1 2).  1 = prime, 2-5 = parallel.  Based on the JWST '
                                 'file naming convention: '
                                 'jw<ppppp><ooo><vvv>_<gg><s><aa>_...  If not provided, files '
                                 'are not filtered by parallel sequence.')
        parser.add_argument('--a', type=str, nargs='+', default=None,
                            help='Select only files matching these JWST activity number(s) '
                                 '(base-36, e.g., 01 02).  Based on the JWST file naming '
                                 'convention: jw<ppppp><ooo><vvv>_<gg><s><aa>_...  If not '
                                 'provided, files are not filtered by activity number.')
        parser.add_argument('--exp', type=str, nargs='+', default=None,
                            help='Select only files matching these JWST exposure number(s) '
                                 '(e.g., 00001 00002).  Based on the JWST file naming '
                                 'convention: jw<ppppp><ooo><vvv>_<gg><s><aa>_<eeeee>_...  '
                                 'If not provided, files are not filtered by exposure number.')

        parser.add_argument('-l', '--list', default=False, action='store_true',
                            help='List the files that would be processed and exit.  No '
                                 'preprocessing is performed.')
        parser.add_argument('--source_type', type=str, default='EXTENDED',
                            choices=['POINT', 'EXTENDED'],
                            help='Source type for the JWST Spec2Pipeline srctype step.'
                                 ' Default is EXTENDED.')
        parser.add_argument('--det1_config', type=str, default=None,
                            help='R|JSON file with parameters for the JWST Detector1Pipeline '
                                 '(Stage 1).  Parameters provided in this file are merged with '
                                 'the defaults, overriding any matching keys.  If not provided, '
                                 'the built-in defaults are used.  Use --save_config to generate '
                                 'template files with the default parameters.  Example JSON '
                                 'format:\n\n'
                                 'F|    {\n'
                                 'F|        "clean_flicker_noise": {\n'
                                 'F|            "skip": false,\n'
                                 'F|            "fit_method": "fft",\n'
                                 'F|            "n_sigma": 2\n'
                                 'F|        }\n'
                                 'F|    }\n')
        parser.add_argument('--spec2_config', type=str, default=None,
                            help='R|JSON file with parameters for the JWST Spec2Pipeline '
                                 '(Stage 2).  Parameters provided in this file are merged with '
                                 'the defaults, overriding any matching keys.  If not provided, '
                                 'the built-in defaults are used.  Use --save_config to generate '
                                 'template files with the default parameters.  Example JSON '
                                 'format:\n\n'
                                 'F|    {\n'
                                 'F|        "extract_1d": {"skip": false},\n'
                                 'F|        "bkg_subtract": {"skip": false}\n'
                                 'F|    }\n')
        parser.add_argument('--save_config', default=False, action='store_true',
                            help='Save the default Detector1Pipeline and Spec2Pipeline '
                                 'parameter dictionaries to JSON files '
                                 '(det1_config.json, spec2_config.json) in the output '
                                 'directory and exit.  No preprocessing is performed when '
                                 'this flag is set.  These files can then be edited and '
                                 'passed back via --det1_config and --spec2_config.  Note '
                                 'that even without this flag, the parameters used for each '
                                 'run are always saved to the output directory.')
        parser.add_argument('--overwrite_stage1', default=False, action='store_true',
                            help='Overwrite existing Stage 1 output files')
        parser.add_argument('--overwrite_stage2', default=False, action='store_true',
                            help='Overwrite existing Stage 2 output files')
        return parser

    @classmethod
    def main(cls, args):
        import json
        import os
        from pathlib import Path

        from pypeit import log
        from pypeit import PypeItError

        # Initialize the log
        cls.init_log(args)

        # ----------------------------------------------------------
        # Handle --save_config: write default JSON files and exit
        # ----------------------------------------------------------
        if args.save_config:
            output_dir = Path(args.output_dir).resolve()
            output_dir.mkdir(parents=True, exist_ok=True)

            det1_file = output_dir / 'det1_config.json'
            spec2_file = output_dir / 'spec2_config.json'

            with open(det1_file, 'w') as f:
                json.dump(cls.default_det1_params(), f, indent=4)
            log.info(f'Saved default Detector1Pipeline parameters to {det1_file}')

            with open(spec2_file, 'w') as f:
                json.dump(cls.default_spec2_params(source_type=args.source_type), f, indent=4)
            log.info(f'Saved default Spec2Pipeline parameters to {spec2_file}')

            log.info('Edit these files and pass them back via --det1_config and --spec2_config.')
            return

        # Set up the CRDS path
        if args.crds_path is not None:
            os.environ['CRDS_PATH'] = args.crds_path
        elif 'CRDS_PATH' not in os.environ:
            raise PypeItError('CRDS_PATH environment variable is not set and --crds_path was not '
                              'provided. Please set one of them.')

        # Import JWST pipeline modules (after setting CRDS_PATH)
        from jwst.pipeline import Detector1Pipeline
        from jwst.pipeline import Spec2Pipeline

        # Set up directories
        raw_dir = Path(args.raw_dir).resolve()
        if not raw_dir.exists():
            raise PypeItError(f'Raw directory {raw_dir} does not exist. Please download the raw '
                              'data first.')

        output_dir = Path(args.output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get the list of uncalibrated files
        uncal_list = sorted(raw_dir.glob('*_uncal.fits'))
        if len(uncal_list) == 0:
            raise PypeItError(f'No uncalibrated files (*_uncal.fits) found in {raw_dir}.')

        # Filter by JWST file naming convention fields
        _any_filter = any(x is not None for x in [args.pid, args.obs, args.vis,
                                                   args.vg, args.s, args.a,
                                                   args.exp])
        if _any_filter:
            uncal_list = cls._filter_jwst_files(uncal_list, program=args.pid,
                                                observation=args.obs, visit=args.vis,
                                                visit_group=args.vg,
                                                parallel_seq=args.s,
                                                activity=args.a,
                                                exposure=args.exp)
            if len(uncal_list) == 0:
                raise PypeItError('No uncalibrated files match the specified filters.')

        log.info(f'Found {len(uncal_list)} uncalibrated files to process.')

        # If --list is set, print the file list and exit
        if args.list:
            log.info('Files to be processed:')
            for f in uncal_list:
                log.info(f'  {f.name}')
            return

        # ---------------------------------------------------------------
        # Build Stage 1 parameters: defaults + user overrides
        # ---------------------------------------------------------------
        parameter_dict_det1 = cls.default_det1_params()
        if args.det1_config is not None:
            det1_config_file = Path(args.det1_config)
            if not det1_config_file.exists():
                raise PypeItError(f'Detector1Pipeline config file {det1_config_file} not found.')
            with open(det1_config_file) as f:
                user_det1 = json.load(f)
            parameter_dict_det1 = cls._deep_merge(parameter_dict_det1, user_det1)
            log.info(f'Merged user Detector1Pipeline parameters from {det1_config_file}.')
        log.info(f'Detector1Pipeline parameters: {parameter_dict_det1}')

        # Save the parameters used for this run
        det1_config_out = output_dir / 'det1_config.json'
        with open(det1_config_out, 'w') as f:
            json.dump(parameter_dict_det1, f, indent=4)
        log.info(f'Saved Detector1Pipeline parameters to {det1_config_out}')

        # ---------------------------------------------------------------
        # Stage 1: Run the Detector1Pipeline on each uncalibrated file
        # ---------------------------------------------------------------
        for uncal_file in uncal_list:
            rate_file = output_dir / uncal_file.name.replace('_uncal.fits', '_rate.fits')
            if rate_file.exists() and not args.overwrite_stage1:
                log.info(f'Stage 1 output {rate_file} already exists. Skipping stage 1 '
                         'processing.')
                continue
            log.info(f'Running Detector1Pipeline on {uncal_file}')
            Detector1Pipeline.call(str(uncal_file), save_results=True,
                                   output_dir=str(output_dir), steps=parameter_dict_det1)

        # ---------------------------------------------------------------
        # Build Stage 2 parameters: defaults + user overrides
        # ---------------------------------------------------------------
        param_dict_spec2 = cls.default_spec2_params(source_type=args.source_type)
        if args.spec2_config is not None:
            spec2_config_file = Path(args.spec2_config)
            if not spec2_config_file.exists():
                raise PypeItError(f'Spec2Pipeline config file {spec2_config_file} not found.')
            with open(spec2_config_file) as f:
                user_spec2 = json.load(f)
            param_dict_spec2 = cls._deep_merge(param_dict_spec2, user_spec2)
            log.info(f'Merged user Spec2Pipeline parameters from {spec2_config_file}.')
        log.info(f'Spec2Pipeline parameters: {param_dict_spec2}')

        # Save the parameters used for this run
        spec2_config_out = output_dir / 'spec2_config.json'
        with open(spec2_config_out, 'w') as f:
            json.dump(param_dict_spec2, f, indent=4)
        log.info(f'Saved Spec2Pipeline parameters to {spec2_config_out}')

        # ---------------------------------------------------------------
        # Stage 2: Run the Spec2Pipeline on each rate file
        # ---------------------------------------------------------------
        rate_list = sorted(output_dir.glob('*_rate.fits'))
        if len(rate_list) == 0:
            raise PypeItError(f'No rate files (*_rate.fits) found in {output_dir}.')

        # Apply the same filtering to rate files
        if _any_filter:
            rate_list = cls._filter_jwst_files(rate_list, program=args.pid,
                                               observation=args.obs, visit=args.vis,
                                               visit_group=args.vg,
                                               parallel_seq=args.s,
                                               activity=args.a,
                                               exposure=args.exp)
            if len(rate_list) == 0:
                raise PypeItError('No rate files match the specified filters.')

        log.info(f'Found {len(rate_list)} rate files for Stage 2 processing.')

        for rate_file in rate_list:
            assign_wcs_file = output_dir / rate_file.name.replace('_rate.fits',
                                                                  '_msa._assign_wcs')
            if assign_wcs_file.exists() and not args.overwrite_stage2:
                log.info(f'Stage 2 output {assign_wcs_file} already exists. Skipping stage 2 '
                         'processing.')
                continue
            log.info(f'Running Spec2Pipeline on {rate_file}')
            Spec2Pipeline.call(rate_file, save_results=True, output_dir=str(output_dir),
                               steps=param_dict_spec2, input_dir=str(raw_dir))

