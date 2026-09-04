"""
Prepare datacube coadd and extraction setup files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pathlib import Path

from astropy.table import Table

from pypeit import log
from pypeit import PypeItError
from pypeit import inputfiles
from pypeit import outputfiles
from pypeit.par.pypeitpar import CubeExtractionPar
from pypeit.scripts import scriptbase


def _parse_whitelight_range(value):
    """
    Argparse ``type`` callable: parse a comma-separated white-light wavelength range
    into a two-element list, with either entry allowed to be ``none``
    (case-insensitive) for an unset bound. Raising ``ValueError`` here lets argparse
    format a clean CLI error message; the resulting list's length/content is validated
    for real by :class:`~pypeit.par.pypeitpar.CubeExtractionPar`.
    """
    entries = [v.strip() for v in value.split(',')]
    if len(entries) != 2:
        raise ValueError(f'whitelight_range must be provided as min,max; e.g., 9400,10000. Got: {value}')
    return [None if e.lower() == 'none' else float(e) for e in entries]


class SetupDataCube(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Prepare .coadd3d and .extract files for point-source datacube work.',
            width=width, default_log_file=True
        )
        parser.add_argument('pypeit_file', type=str, help='PypeIt reduction file.')
        parser.add_argument('target', type=str, help='Target name, e.g. J0750+6927.')
        parser.add_argument(
            '--sensfile', type=str, default=None,
            help='Optional sensitivity function file. If omitted, the generated .coadd3d '
                 'file will produce unfluxed cubes.'
        )
        parser.add_argument(
            '--whitelight_range', '--wl_range', dest='whitelight_range',
            type=_parse_whitelight_range, default='None,None',
            help='White-light wavelength range, e.g. 9400,10000. Defaults to None,None.'
        )
        parser.add_argument(
            '--manual', type=str, default=None,
            help='Manual extraction position in colon-separated x:y coordinates, e.g. '
                 '10.0:14.0. Do not use commas. If provided, the extract file uses '
                 'opt_prof_method=user_gauss.'
        )
        parser.add_argument('--fwhm', type=float, default=1.1, help='Extraction FWHM in arcsec.')
        parser.add_argument(
            '--snr_thresh', type=float, default=4.0,
            help='S/N threshold for automatic source finding in the extraction file.'
        )
        parser.add_argument(
            '--spatial_delta', type=float, default=None,
            help='Output cube spatial sampling in arcsec. Defaults to 0.678924 for keck_kcrm; '
                 'otherwise omitted.'
        )
        parser.add_argument('--det', type=int, default=1, help='Detector number.')
        parser.add_argument(
            '-o', '--overwrite', default=False, action='store_true',
            help='Overwrite an existing .extract file. The .coadd3d file is always refreshed.'
        )
        parser.add_argument(
            '--append', default=False, action='store_true',
            help='Append newly reduced spec2d files to an existing .coadd3d file without changing '
                 'any other lines. The .extract file is left unchanged.'
        )
        return parser

    @classmethod
    def main(cls, args):
        cls.init_log(args)

        pypeit_path = Path(args.pypeit_file).absolute()
        pypeit_file = inputfiles.PypeItFile.from_file(str(pypeit_path))
        if 'rdx' not in pypeit_file.config or 'spectrograph' not in pypeit_file.config['rdx']:
            raise PypeItError('The PypeIt file must define [rdx] spectrograph.')
        spectrograph = pypeit_file.config['rdx']['spectrograph']
        # require_rawfile=False: this script only needs the parameter set (to locate the
        # Science directory below); the raw data referenced by the .pypeit file may no
        # longer be available, and get_meta_value() already prefers the .pypeit file's
        # own table columns over reading the raw file in the first place.
        spec, par, _ = pypeit_file.get_pypeitpar(require_rawfile=False)

        sensfile = None
        if args.sensfile is not None:
            sensfile = Path(args.sensfile).expanduser().absolute()
            if not sensfile.is_file():
                raise PypeItError(f'Sensitivity function does not exist: {sensfile}')

        if args.manual is not None and ',' in args.manual:
            raise PypeItError(
                'Manual datacube extraction positions must use colon-separated x:y '
                f'syntax, not commas: {args.manual}. For example, use --manual 10.0:14.0.'
            )
        try:
            # Structural validation (whitelight_range length; manual's x:y format) --
            # the same checks CubePar/CubeExtractionPar apply when the written files are
            # later read back in.
            CubeExtractionPar(whitelight_range=args.whitelight_range, manual=args.manual)
        except ValueError as exc:
            raise PypeItError(str(exc)) from exc
        if args.manual is not None:
            x_str, y_str = args.manual.split(':')
            try:
                float(x_str)
                float(y_str)
            except ValueError as exc:
                raise PypeItError(
                    f'Manual datacube extraction positions must be numeric x:y values. '
                    f'Invalid value: {args.manual}.'
                ) from exc

        # NOTE: When [rdx][redux_path] isn't defined in the .pypeit file, get_pypeitpar()
        # defaults it to the current working directory, not the .pypeit file's own
        # directory. So first try the parameter value, then fall back to the .pypeit
        # file's parent directory -- the same two-step approach setup_coadd2d.py uses.
        sci_dir = outputfiles.science_path(par)
        if not sci_dir.is_dir():
            sci_dir = pypeit_path.parent / par['rdx']['scidir']
        if not sci_dir.is_dir():
            raise PypeItError(f'Expected Science directory does not exist: {sci_dir}')

        spec2d_files, missing, target_name = outputfiles.existing_spec2d_files(
            pypeit_file, args.target, sci_dir, spec
        )
        for raw_stem in missing:
            log.warning(f'Expected spec2d product for {raw_stem} not found yet; skipping for now.')
        if len(spec2d_files) == 0:
            raise PypeItError(
                f'No reduced spec2d files found for target={args.target} in {sci_dir}.'
            )

        target_stub = target_name
        source_dir = pypeit_path.parent / 'sources' / target_stub
        source_dir.mkdir(parents=True, exist_ok=True)

        coadd3d_file = source_dir / f'{target_stub}.coadd3d'
        extract_file = source_dir / f'{target_stub}.extract'

        if args.append:
            if not coadd3d_file.is_file():
                raise PypeItError(f'Cannot append to missing .coadd3d file: {coadd3d_file}')
            existing = inputfiles.Coadd3DFile.from_file(str(coadd3d_file), preserve_comments=True)
            existing_names = {Path(f).name for f in existing.data['filename']}
            new_files = [f for f in spec2d_files if f.name not in existing_names]
            if len(new_files) == 0:
                log.info(f'No new spec2d files to append to {coadd3d_file}.')
            else:
                for f in new_files:
                    existing.data.add_row([f.name])
                existing.write(str(coadd3d_file))
                log.info(f'Appended {len(new_files)} spec2d file(s) to {coadd3d_file}:')
                for f in new_files:
                    log.info(f'  {f.name}')
            if extract_file.exists():
                log.info(f'Leaving existing extract file unchanged: {extract_file}')
            return

        coadd3d_cfg = {
            'rdx': {'spectrograph': spectrograph, 'detnum': args.det},
            'reduce': {'cube': {
                'whitelight_range': args.whitelight_range,
                'output_filename': target_stub,
                'combine': True,
                'alignment_method': 'none',
                'method': 'ngp',
                'spat_subpixel': 1,
                'slice_subpixel': 1,
                'spec_subpixel': 1,
                'astrometric': False,
                'save_whitelight': True,
                'weight_method': 'uniform',
            }}
        }
        if sensfile is not None:
            coadd3d_cfg['reduce']['cube']['sensfile'] = str(sensfile)
        if args.spatial_delta is not None:
            coadd3d_cfg['reduce']['cube']['spatial_delta'] = args.spatial_delta

        tbl = Table()
        tbl['filename'] = [f.name for f in spec2d_files]
        coadd3d = inputfiles.Coadd3DFile(
            config=coadd3d_cfg, file_paths=[str(sci_dir)], data_table=tbl
        )
        # Preserve the same helpful, commented-out example this file used to include.
        coadd3d.config['reduce']['cube'].comments['weight_method'] = [
            '# Optional initial object position for relative/auto weighting, in x:y format.',
            '# Use this if the automatic weight-position finder locks onto the wrong source.',
            '# The x,y values should match the white-light image coordinates read from Ginga/DS9.',
            '# weights_init_obj_pos = x:y',
        ]
        coadd3d.write(str(coadd3d_file))
        log.info(f'Wrote {coadd3d_file}')

        if extract_file.exists() and not args.overwrite:
            log.warning(f'{extract_file} exists; leaving it unchanged. Use -o/--overwrite to replace it.')
        else:
            extract_cfg = {
                'reduce': {'cube': {'extraction': {
                    'output_filename': f'{target_stub}_extract',
                    'whitelight_range': args.whitelight_range,
                    'snr_thresh': args.snr_thresh,
                    'opt_prof_method': 'fit_gauss' if args.manual is None else 'user_gauss',
                    'fwhm': args.fwhm,
                }}}
            }
            if args.manual is not None:
                extract_cfg['reduce']['cube']['extraction']['manual'] = args.manual
            inputfiles.ExtractFile(config=extract_cfg).write(str(extract_file))
            log.info(f'Wrote {extract_file}')
