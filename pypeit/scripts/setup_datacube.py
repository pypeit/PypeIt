"""
Prepare datacube coadd and extraction setup files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pathlib import Path

from astropy.io import fits

from pypeit import log
from pypeit import PypeItError
from pypeit import inputfiles
from pypeit.scripts import scriptbase


def target_match_key(target):
    """
    Normalize a target name for permissive matching.

    Parameters
    ----------
    target : :obj:`str`
        Target name, e.g. ``J0750+6927``.

    Returns
    -------
    :obj:`str`
        Normalized target name, e.g. ``J0750p6927``.
    """
    return str(target).strip().replace(' ', '').replace('+', 'p').replace('-', 'm')


def target_matches(value, target):
    """
    Compare target names, allowing either the literal or sanitized form.

    Parameters
    ----------
    value : :obj:`str`
        Target name to test, e.g. as read from a PypeIt file data table.
    target : :obj:`str`
        Target name to compare against, e.g. as provided on the command
        line.

    Returns
    -------
    :obj:`bool`
        True if ``value`` and ``target`` are identical, or if their
        :func:`target_match_key`-normalized forms are identical.
    """
    _value = str(value).strip()
    _target = str(target).strip()
    return _value == _target or target_match_key(_value) == target_match_key(_target)


def validate_whitelight_range(value):
    """
    Validate and normalize a command-line white-light wavelength range.

    Parameters
    ----------
    value : :obj:`str`
        Comma-separated wavelength range, e.g. ``9400,10000``. Either entry
        may instead be ``none`` (case-insensitive) to leave that bound
        unset.

    Returns
    -------
    :obj:`str`
        The comma-separated range, with each entry stripped of surrounding
        whitespace.

    Raises
    ------
    PypeItError
        If ``value`` does not contain exactly two comma-separated entries,
        or if a non-``none`` entry cannot be parsed as a float.
    """
    entries = [v.strip() for v in value.split(',')]
    if len(entries) != 2:
        raise PypeItError('whitelight_range must be provided as min,max; e.g., 9400,10000.')
    for entry in entries:
        if entry.lower() == 'none':
            continue
        try:
            float(entry)
        except ValueError as exc:
            raise PypeItError(
                'whitelight_range entries must be numeric or None; '
                f'could not parse {entry}.'
            ) from exc
    return ','.join(entries)


def spec2d_target(spec2d_file):
    """
    Read the target name from a spec2d primary header.

    Parameters
    ----------
    spec2d_file : :obj:`str`, `Path`_
        Path to a spec2d FITS file.

    Returns
    -------
    :obj:`str`
        The first non-empty value found among the ``TARGET``, ``TARGNAME``,
        and ``OBJECT`` header keywords, in that order. Returns None if none
        of those keywords are present or non-empty.
    """
    header = fits.getheader(spec2d_file, 0)
    for key in ('TARGET', 'TARGNAME', 'OBJECT'):
        if key in header and str(header[key]).strip() != '':
            return str(header[key]).strip()
    return None


def science_directory(pypeit_file):
    """
    Determine the default Science directory for a PypeIt reduction file.

    Parameters
    ----------
    pypeit_file : :obj:`str`, `Path`_
        Path to a PypeIt reduction file.

    Returns
    -------
    `Path`_
        Absolute path to the ``Science`` directory expected alongside
        ``pypeit_file``.
    """
    return Path(pypeit_file).absolute().parent / 'Science'


def matching_science_rows(pypeit_file, target):
    """
    Select science rows matching a target from a PypeIt input file.

    Parameters
    ----------
    pypeit_file : :class:`~pypeit.inputfiles.PypeItFile`
        The parsed PypeIt reduction file.
    target : :obj:`str`
        Target name to match against the data block's ``target`` column;
        see :func:`target_matches`.

    Returns
    -------
    :obj:`list`
        List of data-table rows with frame type ``science`` whose target
        matches ``target``.

    Raises
    ------
    PypeItError
        If ``pypeit_file`` has no data block, the data block has no
        ``target`` column, or no matching science rows are found.
    """
    if pypeit_file.data is None:
        raise PypeItError('The PypeIt file has no data block.')
    if 'target' not in pypeit_file.data.keys():
        raise PypeItError('The PypeIt file data block has no target column.')

    rows = [
        row for row in pypeit_file.data
        if any(ft.strip() == 'science' for ft in str(row['frametype']).split(','))
        and target_matches(row['target'], target)
    ]
    if len(rows) == 0:
        unique_targets = sorted({str(row['target']).strip() for row in pypeit_file.data})
        raise PypeItError(
            f'No science rows found for target={target}. Available targets are: '
            f'{", ".join(unique_targets)}.'
        )
    return rows


def group_science_rows(rows):
    """
    Group science rows by comb_id, using the first row in each group as the expected spec2d stem.

    Parameters
    ----------
    rows : :obj:`list`
        Science rows to group, e.g. as returned by :func:`matching_science_rows`.

    Returns
    -------
    :obj:`list`
        List of groups, where each group is a :obj:`list` of rows sharing
        the same ``comb_id`` (or, if a row has no usable ``comb_id``, a
        single-row group keyed by its own ``filename``).
    """
    def group_key(row):
        if 'comb_id' not in row.colnames:
            return str(row['filename'])
        comb_id = str(row['comb_id']).strip()
        return comb_id if comb_id.lower() not in ('', 'none', '-1') else str(row['filename'])

    groups = []
    seen = set()
    for row in rows:
        key = group_key(row)
        if key in seen:
            continue
        seen.add(key)
        group = [r for r in rows if group_key(r) == key]
        groups.append(group)
    return groups


def find_reduced_spec2d(science_dir, raw_filename, target):
    """
    Find the existing spec2d product for a raw file and target.

    Parameters
    ----------
    science_dir : :obj:`str`, `Path`_
        Directory to search for spec2d files, e.g. as returned by
        :func:`science_directory`.
    raw_filename : :obj:`str`
        Name of the raw science file whose reduced spec2d product is sought.
    target : :obj:`str`
        Target name that the spec2d file's header must match; see
        :func:`target_matches`.

    Returns
    -------
    `Path`_
        Path to the matching spec2d file. If more than one candidate
        matches, the first (alphabetically sorted) is used and a warning is
        logged. Returns None if no candidate matches.
    """
    raw_stem = Path(str(raw_filename).strip()).stem
    candidates = sorted(Path(science_dir).glob(f'spec2d_{raw_stem}-*.fits'))
    matches = []
    for candidate in candidates:
        try:
            hdr_target = spec2d_target(candidate)
        except Exception as exc:
            log.warning(f'Could not read target from {candidate}; skipping. Original error: {exc}')
            continue
        if hdr_target is not None and target_matches(hdr_target, target):
            matches.append(candidate)

    if len(matches) > 1:
        log.warning(
            f'Multiple matching spec2d files found for {raw_stem}; using {matches[0].name}.'
        )
    return None if len(matches) == 0 else matches[0]


def existing_spec2d_files(pypeit_file, target, science_dir):
    """
    Find the current reduced spec2d products expected for a target.

    Parameters
    ----------
    pypeit_file : :class:`~pypeit.inputfiles.PypeItFile`
        The parsed PypeIt reduction file.
    target : :obj:`str`
        Target name to match; see :func:`target_matches`.
    science_dir : :obj:`str`, `Path`_
        Directory to search for reduced spec2d files.

    Returns
    -------
    files : :obj:`list`
        Paths to the reduced spec2d files found for the target, one per
        combination group.
    missing : :obj:`list`
        Raw-file stems (see :func:`group_science_rows`) for combination
        groups whose spec2d product has not yet been reduced.
    target_name : :obj:`str`
        The literal target name as recorded in the PypeIt file's data
        block.
    """
    rows = matching_science_rows(pypeit_file, target)
    files = []
    missing = []
    for group in group_science_rows(rows):
        raw_filename = group[0]['filename']
        spec2d = find_reduced_spec2d(science_dir, raw_filename, target)
        if spec2d is None:
            missing.append(Path(str(raw_filename).strip()).stem)
            continue
        files.append(spec2d)
    return files, missing, str(rows[0]['target']).strip()


def write_coadd3d_file(
        coadd3d_file, spectrograph, target_stub, science_dir, spec2d_files, whitelight_range,
        sensfile, spatial_delta=None, det=1):
    """
    Write a .coadd3d file for datacube construction.

    Parameters
    ----------
    coadd3d_file : :obj:`str`, `Path`_
        Output path for the ``.coadd3d`` file. Any existing file at this
        path is overwritten.
    spectrograph : :obj:`str`
        PypeIt name of the spectrograph, written to ``[rdx] spectrograph``.
    target_stub : :obj:`str`
        Target name used to set the ``output_filename`` cube parameter.
    science_dir : :obj:`str`, `Path`_
        Directory containing the reduced spec2d files, written as the
        ``spec2d read`` block's ``path``.
    spec2d_files : :obj:`list`
        Paths to the reduced spec2d files to list in the ``spec2d read``
        block.
    whitelight_range : :obj:`str`
        Comma-separated white-light wavelength range; see
        :func:`validate_whitelight_range`.
    sensfile : :obj:`str`, `Path`_, None
        Sensitivity function file to set as ``sensfile``. If None, no
        ``sensfile`` line is written and the cube is left unfluxed.
    spatial_delta : :obj:`float`, optional
        Output cube spatial sampling in arcsec. If None, no
        ``spatial_delta`` line is written and the code default is used.
    det : :obj:`int`, optional
        Detector number, written to ``[rdx] detnum``.
    """
    spatial_delta_line = (
        '' if spatial_delta is None else f'        spatial_delta = {spatial_delta}\n'
    )
    sensfile_line = '' if sensfile is None else f'        sensfile = {sensfile}\n'
    lines = [
        '# User-defined execution parameters\n',
        '[rdx]\n',
        f'    spectrograph = {spectrograph}\n',
        f'    detnum = {det}\n',
        '[reduce]\n',
        '    [[cube]]\n',
        f'        whitelight_range = {whitelight_range}\n',
        f'        output_filename = {target_stub}\n',
        '        combine = True\n',
        '        alignment_method = none\n',
        '        method = ngp\n',
        '        spat_subpixel = 1\n',
        '        slice_subpixel = 1\n',
        '        spec_subpixel = 1\n',
        '        astrometric = False\n',
        sensfile_line,
        '        save_whitelight = True\n',
        spatial_delta_line,
        '        # Optional initial object position for relative/auto weighting, in x:y format.\n',
        '        # Use this if the automatic weight-position finder locks onto the wrong source.\n',
        '        # The x,y values should match the white-light image coordinates read from Ginga/DS9.\n',
        '        # weights_init_obj_pos = x:y\n',
        '        weight_method = uniform\n',
        '\n',
        '# Read in the data\n',
        'spec2d read\n',
        f'  path {Path(science_dir).absolute()}\n',
        'filename\n',
    ]
    lines += [f'{f.name}\n' for f in spec2d_files]
    lines += ['spec2d end\n']
    Path(coadd3d_file).write_text(''.join(lines))


def append_spec2d_files(coadd3d_file, spec2d_files):
    """
    Append new spec2d filenames to an existing .coadd3d file.

    This preserves all existing file contents except for inserting new filenames
    immediately before the ``spec2d end`` line.

    Parameters
    ----------
    coadd3d_file : :obj:`str`, `Path`_
        Path to an existing ``.coadd3d`` file to modify in place.
    spec2d_files : :obj:`list`
        Paths to spec2d files to append. Any whose filename already appears
        in the file's ``spec2d read`` block is skipped.

    Returns
    -------
    :obj:`list`
        Filenames (not full paths) that were newly appended. Empty if all
        of ``spec2d_files`` were already present.

    Raises
    ------
    PypeItError
        If ``coadd3d_file`` does not exist, or if it has no ``spec2d
        read``/``spec2d end`` block.
    """
    coadd3d_path = Path(coadd3d_file)
    if not coadd3d_path.is_file():
        raise PypeItError(f'Cannot append to missing .coadd3d file: {coadd3d_path}')

    lines = coadd3d_path.read_text().splitlines(keepends=True)
    block_start = None
    block_end = None
    for idx, line in enumerate(lines):
        stripped = line.split('#')[0].strip()
        if block_start is None and stripped == 'spec2d read':
            block_start = idx
            continue
        if block_start is not None and stripped == 'spec2d end':
            block_end = idx
            break
    if block_start is None or block_end is None:
        raise PypeItError(f'Could not find spec2d read/end block in {coadd3d_path}')

    existing = set()
    for line in lines[block_start + 1:block_end]:
        stripped = line.split('#')[0].strip()
        if stripped == '' or stripped.startswith('path ') or stripped == 'filename':
            continue
        existing.add(Path(stripped).name)

    new_names = [f.name for f in spec2d_files if f.name not in existing]
    if len(new_names) == 0:
        return []

    insert_lines = [f'{name}\n' for name in new_names]
    lines[block_end:block_end] = insert_lines
    coadd3d_path.write_text(''.join(lines))
    return new_names


def validate_manual(manual):
    """
    Validate and normalize the manual datacube extraction position.

    The PypeIt configuration parser treats comma-separated values as lists, but
    the datacube extraction parser expects the manual position as one string in
    x:y format.

    Parameters
    ----------
    manual : :obj:`str`, None
        Manual extraction position(s), as one or more semicolon-separated
        ``x:y`` pairs, e.g. ``10.0:14.0`` or ``10.0:14.0;20.0:25.0``. If
        None, no validation is performed.

    Returns
    -------
    :obj:`str`, None
        The stripped, validated position string, or None if ``manual`` is
        None.

    Raises
    ------
    PypeItError
        If ``manual`` contains a comma, if any ``;``-separated entry is not
        a single ``:``-separated pair, or if either field of a pair is not
        numeric.
    """
    if manual is None:
        return None

    _manual = manual.strip()
    if ',' in _manual:
        raise PypeItError(
            'Manual datacube extraction positions must use colon-separated x:y '
            f'syntax, not commas: {_manual}. For example, use --manual 10.0:14.0.'
        )

    positions = _manual.split(';')
    for pos in positions:
        fields = pos.split(':')
        if len(fields) != 2:
            raise PypeItError(
                'Manual datacube extraction positions must be provided as x:y. '
                f'Invalid value: {pos}.'
            )
        try:
            float(fields[0])
            float(fields[1])
        except ValueError:
            raise PypeItError(
                'Manual datacube extraction positions must be numeric x:y values. '
                f'Invalid value: {pos}.'
            )

    return _manual


def write_extract_file(extract_file, target_stub, whitelight_range, manual=None, fwhm=1.1,
                       snr_thresh=4.0):
    """
    Write a .extract file for point-source datacube extraction.

    Parameters
    ----------
    extract_file : :obj:`str`, `Path`_
        Output path for the ``.extract`` file. Any existing file at this
        path is overwritten.
    target_stub : :obj:`str`
        Target name used to set the ``output_filename`` extraction
        parameter (with ``_extract`` appended).
    whitelight_range : :obj:`str`
        Comma-separated white-light wavelength range; see
        :func:`validate_whitelight_range`.
    manual : :obj:`str`, optional
        Manual extraction position(s) in ``x:y`` format; see
        :func:`validate_manual`. If provided, ``opt_prof_method`` is set to
        ``user_gauss`` and a ``manual`` line is written; otherwise
        ``opt_prof_method`` is ``fit_gauss``.
    fwhm : :obj:`float`, optional
        Extraction FWHM in arcsec.
    snr_thresh : :obj:`float`, optional
        S/N threshold for automatic source finding.
    """
    manual = validate_manual(manual)
    opt_prof_method = 'fit_gauss' if manual is None else 'user_gauss'
    lines = [
        '# User-defined execution parameters\n',
        '[reduce]\n',
        '    [[cube]]\n',
        '     [[[extraction]]]\n',
        f'       output_filename = {target_stub}_extract\n',
        f'       whitelight_range = {whitelight_range}\n',
        f'       snr_thresh = {snr_thresh}\n',
        f'       opt_prof_method = {opt_prof_method}\n',
        f'       fwhm = {fwhm}\n',
    ]
    if manual is not None:
        lines.append(f'       manual = {manual}\n')
    Path(extract_file).write_text(''.join(lines))


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
            '--whitelight_range', '--wl_range', dest='whitelight_range', type=str,
            default='None,None',
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

        wl_range = validate_whitelight_range(args.whitelight_range)
        sensfile = None
        if args.sensfile is not None:
            sensfile = Path(args.sensfile).expanduser().absolute()
            if not sensfile.is_file():
                raise PypeItError(f'Sensitivity function does not exist: {sensfile}')

        sci_dir = science_directory(pypeit_path)
        if not sci_dir.is_dir():
            raise PypeItError(f'Expected Science directory does not exist: {sci_dir}')

        spec2d_files, missing, target_name = existing_spec2d_files(
            pypeit_file, args.target, sci_dir
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
            new_names = append_spec2d_files(coadd3d_file, spec2d_files)
            if len(new_names) == 0:
                log.info(f'No new spec2d files to append to {coadd3d_file}.')
            else:
                log.info(f'Appended {len(new_names)} spec2d file(s) to {coadd3d_file}:')
                for name in new_names:
                    log.info(f'  {name}')
            if extract_file.exists():
                log.info(f'Leaving existing extract file unchanged: {extract_file}')
            return

        write_coadd3d_file(
            coadd3d_file, spectrograph, target_stub, sci_dir, spec2d_files, wl_range, sensfile,
            spatial_delta=args.spatial_delta, det=args.det
        )
        log.info(f'Wrote {coadd3d_file}')

        if extract_file.exists() and not args.overwrite:
            log.warning(f'{extract_file} exists; leaving it unchanged. Use -o/--overwrite to replace it.')
        else:
            write_extract_file(
                extract_file, target_stub, wl_range, manual=args.manual, fwhm=args.fwhm,
                snr_thresh=args.snr_thresh
            )
            log.info(f'Wrote {extract_file}')
