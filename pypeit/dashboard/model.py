"""
Headless (Qt-free) data/logic layer for the PypeIt Dashboard.

This module holds the dashboard's PypeIt-facing logic so it can be unit
tested with plain ``pytest`` (no display required), keeping the Qt views
thin.  Stage 0 only needs the cheap reduction metadata used by the header
banner; the full reduction-state layer (loading/deriving
``RunPypeItState``) is added in Stage 1.
"""

import json
from pathlib import Path
from dataclasses import dataclass

from pydantic import ValidationError

from pypeit import log
from pypeit.state.run_state import RunPypeItState, same_det, state_table

# Map the spectrograph ``pypeline`` attribute to the reduction-path label
# shown in the header.  Internally PypeIt names the IFU path 'SlicerIFU';
# the dashboard displays the shorter, user-facing 'IFU'.
PYPELINE_DISPLAY = {
    'MultiSlit': 'MultiSlit',
    'Echelle': 'Echelle',
    'SlicerIFU': 'IFU',
}


@dataclass
class HeaderInfo:
    """
    Cheap, no-reduction metadata describing a reduction, used to populate
    the dashboard's shared header banner (design requirement R6).

    Attributes:
        pypeit_file (str):
            Name of the ``.pypeit`` file (basename, for display).
        spectrograph (str):
            The spectrograph name (``PYP_SPEC``), e.g. ``shane_kast_blue``.
        setup (str):
            The setup/configuration identifier, e.g. ``A``.  May be ``None``
            if the ``.pypeit`` file does not define one.
        path (str):
            The reduction path label: ``MultiSlit``, ``Echelle``, or ``IFU``.
        redux_dir (str):
            The reduction directory (absolute path).
    """
    pypeit_file: str
    spectrograph: str
    setup: str
    path: str
    redux_dir: str


def read_header_info(pypeit_file, redux_path=None):
    """
    Parse a ``.pypeit`` file and return the cheap metadata needed by the
    header banner, performing no reduction.

    Args:
        pypeit_file (:obj:`str`, :obj:`pathlib.Path`):
            Path to the ``.pypeit`` reduction file.
        redux_path (:obj:`str`, optional):
            Reduction directory.  If ``None``, the directory containing
            ``pypeit_file`` is used.

    Returns:
        :class:`HeaderInfo`: The parsed header metadata.

    Raises:
        FileNotFoundError:
            If ``pypeit_file`` does not exist (R11 edge case: the named
            ``.pypeit`` file is absent).
    """
    header_info, _ = _parse_pypeit(pypeit_file, redux_path=redux_path)
    return header_info


def _parse_pypeit(pypeit_file, redux_path=None):
    """
    Parse a ``.pypeit`` file once, returning both the header metadata and
    the instantiated spectrograph (the latter is reused for, e.g.,
    detector-name formatting).

    Args:
        pypeit_file (:obj:`str`, :obj:`pathlib.Path`):
            Path to the ``.pypeit`` reduction file.
        redux_path (:obj:`str`, optional):
            Reduction directory.  If ``None``, the directory containing
            ``pypeit_file`` is used.

    Returns:
        tuple: ``(HeaderInfo, Spectrograph)``.

    Raises:
        FileNotFoundError:
            If ``pypeit_file`` does not exist (R11 edge case).
    """
    # Import here to avoid importing the heavy pypeit machinery at module
    # load time (keeps GUI startup snappy and module imports cheap).
    from pypeit.inputfiles import PypeItFile

    path = Path(pypeit_file)
    if not path.is_file():
        # R11: surface a clear, typed error rather than crashing later.
        raise FileNotFoundError(f'PypeIt file not found: {path}')

    # vet=False: we only need the configuration metadata for the header, not
    # a full validation of the data block (that is run_pypeit's job).
    pyfile = PypeItFile.from_file(str(path), vet=False)
    spec = pyfile.get_spectrograph()

    # Map the internal pypeline name to the user-facing path label.
    path_label = PYPELINE_DISPLAY.get(spec.pypeline, spec.pypeline)

    redux_dir = Path(redux_path) if redux_path is not None else path.parent

    header_info = HeaderInfo(pypeit_file=path.name,
                             spectrograph=spec.name,
                             setup=pyfile.setup_name,
                             path=path_label,
                             redux_dir=str(redux_dir.resolve()))
    return header_info, spec


# How the reduction state was obtained / why it is unavailable.  Used by the
# views (Stage 2) to choose what to render (R11).
LOAD_STATE_FILE = 'state_file'        # loaded <root>_state.json (R4)
LOAD_DERIVED = 'derived'              # derived via pypeit_status path (R5)
LOAD_NOT_STARTED = 'not_started'      # no state file and not derived
LOAD_MALFORMED = 'malformed'          # state file present but unreadable
LOAD_FILE_NOT_FOUND = 'file_not_found'  # the .pypeit file is missing
LOAD_ERROR = 'error'                  # deriving the state raised

# The normalized status-table columns the views consume.
STATUS_COLUMNS = ['calibration_group', 'detector', 'step', 'in_pipeline',
                  'required', 'status', 'output_file']

# The per-frame science-status columns (mirrors get_science_status); the
# Science view (Stage 6) consumes these.
SCIENCE_COLUMNS = ['frame', 'detector', 'calib', 'objtype', 'process',
                   'findobj', 'skysub', 'extract', 'nobj', 'spec2d', 'spec1d']

# Sentinel used by RunPypeItState.get_status() for an absent entry, mapped
# here to an explicit "not present" status.
_ABSENT = 'absent'

# Per-``.pypeit`` cache of the *planned* science/standard frames (computed once
# from the metadata; see DashboardModel._planned_science_frames).  Caching
# means the (expensive) metadata read happens once per session, so re-seeding
# planned frames on every state reload (live monitoring, post-(re)build) is
# cheap (Round-4).  Keyed by the resolved ``.pypeit`` path.
_PLANNED_SCIENCE_CACHE = {}


class DashboardModel:
    """
    Headless (Qt-free) model for one PypeIt reduction.

    It acquires a :class:`pypeit.state.run_state.RunPypeItState` by source priority
    (load ``<root>_state.json`` if present, R4; otherwise derive it the
    ``pypeit_status`` way, R5), and exposes a clean, normalized API the Qt
    views consume — a status table, the ``(calib_id, det)`` pairs, the
    path-aware step order, and graceful edge states (R11).  It never raises
    for a missing/empty/malformed state: such cases are reported via
    :attr:`load_status`.

    Args:
        pypeit_file (:obj:`str`, :obj:`pathlib.Path`):
            Path to the ``.pypeit`` reduction file.
        redux_path (:obj:`str`, optional):
            Reduction directory (where the state file lives).  Defaults to
            the directory containing ``pypeit_file``.
        derive (:obj:`bool`, optional):
            If True (default) and no state file is present, derive the state
            via :class:`~pypeit.pypeit.PypeIt` (requires the raw data).  Set
            False to skip deriving (e.g. CI tests without RAW_DATA): the
            model then reports ``not_started`` when no state file exists.

    Attributes:
        pypeit_file (:obj:`pathlib.Path`): The ``.pypeit`` path.
        redux_dir (:obj:`pathlib.Path`): The reduction directory.
        state_path (:obj:`pathlib.Path`): The expected ``<root>_state.json`` path.
        header_info (:class:`HeaderInfo`): Header metadata, or ``None``.
        run_state (:class:`~pypeit.state.run_state.RunPypeItState`): The state, or
            ``None`` when unavailable.
        load_status (str): One of the ``LOAD_*`` constants.
        error (str): An error message when relevant, else ``None``.
    """

    def __init__(self, pypeit_file, redux_path=None, derive=True):
        self.pypeit_file = Path(pypeit_file)
        self.redux_dir = Path(redux_path) if redux_path is not None \
            else self.pypeit_file.parent
        # The state file PypeIt writes: <pypeit_root>_state.json (see
        # RunPypeItState.outfile), located in the reduction directory.
        self.state_path = self.redux_dir / \
            f'{self.pypeit_file.stem}_state.json'

        self.header_info = None
        self.spectrograph = None
        self.run_state = None
        self.error = None

        # Header metadata + spectrograph (cheap; the spectrograph also gives
        # the pipeline and detector names).  A missing .pypeit file is an
        # edge state, not a crash (R11).
        try:
            self.header_info, self.spectrograph = \
                _parse_pypeit(self.pypeit_file, redux_path=redux_path)
        except FileNotFoundError as e:
            self.load_status = LOAD_FILE_NOT_FOUND
            self.error = str(e)
            return

        self._acquire_state(derive=derive)

    def det_name(self, det):
        """
        Return a human-readable name for a detector or mosaic, for display
        in the views (drop-downs, navigator).

        Args:
            det (:obj:`int`, :obj:`tuple`, :obj:`list`):
                The detector (int) or mosaic (tuple/list of ints), as stored
                in the state.

        Returns:
            str: The PypeIt detector name (e.g. ``DET01``, ``MSC01``) when it
            can be resolved, else a plain readable fallback (``Det 1`` /
            ``Mosaic (1, 5)``).
        """
        # Mosaics are stored as lists by pydantic; get_det_name wants a tuple.
        key = tuple(det) if isinstance(det, (list, tuple)) else det
        if self.spectrograph is not None:
            try:
                return self.spectrograph.get_det_name(key)
            except Exception:
                # Fall through to a readable, never-failing label.
                pass
        if isinstance(key, tuple):
            return f'Mosaic {key}'
        return f'Det {key}'

    def _acquire_state(self, derive):
        """
        Acquire the reduction state by source priority (R4 then R5).

        Args:
            derive (:obj:`bool`):
                Whether to derive the state when no state file is present.

        Returns:
            None.  Sets :attr:`run_state` and :attr:`load_status`.
        """
        # R4: a state file present is the fast, authoritative source.
        if self.state_path.is_file():
            try:
                with open(self.state_path, 'rt') as fh:
                    data = json.load(fh)
                self.run_state = RunPypeItState.model_validate(data)
                self.load_status = LOAD_STATE_FILE
                # Re-seed the planned science/standard frames (Round-4): a
                # calibration (re)build rewrites the state file with no science
                # entries, so without this the Science view would empty out
                # after building calibrations.  Cheap once the planned-frame
                # list is cached; a one-time build is allowed only when the
                # caller permits deriving (launch, not a live-monitor reload).
                self._seed_planned_science(allow_build=derive)
            except (json.JSONDecodeError, ValidationError, OSError) as e:
                # Malformed / partial / schema-mismatched state file (R11).
                self.run_state = None
                self.load_status = LOAD_MALFORMED
                self.error = str(e)
                log.warning(f'Could not read reduction state '
                            f'{self.state_path}: {e}')
            return

        # No state file.  Optionally derive it the pypeit_status way (R5).
        if not derive:
            self.load_status = LOAD_NOT_STARTED
            return
        try:
            # Heavy imports kept local so module import stays cheap.
            from pypeit import pypeit
            pypeIt = pypeit.PypeIt(str(self.pypeit_file), reuse_calibs=True,
                                   calib_only=True)
            pypeIt.calib_all(status_only=True, reload_only=True)
            # Also derive science-frame status from the on-disk products
            # (Stage 6, S6-Q7): best-effort — never let it block the
            # calibration status the dashboard primarily needs.
            try:
                from pypeit.state import science_status
                science_status.derive_science_from_disk(
                    pypeIt.run_state, str(self.redux_dir),
                    fitstbl=pypeIt.fitstbl)
            except Exception as e:
                log.warning(f'Could not derive science state for '
                            f'{self.pypeit_file}: {e}')
            self.run_state = pypeIt.run_state
            self.load_status = LOAD_DERIVED
            # Cache the planned science/standard frames from the fitstbl we
            # just built (free here), then seed them so the Science view lists
            # upcoming frames before they are reduced (Round-3 #2), mirroring
            # the planned calibrations.  The cache makes later state reloads
            # cheap (Round-4).
            try:
                from pypeit.state import science_status
                _PLANNED_SCIENCE_CACHE[str(self.pypeit_file)] = \
                    science_status.planned_science_from_fitstbl(pypeIt.fitstbl)
                self._seed_planned_science()
            except Exception as e:
                log.warning(f'Could not seed planned science frames for '
                            f'{self.pypeit_file}: {e}')
        except Exception as e:
            # Deriving needs the raw data and a full PypeIt setup; if that
            # fails, report it rather than crash the GUI (R11).
            self.load_status = LOAD_ERROR
            self.error = str(e)
            log.warning(f'Could not derive reduction state for '
                        f'{self.pypeit_file}: {e}')

    def _planned_science_frames(self, allow_build=False):
        """
        Return the cached list of planned science/standard frames for this
        ``.pypeit`` file (Round-4).  On a cache miss, build the metadata to
        compute it **only** when ``allow_build`` is True (so live-monitoring
        reloads and CI loads never trigger a heavy, possibly-failing build).

        Args:
            allow_build (:obj:`bool`, optional): Permit a one-time metadata
                build on a cache miss (True at launch; False on refresh/CI).

        Returns:
            list: Planned-frame dicts (see
            :func:`~pypeit.state.science_status.planned_science_from_fitstbl`);
            empty when uncached and a build is not permitted/possible.
        """
        key = str(self.pypeit_file)
        if key in _PLANNED_SCIENCE_CACHE:
            return _PLANNED_SCIENCE_CACHE[key]
        if not allow_build:
            return []
        from pypeit.state import science_status
        planned = []
        try:
            # Cold cache + a build is allowed (e.g. relaunch with a state file
            # already present): build the metadata once to learn the planned
            # frames, then cache it for the rest of the session.
            from pypeit import pypeit
            pypeIt = pypeit.PypeIt(str(self.pypeit_file), reuse_calibs=True,
                                   calib_only=True)
            planned = science_status.planned_science_from_fitstbl(
                pypeIt.fitstbl)
        except Exception as e:
            log.warning(f'Could not determine planned science frames for '
                        f'{self.pypeit_file}: {e}')
        _PLANNED_SCIENCE_CACHE[key] = planned
        return planned

    def _seed_planned_science(self, allow_build=False):
        """
        Seed the planned science/standard frames into :attr:`run_state` (from
        the cache, or a one-time build when ``allow_build``), keyed to the
        calibration ``(group, det)`` pairs in the current state (Round-3 #2 /
        Round-4).  Best-effort: never raises into the load/derive path.

        Args:
            allow_build (:obj:`bool`, optional): Permit a one-time metadata
                build on a cache miss.

        Returns:
            None.
        """
        if self.run_state is None or self.header_info is None:
            return
        try:
            from pypeit.state import science_status
            planned = self._planned_science_frames(allow_build=allow_build)
            if not planned:
                return
            group_dets = {}
            for g, d in self.calib_det_pairs():
                group_dets.setdefault(g, []).append(d)
            science_status.seed_planned_science_entries(
                self.run_state, planned, group_dets)
        except Exception as e:
            log.warning(f'Could not seed planned science frames for '
                        f'{self.pypeit_file}: {e}')

    def is_started(self):
        """
        Whether the reduction has any calibration-state entries.

        Returns:
            bool: True if a state with at least one entry is available.
        """
        return self.run_state is not None \
            and self.run_state.get_status() is not None

    def default_steps(self):
        """
        Return the active spectrograph's calibration steps, in pipeline
        order (path-aware: MultiSlit/Echelle vs IFU), including ``bpm``.

        Returns:
            list: The ordered step names, or ``[]`` if the pipeline is
            unknown (e.g. the ``.pypeit`` file was not found).
        """
        if self.header_info is None:
            return []
        # Map the pipeline label to the Calibrations subclass that defines
        # the step order (no PypeIt instance / RAW_DATA needed).
        from pypeit.calibrations import (MultiSlitCalibrations,
                                         IFUCalibrations)
        if self.header_info.path == 'IFU':
            return IFUCalibrations.default_steps()
        return MultiSlitCalibrations.default_steps()

    def step_order(self, include_bpm=False):
        """
        Return the path-aware step order for the calibration button row.

        Args:
            include_bpm (:obj:`bool`, optional):
                Include the ``bpm`` step.  Default False — ``bpm`` runs
                internally but has no standalone output/QA, so it is omitted
                from the button row (design C13).

        Returns:
            list: The ordered step names.
        """
        steps = self.default_steps()
        if include_bpm:
            return steps
        return [s for s in steps if s != 'bpm']

    def status_table(self):
        """
        Return the normalized calibration-status table the views consume.

        Built on :meth:`pypeit.state.run_state.RunPypeItState.get_status` but
        normalized: ``required`` becomes a real :obj:`bool` (or ``None``),
        the ``"--"`` sentinels become ``None``/``"absent"``, an
        ``in_pipeline`` column is added (membership in
        :meth:`default_steps`), and an empty/unavailable state yields an
        empty table.

        Returns:
            `astropy.table.Table`_: Columns :data:`STATUS_COLUMNS`; empty when
            the reduction has not started or the state is unavailable.
        """
        if self.run_state is None:
            return state_table([], names=STATUS_COLUMNS)
        raw = self.run_state.get_status()
        if raw is None:
            return state_table([], names=STATUS_COLUMNS)

        in_pipeline = set(self.default_steps())
        # Map the stringified "required" to a real bool (or None).
        required_map = {'True': True, 'False': False}
        rows = []
        for row in raw:
            rows.append({
                'calibration_group': row['calibration_group'],
                'detector': row['detector'],
                'step': row['steps'],
                'in_pipeline': row['steps'] in in_pipeline,
                'required': required_map.get(str(row['required']), None),
                'status': _ABSENT if row['status'] == '--' else row['status'],
                'output_file': None if row['output_file'] == '--'
                else row['output_file'],
            })
        return state_table(rows, names=STATUS_COLUMNS)

    def calib_det_pairs(self):
        """
        Enumerate the ``(calibration_group, detector)`` pairs present in
        the state, in sorted order, for the scope drop-downs/selectors
        (R16 / C2).

        Returns:
            list: ``(calib_id, det)`` tuples; ``det`` is left in its raw
            form (an :obj:`int` or a mosaic tuple) for the view to format.
        """
        table = self.status_table()
        seen = []
        for row in table:
            pair = (row['calibration_group'], row['detector'])
            if pair not in seen:
                seen.append(pair)
        return sorted(seen, key=lambda p: (p[0], str(p[1])))

    def calibrations_ready(self, group, det):
        """
        Whether the **required** calibrations for one ``(group, det)`` are all
        built successfully — the precondition for a science (Re)Build
        (Round-3 #2): until the calibrations a science frame depends on exist,
        its (Re)Build is disabled.

        Args:
            group: The calibration group ID of the science frame.
            det: The detector (int) or mosaic the science frame is on.

        Returns:
            bool: True only if at least one required, in-pipeline calibration
            step exists for ``(group, det)`` and **all** such steps have a
            ``success``/``complete`` status; False otherwise (including when no
            matching calibration entries are known).
        """
        table = self.status_table()
        # The required, in-pipeline calibration entries for this (group, det).
        req = [row for row in table
               if row['calibration_group'] == group
               and same_det(row['detector'], det)
               and row['required'] is True and row['in_pipeline']]
        if len(req) == 0:
            return False
        return all(row['status'] in ('success', 'complete') for row in req)

    def is_stale(self):
        """
        Whether the loaded ``*_state.json`` looks out of date relative to
        the ``.pypeit`` file or the calibration outputs (R13).

        Returns:
            bool: True if the state file's mtime is older than the
            ``.pypeit`` file or the newest file in ``Calibrations/``; False
            if there is no state file or nothing newer.
        """
        if not self.state_path.is_file():
            return False
        try:
            state_mtime = self.state_path.stat().st_mtime
        except OSError:
            return False
        newest = 0.0
        if self.pypeit_file.is_file():
            newest = max(newest, self.pypeit_file.stat().st_mtime)
        calib_dir = self.redux_dir / 'Calibrations'
        if calib_dir.is_dir():
            for entry in calib_dir.iterdir():
                try:
                    newest = max(newest, entry.stat().st_mtime)
                except OSError:
                    pass
        return newest > state_mtime

    # ------------------------------------------------------------------
    # Stage 3: per-step detail accessors (for the Calibrations view).
    # ------------------------------------------------------------------

    @property
    def calib_dir(self):
        """
        The reduction's ``Calibrations/`` directory (where output files
        live).

        Returns:
            :obj:`pathlib.Path`: ``<redux_dir>/Calibrations``.
        """
        return self.redux_dir / 'Calibrations'

    def step_entry(self, step, group, det):
        """
        Return the raw ``RunPypeItState`` entry for one
        ``(step, group, det)``, or ``None`` if absent.

        Args:
            step (:obj:`str`): Calibration step name (e.g. ``wv_calib``).
            group (:obj:`int`): Calibration group ID.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            The matching pydantic step entry, or ``None``.
        """
        if self.run_state is None or not hasattr(self.run_state, step):
            return None
        for item in getattr(self.run_state, step):
            if item.calib_id == group and same_det(item.det, det):
                return item
        return None

    def output_path(self, step, group, det):
        """
        Full path to a step's processed output file, or ``None``.

        Args:
            step (:obj:`str`): Calibration step name.
            group (:obj:`int`): Calibration group ID.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            :obj:`pathlib.Path` or ``None``: ``Calibrations/<output_file>``.
        """
        entry = self.step_entry(step, group, det)
        if entry is None or not entry.output_file:
            return None
        # output_file is stored as a basename; join under Calibrations/.
        return self.calib_dir / Path(entry.output_file).name

    @property
    def log_path(self):
        """
        The reduction ``.log`` file that ``run_pypeit`` /
        ``pypeit_run_to_calibstep`` write to (used by the single-run lock to
        detect an active run via the log's modification time, design X1).

        Returns:
            :obj:`pathlib.Path`: ``<pypeit_root>.log`` (next to the ``.pypeit``
            file, matching ``pypeit_run_to_calibstep``'s default log path).
        """
        return self.pypeit_file.with_suffix('.log')

    def step_output_files(self, step, group, det):
        """
        The existing ``Calibrations/`` file(s) a (re)build of this step would
        **overwrite** — named explicitly in the clobber confirmation (X2/X3,
        S4-Q4).

        Most steps write a single ``output_file``; ``slits`` writes **both**
        ``Slits_*`` (its ``output_file``) and the companion ``Edges_*``.  Only
        files that currently exist on disk are returned (a step with no output
        yet is a fresh build, nothing to overwrite).

        Args:
            step (:obj:`str`): Calibration step name.
            group (:obj:`int`): Calibration group ID.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            list: Existing :obj:`pathlib.Path` outputs (possibly empty).
        """
        candidates = []
        out = self.output_path(step, group, det)
        if out is not None:
            candidates.append(out)
            if step == 'slits':
                # slits also writes the Edges_* file alongside Slits_*; clear
                # both so the rebuild is genuine (handles *_all_* naming too).
                candidates.append(
                    out.with_name(out.name.replace('Slits', 'Edges', 1)))
        return [f for f in candidates if f.exists()]

    def calib_file_path(self, prefix, group, det, ext='fits'):
        """
        Construct the path to a calibration file by naming convention
        (``<prefix>_<setup>_<group>_<detname>.<ext>``) — used when the
        viewer needs a *different* file than a step's ``output_file``
        (e.g. the ``Edges_*`` file for the ``slits`` step, or the
        ``Slits_*`` file required by ``pypeit_chk_scattlight``).

        Args:
            prefix (:obj:`str`): File-type prefix (e.g. ``Edges``, ``Slits``).
            group (:obj:`int`): Calibration group ID.
            det: Detector (int) or mosaic (tuple/list).
            ext (:obj:`str`, optional): Extension (``fits`` or ``fits.gz``).

        Returns:
            :obj:`pathlib.Path` or ``None``: The constructed path (``None`` if the
            setup is unknown).
        """
        if self.header_info is None or self.header_info.setup is None:
            return None
        name = f'{prefix}_{self.header_info.setup}_{group}_' \
               f'{self.det_name(det)}.{ext}'
        return self.calib_dir / name

    def step_metrics(self, step, group, det):
        """
        Return the entry-level metrics for a step, as a plain dict for the
        detail panel (per-slit metrics come from :meth:`slit_table`).

        Args:
            step (:obj:`str`): Calibration step name.
            group (:obj:`int`): Calibration group ID.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            dict: Metric name → value (empty if the entry is absent or has
            no entry-level metrics). ``bias`` → mean/std; ``slits`` →
            nslits; ``flats`` → corrections/pixelflat_source.
        """
        entry = self.step_entry(step, group, det)
        if entry is None:
            return {}
        metrics = {}
        for field in ('mean', 'std', 'nslits', 'corrections',
                      'pixelflat_source'):
            value = getattr(entry, field, None)
            if value is not None and value != []:
                metrics[field] = value
        return metrics

    def slit_table(self, step, group, det):
        """
        Return the per-slit/order rows for a step that tracks them
        (``slits`` / ``wv_calib`` / ``tilts`` / ``flats``).

        Args:
            step (:obj:`str`): Calibration step name.
            group (:obj:`int`): Calibration group ID.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            list: One dict per slit (sorted by slit id), with keys ``slit``,
            ``status``, and step-specific metrics: ``center``/``slitord_id``
            (slits), ``rms`` (wv_calib/tilts), or ``corrections`` (flats — a
            dict ``name → {'mean':…, 'rms':…}``). Empty if the step has no
            per-slit data.
        """
        entry = self.step_entry(step, group, det)
        slits = getattr(entry, 'slits', None) if entry is not None else None
        if not slits:
            return []
        rows = []
        for slit_id in sorted(slits.keys()):
            slit = slits[slit_id]
            row = {'slit': slit_id, 'status': slit.status}
            for field in ('center', 'slitord_id', 'rms'):
                value = getattr(slit, field, None)
                if value is not None:
                    row[field] = value
            corrections = getattr(slit, 'corrections', None)
            if corrections:
                row['corrections'] = {
                    name: {'mean': m.mean, 'rms': m.rms}
                    for name, m in corrections.items()}
            rows.append(row)
        return rows

    # ------------------------------------------------------------------
    # Stage 6: science-frame accessors (for the Science view).
    # ------------------------------------------------------------------

    def science_table(self):
        """
        Return the normalized per-frame science-status table the Science view
        renders (one row per ``(frame, detector)``).

        Returns:
            `astropy.table.Table`_: Columns :data:`SCIENCE_COLUMNS`; empty when
            there are no science entries.
        """
        if self.run_state is None:
            return state_table([], names=SCIENCE_COLUMNS)
        raw = self.run_state.get_science_status()
        if raw is None or len(raw) == 0:
            return state_table([], names=SCIENCE_COLUMNS)
        return raw

    def has_science(self):
        """
        Whether any science-frame state is available.

        Returns:
            bool: True if the state has at least one science entry.
        """
        return len(self.science_table()) > 0

    def science_frame_entry(self, frame, det):
        """
        Return the raw ``ScienceFrameState`` for one ``(frame, det)``, or
        ``None``.

        Args:
            frame (:obj:`str`): Exposure basename.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            The matching pydantic science entry, or ``None``.
        """
        if self.run_state is None:
            return None
        return self.run_state.science_entry(frame, det)

    def science_slit_table(self, frame, det):
        """
        Return the per-slit science rows for one frame (``ScienceSlit``).

        Args:
            frame (:obj:`str`): Exposure basename.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            list: One dict per slit (sorted), keys ``slit``, ``status``,
            ``nobj``. Empty if none.
        """
        entry = self.science_frame_entry(frame, det)
        slits = getattr(entry, 'slits', None) if entry is not None else None
        if not slits:
            return []
        rows = []
        for slit_id in sorted(slits.keys(), key=lambda k: int(k)):
            slit = slits[slit_id]
            rows.append({'slit': slit_id, 'status': slit.status,
                         'nobj': slit.nobj})
        return rows

    def science_object_table(self, frame, det):
        """
        Return the per-object science rows for one frame (``ScienceObj``).

        Args:
            frame (:obj:`str`): Exposure basename.
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            list: One dict per detected object, keys ``objid``, ``slitid``,
            ``spat_pixpos``, ``fwhm``, ``snr_find``, ``s2n``, ``sign``,
            ``extracted``. Empty if none.
        """
        entry = self.science_frame_entry(frame, det)
        objects = getattr(entry, 'objects', None) if entry is not None else None
        if not objects:
            return []
        rows = []
        for obj in objects:
            rows.append({
                'objid': obj.objid, 'slitid': obj.slitid,
                'spat_pixpos': obj.spat_pixpos, 'fwhm': obj.fwhm,
                'snr_find': obj.snr_find, 's2n': obj.s2n,
                'sign': obj.sign, 'extracted': obj.extracted})
        return rows

    def science_qa_files(self, frame, det):
        """
        Return **all** on-disk QA PNGs for one science ``(frame, det)`` (Stage 6
        Round-1 #3, S6-Q15(a)): the per-object ``obj_prof``/``obj_trace`` figures
        and the frame-level ``spec_flex_*`` flexure figures.

        Derived **at view time** by globbing ``<redux>/QA/PNGs/`` for files whose
        name starts with the frame basename and contains the detector name — so
        no science-state change is needed (mirrors the S6-Q10(a) reconstruct-
        from-disk approach) and it works for both live and derived state.

        Args:
            frame (:obj:`str`): Exposure basename (the science-table key).
            det: Detector (int) or mosaic (tuple/list).

        Returns:
            list: Sorted :obj:`pathlib.Path` of the matching QA PNGs (empty if the
            ``QA/PNGs`` directory or matches are absent).
        """
        qa_dir = self.redux_dir / 'QA' / 'PNGs'
        if not qa_dir.is_dir():
            return []
        detname = self.det_name(det)
        return [p for p in sorted(qa_dir.glob('*.png'))
                if p.name.startswith(str(frame)) and detname in p.name]

    def science_object_qa_files(self, frame, det, slitid):
        """
        Return the per-object QA PNGs (``obj_prof`` / ``obj_trace``) for one
        object's slit, so the per-object table can open them (Stage 6
        Round-1 #3, S6-Q15(c)).  The science QA PNGs encode the slit as
        ``S{spat_id:04d}``, matching ``ScienceObj.slitid``.

        Args:
            frame (:obj:`str`): Exposure basename.
            det: Detector (int) or mosaic (tuple/list).
            slitid (:obj:`int`): The object's slit (SPAT_ID).

        Returns:
            dict: ``{'obj_prof': Path|None, 'obj_trace': Path|None}``.
        """
        out = {'obj_prof': None, 'obj_trace': None}
        if slitid is None:
            return out
        tag = f'S{int(slitid):04d}'
        for path in self.science_qa_files(frame, det):
            if tag not in path.name:
                continue
            if 'obj_prof' in path.name:
                out['obj_prof'] = path
            elif 'obj_trace' in path.name:
                out['obj_trace'] = path
        return out
