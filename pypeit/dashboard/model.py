"""
Headless (Qt-free) data/logic layer for the PypeIt Dashboard.

This module holds the dashboard's PypeIt-facing logic so it can be unit
tested with plain ``pytest`` (no display required), keeping the Qt views
thin.  It provides the cheap reduction metadata used by the header banner
and the full reduction-state layer (loading or deriving a
:class:`~pypeit.state.run_state.RunPypeItState`).

.. include:: ../include/links.rst
"""

import json
from pathlib import Path
from dataclasses import dataclass

from pydantic import ValidationError

from pypeit import log, PypeItError
from pypeit import pypeit
from pypeit.calibrations import Calibrations, MultiSlitCalibrations, IFUCalibrations
from pypeit.inputfiles import PypeItFile
from pypeit.state import science_status
from pypeit.state.run_state import RunPypeItState, same_det, state_table


@dataclass
class HeaderInfo:
    """
    Cheap, no-reduction metadata describing a reduction, used to populate
    the dashboard's shared header banner.

    Attributes
    ----------
    pypeit_file : :obj:`str`
        Name of the ``.pypeit`` file (basename, for display).
    spectrograph : :obj:`str`
        The spectrograph name (``PYP_SPEC``), e.g. ``shane_kast_blue``.
    setup : :obj:`str`
        The setup/configuration identifier, e.g. ``A``.  May be ``None``
        if the ``.pypeit`` file does not define one.
    path : :obj:`str`
        The reduction path: the spectrograph ``pypeline`` value, shown
        verbatim (e.g. ``MultiSlit``, ``Echelle``, ``SlicerIFU``).
    redux_dir : :obj:`str`
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

    Parameters
    ----------
    pypeit_file : :obj:`str`, :obj:`pathlib.Path`
        Path to the ``.pypeit`` reduction file.
    redux_path : :obj:`str`, optional
        Reduction directory.  If ``None``, the directory containing
        ``pypeit_file`` is used.

    Returns
    -------
    :class:`HeaderInfo`
        The parsed header metadata.

    Raises
    ------
    FileNotFoundError
        If ``pypeit_file`` does not exist.
    """
    header_info, _ = _parse_pypeit(pypeit_file, redux_path=redux_path)
    return header_info


def _parse_pypeit(pypeit_file, redux_path=None):
    """
    Parse a ``.pypeit`` file once, returning both the header metadata and
    the instantiated spectrograph (the latter is reused for, e.g.,
    detector-name formatting).

    Parameters
    ----------
    pypeit_file : :obj:`str`, :obj:`pathlib.Path`
        Path to the ``.pypeit`` reduction file.
    redux_path : :obj:`str`, optional
        Reduction directory.  If ``None``, the directory containing
        ``pypeit_file`` is used.

    Returns
    -------
    header_info : :class:`HeaderInfo`
        The parsed header metadata.
    spec : :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
        The instantiated spectrograph.

    Raises
    ------
    FileNotFoundError
        If ``pypeit_file`` does not exist.
    """
    path = Path(pypeit_file)
    if not path.is_file():
        # Surface a clear, typed error rather than crashing later.
        raise FileNotFoundError(f'PypeIt file not found: {path}')

    # vet=False: we only need the configuration metadata for the header, not
    # a full validation of the data block (that is run_pypeit's job).
    pyfile = PypeItFile.from_file(str(path), vet=False)
    spec = pyfile.get_spectrograph()

    redux_dir = Path(redux_path) if redux_path is not None else path.parent

    header_info = HeaderInfo(pypeit_file=path.name,
                             spectrograph=spec.name,
                             setup=pyfile.setup_name,
                             path=spec.pypeline,
                             redux_dir=str(redux_dir.resolve()))
    return header_info, spec


class DashboardModel:
    """
    Headless (Qt-free) model for one PypeIt reduction.

    It acquires a :class:`~pypeit.state.run_state.RunPypeItState` by source
    priority (load ``<root>_state.json`` if present; otherwise derive it
    from the on-disk calibration products, the same way ``pypeit_status``
    does), and exposes a clean, normalized API the Qt views consume — a
    status table, the ``(calib_id, det)`` pairs, the pipeline-aware step
    order, and graceful edge states.  It never raises for a
    missing/empty/malformed state: such cases are reported via
    :attr:`load_status`.

    Parameters
    ----------
    pypeit_file : :obj:`str`, :obj:`pathlib.Path`
        Path to the ``.pypeit`` reduction file.
    redux_path : :obj:`str`, optional
        Reduction directory (where the state file lives).  Defaults to
        the directory containing ``pypeit_file``.
    derive : :obj:`bool`, optional
        If True (default) and no state file is present, derive the state
        via :class:`~pypeit.pypeit.PypeIt` (requires the raw data).  Set
        False to skip deriving (e.g. CI tests without RAW_DATA): the
        model then reports ``not_started`` when no state file exists.

    Attributes
    ----------
    pypeit_file : :obj:`pathlib.Path`
        The ``.pypeit`` path.
    redux_dir : :obj:`pathlib.Path`
        The reduction directory.
    state_path : :obj:`pathlib.Path`
        The expected ``<root>_state.json`` path.
    header_info : :class:`HeaderInfo`
        Header metadata, or ``None``.
    run_state : :class:`~pypeit.state.run_state.RunPypeItState`
        The state, or ``None`` when unavailable.
    load_status : :obj:`str`
        One of the ``LOAD_*`` class attributes.
    error : :obj:`str`
        An error message when relevant, else ``None``.
    """

    # -- How the reduction state was obtained / why it is unavailable.
    # Used by the views to choose what to render.
    #: The state was loaded from an existing ``<root>_state.json`` file.
    LOAD_STATE_FILE = 'state_file'
    #: The state was derived from the on-disk products (no state file).
    LOAD_DERIVED = 'derived'
    #: No state file exists and the state was not derived.
    LOAD_NOT_STARTED = 'not_started'
    #: A state file is present but could not be read/validated.
    LOAD_MALFORMED = 'malformed'
    #: The named ``.pypeit`` file is missing.
    LOAD_FILE_NOT_FOUND = 'file_not_found'
    #: Deriving the state raised an error.
    LOAD_ERROR = 'error'

    #: The normalized status-table columns the views consume.
    STATUS_COLUMNS = ['calibration_group', 'detector', 'step', 'in_pipeline',
                      'required', 'status', 'output_file']

    #: The per-frame science-status columns consumed by the Science view
    #: (mirrors :meth:`~pypeit.state.run_state.RunPypeItState.get_science_status`).
    SCIENCE_COLUMNS = ['frame', 'detector', 'calib', 'objtype', 'process',
                       'findobj', 'skysub', 'extract', 'nobj', 'spec2d',
                       'spec1d']

    #: Status value used for an entry that is not present in the state
    #: (:meth:`~pypeit.state.run_state.RunPypeItState.get_status` reports
    #: these with a ``"--"`` sentinel).
    _ABSENT = 'absent'

    #: The :class:`~pypeit.calibrations.Calibrations` subclass that defines
    #: the calibration-step order, keyed on the spectrograph ``pypeline``
    #: value.  New pipeline types should be added here explicitly; an
    #: unknown pypeline falls back to the full set of known steps (see
    #: :meth:`default_steps`).
    _CALIB_CLASS_BY_PYPELINE = {
        'MultiSlit': MultiSlitCalibrations,
        # Echelle currently shares the MultiSlit step order (there is no
        # dedicated EchelleCalibrations class); revisit if that changes.
        'Echelle': MultiSlitCalibrations,
        'SlicerIFU': IFUCalibrations,
    }

    #: Per-``.pypeit`` cache (class-level, so it spans model reloads within a
    #: session) of the *planned* science/standard frames, keyed by the
    #: ``.pypeit`` path.  Computing the planned frames requires building the
    #: reduction metadata, which is expensive; caching it means the build
    #: happens at most once per session, so re-seeding planned frames on
    #: every state reload (live monitoring, post-(re)build) is cheap.
    _PLANNED_SCIENCE_CACHE = {}

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
        # edge state the views render, not a crash.
        try:
            self.header_info, self.spectrograph = \
                _parse_pypeit(self.pypeit_file, redux_path=redux_path)
        except FileNotFoundError as e:
            self.load_status = self.LOAD_FILE_NOT_FOUND
            self.error = str(e)
            return

        self._acquire_state(derive=derive)

    def det_name(self, det):
        """
        Return a human-readable name for a detector or mosaic, for display
        in the views (drop-downs, navigator).

        Parameters
        ----------
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            The detector (int) or mosaic (tuple/list of ints), as stored
            in the state.

        Returns
        -------
        :obj:`str`
            The PypeIt detector name (e.g. ``DET01``, ``MSC01``) when it
            can be resolved, else a plain readable fallback (``Det 1`` /
            ``Mosaic (1, 5)``).
        """
        # Mosaics are stored as lists by pydantic; get_det_name wants a tuple.
        key = tuple(det) if isinstance(det, (list, tuple)) else det
        if self.spectrograph is not None:
            try:
                return self.spectrograph.get_det_name(key)
            except (PypeItError, TypeError):
                # get_det_name raises PypeItError for a detector/mosaic the
                # spectrograph does not define (TypeError for a malformed
                # det); fall through to a readable label.
                pass
        if isinstance(key, tuple):
            return f'Mosaic {key}'
        return f'Det {key}'

    def _acquire_state(self, derive):
        """
        Acquire the reduction state by source priority: an existing state
        file first, otherwise (optionally) derive it from the on-disk
        products.  Sets :attr:`run_state` and :attr:`load_status`.

        Parameters
        ----------
        derive : :obj:`bool`
            Whether to derive the state when no state file is present.
        """
        if self.state_path.is_file():
            # A state file present is the fast, authoritative source.
            self._load_state_file(derive=derive)
        elif not derive:
            self.load_status = self.LOAD_NOT_STARTED
        else:
            self._derive_state()

    def _load_state_file(self, derive):
        """
        Load and validate the ``<root>_state.json`` state file, reporting a
        malformed/unreadable file via :attr:`load_status` instead of
        raising.

        Parameters
        ----------
        derive : :obj:`bool`
            Whether a one-time metadata build is permitted when seeding
            the planned science frames (see :meth:`_seed_planned_science`).
        """
        try:
            with open(self.state_path, 'rt') as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            self._set_malformed(e)
            return

        try:
            self.run_state = RunPypeItState.model_validate(data)
        except ValidationError as e:
            self._set_malformed(e)
        else:
            self.load_status = self.LOAD_STATE_FILE
            # Re-seed the planned science/standard frames: a calibration
            # (re)build rewrites the state file with no science entries, so
            # without this the Science view would empty out after building
            # calibrations.  Cheap once the planned-frame list is cached; a
            # one-time build is allowed only when the caller permits
            # deriving (launch, not a live-monitor reload).
            self._seed_planned_science(allow_build=derive)

    def _set_malformed(self, error):
        """
        Record a malformed/partial/schema-mismatched state file.

        Parameters
        ----------
        error : Exception
            The error raised while reading/validating the state file.
        """
        self.run_state = None
        self.load_status = self.LOAD_MALFORMED
        self.error = str(error)
        log.warning(f'Could not read reduction state {self.state_path}: '
                    f'{error}')

    def _derive_state(self):
        """
        Derive the reduction state from the on-disk calibration and science
        products (used when no state file exists).  Failures are reported
        via :attr:`load_status`, never raised.
        """
        try:
            pypeIt = pypeit.PypeIt(str(self.pypeit_file), reuse_calibs=True,
                                   calib_only=True)
            pypeIt.calib_all(status_only=True, reload_only=True)
        except Exception as e:
            # Deliberate never-crash-the-GUI guard: deriving needs the raw
            # data and a full PypeIt setup, either of which can fail in many
            # ways; report the failure as an edge state instead.
            self.load_status = self.LOAD_ERROR
            self.error = str(e)
            log.warning(f'Could not derive reduction state for '
                        f'{self.pypeit_file}: {e}')
            return

        # Also derive science-frame status from the on-disk products.
        try:
            science_status.derive_science_from_disk(
                pypeIt.run_state, str(self.redux_dir), fitstbl=pypeIt.fitstbl)
        except Exception as e:
            # Deliberate never-crash-the-GUI guard: the science status is
            # best-effort and must never block the calibration status the
            # dashboard primarily needs.
            log.warning(f'Could not derive science state for '
                        f'{self.pypeit_file}: {e}')

        self.run_state = pypeIt.run_state
        self.load_status = self.LOAD_DERIVED

        # Cache the planned science/standard frames from the metadata we
        # just built (free here), then seed them so the Science view lists
        # upcoming frames before they are reduced, mirroring the planned
        # calibrations.  The cache makes later state reloads cheap.
        try:
            self._PLANNED_SCIENCE_CACHE[str(self.pypeit_file)] = \
                science_status.planned_science_from_fitstbl(pypeIt.fitstbl)
            self._seed_planned_science()
        except Exception as e:
            # Deliberate never-crash-the-GUI guard: planned frames are a
            # display nicety and must not fail the derive.
            log.warning(f'Could not seed planned science frames for '
                        f'{self.pypeit_file}: {e}')

    def _planned_science_frames(self, allow_build=False):
        """
        Return the cached list of planned science/standard frames for this
        ``.pypeit`` file.  On a cache miss, build the metadata to compute
        it **only** when ``allow_build`` is True (so live-monitoring
        reloads and CI loads never trigger a heavy, possibly-failing
        build).

        Parameters
        ----------
        allow_build : :obj:`bool`, optional
            Permit a one-time metadata build on a cache miss (True at
            launch; False on refresh/CI).

        Returns
        -------
        :obj:`list`
            Planned-frame dicts (see
            :func:`~pypeit.state.science_status.planned_science_from_fitstbl`);
            empty when uncached and a build is not permitted/possible.
        """
        key = str(self.pypeit_file)
        if key in self._PLANNED_SCIENCE_CACHE:
            return self._PLANNED_SCIENCE_CACHE[key]
        if not allow_build:
            return []
        planned = []
        try:
            # Cold cache + a build is allowed (e.g. relaunch with a state
            # file already present): build the metadata once to learn the
            # planned frames, then cache it for the rest of the session.
            pypeIt = pypeit.PypeIt(str(self.pypeit_file), reuse_calibs=True,
                                   calib_only=True)
            planned = science_status.planned_science_from_fitstbl(
                pypeIt.fitstbl)
        except Exception as e:
            # Deliberate never-crash-the-GUI guard: the metadata build can
            # fail in many ways (missing raw data, bad headers, ...) and the
            # planned-frame list is a display nicety.
            log.warning(f'Could not determine planned science frames for '
                        f'{self.pypeit_file}: {e}')
        self._PLANNED_SCIENCE_CACHE[key] = planned
        return planned

    def _seed_planned_science(self, allow_build=False):
        """
        Seed the planned science/standard frames into :attr:`run_state`
        (from the cache, or a one-time build when ``allow_build``), keyed
        to the calibration ``(group, det)`` pairs in the current state.
        Best-effort: never raises into the load/derive path.

        Parameters
        ----------
        allow_build : :obj:`bool`, optional
            Permit a one-time metadata build on a cache miss.
        """
        if self.run_state is None or self.header_info is None:
            return
        try:
            planned = self._planned_science_frames(allow_build=allow_build)
            if not planned:
                return
            group_dets = {}
            for g, d in self.calib_det_pairs():
                group_dets.setdefault(g, []).append(d)
            science_status.seed_planned_science_entries(
                self.run_state, planned, group_dets)
        except Exception as e:
            # Deliberate never-crash-the-GUI guard: seeding planned frames
            # must never break loading an otherwise-valid state.
            log.warning(f'Could not seed planned science frames for '
                        f'{self.pypeit_file}: {e}')

    def is_started(self):
        """
        Whether the reduction has any calibration-state entries.

        Returns
        -------
        :obj:`bool`
            True if a state with at least one entry is available.
        """
        return self.run_state is not None \
            and self.run_state.get_status() is not None

    def default_steps(self):
        """
        Return the active spectrograph's calibration steps, in pipeline
        order (specific to the ``pypeline``; see
        :attr:`_CALIB_CLASS_BY_PYPELINE`), including ``bpm``.

        A pypeline without an entry in :attr:`_CALIB_CLASS_BY_PYPELINE`
        (e.g. a newly added pipeline type) falls back to the full set of
        known calibration steps, so the dashboard remains usable until the
        new pipeline's order is registered.

        Returns
        -------
        :obj:`list`
            The ordered step names, or ``[]`` if the pipeline is unknown
            (e.g. the ``.pypeit`` file was not found).
        """
        if self.header_info is None:
            return []
        calib_class = self._CALIB_CLASS_BY_PYPELINE.get(self.header_info.path)
        if calib_class is None:
            log.warning(f'No calibration-step order registered for the '
                        f'{self.header_info.path} pypeline; falling back to '
                        f'the full set of known calibration steps.')
            return list(Calibrations.step_frame_map.keys())
        return calib_class.default_steps()

    def step_order(self, include_bpm=False):
        """
        Return the pipeline-aware step order for the calibration button row.

        Parameters
        ----------
        include_bpm : :obj:`bool`, optional
            Include the ``bpm`` step.  Default False — ``bpm`` runs
            internally but has no standalone output/QA, so it is omitted
            from the button row.

        Returns
        -------
        :obj:`list`
            The ordered step names.
        """
        steps = self.default_steps()
        if include_bpm:
            return steps
        return [s for s in steps if s != 'bpm']

    def status_table(self):
        """
        Return the normalized calibration-status table the views consume.

        Built on :meth:`~pypeit.state.run_state.RunPypeItState.get_status`
        but normalized: ``required`` becomes a real :obj:`bool` (or
        ``None``), the ``"--"`` sentinels become ``None``/``"absent"``, an
        ``in_pipeline`` column is added (membership in
        :meth:`default_steps`), and an empty/unavailable state yields an
        empty table.

        Returns
        -------
        `astropy.table.Table`_
            Columns :attr:`STATUS_COLUMNS`; empty when the reduction has
            not started or the state is unavailable.
        """
        if self.run_state is None:
            return state_table([], names=self.STATUS_COLUMNS)
        raw = self.run_state.get_status()
        if raw is None:
            return state_table([], names=self.STATUS_COLUMNS)

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
                'status': self._ABSENT if row['status'] == '--'
                else row['status'],
                'output_file': None if row['output_file'] == '--'
                else row['output_file'],
            })
        return state_table(rows, names=self.STATUS_COLUMNS)

    def calib_det_pairs(self):
        """
        Enumerate the ``(calibration_group, detector)`` pairs present in
        the state, in sorted order, for the scope drop-downs/selectors.

        Returns
        -------
        :obj:`list`
            ``(calib_id, det)`` tuples; ``det`` is left in its raw form
            (an :obj:`int` or a mosaic tuple) for the view to format.
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
        Whether the **required** calibrations for one ``(group, det)`` are
        all built successfully — the precondition for a science (Re)Build:
        until the calibrations a science frame depends on exist, its
        (Re)Build is disabled.

        Parameters
        ----------
        group : :obj:`int`
            The calibration group ID of the science frame.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            The detector (int) or mosaic the science frame is on.

        Returns
        -------
        :obj:`bool`
            True only if at least one required, in-pipeline calibration
            step exists for ``(group, det)`` and **all** such steps have a
            ``success``/``complete`` status; False otherwise (including
            when no matching calibration entries are known).
        """
        table = self.status_table()
        # The required, in-pipeline calibration entries for this (group, det).
        req = [row for row in table
               if row['calibration_group'] == group
               and same_det(row['detector'], det)
               and row['required'] is True and row['in_pipeline']]
        if len(req) == 0:
            return False
        return all(row['status'] == 'success' for row in req)

    def is_stale(self):
        """
        Whether the loaded ``*_state.json`` looks out of date relative to
        the ``.pypeit`` file or the calibration outputs.

        Returns
        -------
        :obj:`bool`
            True if the state file's mtime is older than the ``.pypeit``
            file or the newest file in ``Calibrations/``; False if there
            is no state file or nothing newer.
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
    # Per-step detail accessors (for the Calibrations view).
    # ------------------------------------------------------------------

    @property
    def calib_dir(self):
        """
        The reduction's ``Calibrations/`` directory (where output files
        live).

        Returns
        -------
        :obj:`pathlib.Path`
            ``<redux_dir>/Calibrations``.
        """
        return self.redux_dir / 'Calibrations'

    def step_entry(self, step, group, det):
        """
        Return the raw ``RunPypeItState`` entry for one
        ``(step, group, det)``, or ``None`` if absent.

        Parameters
        ----------
        step : :obj:`str`
            Calibration step name (e.g. ``wv_calib``).
        group : :obj:`int`
            Calibration group ID.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        object
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

        Parameters
        ----------
        step : :obj:`str`
            Calibration step name.
        group : :obj:`int`
            Calibration group ID.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        :obj:`pathlib.Path` or None
            ``Calibrations/<output_file>``.
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
        detect an active run via the log's modification time).

        Returns
        -------
        :obj:`pathlib.Path`
            ``<pypeit_root>.log`` (next to the ``.pypeit`` file, matching
            ``pypeit_run_to_calibstep``'s default log path).
        """
        return self.pypeit_file.with_suffix('.log')

    def step_output_files(self, step, group, det):
        """
        The existing ``Calibrations/`` file(s) a (re)build of this step
        would **overwrite** — named explicitly in the clobber confirmation
        dialog.

        Most steps write a single ``output_file``; ``slits`` writes
        **both** ``Slits_*`` (its ``output_file``) and the companion
        ``Edges_*``.  Only files that currently exist on disk are returned
        (a step with no output yet is a fresh build, nothing to overwrite).

        Parameters
        ----------
        step : :obj:`str`
            Calibration step name.
        group : :obj:`int`
            Calibration group ID.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        :obj:`list`
            Existing :obj:`pathlib.Path` outputs (possibly empty).
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

        Parameters
        ----------
        prefix : :obj:`str`
            File-type prefix (e.g. ``Edges``, ``Slits``).
        group : :obj:`int`
            Calibration group ID.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).
        ext : :obj:`str`, optional
            Extension (``fits`` or ``fits.gz``).

        Returns
        -------
        :obj:`pathlib.Path` or None
            The constructed path (``None`` if the setup is unknown).
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

        Parameters
        ----------
        step : :obj:`str`
            Calibration step name.
        group : :obj:`int`
            Calibration group ID.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        :obj:`dict`
            Metric name → value (empty if the entry is absent or has no
            entry-level metrics). ``bias`` → mean/std; ``slits`` →
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

        Parameters
        ----------
        step : :obj:`str`
            Calibration step name.
        group : :obj:`int`
            Calibration group ID.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        :obj:`list`
            One dict per slit (sorted by slit id), with keys ``slit``,
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
    # Science-frame accessors (for the Science view).
    # ------------------------------------------------------------------

    def science_table(self):
        """
        Return the normalized per-frame science-status table the Science
        view renders (one row per ``(frame, detector)``).

        Returns
        -------
        `astropy.table.Table`_
            Columns :attr:`SCIENCE_COLUMNS`; empty when there are no
            science entries.
        """
        if self.run_state is None:
            return state_table([], names=self.SCIENCE_COLUMNS)
        raw = self.run_state.get_science_status()
        if raw is None or len(raw) == 0:
            return state_table([], names=self.SCIENCE_COLUMNS)
        return raw

    def has_science(self):
        """
        Whether any science-frame state is available.

        Returns
        -------
        :obj:`bool`
            True if the state has at least one science entry.
        """
        return len(self.science_table()) > 0

    def science_frame_entry(self, frame, det):
        """
        Return the raw ``ScienceFrameState`` for one ``(frame, det)``, or
        ``None``.

        Parameters
        ----------
        frame : :obj:`str`
            Exposure basename.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        object
            The matching pydantic science entry, or ``None``.
        """
        if self.run_state is None:
            return None
        return self.run_state.science_entry(frame, det)

    def science_slit_table(self, frame, det):
        """
        Return the per-slit science rows for one frame (``ScienceSlit``).

        Parameters
        ----------
        frame : :obj:`str`
            Exposure basename.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        :obj:`list`
            One dict per slit (sorted), keys ``slit``, ``status``,
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

        Parameters
        ----------
        frame : :obj:`str`
            Exposure basename.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        :obj:`list`
            One dict per detected object, keys ``objid``, ``slitid``,
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
        Return **all** on-disk QA PNGs for one science ``(frame, det)``:
        the per-object ``obj_prof``/``obj_trace`` figures and the
        frame-level ``spec_flex_*`` flexure figures.

        Derived **at view time** by globbing ``<redux>/QA/PNGs/`` for files
        whose name starts with the frame basename and contains the detector
        name — so no science-state change is needed and it works for both
        live and derived state.

        Parameters
        ----------
        frame : :obj:`str`
            Exposure basename (the science-table key).
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).

        Returns
        -------
        :obj:`list`
            Sorted :obj:`pathlib.Path` of the matching QA PNGs (empty if
            the ``QA/PNGs`` directory or matches are absent).
        """
        qa_dir = self.redux_dir / 'QA' / 'PNGs'
        if not qa_dir.is_dir():
            return []
        detname = self.det_name(det)
        return [p for p in sorted(qa_dir.glob('*.png'))
                if p.name.startswith(str(frame)) and detname in p.name]

    def science_object_qa_files(self, frame, det, slitid):
        """
        Return **all** per-object QA PNGs for one object's slit, so the
        per-object table can open them.  The science QA PNGs encode the
        slit as ``S{spat_id:04d}``, matching ``ScienceObj.slitid``.

        Discovery is purely by the QA file-naming convention (no fixed
        list of QA types), so newly added per-object QA figures are picked
        up with no dashboard change.

        Parameters
        ----------
        frame : :obj:`str`
            Exposure basename.
        det : :obj:`int`, :obj:`tuple`, :obj:`list`
            Detector (int) or mosaic (tuple/list).
        slitid : :obj:`int`
            The object's slit (SPAT_ID).

        Returns
        -------
        :obj:`list`
            Sorted :obj:`pathlib.Path` of the object's QA PNGs (empty if
            none, or if ``slitid`` is None).
        """
        if slitid is None:
            return []
        tag = f'S{int(slitid):04d}'
        return [path for path in self.science_qa_files(frame, det)
                if tag in path.name]


# Module-level aliases for the load-status values, kept so the view modules
# and tests can reference them as ``model.LOAD_*`` without holding a
# DashboardModel instance.  The class attributes are authoritative.
LOAD_STATE_FILE = DashboardModel.LOAD_STATE_FILE
LOAD_DERIVED = DashboardModel.LOAD_DERIVED
LOAD_NOT_STARTED = DashboardModel.LOAD_NOT_STARTED
LOAD_MALFORMED = DashboardModel.LOAD_MALFORMED
LOAD_FILE_NOT_FOUND = DashboardModel.LOAD_FILE_NOT_FOUND
LOAD_ERROR = DashboardModel.LOAD_ERROR
