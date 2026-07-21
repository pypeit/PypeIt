"""
Data model and I/O for the **PypeIt reduction state**.

While ``run_pypeit`` reduces data it records, per calibration group and per
detector/mosaic, the status of each calibration step (plus step-specific
metrics and per-slit/order detail).  This state is held by the pydantic
:class:`RunPypeItState` model and serialized to ``<pypeit_root>_state.json``
in the reduction directory.  It is a *live* record — updated at each step
transition — and is consumed by tools rather than edited by hand: the PypeIt
Dashboard (:ref:`dashboard`) and ``pypeit_status`` both read it (see
:ref:`state` in the documentation).

State I/O is deliberately non-essential to the reduction: it is written
through :meth:`RunPypeItState.safe_write` / :meth:`safe_update_calib`, which
log and swallow any error so that state bookkeeping can never abort a run.

.. include:: ../include/links.rst
"""
import json
from pathlib import Path

from pydantic import BaseModel, Field, ValidationError
from typing import List, Optional, Dict, Literal

import numpy as np
from astropy import table

from pypeit import log


def same_det(det1, det2):
    """
    Compare two detector identifiers for equality, treating a detector
    mosaic the same whether it is stored as a tuple or a list.

    A single detector is an :obj:`int`; a detector mosaic is a tuple of
    ints (e.g. ``(1, 5)``).  Because pydantic coerces the stored ``det``
    field to a :obj:`list`, a direct ``==`` comparison between a stored
    mosaic (``[1, 5]``) and an incoming mosaic (``(1, 5)``) is always
    False.  This helper normalizes both sides before comparing.

    Args:
        det1 (:obj:`int`, :obj:`list`, :obj:`tuple`):
            First detector identifier (single detector or mosaic).
        det2 (:obj:`int`, :obj:`list`, :obj:`tuple`):
            Second detector identifier (single detector or mosaic).

    Returns:
        :obj:`bool`: True if the two identifiers refer to the same
        detector or detector mosaic.
    """
    # Normalize mosaics (list/tuple) to a tuple; leave scalars as-is
    norm1 = tuple(det1) if isinstance(det1, (list, tuple)) else det1
    norm2 = tuple(det2) if isinstance(det2, (list, tuple)) else det2
    return norm1 == norm2


def state_table(rows, names=None):
    """
    Construct an `astropy.table.Table`_ from a list of row dictionaries,
    forcing every column to ``object`` dtype.

    The ``object`` dtype preserves heterogeneous values within one column —
    e.g. a ``detector`` column holding both single-detector integers and
    mosaic tuples, or entries that are ``None`` — which a plain numpy
    string/integer column cannot represent.

    Args:
        rows (:obj:`list`):
            One :obj:`dict` per table row, all sharing the same keys.
        names (:obj:`list`, optional):
            Column names.  Required to build an empty (zero-row) table when
            ``rows`` is empty; otherwise defaults to the keys of the first
            row.

    Returns:
        `astropy.table.Table`_: The assembled table, with all columns of
        ``object`` dtype; zero-length if ``rows`` is empty.
    """
    if names is None:
        names = list(rows[0].keys()) if rows else []
    tbl = table.Table()
    for name in names:
        # Fill through an object array so a tuple-valued cell (a detector
        # mosaic) is stored as one cell, not broadcast into a 2D column.
        col = np.empty(len(rows), dtype=object)
        col[:] = [row[name] for row in rows]
        tbl[name] = col
    return tbl

# Calibration state
class BaseCalibState(BaseModel):
    calib_id: int # Calibration ID
    det: int | List[int]  # Detector number or mosaic tuple
    step: str
    required: bool = False
    input_files: Optional[List[str]] = None
    output_file: Optional[str] = None
    qa_files: Optional[List[str]] = None
    status: Literal["complete", "fail", "undone", "running", "success"] = "undone"

class BiasCalibState(BaseCalibState):
    step: Literal["bias"] = "bias"
    # Metrics
    mean: Optional[float] = None
    std: Optional[float] = None

class WvCalibSlit(BaseModel):
    status: Literal["success", "fail", "undone", ] = "undone"
    # Metrics
    rms: Optional[float] = None

class WvCalibState(BaseCalibState):
    step: Literal["wv_calib"] = "wv_calib"
    slits: Optional[Dict[int, WvCalibSlit]] = Field(default_factory=dict)

class SlitEdges(BaseModel):
    status: Literal["success", "fail", "undone", ] = "undone"
    # Metrics
    center: Optional[float] = None
    slitord_id: Optional[int] = None

class SlitEdgesState(BaseCalibState):
    step: Literal["slits"] = "slits"
    nslits: Optional[int] = None
    slits: Optional[Dict[int, SlitEdges]] = Field(default_factory=dict)

class TiltsSlit(BaseModel):
    status: Literal["success", "fail", "undone", ] = "undone"
    # Metrics
    rms: Optional[float] = None

class TiltsState(BaseCalibState):
    step: Literal["tilts"] = "tilts"
    slits: Optional[Dict[int, TiltsSlit]] = Field(default_factory=dict)

class FlatCorrectionMetric(BaseModel):
    """
    Mean and scatter of one applied flat-field correction over a single
    slit's good, on-slit pixels.
    """
    # Mean of the correction image (these corrections hover about 1.0)
    mean: Optional[float] = None
    # Scatter about the mean (np.std) of the correction image
    rms: Optional[float] = None

class FlatsSlit(BaseModel):
    """
    Per-slit flat-field state: an overall status plus the per-correction
    mean/RMS metrics.
    """
    # 'skip' is flats-specific (SKIPFLATCALIB): the slit was intentionally
    # skipped, distinct from a generation failure ('fail').
    status: Literal["success", "fail", "skip", "undone"] = "undone"
    # Per-correction metrics, keyed by correction name
    # ('pixelflat', 'spat_illum', 'spec_illum')
    corrections: Optional[Dict[str, FlatCorrectionMetric]] = \
        Field(default_factory=dict)

class FlatsState(BaseCalibState):
    step: Literal["flats"] = "flats"
    # Names of the corrections present in the merged FlatImages product:
    # subset of 'pixelflat' (pixel-to-pixel), 'spat_illum' (slit
    # illumination), 'spec_illum' (spectral illumination)
    corrections: Optional[List[str]] = Field(default_factory=list)
    # Provenance of the pixel flat: 'raw', 'user_file', or 'slitless'
    pixelflat_source: Optional[str] = None
    # Raw input frames grouped by role (each may be absent)
    pixelflat_files: Optional[List[str]] = None
    illumflat_files: Optional[List[str]] = None
    lampoff_files: Optional[List[str]] = None
    # Per-slit status + metrics
    slits: Optional[Dict[int, FlatsSlit]] = Field(default_factory=dict)

class DarkCalibState(BaseCalibState):
    step: Literal["dark"] = "dark"

class ArcCalibState(BaseCalibState):
    step: Literal["arc"] = "arc"

class TiltImgCalibState(BaseCalibState):
    step: Literal["tiltimg"] = "tiltimg"

class ScattLightCalibState(BaseCalibState):
    step: Literal["scattlight"] = "scattlight"

class AlignCalibState(BaseCalibState):
    step: Literal["align"] = "align"


# ---------------------------------------------------------------------------
# Science-frame state
#
# Unlike calibrations (one entry per step), a science exposure is modeled as a
# single entry per (frame/basename, detector) carrying the four macro-steps of
# the science reduction (process -> findobj -> skysub -> extract), the data
# products, and per-slit / per-object detail.
# ---------------------------------------------------------------------------

# The ordered macro-steps of a science-frame reduction
science_steps = ['process', 'findobj', 'skysub', 'extract']

class ScienceObj(BaseModel):
    """
    One detected (and possibly extracted) object on a science exposure.
    """
    objid: int
    slitid: Optional[int] = None
    spat_pixpos: Optional[float] = None     # spatial pixel position
    fwhm: Optional[float] = None
    snr_find: Optional[float] = None        # smash_snr from object finding
    s2n: Optional[float] = None             # median S/N from extraction
    sign: Optional[int] = None              # +1, or -1 for negative traces
    extracted: bool = False                 # True once a 1D spectrum exists

class ScienceSlit(BaseModel):
    """
    Per-slit science status (from the slit bitmask).
    """
    # 'fail' = BADSKYSUB/BADEXTRACT flagged; 'skip' reserved for skipped slits
    status: Literal["success", "fail", "skip", "undone"] = "undone"
    nobj: Optional[int] = None              # objects found on this slit

class ScienceStep(BaseModel):
    """
    Status of one macro-step (process/findobj/skysub/extract) of a science
    exposure.
    """
    status: Literal["undone", "running", "success", "fail"] = "undone"

class ScienceFrameState(BaseModel):
    """
    State of one reduced science/standard exposure on one detector/mosaic.
    """
    frame: str                              # combined basename of the exposure
    det: int | List[int]
    calib_id: int
    objtype: str                            # 'science' or 'standard'
    comb_id: Optional[int] = None
    bkg_id: Optional[int] = None
    # Contributing raw frame filename(s) — recorded so a tool (e.g. the
    # Dashboard's (Re)Build) can re-run this exposure via pypeit_reduce_by_step.
    raw_files: Optional[List[str]] = Field(default_factory=list)
    # Macro-step statuses
    process: ScienceStep = Field(default_factory=ScienceStep)
    findobj: ScienceStep = Field(default_factory=ScienceStep)
    skysub: ScienceStep = Field(default_factory=ScienceStep)
    extract: ScienceStep = Field(default_factory=ScienceStep)
    # Data products
    spec2d_file: Optional[str] = None
    spec1d_file: Optional[str] = None
    # Exposure-level + per-slit / per-object detail
    nobj: Optional[int] = None
    slits: Optional[Dict[int, ScienceSlit]] = Field(default_factory=dict)
    objects: Optional[List[ScienceObj]] = Field(default_factory=list)


# NOTE: the key order here is authoritative for the order in which steps
# are iterated/displayed (e.g. get_status, print_status).  It is kept in
# the pipeline's processing order, matching
# MultiSlitCalibrations.default_steps() (bpm is not tracked here); align is
# IFU-only and appended last.
calib_classes = {
    'bias': BiasCalibState,
    'dark': DarkCalibState,
    'slits': SlitEdgesState,
    'arc': ArcCalibState,
    'tiltimg': TiltImgCalibState,
    'wv_calib': WvCalibState,
    'tilts': TiltsState,
    'scattlight': ScattLightCalibState,
    'flats': FlatsState,
    'align': AlignCalibState,
}

slit_classes = {
    'wv_calib': WvCalibSlit,
    'tilts': TiltsSlit,
    'slits': SlitEdges,
    'flats': FlatsSlit,
}

class RunPypeItState(BaseModel):
    """
    The state of a PypeIt run.

    This pydantic model records, for one ``.pypeit`` reduction, the status of
    every calibration step **per calibration group and per detector/mosaic**.
    Each calibration field (``bias``, ``dark``, ``slits``, ``arc``,
    ``tiltimg``, ``wv_calib``, ``tilts``, ``scattlight``, ``flats``, ``align``)
    is a list of per-``(calib_id, det)`` entries (subclasses of
    :class:`BaseCalibState`) carrying that step's ``status``, ``required``
    flag, input/output/QA files, and step-specific metrics (and, for
    ``slits``/``wv_calib``/``tilts``/``flats``, per-slit/order detail).

    It is written to ``<pypeit_root>_state.json`` (see :attr:`outfile`) and
    is updated incrementally during a run; use :meth:`safe_write` /
    :meth:`safe_update_calib` so state I/O can never abort the reduction, and
    :meth:`get_status` for a tabular summary.  It is read by the PypeIt
    Dashboard and ``pypeit_status``; **it is not meant to be edited by hand**.

    Attributes:
        pypeit_file (str): The ``.pypeit`` file this state belongs to.
        current_step (str): The step most recently updated.
        current_det (int): The detector most recently updated.
        current_calibID (int): The calibration group most recently updated.
        previous_step (str): The step updated before ``current_step``.
        path (str): Optional explicit path for the state JSON file; if
            ``None``, :attr:`outfile` is derived from ``pypeit_file``.
    """

    # Required
    pypeit_file: str
    current_step: str
    current_det: int
    current_calibID: int

    # Optional
    previous_step: str = 'none'

    # Calibrations (listed in pipeline processing order, matching
    # calib_classes / default_steps(); align is IFU-only, appended last)
    bias: Optional[List[BiasCalibState]] = Field(default_factory=list)
    dark: Optional[List[DarkCalibState]] = Field(default_factory=list)
    slits: Optional[List[SlitEdgesState]] = Field(default_factory=list)
    arc: Optional[List[ArcCalibState]] = Field(default_factory=list)
    tiltimg: Optional[List[TiltImgCalibState]] = Field(default_factory=list)
    wv_calib: Optional[List[WvCalibState]] = Field(default_factory=list)
    tilts: Optional[List[TiltsState]] = Field(default_factory=list)
    scattlight: Optional[List[ScattLightCalibState]] = Field(default_factory=list)
    flats: Optional[List[FlatsState]] = Field(default_factory=list)
    align: Optional[List[AlignCalibState]] = Field(default_factory=list)
    # Science/standard exposures: one entry per (frame/basename, detector)
    science: Optional[List[ScienceFrameState]] = Field(default_factory=list)
    path: Optional[str] = None

    @property
    def outfile(self):
        """
        Path to the state JSON file.

        Returns:
            str: :attr:`path` if set, else the ``.pypeit`` file name with the
            ``.pypeit`` extension replaced by ``_state.json``.
        """
        if self.path is not None:
            return self.path
        return self.pypeit_file.replace('.pypeit', '_state.json')

    def load(self, path:str=None):
        """
        Load the state from :attr:`outfile`, if it exists.

        Args:
            path (:obj:`str`, optional):
                Unused; retained for backward compatibility (the file read is
                :attr:`outfile`).

        Returns:
            :class:`RunPypeItState`: The state validated from the JSON file,
            or ``self`` unchanged if no state file is present.
        """
        if not Path(self.outfile).is_file():
            return self
        log.info(f'Loading existing state from {self.outfile}')
        with open(self.outfile, 'rt') as fh:
            update_dict = json.load(fh)
        # Return
        return RunPypeItState.model_validate(update_dict)

    def merge_from_disk(self):
        """
        Overlay the calibration and science statuses from the existing on-disk
        state file onto this state (which is updated in place), matching
        entries by ``(calib_id, det)`` and ``(frame, det)``.

        A step-runner script (``pypeit_run_to_calibstep`` /
        ``pypeit_reduce_by_step``) starts with a *fresh* state — its
        calibrations are the required (``undone``) set and it has no science —
        so writing it would **reset the other portion** of the shared
        ``*_state.json`` (e.g. a science step-build would blank out the
        calibration statuses a prior calibration build wrote, and vice versa).
        Calling this first preserves whatever the other runner already recorded.
        This instance's :attr:`outfile` (and ``pypeit_file``) is left untouched,
        so the subsequent write goes to the same file.  Best-effort: a missing
        or unreadable state file leaves ``self`` unchanged.

        Returns
        -------
        None
        """
        if not Path(self.outfile).is_file():
            return
        try:
            with open(self.outfile, 'rt') as fh:
                prev = RunPypeItState.model_validate(json.load(fh))
        except (OSError, json.JSONDecodeError, ValidationError) as e:
            log.warning(f"Could not merge existing state {self.outfile}: {e}")
            return
        # Overlay each calibration step's prior entries (replace a matching
        # (calib_id, det), else append) so prior statuses survive.
        for step in calib_classes:
            current = getattr(self, step)
            for prior in getattr(prev, step):
                for i, item in enumerate(current):
                    if item.calib_id == prior.calib_id \
                            and same_det(item.det, prior.det):
                        current[i] = prior
                        break
                else:
                    current.append(prior)
        # Carry forward prior science entries this run does not already hold
        # (this run updates its own entry on top).
        for prior in prev.science:
            if self.science_entry(prior.frame, prior.det) is None:
                self.science.append(prior)



    def update_calib(self, step:str, calib_id: int, det, key:str, value,
                     slit:str=None):
        """
        Create or update the state entry for a single calibration step.

        Looks up the entry matching ``(calib_id, det)`` for the given
        ``step``; if none exists it is created.  The requested ``key`` is
        then set to ``value`` (appended if the field is a list), either on
        the step entry or, when ``slit`` is given, on the per-slit/order
        sub-entry.

        Args:
            step (:obj:`str`):
                Calibration step name (must be a key of ``calib_classes``).
            calib_id (:obj:`int`):
                Calibration group ID.
            det (:obj:`int`, :obj:`tuple`, :obj:`list`):
                Detector (int) or detector mosaic (tuple/list of ints).
            key (:obj:`str`):
                Name of the field to set on the step (or per-slit) entry.
            value:
                Value to assign.  A list/tuple replaces the field with a
                flat list (e.g. ``input_files``); a scalar appended to a
                field that already holds a list accumulates (e.g. flats
                ``types``); otherwise the scalar is assigned.
            slit (:obj:`int`, optional):
                Slit/order ID.  If provided, ``key`` is set on the
                per-slit sub-entry rather than the step entry itself.
        """
        # Current step
        if self.current_step != step:
            self.previous_step = self.current_step
        self.current_step = step
        # Select items so far
        if step not in calib_classes:
            return
        # Grab the entry
        self_items = getattr(self, step)
        found_it = False
        # Grab the item; use same_det so a mosaic stored as a list still
        # matches the same mosaic passed in as a tuple (see same_det)
        for index, item in enumerate(self_items):
            if item.calib_id == calib_id and same_det(item.det, det):
                found_it = True
                break

        # Create it?
        if not found_it:
            item = calib_classes[step](calib_id=calib_id, det=det)
            self_items.append(item)
            index = -1

        # Set
        if slit is None:
            current = getattr(self_items[index], key)
            if isinstance(value, (list, tuple)):
                # A list/tuple value (e.g. input_files=self.raw_files) is
                # stored wholesale as a flat list.  Doing this explicitly
                # (rather than appending) keeps repeat calls idempotent: it
                # avoids nesting the new list inside the old one, which
                # would corrupt the field and break load() validation.
                setattr(self_items[index], key, list(value))
            elif isinstance(current, list):
                # A scalar added to a field that is already a list (e.g. the
                # entry-level flats 'corrections') is appended.
                current.append(value)
            else:
                setattr(self_items[index], key, value)
        else:
            if slit not in self_items[index].slits.keys():
                self_items[index].slits[slit] = slit_classes[step]()
            setattr(self_items[index].slits[slit], key, value)

    def write(self):
        """
        Serialize the state to its JSON file (``self.outfile``).

        This may raise on I/O or serialization errors; use
        :func:`safe_write` from within a reduction so a failure here can
        never abort the run.
        """
        json_string = self.model_dump_json(exclude_none=True, indent=4, round_trip=True)
        # Write
        with open(self.outfile, 'w', encoding='utf-8') as f:
            f.write(json_string)

    def safe_write(self):
        """
        Write the state, catching and logging any error so that state I/O
        can never crash a PypeIt reduction.

        Returns:
            :obj:`bool`: True if the write succeeded, False if it failed
            (in which case a warning is logged and the run continues).
        """
        try:
            self.write()
        except Exception as e:
            # Deliberately broad: state I/O is non-essential bookkeeping and
            # must never abort the reduction.  Log and carry on.
            log.warning(f"Failed to write reduction state to {self.outfile}: "
                        f"{e}")
            return False
        else:
            return True

    def safe_update_calib(self, step:str, calib_id:int, det, key:str, value,
                          slit:str=None):
        """
        Wrapper on :func:`update_calib` that catches and logs any error so
        that state bookkeeping can never crash a PypeIt reduction.

        Args:
            step (:obj:`str`):
                Calibration step name.
            calib_id (:obj:`int`):
                Calibration group ID.
            det (:obj:`int`, :obj:`tuple`, :obj:`list`):
                Detector or detector mosaic.
            key (:obj:`str`):
                Field to set on the step (or per-slit) entry.
            value:
                Value to assign (or append, for list-valued fields).
            slit (:obj:`int`, optional):
                Slit/order ID for a per-slit update.

        Returns:
            :obj:`bool`: True if the update succeeded, False otherwise.
        """
        try:
            self.update_calib(step, calib_id, det, key, value, slit=slit)
        except Exception as e:
            # Deliberately broad: state bookkeeping is non-essential and must
            # never abort the reduction.  Log and carry on.
            log.warning(f"Failed to update reduction state "
                        f"(step={step}, calib_id={calib_id}, det={det}, "
                        f"key={key}): {e}")
            return False
        else:
            return True

    # -------------------------------------------------------------------
    # Science-frame state
    # -------------------------------------------------------------------
    def science_entry(self, frame:str, det):
        """
        Return the science entry matching ``(frame, det)``, or None.

        Args:
            frame (:obj:`str`): The exposure basename.
            det (:obj:`int`, :obj:`tuple`, :obj:`list`): Detector or mosaic.

        Returns:
            :class:`ScienceFrameState` or None.
        """
        for item in self.science:
            if item.frame == frame and same_det(item.det, det):
                return item
        return None

    def add_or_get_science(self, frame:str, det, calib_id:int=-1,
                           objtype:str='science', comb_id:int=None,
                           bkg_id:int=None):
        """
        Find the science entry for ``(frame, det)`` or create it.

        Args:
            frame (:obj:`str`): The exposure basename.
            det (:obj:`int`, :obj:`tuple`, :obj:`list`): Detector or mosaic.
            calib_id (:obj:`int`, optional): Calibration group ID.
            objtype (:obj:`str`, optional): 'science' or 'standard'.
            comb_id (:obj:`int`, optional): Combination group ID.
            bkg_id (:obj:`int`, optional): Background (A-B) group ID.

        Returns:
            :class:`ScienceFrameState`: The existing or newly-created entry.
        """
        entry = self.science_entry(frame, det)
        if entry is None:
            entry = ScienceFrameState(frame=frame, det=det, calib_id=calib_id,
                                      objtype=objtype, comb_id=comb_id,
                                      bkg_id=bkg_id)
            self.science.append(entry)
        return entry

    def safe_update_science(self, frame:str, det, step:str=None,
                            status:str=None, calib_id:int=-1,
                            objtype:str='science', **fields):
        """
        Get-or-create a science entry and update a step status and/or
        top-level fields, catching and logging any error so science
        bookkeeping can never crash a reduction.

        Args:
            frame (:obj:`str`): The exposure basename.
            det (:obj:`int`, :obj:`tuple`, :obj:`list`): Detector or mosaic.
            step (:obj:`str`, optional): One of ``science_steps``; if given,
                its ``status`` is set to ``status``.
            status (:obj:`str`, optional): Status to assign to ``step``.
            calib_id (:obj:`int`, optional): Calibration group (on create).
            objtype (:obj:`str`, optional): 'science'/'standard' (on create).
            **fields: Other top-level :class:`ScienceFrameState` fields to set
                (e.g. ``nobj``, ``spec1d_file``, ``spec2d_file``).

        Returns:
            :class:`ScienceFrameState` or None: The entry (None on failure).
        """
        try:
            entry = self.add_or_get_science(frame, det, calib_id=calib_id,
                                            objtype=objtype)
            if step is not None and status is not None:
                getattr(entry, step).status = status
            for key, value in fields.items():
                setattr(entry, key, value)
        except Exception as e:
            # Deliberately broad: science-state bookkeeping is non-essential
            # and must never abort the reduction.  Log and carry on.
            log.warning(f"Failed to update science state "
                        f"(frame={frame}, det={det}, step={step}): {e}")
            return None
        else:
            return entry

    def get_science_status(self):
        """
        Summarize the science-frame state as a table: one row per
        ``(frame, det)`` with the four macro-step statuses, object count,
        and product presence.

        Returns:
            `astropy.table.Table`_ or None: None if there are no science
            entries.
        """
        if self.science is None or len(self.science) == 0:
            return None
        rows = []
        for item in self.science:
            det = tuple(item.det) if isinstance(item.det, list) else item.det
            rows.append({
                "frame": item.frame,
                "detector": det,
                "calib": item.calib_id,
                "objtype": item.objtype,
                "process": item.process.status,
                "findobj": item.findobj.status,
                "skysub": item.skysub.status,
                "extract": item.extract.status,
                "nobj": item.nobj if item.nobj is not None else "--",
                "spec2d": "yes" if item.spec2d_file is not None else "--",
                "spec1d": "yes" if item.spec1d_file is not None else "--",
            })
        return state_table(rows)

    def get_status(self):
        """
        Summarize the state as a tabular, per-step status overview.

        This is intentionally less detailed than the full serialized state
        (:meth:`write`): it omits the per-slit detail and metrics, giving a
        scannable health overview (used by :meth:`print_status` and the
        dashboard's Status view).

        Returns:
            `astropy.table.Table`_ or None: One row per
            ``(calibration_group, detector, step)`` with columns
            ``calibration_group``, ``detector``, ``steps``, ``required``,
            ``status``, and ``output_file``; missing entries are filled with
            ``"--"``.  Returns ``None`` if no calibration entries exist.
        """
        # Collect all unique (calib_id, det) pairs across all steps
        pairs = {
                (item.calib_id, tuple(item.det) if isinstance(item.det, list) else item.det)
                for step in calib_classes for item in getattr(self, step)
                }

        if len(pairs) == 0:
            return None

        # Build all rows in one pass.  Sort with the detector stringified so
        # a run mixing single detectors (int) and mosaics (tuple) still sorts
        # (int and tuple are not comparable).
        rows = []
        for calib_id, det in sorted(pairs, key=lambda p: (p[0], str(p[1]))):
            for step_name, step_class in calib_classes.items():
                # Find the matching entry: next() returns the first entry in
                # this step's list matching (calib_id, det), or None if the
                # step has no entry for this pair.
                items = getattr(self, step_name)
                entry = next(
                        (item for item in items if item.calib_id == calib_id and
                        same_det(item.det, det)), None
                        )
                rows.append({
                    "calibration_group": calib_id,
                    "detector": det,
                    "steps": step_name,
                    "required": "--" if entry is None else str(entry.required),
                    "status": "--" if entry is None else entry.status,
                    "output_file": Path(entry.output_file).name
                        if entry is not None and entry.output_file is not None else "--"
                })

        # Make a single table
        return state_table(rows)

    def print_status(self):
        """
        Pretty-print the reduction state to the terminal: a per-(group,
        detector) calibration table (from :meth:`get_status`) followed by the
        science-frame status (:meth:`_print_science_status`).  Used by
        :ref:`pypeit_status`.

        Returns:
            None.
        """
        status_tbl = self.get_status()

        print(f'PypeIt Reduction Status: {Path(self.pypeit_file).name}')
        print('=' * 70)

        if status_tbl is None or len(status_tbl) == 0:
            print('  No calibration state entries found.')
        else:
            # Column widths for formatting
            widths = {'steps': 14, 'required': 10, 'status': 10,
                      'output_file': 20}
            col_header = (f"  {'Step':<13} {'Required':<9} {'Status':<9} "
                          f"{'Output File'}")
            separator = ('  ' + '-' * (widths['steps'] - 1)
                         + ' ' + '-' * (widths['required'] - 1)
                         + ' ' + '-' * (widths['status'] - 1)
                         + ' ' + '-' * widths['output_file'])
            # get_status() emits the rows already grouped by (calibration
            # group, detector), so print a new group header whenever the
            # pair changes.
            current = None
            for row in status_tbl:
                pair = (row['calibration_group'], row['detector'])
                if pair != current:
                    current = pair
                    print(f'\n  Calibration Group: {pair[0]}, '
                          f'Detector: {pair[1]}')
                    print(col_header)
                    print(separator)
                print(f"  {row['steps']:<{widths['steps']}}"
                      f"{str(row['required']):<{widths['required']}}"
                      f"{row['status']:<{widths['status']}}"
                      f"{row['output_file']}")

        # Science-frame status (if any)
        self._print_science_status()

    def _print_science_status(self):
        """
        Print the per-exposure science-frame status table to stdout (the
        four macro-steps, object count, and product presence).
        """
        sci_tbl = self.get_science_status()
        if sci_tbl is None or len(sci_tbl) == 0:
            return
        print('\n  Science Frames')
        col_header = (f"  {'Frame':<28} {'Type':<9} {'process':<8} "
                      f"{'findobj':<8} {'skysub':<8} {'extract':<8} "
                      f"{'nobj':<5} {'spec2d':<7} {'spec1d'}")
        print(col_header)
        print('  ' + '-' * (len(col_header) - 2))
        for row in sci_tbl:
            print(f"  {str(row['frame']):<28} {str(row['objtype']):<9} "
                  f"{row['process']:<8} {row['findobj']:<8} "
                  f"{row['skysub']:<8} {row['extract']:<8} "
                  f"{str(row['nobj']):<5} "
                  f"{row['spec2d']:<7} {row['spec1d']}")
