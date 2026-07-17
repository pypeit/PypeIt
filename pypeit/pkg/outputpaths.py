"""
Centralized, singleton-managed output-path resolution for PypeIt reductions.

A single instance, :obj:`pypeit.outputPaths`, is created at import time with
cwd-based defaults (mirroring :obj:`pypeit.log` and :obj:`pypeit.dataPaths`),
but is not yet *configured*. Orchestration-layer code (a ``scripts/`` script
class, or :class:`~pypeit.pypeit.PypeIt`'s ``__init__`` for entry points
where a script can't build the parameter set itself, e.g. ``run_pypeit``)
calls :func:`PypeItOutputPaths.configure` exactly once per CLI execution to
lock in the real, parameter-derived values. Code in ``pypeit/core/`` must
NOT import or query this object; it only ever receives already-resolved
``Path``/``str`` arguments.

Each managed directory (``redux``, ``science``, ``qa``, ``qa_pngs``,
``calibrations``, ``coadd_science``, ``coadd_qa``, ``coadd_qa_pngs``,
``collate``) is tracked internally as a small :class:`_ManagedPath` record
holding ``parent``, ``name``, ``full``, and ``ready``. :func:`configure`
(and ``__init__``) set only ``parent``/``name`` for every record -- pure
path arithmetic, no I/O, no change to ``ready``. The corresponding
``@property`` (e.g. :attr:`~PypeItOutputPaths.science`) is what resolves and
caches ``full`` (on first access, if not already cached) and creates the
directory on disk exactly once, flipping that record's ``ready`` flag --
but only if the instance has been configured and is not in dry-run mode. In
dry-run mode (or before :func:`configure` has ever been called), a property
still resolves and returns ``full``; it just never touches the filesystem
or sets ``ready``.

.. include:: ../include/links.rst
"""
import dataclasses
import logging
from pathlib import Path

from .exceptions import PypeItPathError

# Retrieves the already-configured "pypeit" logger by name from Python's
# logging registry -- does NOT recreate or reconfigure it (that already
# happened once, in pypeit/__init__.py, before this module is imported).
# A relative/absolute `from pypeit import log` is deliberately avoided here:
# pkg/ modules do not import back from the partially-initialized top-level
# pypeit package (see the circular-import warning atop pkg/logger.py).
log = logging.getLogger('pypeit')


@dataclasses.dataclass
class _ManagedPath:
    """Internal record for one managed output directory."""
    parent: Path
    """Parent directory."""
    name: str
    """Directory name, relative to `parent`."""
    full: Path = None
    """Full, absolute path (`parent`/`name`); resolved and cached on first
    access."""
    ready: bool = False
    """Whether the directory has been created on disk (or was already
    found to exist) and is available to receive files."""


class PypeItOutputPaths:
    """
    Single source of truth for PypeIt's output directory structure.

    See the module docstring for the configure-once / lazy-create-on-access
    lifecycle. Parameters mirror the ``PypeItPar`` keys that define the same
    concepts (:class:`~pypeit.par.pypeitpar.ReduxPar`,
    :class:`~pypeit.par.pypeitpar.CalibrationsPar`).

    Parameters
    ----------
    redux_path : :obj:`str`, :obj:`~pathlib.Path`, optional
        Root reduction directory. Defaults to the current working
        directory.
    scidir : :obj:`str`, optional
        Name of the science-output subdirectory.
    qadir : :obj:`str`, optional
        Name of the quality-assessment subdirectory.
    calib_dir : :obj:`str`, optional
        Name of the processed-calibrations subdirectory.
    coadd_suffix : :obj:`str`, optional
        Suffix appended to ``scidir``/``qadir`` for the 2D-coadd output
        directories.
    collate_outdir : :obj:`str`, :obj:`~pathlib.Path`, optional
        Output directory for :ref:`pypeit_collate_1d`. Defaults to
        ``redux_path`` if not given.
    dryrun : :obj:`bool`, optional
        If True, properties resolve ``full`` but never create directories
        or set ``ready``, and :func:`configure` may be called more than
        once. Intended for inspection/testing only; always ``False`` in
        production code paths.
    """

    def __init__(self, redux_path=None, scidir='Science', qadir='QA',
                 calib_dir='Calibrations', coadd_suffix='_coadd',
                 collate_outdir=None, dryrun=False):
        self._configured = False
        self._dryrun = dryrun
        self._paths = {}
        self._apply(redux_path, scidir, qadir, calib_dir, coadd_suffix,
                    collate_outdir)

    # ---- internal state management -------------------------------------
    def _apply(self, redux_path, scidir, qadir, calib_dir, coadd_suffix,
               collate_outdir):
        """
        Compute and store ``parent``/``name`` for every managed directory,
        replacing any existing records (so ``full``/``ready`` are always
        reset to unresolved/not-ready here). Pure path arithmetic -- no I/O,
        no change to ``_configured``.
        """
        redux_full = Path(redux_path).absolute() if redux_path is not None \
                        else Path.cwd()
        qa_full = redux_full / qadir
        coadd_qa_name = f'{qadir}{coadd_suffix}'
        coadd_qa_full = redux_full / coadd_qa_name
        collate_full = Path(collate_outdir).absolute() if collate_outdir is not None \
                            else redux_full

        self._paths = {
            'redux':         _ManagedPath(redux_full.parent, redux_full.name),
            'science':       _ManagedPath(redux_full, scidir),
            'qa':            _ManagedPath(redux_full, qadir),
            'qa_pngs':       _ManagedPath(qa_full, 'PNGs'),
            'calibrations':  _ManagedPath(redux_full, calib_dir),
            'coadd_science': _ManagedPath(redux_full, f'{scidir}{coadd_suffix}'),
            'coadd_qa':      _ManagedPath(redux_full, coadd_qa_name),
            'coadd_qa_pngs': _ManagedPath(coadd_qa_full, 'PNGs'),
            'collate':       _ManagedPath(collate_full.parent, collate_full.name),
        }

    def _get(self, key: str) -> Path:
        """Resolve (and, if eligible, create) one managed directory."""
        rec = self._paths[key]
        if rec.full is None:
            rec.full = rec.parent / rec.name
        if not rec.ready and self._configured and not self._dryrun:
            rec.full.mkdir(parents=True, exist_ok=True)
            rec.ready = True
            log.info(f'Output directory ready ({key}): {rec.full}')
        return rec.full

    # ---- resolved-path properties ---------------------------------------
    @property
    def redux(self) -> Path:
        """Root reduction directory."""
        return self._get('redux')

    @property
    def science(self) -> Path:
        """Science-output directory."""
        return self._get('science')

    @property
    def qa(self) -> Path:
        """Quality-assessment directory."""
        return self._get('qa')

    @property
    def qa_pngs(self) -> Path:
        """QA PNG-plot directory."""
        return self._get('qa_pngs')

    @property
    def calibrations(self) -> Path:
        """Processed-calibrations directory."""
        return self._get('calibrations')

    @property
    def coadd_science(self) -> Path:
        """2D-coadd science-output directory."""
        return self._get('coadd_science')

    @property
    def coadd_qa(self) -> Path:
        """2D-coadd quality-assessment directory."""
        return self._get('coadd_qa')

    @property
    def coadd_qa_pngs(self) -> Path:
        """2D-coadd QA PNG-plot directory."""
        return self._get('coadd_qa_pngs')

    @property
    def collate(self) -> Path:
        """Output directory for :ref:`pypeit_collate_1d`."""
        return self._get('collate')

    # ---- object-level state (read-only) ---------------------------------
    @property
    def configured(self) -> bool:
        """Whether :func:`configure` has been called on this instance."""
        return self._configured

    @property
    def dryrun(self) -> bool:
        """Whether this instance is in dry-run (inspection-only) mode."""
        return self._dryrun

    # ---- one-time (re)configuration --------------------------------------
    def configure(self, par=None, redux_path=None, scidir=None, qadir=None,
                  calib_dir=None, coadd_suffix=None, collate_outdir=None,
                  dryrun=None, caller=None):
        """
        Configure the instance from a :class:`~pypeit.par.pypeitpar.PypeItPar`
        and/or explicit overrides, exactly as if it had just been
        constructed with these arguments -- every managed directory's
        ``full``/``ready`` state is discarded and recomputed from scratch.
        May only be called once per instance unless the instance is in
        dry-run mode.

        Parameters
        ----------
        par : :class:`~pypeit.par.pypeitpar.PypeItPar`, optional
            Full parameter set to draw ``redux_path``/``scidir``/
            ``qadir``/``calib_dir``/``collate_outdir`` from. Explicit
            keyword arguments take precedence over the values in ``par``.
        redux_path : :obj:`str`, :obj:`~pathlib.Path`, optional
            Explicit override for the reduction root.
        scidir : :obj:`str`, optional
            Explicit override for the science-output subdirectory name.
        qadir : :obj:`str`, optional
            Explicit override for the QA subdirectory name.
        calib_dir : :obj:`str`, optional
            Explicit override for the calibrations subdirectory name.
        coadd_suffix : :obj:`str`, optional
            Explicit override for the 2D-coadd directory suffix.
        collate_outdir : :obj:`str`, :obj:`~pathlib.Path`, optional
            Explicit override for the collate1d output directory.
        dryrun : :obj:`bool`, optional
            If not None, set the instance's dry-run mode.
        caller : :obj:`str`, optional
            Free-form string identifying the calling code, included in the
            log message for easier debugging.

        Raises
        ------
        PypeItPathError
            If the instance has already been configured and is not in
            dry-run mode.
        """
        if self._configured and not self._dryrun:
            raise PypeItPathError(
                'PypeItOutputPaths has already been configured for this '
                'execution and cannot be reconfigured (only permitted when '
                'the instance is in dry-run mode, for inspection/testing).')

        if par is not None:
            rdx, cal = par['rdx'], par['calibrations']
            redux_path = redux_path if redux_path is not None else rdx['redux_path']
            scidir = scidir if scidir is not None else rdx['scidir']
            qadir = qadir if qadir is not None else rdx['qadir']
            calib_dir = calib_dir if calib_dir is not None else cal['calib_dir']
            if collate_outdir is None and 'collate1d' in par.keys() \
                    and par['collate1d']['outdir'] is not None:
                collate_outdir = par['collate1d']['outdir']

        if dryrun is not None:
            self._dryrun = dryrun

        self._apply(redux_path,
                    scidir if scidir is not None else 'Science',
                    qadir if qadir is not None else 'QA',
                    calib_dir if calib_dir is not None else 'Calibrations',
                    coadd_suffix if coadd_suffix is not None else '_coadd',
                    collate_outdir)
        self._configured = True

        ctx = f' [caller={caller}]' if caller else ''
        dry = ' (dry run)' if self._dryrun else ''
        log.info(f'Output paths configured{ctx}{dry}: '
                 f'redux_path={self._paths["redux"].parent / self._paths["redux"].name}')

    # ---- held in reserve — stub only, not wired to any caller ------------
    def derive(self, **overrides):
        """
        Return an independent :class:`PypeItOutputPaths` seeded from current
        state with the given overrides applied.

        Reserved for future work; not called anywhere in the current
        implementation. Its original motivation (letting coadd/collate
        compute derived paths without touching the shared singleton) is
        largely superseded by treating `coadd_science`/`coadd_qa`/`collate`
        as first-class properties on this same object -- kept as a stub
        rather than removed, pending a decision on whether it's still
        needed for some other future use.

        Raises
        ------
        NotImplementedError
            Always -- this method is not yet implemented.
        """
        raise NotImplementedError(
            'PypeItOutputPaths.derive() is reserved for future work and is '
            'not yet implemented.')

    def __repr__(self):
        redux = self._paths['redux']
        return (f'<{self.__class__.__name__}: '
                f'redux_path={redux.parent / redux.name}, '
                f'configured={self._configured}, dryrun={self._dryrun}>')
