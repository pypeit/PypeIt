# PypeIt

PypeIt is a pure Python package for processing raw spectroscopic data from
astronomical telescopes into calibrated spectra for scientific analysis.

## Development Setup

- Create and activate a fresh python environment

- Clone your fork of the main pypeit repo, hosted at
  [GitHub](https://github.com/pypeit/PypeIt), and add a remote git connection to
  the main repo called `upstream`.

- Install pypeit using the zshell (or bash equivalent) command `pip install -e ".[dev]"`.

- Create a new git branch, where simple bug fixes branch from `release` and all
  other development branches from `develop`.

- Frequently fetch changes to the upstream base branch (`release` or `develop`)
  and rebase, as necessary.

- Use `pathlib` (e.g. `Path(...).name`, `.stem`, `.with_name()`, `.exists()`)
  for all new filesystem-path handling; prefer it over `os.path` in new code.

- Additional development guidelines is provided in `doc/dev/development.rst`.

## Key Architectural Components

- PypeIt implements four processing paths (pypelines), depending on the format
  of the spectrograph used to collect the raw data:

    - **MultiSlit**: Standard long-slit and multi-slit spectrographs
    - **Echelle**: Cross-dispersed echelle spectrographs
    - **SlicerIFU**: Slicer-based integral-field units
    - **Fiber**: Fiber-fed spectrographs (e.g., IFU fiber bundles)

- Instrument specifications and data-processing considerations that are specific
  to each spectrograph (including the processing path that should be used) are
  isolated in their respective spectrograph classes, all of which inherit from
  `pypeit.spectrographs.spectrograph.Spectrograph`.

- Users primarily modify the code performance via user-level parameters held by
  the `pypeit.par.pypeitpar.PypeItPar` class, which packages the hierarchy of
  parameter sets used throughout the code.

- Core processing modules, particularly for calibrations, produce FITS files
  that are saved to disk and reused as necessary.  These modules commonly use
  `pypeit.datamodel.DataContainer` as a base class to enforce strict adherence
  to a well-defined datamodel and to provide a common IO interface.

- Significant portions of the `pypeit/data` directory are not included in the
  package distribution of the code, but rely on a cache system that downloads
  files as needed for processing.

- The pypeline dispatch system uses class name conventions: setting
  `pypeline = 'Fiber'` causes `FindObjects.get_instance()` and
  `Extract.get_instance()` to find `FiberFindObjects` and `FiberExtract`
  respectively.  When adding a new pypeline, the string must be added to
  validation lists in `specobj.py`, `specobjs.py`, `slittrace.py`,
  `show_2dspec.py`, `spectrograph.py`, and `pypeit_steps.py`.

- Calibrations: The `Fiber` pypeline falls through to `IFUCalibrations` in
  `calibrations.py` (not in `['MultiSlit', 'Echelle']`), sharing the
  calibration flow with `SlicerIFU`.

### Fiber Pypeline

The `Fiber` pypeline handles fiber-fed spectrographs where each fiber IS the
object — no peak-detection object finding is needed.

- **FiberFindObjects** (`find_objects.py`): Creates one `SpecObj` per fiber
  with the trace at the fiber center (midpoint of slit edges). Uses joint sky
  subtraction inherited from `SlicerIFUFindObjects`. Must reset `reduce_bpm`
  after `global_skysub` because per-slit sky fitting rejects narrow fibers
  before the joint fit runs.

- **FiberExtract** (`extraction.py`): Performs boxcar and Horne (1986) optimal
  extraction using the global sky model directly (no local sky subtraction).
  Builds empirical spatial profiles from the flat field.

- **Fiber metadata** (`spectrograph.py`): The `get_fiber_metadata()` method
  allows spectrographs to map detected trace positions to instrument-defined
  fiber IDs, names, and types. This populates `MASKDEF_ID` and
  `MASKDEF_OBJNAME` on each `SpecObj`.

### Binospec IFU Notes

- MMT Binospec IFU (`mmt_binospec.py`) uses `pypeline = 'Fiber'`.
- Has 360 fibers on side A (DET01) and 356 on side B (DET02), including 40
  dedicated sky fibers per side.
- Spectral flexure correction is disabled (`spec_method = 'skip'`) because
  Binospec has active flexure control.
- Fiber identification uses cross-correlation against a reference profile
  (`fiber_ref_profile.fits`). Sky fibers are identified by `FIB_NAME`
  starting with `'SKY'` (the `FIB_TYPE` field in the reference file is
  unreliable — all fibers are marked as `SKY`).
- The `sky_fiber_indices_0based` class attribute contains *array position
  indices*, not physical `FIB_ID` values. Do not use `sky_fiber_ids`
  (= indices + 1) to match against `FIB_ID` from the reference profile.
- Header metadata cards: temperature is `TEMP`, humidity is `HUMID`,
  parallactic angle is `PA` (in `headarr[1]`).
- **Block-slit architecture**: Fibers are grouped into 21 block-slits per
  detector (42 total), each containing 8–20 fibers. Individual fibers are
  objects within their block-slit, not separate slits. Edge detection uses
  a high Sobel threshold (100σ) to find block boundaries at ~70 px
  inter-block gaps.
- Spatial illumination correction is not applied (no `IllumFlat`); fiber
  throughput variations are corrected post-extraction using per-fiber
  throughput weights derived from the flat field.
- Wavelength calibration runs per block-slit (42 total vs 720 individual
  fibers previously), providing significant performance improvement.
- Throughput corrections are applied post-extraction: sky fibers are
  extracted and throughput-corrected before building the 2D sky model.

## Testing

- All tests are collected in the `pypeit/tests` directory.

- Tests are written as plain module-level functions, not wrapped in `TestX`
  classes with `setup_class`. Use module-level constants or pytest fixtures
  for shared setup.

- Tests in `pypeit/tests` should be limited to unit tests that do not require
  the use of large data files.

- Tests should be deterministic; i.e., all random-number generators should use a
  fixed seed.

- Test coverage is supplemented by the PypeIt development suite, hosted at
  https://github.com/pypeit/PypeIt-development-suite, which requires data files
  hosted on Google Drive; see README.rst for the link.

## Documentation

- All functions and classes, except for tests in `pypeit/tests` should include
  docstrings that explain their purpose, input arguments, and output objects.

- The docstring style is currently not consistent within the repository, but
  Numpy style docstrings are preferred; see
  https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html

- Package documentation for users is held in the `doc` directory, which is built
  using Sphinx and hosted on ReadTheDocs at https://pypeit.readthedocs.io/.

- A complete rebuild of the documentation is performed by executing the bash
  command `cd doc; make clean ; make html`.  This requires access to the
  internet and the `PYPEIT_DEV` environmental variable, which points to the
  directory containing containing the `RAW_DATA` directory copied from the
  PypeIt development suite Google Drive (see "Testing" above).  If these
  requirements are not met, a limited rebuild of the documentation can be
  achieved by executing `cd doc ; make htmlonly`.

- Documentation should be updated with each GitHub pull request.

## Usage

- Users interact with the code base via execution of the command-line scripts
  found in `pypeit/scripts`.

- The primary data-processing script is `run_pypeit.py`, which primarily
  instantiates the `pypeit.pypeit.PypeIt` class and runs its methods.

## External Resources

- **Documentation**: https://pypeit.readthedocs.io/
- **Development Suite**: https://github.com/pypeit/PypeIt-development-suite

