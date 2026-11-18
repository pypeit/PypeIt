# PypeIt

PypeIt is a pure Python package for processing raw spectroscopic data from
astronomical telescopes into calibrated spectra for scientific analysis.

## Development Workflow

- Install in editable mode with development dependencies: `pip install -e ".[dev]"`.
  Supported Python versions are `>=3.11,<3.15` (see `pyproject.toml`).

- Branching: simple bug fixes branch from `release`; all other development
  branches from `develop`.  The main remote is `upstream`
  (https://github.com/pypeit/PypeIt).

- Logging and errors: use `from pypeit import log` (`log.info`, `log.warning`,
  ...) and raise `PypeItError` (`from pypeit import PypeItError`).  The legacy
  `msgs` interface (`msgs.error`, `msgs.warn`, ...) has been removed; do not
  use it.

- Additional development guidelines are provided in `doc/dev/development.rst`.

### Common commands

- Run unit tests: `pytest pypeit/tests` (or a single file, e.g.
  `pytest pypeit/tests/test_metadata.py`); `tox` runs the configured test
  environments.
- Code-style check: `tox -e codestyle` (currently only checks syntax errors,
  `E9`).
- Quick documentation build (no network or dev-suite data needed):
  `cd doc ; make htmlonly`.

### When changing X, also update Y

- **New or changed user-level parameter** (`pypeit/par/pypeitpar.py`): update
  the parameter docstrings/options and rebuild the docs so the auto-generated
  parameter tables are regenerated.
- **New or changed datamodel** (a `DataContainer` subclass): bump its version
  as appropriate and rebuild the docs so the datamodel tables are regenerated.
- **New command-line script** (`pypeit/scripts/`): follow the `ScriptBase`
  pattern and register the entry point under `[project.scripts]` in
  `pyproject.toml`.
- **New pypeline**: add the pypeline string to the validation lists in
  `specobj.py`, `specobjs.py`, `slittrace.py`, `scripts/show_2dspec.py`,
  `spectrographs/spectrograph.py`, and `pypeit_steps.py`, and implement the
  `<Pypeline>FindObjects` and `<Pypeline>Extract` classes (see Architecture).
- **Any user-facing change**: add an entry to the in-development release notes
  in `doc/releases/` (e.g. `3.1.0dev.rst`).  `CHANGES.rst` is deprecated.
- **Any pull request**: update the relevant documentation in `doc/`.

## Coding Conventions

**These conventions apply to new code.**  Much of the existing code base does
not yet follow them; it is migrated gradually, and only when a developer asks
for it.  Do not reformat or restyle code that you are not otherwise changing.

### Formatting

- **Lines are at most 99 characters.**  This is not yet enforced by
  `tox -e codestyle`; it will be once enough of the code base has migrated.
- **Anything that wraps onto multiple lines** (a signature, a call, an `if`
  condition, a collection literal) puts its contents on the line *below* the
  opening construct, indented, with the closing bracket (and any trailing
  colon) alone on a line aligned with the line that opened it:

  ```python
  def fit_spatial_profile(
      self, image, ivar, slit_index, gpm=None, nbreak=10, maxiter=20, debug=False
  ):

  profile, flag = self.fit_spatial_profile(
      image=self.sciImg.image, ivar=self.sciImg.ivar, slit_index=i, gpm=gpm,
      maxiter=maxiter
  )

  if (
      self.par['reduce']['skysub']['joint_fit'] and self.spectrograph.pypeline == 'Fiber'
  ):
  ```

  Pack as many arguments or elements onto each line as fit within the line
  limit; do not put one per line.  Implicitly concatenated strings continue on
  the following lines, aligned under the first string.  Do not add a trailing
  comma after the last element:

  ```python
  model = fit_profile_model(
      xdata=xdata, ydata=ydata, ivar=ivar, basis=profile_basis, gpm=gpm, upper=upper,
      lower=lower, maxiter=maxiter
  )

  parser.add_argument(
      '-r', '--rdxdir', type=str,
      help='Path for the reduced data.  If not provided, set to the parent '
           'directory of rawdir.'
  )
  ```
- **Do not reformat code with `black`** or a similar tool.  Formatting changes
  are limited to line breaks and indentation; keep the existing spacing and
  literals, e.g., compact slices like `[:,i]` (no space after the comma) and
  floats like `3.` (not `3.0`).

### Logic and numerics

- **Always pass arguments by keyword** in any call with more than one
  argument, e.g. `np.sort(arr, kind='stable')`, not `np.sort(arr, 'stable')`.
  This guards against signature changes and makes call sites self-documenting.
- **Always use explicit condition checks.**  Never rely on the truthiness or
  falsiness of a non-boolean object: write `if x is None`, `if len(x) == 0`,
  `if x != 0`, not `if not x` or `if x`.  Truthiness on an actual `bool` is
  fine.
- **In a conditional expression with a `None` branch, put `None` first:**
  `None if x is None else f(x)`, not `f(x) if x is not None else None`.
- **Never use NaN intentionally.**  A value that is invalid, missing or to be
  ignored is flagged with a mask, never by setting it to NaN.  A NaN always
  signifies a computational error.  Do not write NaN as a placeholder, a fill
  value, an "undefined" result, or a gap in a plotted line.  The one
  exception: a masked float column written to a FITS binary table is stored as
  NaN (the FITS null for floats) by astropy's writer, and `Table.read` masks it
  again; the code still uses masked columns and never writes NaN itself.
- **Use `pathlib`, not `os.path`,** for paths and file-system operations
  (e.g. `Path(...).name`, `.stem`, `.with_name()`, `.exists()`).
- **Always sort directory listings.**  The order returned by `Path.glob`,
  `Path.rglob`, and `Path.iterdir` depends on the platform and file system, so
  wrap them in `sorted(...)`.  The only exception is when the order cannot
  matter, e.g., counting the files or deleting all of them.

### Structure and imports

- **No nested function definitions.**  A helper belongs at module level,
  private by a leading underscore if it is not part of the interface; pass the
  values it needs explicitly (e.g., via `args=` to `scipy.optimize.minimize`).
  The exception is test code, where a closure over the fixture under test is
  often the clearest way to express a check.
- **Keep module-level globals to a minimum.**  Do not define a global constant
  unless (1) it is used in many different functions, or (2) it is in test
  code.  Otherwise compute the value where it is used and explain it in that
  function's docstring.
- **Imports go at the top of the file.**  The one exception is a script: the
  imports a `ScriptBase` subclass needs belong inside its `main`, because every
  registered entry point is imported when the console scripts are resolved, so
  a module-level import of the heavy machinery makes `--help` pay for it.
- **Never import private functions into other modules.**  A function (or
  class, or variable) with a leading underscore is private to its module: do
  not import it elsewhere, or call it as `module._name`.  If another module
  needs it, make it public by removing the underscore (and document it as part
  of the interface).  The exception is test code, which may exercise private
  helpers directly.
- **Avoid `import ... as`.**  Industry-standard aliases (`import numpy as np`)
  are fine; otherwise import the module or name directly.
- **Adapted code must have a license compatible with PypeIt's BSD-3 license**
  and carries an attribution note in its module docstring naming the source
  project and file.  Never copy code whose license prohibits redistribution;
  import it instead.

### Docstrings and language

- Docstrings are **NumPy style**, rendered by Sphinx's `napoleon` extension;
  see https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html
- **The opening `"""` is alone on its line**; the docstring text starts on the
  next line:

  ```python
  def foo():
      """
      Do the thing.
      """
  ```
- **Docstring text wraps at 79 characters** (stricter than the 99-character
  code limit).
- **Any docstring containing a backslash is an r-string** (`r"""`), so that
  LaTeX in `:math:` roles reads as written rather than being double-escaped.
- **Cross-references use intersphinx roles**, e.g.
  `` :class:`numpy.ndarray` ``, `` :func:`astropy.stats.sigma_clip` ``,
  `` :class:`~pypeit.specobjs.SpecObjs` ``.  Existing docstrings use the older
  link-target style (`` `numpy.ndarray`_ ``); leave those in place unless the
  developer asks for them to be converted.  Bare package names in prose are
  literals: ``` ``astropy`` ```.
- **`ln` is the natural logarithm and `log` the base-10 logarithm** in
  docstrings, comments, and the documentation (e.g., `:math:` roles).
- **Write in American English**, in names, comments, docstrings, messages, and
  the documentation: `color`, `behavior`, `labeled`, `normalize`, `center`,
  `catalog`.  The exception is quoted text, which keeps its original spelling.
- Parameter, datamodel and bitmask documentation is **generated by reflection**
  so it cannot drift from the code.  Keep that property.

## Architecture

- PypeIt implements four processing paths (pypelines), depending on the format
  of the spectrograph used to collect the raw data:

    - **MultiSlit**: Standard long-slit and multi-slit spectrographs
    - **Echelle**: Cross-dispersed echelle spectrographs
    - **SlicerIFU**: Slicer-based integral-field units
    - **Fiber**: Fiber-fed spectrographs (e.g., IFU fiber bundles)

- Processing flow: `run_pypeit` instantiates `pypeit.pypeit.PypeIt`, which
  runs calibrations (`Calibrations`), object finding (`FindObjects`), and
  extraction (`Extract`), writing `spec2d` and `spec1d` (`SpecObjs`) outputs.
  Post-processing (fluxing, coadding, etc.) is done by separate scripts.
  Lower-level numerical algorithms live in `pypeit/core/`; the top-level
  modules wrap them in classes.

- Instrument specifications and data-processing considerations that are
  specific to each spectrograph (including the pypeline to use) are isolated
  in their respective spectrograph classes, all of which inherit from
  `pypeit.spectrographs.spectrograph.Spectrograph`.  Instrument-specific
  behavior belongs in a spectrograph method or hook, not in branches on the
  spectrograph name in core code.

- Pypeline dispatch uses class-name conventions: `FindObjects.get_instance()`
  and `Extract.get_instance()` select the subclass named
  `<pypeline>FindObjects` / `<pypeline>Extract` (e.g. `FiberFindObjects`).
  `Calibrations.get_instance()` uses `MultiSlitCalibrations` for `MultiSlit`
  and `Echelle`, and `IFUCalibrations` for all other pypelines.

- Users primarily modify the code performance via user-level parameters held
  by the `pypeit.par.pypeitpar.PypeItPar` class, which packages the hierarchy
  of parameter sets used throughout the code.

- Frame metadata is held in `pypeit.metadata.PypeItMetaData` (the `fitstbl`).
  Select and locate frames with its query methods (`keys`, `find_frames`,
  `find_frame_files`, `frame_paths`, `find_calib_group`, etc.) rather than
  hand-rolling operations on the underlying `fitstbl.table`.  Frame types (the
  `framebit`/`frametype` columns) are only assigned once
  `get_frame_types`/`finalize_usr_build` has run, which is *after*
  `PypeItMetaData.__init__` (and therefore after the
  `Spectrograph.cache_metadata` hook); defer type-based lookups until after
  typing.

- Core processing modules, particularly for calibrations, produce FITS files
  that are saved to disk and reused as necessary.  These modules commonly use
  `pypeit.datamodel.DataContainer` as a base class to enforce strict adherence
  to a well-defined datamodel and to provide a common IO interface.

- Significant portions of the `pypeit/data` directory are not included in the
  package distribution, but rely on a cache system (`pypeit/pkg/cache.py`,
  `pypeit/pkg/pypeitdata.py`) that downloads files as needed.  Do not assume a
  data file exists locally just because the code references it.

- Dependency direction: core processing modules and `pypeit/spectrographs`
  must not import from `pypeit.scripts`, `pypeit.gui`, or
  `pypeit.display.qlview` (these depend on the core, not the reverse).  The
  ginga helpers in `pypeit.display.display` are the exception.

## Testing

- All unit tests are collected in the `pypeit/tests` directory.

- Tests are written as plain module-level functions, not wrapped in `TestX`
  classes with `setup_class`.  Use module-level constants or pytest fixtures
  for shared setup.

- Tests in `pypeit/tests` should be limited to unit tests that do not require
  large data files.  Tests that need real data belong in the PypeIt
  development suite.

- Tests should be deterministic; i.e., all random-number generators should use
  a fixed seed.

- **Every `assert` in a test carries a message describing the expected
  behavior.**  A bare assertion reports only that a comparison failed; the
  message is what tells the next reader what was supposed to happen.

- **Test data goes in `pypeit/data/tests/`, never under `pypeit/tests/`.**
  `pypeit/tests/` holds test *code* only.

- Test coverage is supplemented by the PypeIt development suite, hosted at
  https://github.com/pypeit/PypeIt-development-suite, which requires data
  files hosted on Google Drive; see its README.rst for the link.

## Documentation

- All functions and classes, except for tests in `pypeit/tests`, should
  include docstrings that explain their purpose, input arguments, and output
  objects.

- The docstring style is currently not consistent within the repository; new
  docstrings follow the rules under "Coding Conventions" above.

- Package documentation for users is held in the `doc` directory, which is
  built using Sphinx and hosted on ReadTheDocs at
  https://pypeit.readthedocs.io/.

- A complete rebuild of the documentation is performed by executing
  `cd doc ; make clean ; make html`.  This requires internet access and the
  `PYPEIT_DEV` environment variable, which points to the directory containing
  the `RAW_DATA` directory copied from the PypeIt development suite Google
  Drive (see "Testing" above).  If these requirements are not met, a limited
  rebuild can be done with `cd doc ; make htmlonly`.

## Usage

- Users interact with the code base via the command-line scripts in
  `pypeit/scripts`.

- The primary data-processing script is `run_pypeit.py`, which instantiates
  the `pypeit.pypeit.PypeIt` class and runs its methods.

## External Resources

- **Documentation**: https://pypeit.readthedocs.io/
- **Development Suite**: https://github.com/pypeit/PypeIt-development-suite
