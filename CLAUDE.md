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

- Additional development guidelines is provided in `doc/dev/development.rst`.

## Key Architectural Components

- PypeIt implements three primary processing paths, depending on the format of
  the spectrograph used to collect the raw data:

    - **MultiSlit**: Standard long-slit and multi-slit spectrographs
    - **Echelle**: Cross-dispersed echelle spectrographs
    - **SlicerIFU**: Slicer-based integral-field units

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

## Testing

- All tests are collected in the `pypeit/tests` directory.

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

