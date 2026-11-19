# `pypeit/wavemodel.py` dependency audit

Goal: identify which functions in `pypeit/wavemodel.py` can move into
`pypeit/core`, where a module may depend only on:

- the Python standard library,
- popular external packages (numpy, scipy, astropy, matplotlib),
- other `pypeit.core` modules.

The companion file `wavemodel_dependencies.csv` has the same information with
one row per function, ready to import into a spreadsheet.

## Module-level imports

| Import | Category | Used by |
|---|---|---|
| `os` | stdlib | **unused** |
| `re` | stdlib | `iraf_datareader` |
| `astropy` (bare; uses `astropy.constants`) | external | `blackbody` |
| `scipy` (bare; uses `scipy.interpolate`) | external | `transparency`, `nearIR_modelsky`, `optical_modelThAr`, `create_linelist` |
| `numpy` | external | all |
| `matplotlib.pyplot` | external | debug plots in 6 functions |
| `astropy.io.fits` | external | `nearIR_modelsky`, `optical_modelThAr` (optional output file) |
| `astropy.convolution` | external | `conv2res` |
| `astropy.table.Table` | external | `create_linelist` |
| `astropy.units` | external | `create_linelist` |
| `pypeit.core.arc` | core | `create_linelist` (`detect_lines`) |
| `pypeit.core.plot` | core | debug plots (`pyplot_rcparams*`) |
| `pypeit.core.wave` | core | `create_linelist` (`airtovac`) |
| `IPython.embed` | external | **unused** |
| `pypeit.log` | **non-core** | all 13 functions |
| `pypeit.dataPaths` | **non-core** | `oh_lines`, `transparency`, `h2o_lines` |
| `pypeit.io` | **non-core** | `thar_lines` (`load_thar_spec`) |

`import astropy` / `import scipy` followed by `astropy.constants.*` /
`scipy.interpolate.*` only works because some other module has already imported
those subpackages. Any move should import them explicitly.

## Per-function dependencies

Legend: **A** = can move to core after the log fix only;
**B** = needs a split or refactor (it reads PypeIt data files or calls a function
that does).

| # | Function | Line | External | `pypeit.core` | Non-core pypeit | Calls in this module | Status |
|---|---|---|---|---|---|---|---|
| 1 | `blackbody` | 31 | numpy, astropy.constants, mpl | plot | log | – | **A** |
| 2 | `addlines2spec` | 85 | numpy, mpl | plot | log | – | **A** |
| 3 | `oh_lines` | 145 | numpy | – | log, **dataPaths** | – | **B** |
| 4 | `transparency` | 160 | numpy, scipy.interpolate, mpl | plot | log, **dataPaths** | – | **B** |
| 5 | `h2o_lines` | 222 | numpy | – | log, **dataPaths** | – | **B** |
| 6 | `thar_lines` | 240 | numpy | – | log, **io.load_thar_spec** | – | **B** |
| 7 | `nearIR_modelsky` | 264 | numpy, scipy.interpolate, astropy.io.fits, mpl | plot | log | blackbody, oh_lines, addlines2spec, h2o_lines, conv2res | **B** (indirect dataPaths) |
| 8 | `optical_modelThAr` | 426 | numpy, scipy.interpolate, astropy.io.fits, mpl | plot | log | thar_lines, conv2res | **B** (indirect io) |
| 9 | `conv2res` | 545 | numpy, astropy.convolution, mpl | plot | log | – | **A** |
| 10 | `iraf_datareader` | 612 | re, numpy | – | log | – | **A** |
| 11 | `create_linelist` | 661 | numpy, scipy.interpolate, astropy.units, astropy.table | arc, wave | log | – | **A** |
| 12 | `create_OHlinelist` | 728 | numpy | – | log | nearIR_modelsky, create_linelist | **B** (indirect dataPaths) |
| 13 | `create_ThArlinelist` | 802 | numpy | – | log | optical_modelThAr, create_linelist | **B** (indirect io) |

### Internal call graph

```
create_OHlinelist ──► nearIR_modelsky ──► blackbody
        │                   ├──────────► oh_lines        [dataPaths]
        │                   ├──────────► addlines2spec
        │                   ├──────────► h2o_lines       [dataPaths]
        │                   └──────────► conv2res
        └──────────► create_linelist ──► core.arc.detect_lines, core.wave.airtovac

create_ThArlinelist ──► optical_modelThAr ──► thar_lines [pypeit.io]
        │                       └───────────► conv2res
        └──────────► create_linelist

transparency [dataPaths]        iraf_datareader        (both standalone)
```

### Callers outside the module

| Caller | Uses |
|---|---|
| `pypeit/core/flexure.py` | `wavemodel.nearIR_modelsky` |
| `pypeit/core/flux_calib.py` | `wavemodel.conv2res` |
| `pypeit/tests/test_wavemodel.py` | `wavemodel.transparency` |

Two `pypeit.core` modules already import `pypeit.wavemodel`, so `core` depends
on a non-core module today. Moving `conv2res` and a core version of
`nearIR_modelsky` would remove that.

## First-pass core candidates

**Group A: blocked only by `pypeit.log`** (5 functions):
`blackbody`, `addlines2spec`, `conv2res`, `iraf_datareader`, `create_linelist`.

**Group B: need a loader/compute split** (8 functions):

- **Pure data loaders:** `oh_lines`, `h2o_lines`, `thar_lines`. These stay
  outside core, or become core functions that take a file path or HDU while a
  thin non-core wrapper resolves the path with `dataPaths` / `io`.
- **`transparency`:** split it into a loader (non-core) and an
  interpolate-and-clip step (core) that takes `(wave_mod, tran_mod)`.
- **Model builders:** `nearIR_modelsky`, `optical_modelThAr`. A core version
  takes the line/template arrays as arguments. A non-core wrapper keeps the
  current signature for existing callers such as `core.flexure`.
- **Drivers:** `create_OHlinelist`, `create_ThArlinelist`. These stay outside
  core, or move once the model builders above have core versions.

## Open questions / cross-cutting issues

1. **Policy for `pypeit.log`.** Every function uses it, so this decision
   unblocks Group A. `pypeit.log` is a stdlib `logging.Logger` built by
   `pypeit.pkg.logger.get_logger`. A core module could use
   `logging.getLogger(__name__)` and still pick up PypeIt's handlers if the
   logger hierarchy is set up for that (this needs checking). Note that 34 of
   the 39 existing `pypeit/core/*.py` modules already import `pypeit.log`, so
   this is a repo-wide decision, not one specific to wavemodel.
2. **Transitive purity of core.** `core.arc`, which `create_linelist` needs,
   imports `pypeit.log`, `pypeit.PypeItError` and `pypeit.utils`. `core.wave`
   imports `pypeit.PypeItError`. These are allowed under the "any `pypeit.core`
   import is fine" rule, but the dependency tree as a whole is not yet pure.
3. **Debug plotting.** matplotlib is external, so it is allowed. If core
   should stay free of plotting, the `debug` branches in 6 functions could
   move to a QA/plot module.
4. **Latent bugs to fix during the move:**
   - `iraf_datareader`: `N_lines` is never set if no `features` line is
     found, and the path is built by string concatenation
     (`database_dir+id_file`).
   - The bare `import astropy` / `import scipy` noted above.
5. **Unused imports:** `os`, `IPython.embed`.
