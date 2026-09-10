# fitramp

Up-the-ramp fitting of non-destructive detector reads, with likelihood-based
jump (cosmic-ray) detection.

This is third-party code vendored into PypeIt, **not** developed by the PypeIt
team.

- **Source:** https://github.com/t-brandt/fitramp
- **Author:** Timothy D. Brandt
- **License:** MIT (see `LICENSE` in this directory), Copyright (c) 2023-2024
  Timothy D. Brandt
- **Algorithm:** described in Brandt (2024a), PASP 136, 045005
  (https://arxiv.org/abs/2404.01326) and Brandt (2024b)
  (https://arxiv.org/abs/2309.08753)

## Do not refactor

`fitramp.py` is kept **byte-identical** to the upstream file so that upstream
fixes can be pulled in with a clean diff.  Do not reformat it or adapt it to
PypeIt conventions.  Any bug fixes should be sent upstream first and then
copied back here.

## Usage in PypeIt

Currently used only by the MMT/MMIRS spectrograph
(`pypeit.spectrographs.mmt_mmirs`) for up-the-ramp fitting of raw ramp cubes:

```python
from pypeit.ext.fitramp import fitramp
```
