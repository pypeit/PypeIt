"""
The PypeIt reduction-**state** package.

This package holds the data model and I/O for the reduction state recorded
while ``run_pypeit`` runs (and read by the PypeIt Dashboard and
``pypeit_status``):

- :mod:`pypeit.state.run_state` — the pydantic data model
  (:class:`~pypeit.state.run_state.RunPypeItState` and the per-step / per-frame
  classes) and its JSON I/O.  This is Qt-free and free of product imports.
- :mod:`pypeit.state.science_status` — helpers that translate the science
  reduction products (``SpecObjs`` / ``AllSpec2DObj`` / slit bitmasks), or the
  on-disk products, into the science portion of the state.

Following the package-wide import convention, this ``__init__`` only lists
the submodules; import the classes you need directly from them (e.g.
``from pypeit.state.run_state import RunPypeItState``) to avoid import
cycles.
"""

__all__ = ['run_state', 'science_status']
