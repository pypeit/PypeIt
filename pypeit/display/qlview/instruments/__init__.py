"""
Instrument classes used by the quicklook viewer to configure the frontend.

These classes serve two purposes:

1. Each instrument registered in
:class:`~pypeit.display.qlview.instruments.registry.InstrumentRegistry` is
displayed as an option in the instrument dropdown.

2. They provide configuration information needed to display instrument info,
for example column headers, rendering the raw image, etc.

All instruments subclass
:class:`~pypeit.display.qlview.instruments.base.Instrument`, and each lives in
its own module named after the corresponding PypeIt spectrograph (e.g.
``keck_deimos.py``).  To add an instrument, create a new module here and add
its class to the ``InstrumentRegistry``; classes are not seen by the viewer
until they are registered!
"""

from .base import Instrument
from .registry import InstrumentRegistry

__all__ = ["Instrument", "InstrumentRegistry"]
