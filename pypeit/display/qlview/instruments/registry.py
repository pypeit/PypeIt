"""
Registry of the instruments offered in the quicklook viewer's instrument dropdown.
"""

from __future__ import annotations

from typing import List

from .base import Instrument
from .keck_deimos import DEIMOS
from .keck_mosfire import MOSFIRE


class InstrumentRegistry:
    def __init__(self, logger) -> None:
        """Initialise the registry and register all built-in instrument classes.

        Parameters
        ----------
        logger : logging.Logger
            Ginga application logger, forwarded to each
            :class:`Instrument` instance created via :meth:`create`.
        """
        self.logger = logger
        # Display name -> Instrument class.  The commented-out instruments are
        # untested; to enable one, uncomment it here and import its class
        # from the corresponding keck_*.py module.
        self._registry = {
            "DEIMOS": DEIMOS,
            # "HIRES": HIRES,
            # "LRIS Blue": LRISBlue,
            # "LRIS Red": LRISRed,
            "MOSFIRE": MOSFIRE,
            # "NIRES": NIRES,
            # "NIRSPEC": NIRSPEC,
        }

    def create(self, name: str) -> Instrument:
        """Instantiate and return the :class:`Instrument` for *name*.

        Parameters
        ----------
        name : str
            Display name as it appears in the instrument combo box
            (e.g. ``"DEIMOS"``, ``"LRIS Blue"``).

        Returns
        -------
        Instrument
            A freshly constructed instrument instance.

        Notes
        -----
        Falls back to :class:`DEIMOS` and logs an error when *name* is not
        found in the registry, so callers always receive a usable object.
        """
        cls = self._registry.get(name)
        if cls is None:
            self.logger.error(f"Instrument not recognized: {name}")
            cls = DEIMOS
        return cls(self.logger)

    def names(self) -> List[str]:
        """Return the list of registered instrument display names.

        Returns
        -------
        list of str
            Names in insertion order, matching the order shown in the
            instrument combo box.
        """
        return list(self._registry.keys())

    def instrume_values(self) -> List[tuple]:
        """Return ``(display_name, instrume_value)`` pairs for all registered instruments.

        Reads ``instrume_value`` directly from each class object rather than
        constructing an instance, making this safe to call when you only need
        the FITS keyword value for matching purposes.

        Returns
        -------
        list of (str, str)
            ``(display_name, instrume_value)`` tuples in registration order.
        """
        return [(name, cls.instrume_value) for name, cls in self._registry.items()]
