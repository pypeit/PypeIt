"""
Build a datacube from Binospec IFU spec1d files.

Thin alias of the generic fiber-IFU cube script, :ref:`pypeit_fiber_cube`
(:class:`~pypeit.scripts.fiber_cube.FiberCube`), restricted to the MMT
Binospec IFU and with its historical 0.27 arcsec default spatial scale:

1. Reading extracted 1D fiber spectra from spec1d files (already
   sky-subtracted by the pipeline)
2. Mapping 320 science fibers per side to sky positions
3. Combining both detectors (640 science fibers total)
4. Interpolating scattered fiber positions onto a regular spatial grid

Each input file produces a separate output datacube.

.. include:: ../include/links.rst
"""

from __future__ import annotations

import argparse

from pypeit.scripts.fiber_cube import FiberCube


class BinospecIFUCube(FiberCube):
    """
    Binospec alias of :class:`~pypeit.scripts.fiber_cube.FiberCube`.
    """

    required_spectrograph = 'mmt_binospec_ifu'
    default_spatial_scale = 0.27

    @classmethod
    def get_parser(cls, width: int | None = None) -> argparse.ArgumentParser:
        """
        Construct the command-line argument parser (unchanged options).

        Args:
            width (:obj:`int`, optional):
                Restrict the width of the formatted help output to be no
                longer than this number of characters.

        Returns:
            `argparse.ArgumentParser`_: Command-line interpreter.
        """
        # Bypass FiberCube.get_parser to keep the original Binospec options
        parser = super(FiberCube, cls).get_parser(
            description='Build a datacube from Binospec IFU spec1d files.',
            width=width, default_log_file=True)
        cls.add_common_arguments(parser)
        parser.add_argument('--spatial_scale', type=float, default=0.27,
                            help='Output spatial pixel scale in arcsec '
                                 '(default: 0.27)')
        return parser
