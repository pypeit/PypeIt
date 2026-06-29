"""Interactive 1D fiber spectrum extractor for the MMT Binospec IFU.

Loads a single PypeIt ``spec1d`` file (covering both detector sides),
displays the IFU as a hexagonal flux map, and writes the user's combined
fiber selection as a :class:`~pypeit.onespec.OneSpec` FITS file.

.. include:: ../include/links.rst
"""
from __future__ import annotations

import argparse

from pypeit import PypeItError
from pypeit.scripts import scriptbase


class BinospecIFUExtract(scriptbase.ScriptBase):
    """Interactive fiber-spectrum extractor for the MMT Binospec IFU."""

    @classmethod
    def get_parser(cls, width: int | None = None) -> argparse.ArgumentParser:
        parser = super().get_parser(
            description='Interactive 1D fiber spectrum extractor for the '
                        'MMT Binospec IFU.',
            width=width,
            default_log_file=True)
        parser.add_argument('spec1d_file', type=str,
                            help='PypeIt spec1d FITS file (covering both '
                                 'detector sides)')
        parser.add_argument('-o', '--output', type=str, default=None,
                            help='Output OneSpec FITS filename '
                                 '(default: spec1d_<base>.fits → '
                                 'extract1d_<base>.fits)')
        parser.add_argument('--boxcar', default=False, action='store_true',
                            help='Use boxcar (BOX) extraction columns '
                                 'instead of optimal (OPT)')
        return parser

    @classmethod
    def main(cls, args: argparse.Namespace) -> None:
        from pathlib import Path

        import numpy as np
        from astropy.io import fits

        from pypeit import log
        from pypeit.core import datacube
        from pypeit.gui.binospec_ifu_extract import BinospecIFUExtractGUI
        from pypeit.specobjs import SpecObjs
        from pypeit.spectrographs.util import load_spectrograph

        cls.init_log(args)

        if not Path(args.spec1d_file).name.startswith('spec1d_'):
            raise PypeItError(
                f"Only spec1d files are supported; got "
                f"{Path(args.spec1d_file).name}")

        log.info(f"Loading {args.spec1d_file}")
        sobjs = SpecObjs.from_fitsfile(args.spec1d_file)
        if sobjs.nobj == 0:
            raise PypeItError(f"No objects in {args.spec1d_file}")

        raw_hdr = fits.getheader(args.spec1d_file)

        spectrograph = load_spectrograph(raw_hdr['PYP_SPEC'])
        if spectrograph.name != 'mmt_binospec_ifu':
            raise PypeItError(
                f"This script requires an mmt_binospec_ifu spec1d; "
                f"got PYP_SPEC={spectrograph.name!r}")
        targetx, targety = spectrograph.load_sky_layout()

        prefix = 'BOX' if args.boxcar else 'OPT'
        log.info(f"  Using {prefix} extraction")

        fibers = datacube.load_fibers(sobjs, spectrograph, targetx, targety,
                                      prefix=prefix)
        log.info(f"  Loaded {len(fibers)} science fibers")

        # Project each fiber's instrument-frame x,y to sky.
        x = np.array([f['x'] for f in fibers])
        y = np.array([f['y'] for f in fibers])
        ra, dec = datacube.project_to_sky(x, y, raw_hdr, spectrograph)
        for f, r, d in zip(fibers, ra, dec):
            f['ra'] = float(r)
            f['dec'] = float(d)

        # Output filename
        if args.output is not None:
            outfile = args.output
        else:
            base = Path(args.spec1d_file).stem
            if 'spec1d_' in base:
                base = base.replace('spec1d_', 'extract1d_')
            outfile = base + '.fits'

        gui = BinospecIFUExtractGUI(fibers, raw_hdr, spectrograph.name, outfile)
        gui.show()
