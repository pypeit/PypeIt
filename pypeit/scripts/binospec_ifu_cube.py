"""
Build a datacube from Binospec IFU spec1d files.

Unlike the general ``pypeit_coadd_datacube`` script designed for slicer-based
IFUs, this script handles the fiber-fed Binospec IFU by:

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

from pypeit.scripts import scriptbase


class BinospecIFUCube(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width: int | None = None) -> argparse.ArgumentParser:
        parser = super().get_parser(
            description='Build a datacube from Binospec IFU spec1d files.',
            width=width,
            default_log_file=True)
        parser.add_argument('files', type=str, nargs='+',
                            help='One or more PypeIt spec1d files, or a text '
                                 'file listing them (one per line)')
        parser.add_argument('-o', '--output', type=str, default=None,
                            help='Output FITS filename (only valid for a single '
                                 'input file; default: auto-generated)')
        parser.add_argument('--spatial_scale', type=float, default=0.27,
                            help='Output spatial pixel scale in arcsec (default: 0.27)')
        parser.add_argument('--boxcar', default=False, action='store_true',
                            help='Use boxcar (BOX) extraction columns instead '
                                 'of the default optimal (OPT) columns')
        parser.add_argument('--method', type=str, default='linear',
                            choices=['nearest', 'linear', 'cubic'],
                            help='Spatial interpolation method (default: linear)')
        return parser

    @classmethod
    def main(cls, args: argparse.Namespace) -> None:
        from pathlib import Path

        from astropy.io import fits

        from pypeit import log, PypeItError
        from pypeit.core import datacube
        from pypeit.spectrographs.util import load_spectrograph

        cls.init_log(args)

        # ----------------------------------------------------------------
        # Parse input files
        # ----------------------------------------------------------------
        input_files: list[str] = []
        for f in args.files:
            if f.endswith('.txt'):
                with open(f, 'r') as fh:
                    for line in fh:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            input_files.append(line)
            else:
                input_files.append(f)

        if len(input_files) == 0:
            raise PypeItError("No input files provided.")

        if args.output is not None and len(input_files) > 1:
            raise PypeItError("--output can only be used with a single input file.")

        # Only spec1d files are supported
        for f in input_files:
            if not Path(f).name.startswith('spec1d'):
                raise PypeItError(
                    f"Only spec1d files are supported; got {Path(f).name}")

        log.info(f"Processing {len(input_files)} spec1d file(s)")

        # Load spectrograph and fiber layout (shared across all files)
        with fits.open(input_files[0]) as hdu:
            spectrograph = load_spectrograph(hdu[0].header['PYP_SPEC'])
        if spectrograph.name != 'mmt_binospec_ifu':
            raise PypeItError(
                f"This script requires an mmt_binospec_ifu spec1d; "
                f"got PYP_SPEC={spectrograph.name!r}")

        targetx, targety = spectrograph.load_sky_layout()

        # ----------------------------------------------------------------
        # Process each file into a separate datacube
        # ----------------------------------------------------------------
        for input_file in input_files:
            log.info(f"Building datacube for {Path(input_file).name}")
            datacube.build_cube_from_spec1d(
                input_file, spectrograph, targetx, targety,
                boxcar=args.boxcar, spatial_scale=args.spatial_scale,
                method=args.method, output=args.output)
