"""
Preprocess up-the-ramp cubes into 2D count-rate images.

Each multi-read raw cube is fit up the ramp (with jump detection; see
:mod:`pypeit.ext.fitramp.fitramp`) and the result is written to the ``RampFit``
directory inside the reduction directory (``--odir``, defaulting to the
current directory), with the same file name as the raw cube.  The reduction
finds and reuses these files automatically — and creates them itself when
missing — so running this script is optional: it lets users inspect the
fitted images (units of e-/s) before a full reduction and front-loads the
fitting cost.

Up-the-ramp fitting is currently only implemented for MMT/MMIRS, so the
spectrograph must be one of :data:`RAMP_SPECTROGRAPHS`.

.. include:: ../include/links.rst
"""

from __future__ import annotations

import argparse

from pypeit.scripts import scriptbase

#: Spectrographs for which up-the-ramp preprocessing is implemented.  Adding a
#: new instrument here is the only change needed to enable it for this script,
#: provided the spectrograph implements the ramp-fitting interface used below.
RAMP_SPECTROGRAPHS = ('mmt_mmirs',)


class FitRamp(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width: int | None = None) -> argparse.ArgumentParser:
        from pypeit.spectrographs.util import available_spectrographs

        parser = super().get_parser(
            description='Preprocess up-the-ramp cubes into 2D count-rate '
                        'images (e-/s), written to the RampFit directory '
                        'inside the reduction directory.  Currently only '
                        f'supports: {", ".join(RAMP_SPECTROGRAPHS)}.',
            width=width,
            default_log_file=True)
        parser.add_argument('spectrograph', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))
        parser.add_argument('files', type=str, nargs='+',
                            help='One or more raw multi-read cubes')
        parser.add_argument('--odir', type=str, default='.',
                            help='Reduction directory in which the RampFit '
                                 'output directory is created (default: '
                                 'current directory)')
        parser.add_argument('--sig', type=float, default=None,
                            help='Force this single-read noise (e-) instead '
                                 'of calibrating it')
        parser.add_argument('--dark', type=str, default=None,
                            help='Raw dark cube used to calibrate the '
                                 'single-read noise (calibrated once, used '
                                 'for all files). Ignored if --sig is given.')
        parser.add_argument('--force', default=False, action='store_true',
                            help='Re-fit and overwrite existing up-to-date '
                                 'preprocessed images')
        return parser

    @classmethod
    def main(cls, args: argparse.Namespace) -> None:
        from pathlib import Path

        from pypeit import io, log, PypeItError
        from pypeit.spectrographs import mmt_mmirs
        from pypeit.spectrographs.util import load_spectrograph

        cls.init_log(args)

        spec = load_spectrograph(args.spectrograph)
        if spec.name not in RAMP_SPECTROGRAPHS:
            raise PypeItError(
                f'Up-the-ramp fitting is not implemented for {spec.name}; '
                f'supported spectrographs: {", ".join(RAMP_SPECTROGRAPHS)}.')

        if args.sig is not None:
            # Seeds the sigma cache: used for every frame
            spec._ramp_sigma = float(args.sig)
            log.info(f'Using forced single-read noise: {args.sig:.2f} e-')
        elif args.dark is not None:
            dark = Path(args.dark)
            with io.fits_open(dark) as dhdu:
                n_dark_reads = mmt_mmirs.mmirs_count_reads(dhdu)
            if n_dark_reads < spec.ramp_min_cal_groups:
                log.warning(f'{dark.name} has only {n_dark_reads} read(s), '
                            f'fewer than the {spec.ramp_min_cal_groups} '
                            'required for read-noise calibration; '
                            'self-calibration will be used instead')
            # Reuse the reduction's dark-calibration path (calibrated on
            # first use, then cached)
            spec._ramp_dark_files = [dark]

        for f in args.files:
            raw = Path(f)
            rampfit_file = mmt_mmirs.mmirs_rampfit_path(raw, args.odir)
            if not args.force and mmt_mmirs.mmirs_rampfit_fresh(rampfit_file,
                                                                raw):
                log.info(f'{raw.name}: up-to-date preprocessed image exists; '
                         'skipping (use --force to re-fit)')
                continue
            with io.fits_open(raw) as hdu:
                if hdu[0].header.get('RAMPFIT') is not None:
                    log.warning(f'{raw.name} is already a preprocessed '
                                'image; skipping')
                    continue
                n_reads = mmt_mmirs.mmirs_count_reads(hdu)
                if n_reads < spec.ramp_min_reads:
                    log.warning(f'{raw.name}: only {n_reads} read(s); '
                                f'up-the-ramp fitting requires at least '
                                f'{spec.ramp_min_reads}. Skipping.')
                    continue
                log.info(f'{raw.name}: fitting {n_reads} reads')
                detector_par = spec.get_detector_par(1, hdu=hdu)
                rate, sig, eff_ronoise = spec._ramp_fit_image(hdu,
                                                              detector_par)
                try:
                    mmt_mmirs.mmirs_write_rampfit(rampfit_file, rate, hdu,
                                                  sig, eff_ronoise,
                                                  raw.stat().st_mtime)
                except OSError as e:
                    log.error(f'{raw.name}: could not write {rampfit_file} '
                              f'({e}); skipping')
                    continue
            log.info(f'{raw.name}: single-read noise {sig:.2f} e-, effective '
                     f'read noise {eff_ronoise:.2f} e- -> {rampfit_file}')
