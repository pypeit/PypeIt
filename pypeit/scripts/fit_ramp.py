"""
Preprocess up-the-ramp cubes into 2D count-rate images.

Given a :ref:`pypeit_file`, each multi-read raw cube it lists is fit up the ramp
(with jump detection; see :mod:`pypeit.ext.fitramp.fitramp`) and the result is
written to the ramp-fit directory inside the reduction directory (the ``[rdx]``
``redux_path`` and ``rampfit_dir`` parameters), with the same file name as the
raw cube.  ``run_pypeit`` finds and reuses these files automatically — and
creates them itself when missing — so running this script is optional; it lets
users inspect the fitted images (units of e-/s) before a full reduction and
front-loads the fitting cost.  The intended workflow is::

    pypeit_setup ...
    pypeit_fit_ramp blah.pypeit   # optional
    run_pypeit blah.pypeit

Because the script reads the same pypeit file as the reduction, it uses the same
reduction directory, ramp-fit directory, and dark frames (for the per-read noise
calibration) automatically.

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
        parser = super().get_parser(
            description='Preprocess up-the-ramp cubes into 2D count-rate '
                        'images (e-/s) ahead of a reduction, using the same '
                        'pypeit file that run_pypeit will use.  This step is '
                        'optional: run_pypeit fits any ramp it does not find '
                        'already preprocessed.  Currently only supports: '
                        f'{", ".join(RAMP_SPECTROGRAPHS)}.',
            width=width,
            default_log_file=True)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (see pypeit_setup).  '
                                 'The raw frames it lists are fit up the ramp '
                                 'and written to the reduction directory (the '
                                 '[rdx] redux_path and rampfit_dir), where '
                                 'run_pypeit reuses them.')
        parser.add_argument('--force', default=False, action='store_true',
                            help='Re-fit and overwrite existing up-to-date '
                                 'preprocessed images')
        return parser

    @classmethod
    def main(cls, args: argparse.Namespace) -> None:
        from pathlib import Path

        import numpy as np

        from pypeit import inputfiles, PypeItError
        from pypeit.metadata import PypeItMetaData

        cls.init_log(args)

        # Read the pypeit file and build the metadata table.  Building the
        # metadata calls the spectrograph's cache_metadata(), which configures
        # the ramp-fit output directory ([rdx] redux_path + rampfit_dir), the
        # dark frames used to calibrate the per-read noise, and the threading
        # parameters -- exactly as run_pypeit does, so the preprocessed images
        # land where the reduction will look for them.
        pypeitFile = inputfiles.PypeItFile.from_file(args.pypeit_file)
        # Guard on the spectrograph name before reading any frame.
        if pypeitFile.get_spectrograph().name not in RAMP_SPECTROGRAPHS:
            raise PypeItError(
                'Up-the-ramp fitting is not implemented for '
                f'{pypeitFile.get_spectrograph().name}; supported '
                f'spectrographs: {", ".join(RAMP_SPECTROGRAPHS)}.')

        spec, par, _ = pypeitFile.get_pypeitpar()
        fitstbl = PypeItMetaData(spec, par, files=pypeitFile.filenames,
                                 usrdata=pypeitFile.data, strict=True)
        fitstbl.finalize_usr_build(pypeitFile.frametypes, pypeitFile.setup_name)

        redux_path = Path(par['rdx']['redux_path'])
        rampfit_dir = par['rdx']['rampfit_dir']

        # Fit each listed frame through the spectrograph interface, so the
        # script stays agnostic of the instrument-specific ramp code.
        for raw in map(Path, fitstbl.frame_paths(np.arange(len(fitstbl)))):
            spec.preprocess_ramp_file(raw, redux_path, rampfit_dir,
                                      force=args.force)
