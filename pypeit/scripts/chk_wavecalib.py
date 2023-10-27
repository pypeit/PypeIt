"""
This script displays the wavelength calibration diagnostics.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ChkWaveCalib(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Print QA on Wavelength Calib to the screen',
                                    width=width)

        parser.add_argument('input_file', type=str, nargs='+',
                            help='One or more PypeIt WaveCalib file [e.g. WaveCalib_A_1_DET01.fits] or '
                                 'spec2d file')
        return parser

    @staticmethod
    def main(args):

        from IPython import embed
        from pypeit import wavecalib, spec2dobj
        from pypeit.pypmsgs import PypeItError

        # Loop over the input files
        for in_file in args.input_file:
            try:
                # Load
                waveCalib = wavecalib.WaveCalib.from_file(in_file)
                                                        #, chk_version=(not args.try_old))
            except FileNotFoundError:
                pass
            else:
                waveCalib.wave_diagnostics(print_diag=True)
                continue

            try:
                allspec2D = spec2dobj.AllSpec2DObj.from_fits(in_file, chk_version=False)
            except PypeItError:
                pass
            else:
                for det in allspec2D.detectors:
                    print('='*44 + f'{det:^7}' + '='*44)
                    wave_diag = allspec2D[det].wavesol
                    for colname in ['minWave', 'Wave_cen', 'maxWave', 'IDs_Wave_cov(%)']:
                        wave_diag[colname].format = '0.1f'
                    for colname in ['dWave', 'RMS']:
                        wave_diag[colname].format = '0.3f'
                    print(wave_diag)
                continue

            # Should not get here unless it can't read either file type
            raise IOError("Unrecognized file type. Must be a WaveCalib or spec2d file.")

