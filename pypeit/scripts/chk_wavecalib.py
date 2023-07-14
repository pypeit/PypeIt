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

        parser.add_argument('input_file', type=str,
                            help='PypeIt WaveCalib file [e.g. WaveCalib_A_1_DET01.fits] or '
                                 'spec2d file')
        return parser

    @staticmethod
    def main(args):

        from astropy.io import fits
        from IPython import embed
        from pypeit import wavecalib, spec2dobj

        try:
            # Load
            waveCalib = wavecalib.WaveCalib.from_file(args.input_file)
                                                    #, chk_version=(not args.try_old))
        # TODO: Should this specify the type of exception to pass?
        except:
            pass
        else:
            waveCalib.wave_diagnostics(print_diag=True)
            return

        try:
            allspec2D = spec2dobj.AllSpec2DObj.from_fits(args.input_file, chk_version=False)
        # TODO: Should this specify the type of exception to pass?
        except:
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
            return

        # Should not get here unless it can't read either file type
        raise IOError("Unrecognized file type. Must be a WaveCalib or spec2d file.")

