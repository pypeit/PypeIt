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

        # What kind of file are we?? Similar to pypeit_parse_slits
        hdul = fits.open(args.input_file)
        head0 = hdul[0].header
        if 'CALIBTYP' in head0.keys() and head0['CALIBTYP'].strip() == 'WaveCalib':
            file_type = 'WaveCalib'
        elif 'PYP_CLS' in head0.keys() and head0['PYP_CLS'].strip() == 'AllSpec2DObj':
            file_type = 'AllSpec2D'
        else:
            raise IOError("Unrecognized file type. Must be a WaveCalib or spec2d file.")

        if file_type == 'WaveCalib':
            # Load
            waveCalib = wavecalib.WaveCalib.from_file(args.input_file)
                                                    #, chk_version=(not args.try_old))
            # print
            waveCalib.wave_diagnostics(print_diag=True)

        elif file_type == 'AllSpec2D':
            # Load
            allspec2D = spec2dobj.AllSpec2DObj.from_fits(args.input_file, chk_version=False)
            for det in allspec2D.detectors:
                print('='*44 + f'{det:^7}' + '='*44)
                wave_diag = allspec2D[det].wavesol
                for colname in ['minWave', 'Wave_cen', 'maxWave', 'IDs_Wave_cov(%)']:
                    wave_diag[colname].format = '0.1f'
                for colname in ['dWave', 'RMS']:
                    wave_diag[colname].format = '0.3f'
                print(wave_diag)


