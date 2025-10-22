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
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        from IPython import embed
        from astropy.io import fits
        from pypeit import wavecalib, spec2dobj, msgs

        chk_version = not args.try_old

        # Loop over the input files
        for in_file in args.input_file:

            # What kind of file are we??
            hdul = fits.open(in_file)
            head0 = hdul[0].header
            head1 = hdul[1].header
            file_type = None
            if 'CALIBTYP' in head1.keys() and head1['CALIBTYP'].strip() == 'WaveCalib':
                file_type = 'WaveCalib'
            elif 'PYP_CLS' in head0.keys() and head0['PYP_CLS'].strip() == 'AllSpec2DObj':
                file_type = 'AllSpec2D'
            else:
                msgs.error("Bad file type input!")

            if file_type == 'WaveCalib':
                waveCalib = wavecalib.WaveCalib.from_file(in_file, chk_version=chk_version)
                waveCalib.wave_diagnostics(print_diag=True)
                continue

            elif file_type == 'AllSpec2D':
                allspec2D = spec2dobj.AllSpec2DObj.from_fits(in_file, chk_version=chk_version)
                for det in allspec2D.detectors:
                    print('')
                    print('='*50 + f'{det:^7}' + '='*51)
                    wave_diag = allspec2D[det].wavesol
                    for colname in ['minWave', 'Wave_cen', 'maxWave', 'IDs_Wave_cov(%)',
                                    'measured_fwhm']:
                        wave_diag[colname].format = '0.1f'
                    for colname in ['dWave', 'RMS']:
                        wave_diag[colname].format = '0.3f'
                    wave_diag.pprint_all()
                continue
            else:
                # Should not get here unless it can't read either file type
                msgs.error("Unrecognized file type. Must be a WaveCalib or spec2d file.")

