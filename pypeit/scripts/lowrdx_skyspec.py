"""
This script generates a sky spectrum from a LowRedux IDL save file

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class LowRDXSkySpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Read an IDL save file with a LowRedux sky '
                                                'spectrum and convert it into a pypeit file.',
                                    width=width)
        parser.add_argument('lowrdx_sky', type=str, default=None,
                            help = 'LowRedux Sky Spectrum (IDL save file)')
        parser.add_argument('new_file', type=str, default=None, help='PYPIT FITS sky spectrum')
        return parser

    @staticmethod
    def main(args):
        from scipy.io.idl import readsav
        from linetools.spectra.xspectrum1d import XSpectrum1D

        # Read
        lrdx_sky = readsav(args.lowrdx_sky)
        # Generate
        xspec = XSpectrum1D.from_tuple((lrdx_sky['wave_calib'], lrdx_sky['sky_calib']))
        # Write
        xspec.write_to_fits(args.new_file)


