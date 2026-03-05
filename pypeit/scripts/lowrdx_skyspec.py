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

    @classmethod
    def main(cls, args):
        from scipy.io.idl import readsav
        from astropy.io import fits

        # Initialize the log
        cls.init_log(args)

        # Read
        lrdx_sky = readsav(args.lowrdx_sky)
        wave = lrdx_sky['wave_calib']
        sky = lrdx_sky['sky_calib']

#        # Write
#        prihdu = fits.PrimaryHDU(sky)
#        prihdu.name = 'FLUX'
#        hdul = fits.HDUList([prihdu])
#        wvhdu = fits.ImageHDU(wave)
#        hdul.append(wvhdu)
#        
#        prihdu.header['NSPEC'] = 1
#        prihdu.header['NPIX'] = len(sky)
#        hdul.writeto(args.new_file, overwrite=True)

        # Write
        hdr = fits.Header()
        hdr['NSPEC'] = (1, 'Number of spectra')
        hdr['NPIX'] = (len(sky), 'Number of pixels per spectrum')
        hdu = fits.HDUList([
            fits.PrimaryHDU(data=sky, header=hdr),
            fits.ImageHDU(data=wave, name='WAVE')
        ])
        hdu[0].name = 'FLUX'
        hdu.writeto(args.new_file, overwrite=True)


