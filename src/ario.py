from astropy.io import fits
from astropy.units import Quantity

import armsgs

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

def write_1d_spectra(slf, clobber=True):
    """ Write 1D spectra to a multi-extension FITS file

    Parameters
    ----------
    slf
    clobber : bool, optional

    Returns
    -------
    """
    # Primary header
    prihdu = fits.PrimaryHDU()
    hdus = [prihdu]

    # Loop on spectra
    ext = 0
    for kk in xrange(slf._spect['mosaic']['ndet']):
        det = kk+1
        # Loop on spectra
        for specobj in slf._specobjs[det-1]:
            ext += 1
            # Add header keyword
            keywd = 'EXT{:04d}'.format(ext)
            prihdu.header[keywd] = specobj.idx
            # Add Spectrum Table
            cols = []
            for key in specobj.boxcar.keys():
                if isinstance(specobj.boxcar[key], Quantity):
                    cols += [fits.Column(array=specobj.boxcar[key].value,
                                         name=key, format=specobj.boxcar[key].value.dtype)]
                else:
                    cols += [fits.Column(array=specobj.boxcar[key],
                                         name=key, format=specobj.boxcar[key].dtype)]
            coldefs = fits.ColDefs(cols)
            tbhdu = fits.BinTableHDU.from_columns(coldefs)
            hdus += [tbhdu]
    # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto('Science/spec1d_{:s}.fits'.format(slf._basename), clobber=clobber)

#def write_sensitivity():
    #sensfunc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.yaml".format(os.getcwd(), slf._argflag['run']['masterdir'], slf._fitsdict['target'][scidx[0]], 0, "sensfunc")
    #msgs.info("Writing sensfunc: {:s}".format(sensfunc_name))
    #with open(sensfunc_name, 'w') as yamlf:
    #    yamlf.write( yaml.dump(slf._sensfunc))
    #with io.open(sensfunc_name, 'w', encoding='utf-8') as f:
    #    f.write(unicode(json.dumps(slf._sensfunc, sort_keys=True, indent=4, separators=(',', ': '))))
