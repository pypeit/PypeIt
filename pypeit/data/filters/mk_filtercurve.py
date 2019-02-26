from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt

'''
Code for generating a fits file of filter transmission curves.
Wavelength units are Angstrom
Data are collected by Feige Wang in 2019 Feb.
The orginal transmission curves from each instruments/telescopes are in the dev suite. 
'''
pri_hdu = fits.PrimaryHDU()
hdulist = fits.HDUList(pri_hdu)

def tohdu(wave,trans,extname):
    col1 = fits.Column(name='lam',format='1D',array=wave)
    col2 = fits.Column(name='Rlam',format='1D',array=trans)
    cols = fits.ColDefs([col1,col2])
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header.update({"EXTNAME":extname})
    return hdu

## SDSS filter curves
# Download from https://www.sdss3.org/instruments/camera.php#Filters
# filter_curves.fits
par = fits.open('filter_curves_sdss.fits')
for i in ['SDSS-U','SDSS-G','SDSS-I','SDSS-G','SDSS-Z']:
    wave, trans = par[i[-1]].data['wavelength'], par[i[-1]].data['respt']
    hdu = tohdu(wave,trans,i)
    hdulist.append(hdu)

## PS1 filter curves
ps1_tbl = ascii.read('panstarrs_Tonry2012.txt')
wave = ps1_tbl['col1'] * 10.0
for i,ifilter in enumerate(['PS1-G','PS1-R','PS1-I','PS1-Z','PS1-Y']):
    trans = ps1_tbl['col{:d}'.format(i+3)]
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## DECAM filter curves
des_u = ascii.read('DECAM_ufilter.dat')
wave_u,trans_u = des_u['LAMBDA'],des_u['u']
hdu = tohdu(wave_u, trans_u, 'DECAM-U')
hdulist.append(hdu)

des_tbl = ascii.read('STD_BANDPASSES_DR1.dat')
wave = des_tbl['LAMBDA']
for i,ifilter in enumerate(['DECAM-G','DECAM-R','DECAM-I','DECAM-Z','DECAM-Y']):
    trans = des_tbl[ifilter[-1]]
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

##BASS+MZLS filter curves
# dowload from http://legacysurvey.org/dr6/description/
for i,ifilter in enumerate(['BASS-G','BASS-R']):
    tbl = ascii.read('BASS_{:s}_corr.bp'.format(ifilter[-1:]))
    wave = tbl['col1']
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)
tbl = ascii.read('kpzdccdcorr3.txt')
wave = tbl['col1']
trans = tbl['col2']
hdu = tohdu(wave, trans, 'MZLS-Z')
hdulist.append(hdu)

## HSC SSP survey filter curves
# ToDo: I don't think they include atmosphere transmission though.
# download from https://www.naoj.org/Projects/HSC/forobservers.html
for i,ifilter in enumerate(['HSC-G','HSC-R','HSC-I','HSC-Z','HSC-Y']):
    hsc_tbl = ascii.read(ifilter+'.dat')
    wave = hsc_tbl['col1']
    trans = hsc_tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## LSST filter curves
# download from https://raw.githubusercontent.com/lsst-pst/syseng_throughputs/master/components/camera/filters/y_band_Response.dat
for i,ifilter in enumerate(['LSST-U','LSST-G','LSST-R','LSST-I','LSST-Z','LSST-Y']):
    tbl = ascii.read('LSST_{:s}_band_Response.dat'.format(ifilter[-1:]))
    wave = tbl['col1']*10.
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## CFHT filter curves
# download from http://www.cfht.hawaii.edu/Instruments/Filters/megaprimenew.html
# u = cfh9302, g = cfh9402, r = cfh9602, i = cfh9703, z = cfh9901
for i,ifilter in enumerate(['CFHT-U','CFHT-G','CFHT-R','CFHT-I','CFHT-Z']):
    tbl = ascii.read('cfh_{:s}.dat'.format(ifilter[-1:]))
    wave = tbl['col1']*10.
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## UKIRT filter curves
for i,ifilter in enumerate(['UKIRT-Z','UKIRT-Y','UKIRT-J','UKIRT-H','UKIRT-K']):
    tbl = ascii.read('ukidss{:s}.txt'.format(ifilter[-1:]))
    wave = tbl['col1']*1e4
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## VISTA filter curves
#download from http://www.eso.org/sci/facilities/paranal/instruments/vircam/inst/Filters_QE_Atm_curves.tar.gz
for i,ifilter in enumerate(['VISTA-Z','VISTA-Y','VISTA-J','VISTA-H','VISTA-K']):
    tbl = ascii.read('VISTA_Filters_at80K_forETC_{:s}.txt'.format(ifilter[-1:]))
    wave = tbl['col1']*10.
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

##TWO MASS filter curves
#download from https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html#rsr
for i,ifilter in enumerate(['TMASS-J','TMASS-H','TMASS-K']):
    tbl = ascii.read('2mass{:s}.txt'.format(ifilter[-1:]))
    wave = tbl['col1']*1e4
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## WISE filter curves
for i,ifilter in enumerate(['WISE-W1','WISE-W2','WISE-W3','WISE-W4']):
    tbl = ascii.read('RSR-{:s}.txt'.format(ifilter[-2:]))
    wave = tbl['col1']*1e4
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## GAIA DR2 revised filter curves
#download from https://www.cosmos.esa.int/web/gaia/iow_20180316/
gaia_tbl = ascii.read('GaiaDR2_RevisedPassbands.dat')
wave = gaia_tbl['wave'] * 10.0
for i,ifilter in enumerate(['GAIA-G','GAIA-B','GAIA-R']):
    trans = gaia_tbl[ifilter[-1]]
    gdp = trans<90. # bad values in the file was set to 99.0
    hdu = tohdu(wave[gdp],trans[gdp],ifilter)
    hdulist.append(hdu)

## GALEX filter curves
# download from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=GALEX
for i,ifilter in enumerate(['GALEX-N','GALEX-F']):
    tbl = ascii.read('GALEX_GALEX.{:s}UV.dat'.format(ifilter[-1]))
    wave = tbl['col1']
    trans = tbl['col2']
    hdu = tohdu(wave,trans,ifilter)
    hdulist.append(hdu)

## Write out the final filter curve table
hdulist.writeto('filtercurves_20190215.fits',overwrite=True)
