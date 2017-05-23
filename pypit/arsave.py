""" Output for PYPIT
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)
import numpy as np
import os

import astropy.io.fits as pyfits
from astropy.units import Quantity
from astropy import units as u
from astropy.table import Table

import h5py

from pypit import armsgs
from pypit import arutils
from pypit import arparse as settings

from pypit import ardebug as debugger

try: input = raw_input
except NameError: pass

# Logging
msgs = armsgs.get_logger()


def save_arcids(fname, pixels):
    # Setup the HDU
    hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([hdu]) # Insert the primary HDU (input model)
    for o in range(len(pixels)):
        hdulist.append(pyfits.ImageHDU(pixels[o])) # Add a new Image HDU
    ans = 'y'
    if os.path.exists(fname):
        if settings.argflag['output']['overwrite']:
            os.remove(fname)
        else:
            ans = ''
            while ans != 'y' and ans != 'n' and ans != 'r':
                msgs.warn("File %s exists!" % (fname), verbose=settings.argflag['output']['verbosity'])
                ans = raw_input(msgs.input()+"Overwrite? (y/n)")
            if ans == 'y': os.remove(fname)
    if ans == 'y':
        msgs.info("Arc IDs saved successfully to file:"+msgs.newline()+fname)
        hdulist.writeto(fname)
    return


def save_extraction(slf, sciext, scidx, scierr=None, filename="temp.fits", frametype='Extraction', wave=None, sky=None, skyerr=None, extprops=None):
    msgs.info("Saving {0:s} frame as:".format(frametype)+msgs.newline()+filename)
    # Setup the HDU
    if sky is None:
        savdat = sciext
    else:
        if sciext.ndim != 2 or sciext.shape != sky.shape:
            msgs.error("Could not save extraction"+msgs.newline() +
                       "science and sky frames have different dimensions or shape")
        tsavdat = sciext[:,:,np.newaxis]
        tsavdat = np.append(tsavdat,sky[:,:,np.newaxis],axis=2)
        if scierr is not None:
            if skyerr is None:
                msgs.error("An error frame is missing for the sky")
            savdat = tsavdat[:,:,:,np.newaxis]
            tsaverr = scierr[:,:,np.newaxis]
            tsaverr = np.append(tsaverr,skyerr[:,:,np.newaxis],axis=2)
            savdat = np.append(savdat,tsaverr[:,:,:,np.newaxis],axis=3)
        else:
            savdat = tsavdat

    hdu = pyfits.PrimaryHDU(savdat)
    hdulist = pyfits.HDUList([hdu])
    # Write some information to the header
    msgs.info("Writing header information")
    hdrname = "FRAMEEXT"
    hdulist[0].header[hdrname] = (slf._fitsdict['filename'][scidx[0]], 'ARMED: Name of file that was extracted'.format(frametype))
    hdulist[0].header["FRAMETYP"] = (frametype, 'ARMED: extraction frame')
    hdulist[0].header["NUMORDS"] = (sciext.shape[1], 'ARMED: Number of orders extracted')
    hdulist[0].header["PIXSIZE"] = (settings.argflag['reduce']['pixelsize'], 'ARMED: The size of each sampled pixel (km/s)')
    # Loop through all orders and write the wavelength into the header
    if wave is not None:
        for i in range(sciext.shape[1]):
            hdrname = "CDELT{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (np.log10(1.0 + settings.argflag['reduce']['pixelsize']/299792.458), 'ARMED: log10(1+pixsize/c)'.format(frametype))
            hdrname = "CRVAL{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (np.log10(wave[0,i]), 'ARMED: log10(lambda_0)'.format(frametype))
            hdrname = "CLINV{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (wave[0,i], 'ARMED: lambda_0'.format(frametype))
            hdrname = "CRPIX{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (0.0, 'ARMED: Offset=0.0'.format(frametype))
            hdrname = "CNPIX{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (np.size(np.where(wave[:,i]!=-999999.9)[0]), 'ARMED: Offset=0.0'.format(frametype))
    if extprops is not None:
        kys = extprops.keys()
        for j in range(len(kys)):
            hkey = kys[j][:5].upper()
            if np.ndim(extprops[kys[j]])==1 and np.size(extprops[kys[j]]==sciext.shape[1]):
                for i in range(sciext.shape[1]):
                    hdrname = "{0:s}{1:03d}".format(hkey,i+1)
                    hdulist[0].header[hdrname] = (extprops[kys[j]][i], 'ARMED: {0:s} for order {1:d}'.format(kys[j],i+1))
    # Write the file to disk
    if os.path.exists(filename):
        if settings.argflag['output']['overwrite'] == True:
            msgs.warn("Overwriting file:"+msgs.newline()+filename)
            os.remove(filename)
            hdulist.writeto(filename)
            msgs.info("{0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
        else:
            msgs.warn("This file already exists")
            rmfil=''
            while rmfil != 'n' and rmfil != 'y' and rmfil != 'a':
                rmfil=raw_input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
            if rmfil == 'n':
                msgs.warn("Not saving {0:s} frame:".format(frametype)+msgs.newline()+filename)
            else:
                os.remove(filename)
                if rmfil == 'a': settings.argflag['output']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("{0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    else:
        hdulist.writeto(filename)
        msgs.info("{0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    return


def save_master(slf, data, filename="temp.fits", frametype="<None>", ind=[],
                extensions=None, keywds=None, names=None):
    """ Write a MasterFrame
    Parameters
    ----------
    slf
    data
    filename
    frametype
    ind
    extensions : list, optional
      Additional data images to write
    names : list, optional
      Names of the extensions
    keywds : Additional keywords for the Header
    Returns
    -------
    """
    msgs.info("Saving master {0:s} frame as:".format(frametype)+msgs.newline()+filename)
    hdu = pyfits.PrimaryHDU(data)
    hlist = [hdu]
    # Extensions
    if extensions is not None:
        for kk,exten in enumerate(extensions):
            hdu = pyfits.ImageHDU(exten)
            if names is not None:
                hdu.name = names[kk]
            hlist.append(hdu)
    # HDU list
    hdulist = pyfits.HDUList(hlist)
    # Header
    msgs.info("Writing header information")
    for i in range(len(ind)):
        hdrname = "FRAME{0:03d}".format(i+1)
        hdulist[0].header[hdrname] = (slf._fitsdict['filename'][ind[i]], 'PYPIT: File used to generate Master {0:s}'.format(frametype))
    hdulist[0].header["FRAMETYP"] = (frametype, 'PYPIT: Master calibration frame type')
    if keywds is not None:
        for key in keywds.keys():
            hdulist[0].header[key] = keywds[key]
    # Write the file to disk
    if os.path.exists(filename):
        if settings.argflag['output']['overwrite'] is True:
            msgs.warn("Overwriting file:"+msgs.newline()+filename)
            os.remove(filename)
            hdulist.writeto(filename)
            msgs.info("Master {0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
        else:
            msgs.warn("This file already exists")
            rmfil = ''
            while rmfil != 'n' and rmfil != 'y' and rmfil != 'a':
                rmfil = raw_input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
            if rmfil == 'n':
                msgs.warn("Not saving master {0:s} frame:".format(frametype)+msgs.newline()+filename)
            else:
                os.remove(filename)
                if rmfil == 'a': settings.argflag['output']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("Master {0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    else:
        hdulist.writeto(filename)
        msgs.info("Master {0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    return


def save_ordloc(slf, fname):
    # Derive a suitable name
    mstrace_bname, mstrace_bext = os.path.splitext(fname)
    # Save the left order locations
    hdu = pyfits.PrimaryHDU(slf._lordloc)
    hdulist = pyfits.HDUList([hdu])
    # Write the file to disk
    filename = mstrace_bname+"_ltrace"+mstrace_bext
    if os.path.exists(filename):
        if settings.argflag['output']['overwrite'] is True:
            msgs.warn("Overwriting file:"+msgs.newline()+filename)
            os.remove(filename)
            hdulist.writeto(filename)
            msgs.info("Saved left order locations for frame:"+msgs.newline()+fname)
        else:
            msgs.warn("This file already exists:"+msgs.newline()+filename)
            rmfil=''
            while rmfil != 'n' and rmfil != 'y' and rmfil != 'a':
                rmfil=raw_input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
            if rmfil == 'n':
                msgs.warn("Not saving left order traces for file:"+msgs.newline()+fname)
            else:
                os.remove(filename)
                if rmfil == 'a': settings.argflag['output']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("Saved left order locations for frame:"+msgs.newline()+fname)
    else:
        hdulist.writeto(filename)
        msgs.info("Saved left order locations for frame:"+msgs.newline()+fname)
    # Save the right order locations
    hdu = pyfits.PrimaryHDU(slf._rordloc)
    hdulist = pyfits.HDUList([hdu])
    filename = mstrace_bname+"_rtrace"+mstrace_bext
    if os.path.exists(filename):
        if settings.argflag['output']['overwrite'] is True:
            msgs.warn("Overwriting file:"+msgs.newline()+filename)
            os.remove(filename)
            hdulist.writeto(filename)
            msgs.info("Saved right order locations for frame:"+msgs.newline()+fname)
        else:
            msgs.warn("This file already exists:"+msgs.newline()+filename)
            rmfil=''
            while rmfil != 'n' and rmfil != 'y' and rmfil != 'a':
                rmfil=raw_input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
            if rmfil == 'n':
                msgs.warn("Not saving right order traces for file:"+msgs.newline()+fname)
            else:
                os.remove(filename)
                if rmfil == 'a': settings.argflag['output']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("Saved right order locations for frame:"+msgs.newline()+fname)
    else:
        hdulist.writeto(filename)
        msgs.info("Saved right order locations for frame:"+msgs.newline()+fname)
    return


def save_tilts(slf, fname):
    # Derive a suitable name
    msarc_bname, msarc_bext = os.path.splitext(fname)
    # Save the tilts
    hdu = pyfits.PrimaryHDU(slf._tilts)
    hdulist = pyfits.HDUList([hdu])
    # Write the file to disk
    filename = msarc_bname+"_tilts"+msarc_bext
    if os.path.exists(filename):
        if settings.argflag['output']['overwrite'] == True:
            msgs.warn("Overwriting file:"+msgs.newline()+filename)
            os.remove(filename)
            hdulist.writeto(filename)
            msgs.info("Saved order tilts for frame:"+msgs.newline()+fname)
        else:
            msgs.warn("This file already exists:"+msgs.newline()+filename)
            rmfil=''
            while rmfil != 'n' and rmfil != 'y' and rmfil != 'a':
                rmfil=raw_input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
            if rmfil == 'n':
                msgs.warn("Not saving order tilts for file:"+msgs.newline()+fname)
            else:
                os.remove(filename)
                if rmfil == 'a': settings.argflag['output']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("Saved order tilts for frame:"+msgs.newline()+fname)
    else:
        hdulist.writeto(filename)
        msgs.info("Saved order tilts for frame:"+msgs.newline()+fname)
    # Save the saturation mask
    hdu = pyfits.PrimaryHDU(slf._satmask)
    hdulist = pyfits.HDUList([hdu])
    filename = msarc_bname+"_satmask"+msarc_bext
    if os.path.exists(filename):
        if settings.argflag['output']['overwrite'] == True:
            msgs.warn("Overwriting file:"+msgs.newline()+filename)
            os.remove(filename)
            hdulist.writeto(filename)
            msgs.info("Saved saturation mask for frame:"+msgs.newline()+fname)
        else:
            msgs.warn("This file already exists:"+msgs.newline()+filename)
            rmfil=''
            while rmfil != 'n' and rmfil != 'y' and rmfil != 'a':
                rmfil=raw_input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
            if rmfil == 'n':
                msgs.warn("Not saving saturation mask for file:"+msgs.newline()+fname)
            else:
                os.remove(filename)
                if rmfil == 'a': settings.argflag['output']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("Saved saturation mask for frame:"+msgs.newline()+fname)
    else:
        hdulist.writeto(filename)
        msgs.info("Saved saturation mask for frame:"+msgs.newline()+fname)

    return


def save_1d_spectra_hdf5(slf, fitsdict, clobber=True):
    """ Write 1D spectra to an HDF5 file

    Parameters
    ----------
    slf
    clobber

    Returns
    -------

    """
    from linetools import utils as ltu
    import json
    if clobber is False:
        msgs.error("NOT IMPLEMENTED")
    # Open file
    outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.hdf5'.format(slf._basename)
    hdf = h5py.File(outfile,'w')

    # Meta Table
    idict = dict(RA=0., DEC=0.,  # J2000
                 objid=0, slitid=0, det=0, scidx=0,  # specobj IDs
                 FWHM=0.,  # Spatial resolution in arcsec
                 R=0.,     # Spectral resolution (FWHM) in lambda/Dlambda
                 xslit=(0.,0.), nslit=0)
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    meta = Table(lst, names=tkeys)

    # Calculate number of objects and totalpix
    nspec, totpix = 0, 0
    for kk in range(settings.spect['mosaic']['ndet']):
        det = kk+1
        # Loop on slits
        for sl in range(len(slf._specobjs[det-1])):
            nspec += len(slf._specobjs[det-1][sl])
            # Loop on objects
            for specobj in slf._specobjs[det-1][sl]:
                # Calculate max pixels
                totpix = max(totpix, specobj.trace.size)
                # Update meta
                tdict = dict(RA=0., DEC=0.,  # J2000
                             objid=specobj.objid, slitid=specobj.slitid, det=det, scidx=specobj.scidx,  # specobj IDs
                             FWHM=0.,  # Spatial resolution in arcsec
                             R=0.,     # Spectral resolution (FWHM) in lambda/Dlambda
                             xslit=specobj.xslit, nslit=sl+1)  # Slit position and number
                meta.add_row(tdict)
    # Remove dummy row and write
    meta = meta[1:]
    hdf['meta'] = meta

    # Make a Header from fitsdict
    hdict = {}
    for key in fitsdict.keys():
        hdict[key] = fitsdict[key][slf._specobjs[0][0][0].scidx]  # Hopefully this is the right index
    d = ltu.jsonify(hdict)
    hdf['header'] = json.dumps(d)

    # Loop on extraction methods
    for ex_method in ['boxcar', 'optimal']:
        # Check for extraction type
        if not hasattr(slf._specobjs[0][0][0], ex_method):
            continue
        method_grp = hdf.create_group(ex_method)

        # Data arrays are always MaskedArray
        dtypes = []
        for key in getattr(slf._specobjs[0][0][0], ex_method).keys():
            dtype = 'float64' if key == 'wave' else 'float32'
            dtypes.append((str(key), dtype, (totpix)))
        dtypes.append((str('obj_trace'), 'float32', (totpix)))
        data = np.ma.empty((1,), dtype=dtypes)
        # Setup in hdf5
        spec_set = hdf[str(ex_method)].create_dataset('spec', data=data, chunks=True,
                                                      maxshape=(None,), compression='gzip')
        spec_set.resize((nspec,))
        # Fill (and make meta)
        count = 0
        for kk in range(settings.spect['mosaic']['ndet']):
            det = kk+1
            # Loop on slits
            for sl in range(len(slf._specobjs[det - 1])):
                nspec += len(slf._specobjs[det - 1][sl])
                # Loop on spectra
                for specobj in slf._specobjs[det-1][sl]:
                    # Check meta
                    assert meta['objid'][count] == specobj.objid
                    # Trace
                    data['obj_trace'][0][:len(specobj.trace)] = specobj.trace
                    # Rest
                    sdict = getattr(specobj, ex_method)
                    for key in sdict.keys():
                        npix = len(sdict[key])
                        try:
                            data[key][0][:npix] = sdict[key].value
                        except AttributeError:
                            data[key][0][:npix] = sdict[key]
                    # Write
                    spec_set[count] = data
                    count += 1
    #
    hdf.close()

    # Dump into a linetools.spectra.xspectrum1d.XSpectrum1D


def save_1d_spectra_fits(slf, clobber=True):
    """ Write 1D spectra to a multi-extension FITS file

    Parameters
    ----------
    slf
    clobber : bool, optional

    Returns
    -------
    """
    # Primary header
    prihdu = pyfits.PrimaryHDU()
    hdus = [prihdu]

    # Loop on detectors
    ext = 0
    for kk in range(settings.spect['mosaic']['ndet']):
        det = kk+1
        # Loop on slits
        for sl in range(len(slf._specobjs[det-1])):
            # Loop on spectra
            for specobj in slf._specobjs[det-1][sl]:
                ext += 1
                # Add header keyword
                keywd = 'EXT{:04d}'.format(ext)
                prihdu.header[keywd] = specobj.idx
                # Add Spectrum Table
                cols = []
                # Trace
                cols += [pyfits.Column(array=specobj.trace, name=str('obj_trace'), format=specobj.trace.dtype)]
                # Boxcar
                for key in specobj.boxcar.keys():
                    # Skip some
                    if key in ['size']:
                        continue
                    if isinstance(specobj.boxcar[key], Quantity):
                        cols += [pyfits.Column(array=specobj.boxcar[key].value,
                                             name=str('box_'+key), format=specobj.boxcar[key].value.dtype)]
                    else:
                        cols += [pyfits.Column(array=specobj.boxcar[key],
                                             name=str('box_'+key), format=specobj.boxcar[key].dtype)]
                # Optimal
                for key in specobj.optimal.keys():
                    # Skip some
                    if key in ['fwhm']:
                        continue
                    # Generate column
                    if isinstance(specobj.optimal[key], Quantity):
                        cols += [pyfits.Column(array=specobj.optimal[key].value,
                                               name=str('opt_'+key), format=specobj.optimal[key].value.dtype)]
                    else:
                        cols += [pyfits.Column(array=specobj.optimal[key],
                                               name=str('opt_'+key), format=specobj.optimal[key].dtype)]
                # Finish
                coldefs = pyfits.ColDefs(cols)
                tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
                tbhdu.name = specobj.idx
                hdus += [tbhdu]
    # Finish
    hdulist = pyfits.HDUList(hdus)
    hdulist.writeto(settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(slf._basename), overwrite=clobber)


#def write_sensitivity():
    #sensfunc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.yaml".format(os.getcwd(), settings.argflag['run']['directory']['master'], slf._fitsdict['target'][scidx[0]], 0, "sensfunc")
    #msgs.info("Writing sensfunc: {:s}".format(sensfunc_name))
    #with open(sensfunc_name, 'w') as yamlf:
    #    yamlf.write( yaml.dump(slf._sensfunc))
    #with io.open(sensfunc_name, 'w', encoding='utf-8') as f:
    #    f.write(unicode(json.dumps(slf._sensfunc, sort_keys=True, indent=4, separators=(',', ': '))))


def save_obj_info(slf, fitsdict, clobber=True):
    # Lists for a Table
    slits, names, boxsize, opt_fwhm, s2n = [], [], [], [], []
    # Loop on detectors
    for kk in range(settings.spect['mosaic']['ndet']):
        det = kk+1
        dnum = settings.get_dnum(det)
        # Loop on slits
        for sl in range(len(slf._specobjs[det-1])):
            # Loop on spectra
            for specobj in slf._specobjs[det-1][sl]:
                # Append
                names.append(specobj.idx)
                slits.append(sl+1)
                # Boxcar width
                if 'size' in specobj.boxcar.keys():
                    slit_pix = specobj.boxcar['size']
                    # Convert to arcsec
                    binspatial, binspectral = settings.parse_binning(fitsdict['binning'][specobj.scidx])
                    boxsize.append(slit_pix*binspatial*settings.spect[dnum]['platescale'])
                else:
                    boxsize.append(0.)
                # Optimal profile (FWHM)
                if 'fwhm' in specobj.optimal.keys():
                    binspatial, binspectral = settings.parse_binning(fitsdict['binning'][specobj.scidx])
                    opt_fwhm.append(specobj.optimal['fwhm']*binspatial*settings.spect[dnum]['platescale'])
                else:
                    opt_fwhm.append(0.)
                # S2N -- default to boxcar
                sext = (specobj.boxcar if (len(specobj.boxcar) > 0) else specobj.optimal)
                ivar = arutils.calc_ivar(sext['var'])
                is2n = np.median(sext['counts']*np.sqrt(ivar))
                s2n.append(is2n)

    # Generate the table, if we have at least one source
    if len(names) > 0:
        obj_tbl = Table()
        obj_tbl['slit'] = slits
        obj_tbl['slit'].format = 'd'
        obj_tbl['name'] = names
        obj_tbl['box_width'] = boxsize
        obj_tbl['box_width'].format = '.2f'
        obj_tbl['box_width'].unit = u.arcsec
        obj_tbl['opt_fwhm'] = opt_fwhm
        obj_tbl['opt_fwhm'].format = '.3f'
        obj_tbl['opt_fwhm'].unit = u.arcsec
        obj_tbl['s2n'] = s2n
        obj_tbl['s2n'].format = '.2f'
        # Write
        obj_tbl.write(settings.argflag['run']['directory']['science']+'/objinfo_{:s}.txt'.format(slf._basename), format='ascii.fixed_width', overwrite=True)


def save_2d_images(slf, fitsdict, clobber=True):
    """ Write 2D images to the hard drive
    Parameters
    ----------
    slf
    clobber

    Returns
    -------

    """
    import datetime
    # Original header
    scidx = slf._idx_sci[0]
    path = fitsdict['directory'][scidx]
    ifile = fitsdict['filename'][scidx]
    headarr = [None for k in range(settings.spect['fits']['numhead'])]
    for k in range(settings.spect['fits']['numhead']):
        headarr[k] = pyfits.getheader(path+ifile, ext=settings.spect['fits']['headext{0:02d}'.format(k+1)])

    # Primary header
    prihdu = pyfits.PrimaryHDU()
    # Update with original header, skipping a few keywords
    hdus = [prihdu]
    hdukeys = ['BUNIT', 'COMMENT', '', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
               'HISTORY', 'EXTEND', 'DATASEC']
    for key in headarr[0].keys():
        # Use new ones
        if key in hdukeys:
            continue
        # Update unused ones
        prihdu.header[key] = headarr[0][key]
    # History
    if 'HISTORY' in headarr[0].keys():
        # Strip \n
        tmp = str(headarr[0]['HISTORY']).replace('\n', ' ')
        prihdu.header.add_history(str(tmp))

    # PYPIT
    prihdu.header['PIPELINE'] = str('PYPIT')
    prihdu.header['DATE-RDX'] = str(datetime.date.today().strftime('%Y-%b-%d'))
    setup = settings.argflag['reduce']['masters']['setup'].split('_')
    prihdu.header['PYPCNFIG'] = str(setup[0])
    prihdu.header['PYPCALIB'] = str(setup[2])
    prihdu.header['PYPMFDIR'] = str(settings.argflag['run']['directory']['master']+'_'+settings.argflag['run']['spectrograph'])

    ext = 0
    for kk in range(settings.spect['mosaic']['ndet']):
        det = kk+1

        # Processed frame
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = 'DET{:d}-Processed'.format(det)
        hdu = pyfits.ImageHDU(slf._sciframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

        # Variance
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = 'DET{:d}-Var'.format(det)
        hdu = pyfits.ImageHDU(slf._modelvarframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

        # Background subtracted
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = 'DET{:d}-Skysub'.format(det)
        hdu = pyfits.ImageHDU(slf._sciframe[det-1]-slf._bgframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

    # Finish
    hdulist = pyfits.HDUList(hdus)
    hdulist.writeto(settings.argflag['run']['directory']['science']+'/spec2d_{:s}.fits'.format(slf._basename), overwrite=clobber)
