""" Output for PYPEIT
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import datetime

import numpy as np

from astropy import units
from astropy.io import fits
from astropy.table import Table


from pypeit import msgs
from pypeit import specobjs
from pypeit.core import parse
from pypeit import debugger


'''
def save_arcids(fname, pixels):
    # Setup the HDU
    hdu = fits.PrimaryHDU()
    hdulist = fits.HDUList([hdu]) # Insert the primary HDU (input model)
    for o in range(len(pixels)):
        hdulist.append(fits.ImageHDU(pixels[o])) # Add a new Image HDU
    ans = 'y'
    if os.path.exists(fname):
        if settings.argflag['output']['overwrite']:
            os.remove(fname)
        else:
            ans = ''
            while ans != 'y' and ans != 'n' and ans != 'r':
                msgs.warn("File %s exists!" % (fname), verbose=settings.argflag['output']['verbosity'])
                ans = input(msgs.input()+"Overwrite? (y/n)")
            if ans == 'y': os.remove(fname)
    if ans == 'y':
        msgs.info("Arc IDs saved successfully to file:"+msgs.newline()+fname)
        hdulist.writeto(fname)
    return
'''

'''
def save_extraction(slf, sciext, scidx, scierr=None, filename="temp.fits", frametype='Extraction', wave=None, sky=None, skyerr=None, extprops=None):
    msgs.info("Saving {0:s} frame as:".format(frametype)+msgs.newline()+filename)
    # Setup the HDU
    if sky is None:
        savdat = sciext
    else:
        if sciext.ndim != 2 or sciext.shape != sky.shape:
            msgs.error("Could not save extraction"+msgs.newline() +
                       "science and sky frames have different dimensions or shape")
        tsavdat = sciext[:, :, np.newaxis]
        tsavdat = np.append(tsavdat, sky[:, :, np.newaxis], axis=2)
        if scierr is not None:
            if skyerr is None:
                msgs.error("An error frame is missing for the sky")
            savdat = tsavdat[:, :, :, np.newaxis]
            tsaverr = scierr[:, :, np.newaxis]
            tsaverr = np.append(tsaverr, skyerr[:, :, np.newaxis], axis=2)
            savdat = np.append(savdat, tsaverr[:, :, :, np.newaxis], axis=3)
        else:
            savdat = tsavdat

    hdu = fits.PrimaryHDU(savdat)
    hdulist = fits.HDUList([hdu])
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
            hdulist[0].header[hdrname] = (np.log10(wave[0, i]), 'ARMED: log10(lambda_0)'.format(frametype))
            hdrname = "CLINV{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (wave[0, i], 'ARMED: lambda_0'.format(frametype))
            hdrname = "CRPIX{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (0.0, 'ARMED: Offset=0.0'.format(frametype))
            hdrname = "CNPIX{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (np.size(np.where(wave[:, i] != -999999.9)[0]), 'ARMED: Offset=0.0'.format(frametype))
    if extprops is not None:
        kys = extprops.keys()
        for j in range(len(kys)):
            hkey = kys[j][:5].upper()
            if np.ndim(extprops[kys[j]]) == 1 and np.size(extprops[kys[j]] == sciext.shape[1]):
                for i in range(sciext.shape[1]):
                    hdrname = "{0:s}{1:03d}".format(hkey, i+1)
                    hdulist[0].header[hdrname] = (extprops[kys[j]][i], 'ARMED: {0:s} for order {1:d}'.format(kys[j], i+1))
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
                rmfil=input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
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
'''



'''
def save_ordloc(slf, fname):
    # Derive a suitable name
    mstrace_bname, mstrace_bext = os.path.splitext(fname)
    # Save the left order locations
    hdu = fits.PrimaryHDU(slf._lordloc)
    hdulist = fits.HDUList([hdu])
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
                rmfil=input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
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
    hdu = fits.PrimaryHDU(slf._rordloc)
    hdulist = fits.HDUList([hdu])
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
                rmfil=input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
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
    hdu = fits.PrimaryHDU(slf._tilts)
    hdulist = fits.HDUList([hdu])
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
                rmfil=input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
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
    hdu = fits.PrimaryHDU(slf._satmask)
    hdulist = fits.HDUList([hdu])
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
                rmfil=input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
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
'''


'''
def save_1d_spectra_hdf5(slf, fitsdict, clobber=True):
    """ Write 1D spectra to an HDF5 file

    Parameters
    ----------
    slf
    clobber

    Returns
    -------

    """
    debugger.set_trace()  # NEEDS REFACTORING
    if clobber is False:
        msgs.error("NOT IMPLEMENTED")
    # Open file
    outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.hdf5'.format(slf._basename)
    hdf = h5py.File(outfile, 'w')

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
    detref = None
    for kk in range(settings.spect['mosaic']['ndet']):
        det = kk+1
        if slf._specobjs[det-1] is None:
            continue
        if detref is None:
            detref = det-1
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
        hdict[key] = fitsdict[key][slf._specobjs[detref][0][0].scidx]  # Hopefully this is the right index
    d = linetools.utils.jsonify(hdict)
    hdf['header'] = json.dumps(d)

    # Loop on extraction methods
    for ex_method in ['boxcar', 'optimal']:
        # Check for extraction type
        if not hasattr(slf._specobjs[detref][0][0], ex_method):
            continue
        method_grp = hdf.create_group(ex_method)

        # Data arrays are always MaskedArray
        dtypes = []
        for key in getattr(slf._specobjs[detref][0][0], ex_method).keys():
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
            if slf._specobjs[det - 1] is None:
                continue
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
'''

def save_1d_spectra_fits(specObjs, header, pypeline, instrume,
                         outfile, helio_dict=None, telescope=None, overwrite=True,
                         update_det=None):
    """ Write 1D spectra to a multi-extension FITS file

    Parameters
    ----------
    specobjs : SpecObjs object
    header : dict or Row (dict-like)
    pypeline : str
      Name of PypeIt pipeline (e.g. 'MultiSlit')
    instrume : str
      Name of instrument
    outfile : str
    overwrite : bool, optional
    update_det : int or list, optional
      If provided, do not clobber the existing file but only update
      the indicated detectors.  Useful for re-running on a subset of detectors

    Returns
    -------
    outfile : str
    """
    hdus, prihdu = init_hdus(update_det, outfile)
    sobjs_key = specobjs.SpecObj.sobjs_key()
    # Init for spec1d as need be
    if hdus is None:
        prihdu = fits.PrimaryHDU()
        hdus = [prihdu]
        # Add critical data to header
        for key in ['ra', 'dec', 'exptime', 'target', 'airmass', 'filename']:
            # Allow for fitstbl vs. header
            try:
                prihdu.header[key.upper()] = header[key.upper()]
            except KeyError:
                prihdu.header[key.upper()] = header[key]
        try:
            prihdu.header['MJD-OBS'] = header['mjd']  # recorded as 'mjd' in fitstbl
        except KeyError:
            prihdu.header['MJD-OBS'] = header['MJD-OBS']
        prihdu.header['INSTRUME'] = instrume

        # Specify which pipeline created this file
        prihdu.header['PYPELINE'] = pypeline

        # Observatory
        if telescope is not None:
            prihdu.header['LON-OBS'] = telescope['longitude']
            prihdu.header['LAT-OBS'] = telescope['latitude']
            prihdu.header['ALT-OBS'] = telescope['elevation']
        # Helio
        if helio_dict is not None:
            prihdu.header['VEL-TYPE'] = helio_dict['refframe'] # settings.argflag['reduce']['calibrate']['refframe']
            prihdu.header['VEL'] = helio_dict['vel_correction'] # slf.vel_correction

    npix = 0
    ext = len(hdus)-1
    # Loop on specobjs
    for sobj in specObjs.specobjs:
        if sobj is None:
            continue
        ext += 1
        # Add header keyword
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = sobj.idx

        # Add Spectrum Table
        cols = []
        # Trace
        if sobj.trace_spat is not None:
            cols += [fits.Column(array=sobj.trace_spat, name=str('TRACE'), format=sobj.trace_spat.dtype)]
        # FWHM fit from extraction
        if sobj.fwhmfit is not None:
            cols += [fits.Column(array=sobj.fwhmfit, name=str('FWHM'), format=sobj.fwhmfit.dtype)]
        if ext == 1:
            # TODO -- FIX THIS KLUDGE
            try:
                npix = len(sobj['trace'])
            except:  # THIS IS A DUMB KLUDGE
                npix = len(sobj['trace_spat'])
        # Boxcar
        for key in sobj.boxcar.keys():
            # Skip some
            if key in ['BOX_RADIUS']:
                continue
            if isinstance(sobj.boxcar[key], units.Quantity):
                cols += [fits.Column(array=sobj.boxcar[key].value,
                                     name=str('BOX_'+key), format=sobj.boxcar[key].value.dtype)]
            else:
                cols += [fits.Column(array=sobj.boxcar[key],
                                     name=str('BOX_'+key), format=sobj.boxcar[key].dtype)]
        # Optimal
        for key in sobj.optimal.keys():
            # Skip some
            #if key in ['fwhm']:
            #    continue
            # Generate column
            if isinstance(sobj.optimal[key], units.Quantity):
                cols += [fits.Column(array=sobj.optimal[key].value,
                                       name=str('OPT_'+key), format=sobj.optimal[key].value.dtype)]
            else:
                cols += [fits.Column(array=sobj.optimal[key],
                                       name=str('OPT_'+key), format=sobj.optimal[key].dtype)]
        # Finish
        coldefs = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = sobj.idx
        for attr, hdrcard in sobjs_key.items():
            tbhdu.header[hdrcard] = getattr(sobj,attr)
        hdus += [tbhdu]

    # A few more for the header
    prihdu.header['NSPEC'] = len(hdus) - 1
    prihdu.header['NPIX'] = npix
    # If this is echelle write the objid and the orderindx to the header as well


    # Finish
    hdulist = fits.HDUList(hdus)
    #if outfile is None:
    #    outfile = settings.argflag['run']['directory']['science']+'/spec1d_{:s}.fits'.format(slf._basename)
    hdulist.writeto(outfile, overwrite=overwrite)
    msgs.info("Wrote 1D spectra to {:s}".format(outfile))
    return outfile


#def write_sensitivity():
    #sensfunc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.yaml".format(os.getcwd(), settings.argflag['run']['directory']['master'], slf._fitsdict['target'][scidx[0]], 0, "sensfunc")
    #msgs.info("Writing sensfunc: {:s}".format(sensfunc_name))
    #with open(sensfunc_name, 'w') as yamlf:
    #    yamlf.write( yaml.dump(slf._sensfunc))
    #with io.open(sensfunc_name, 'w', encoding='utf-8') as f:
    #    f.write(unicode(json.dumps(slf._sensfunc, sort_keys=True, indent=4, separators=(',', ': '))))

# TODO: (KBW) I don't think core algorithms should take class
# arguments...
def save_obj_info(all_specobjs, spectrograph, basename, science_dir, binning=None):
    """

    Parameters
    ----------
    all_specobjs : list
    fitstbl : Table

    Returns
    -------

    """
    slits, names, spat_pixpos, spat_fracpos, boxsize, opt_fwhm, s2n = [], [], [], [], [], [], []  # Lists for a Table
    binspatial, binspectral = parse.parse_binning(binning)
    for specobj in all_specobjs:
        if specobj is None:
            continue
        # Append
        names.append(specobj.idx)
        slits.append(specobj.slitid)
        spat_pixpos.append(specobj.spat_pixpos)
        if spectrograph.pypeline == 'MultiSlit':
            spat_fracpos.append(specobj.spat_fracpos)
        elif spectrograph.pypeline == 'Echelle':
            spat_fracpos.append(specobj.ech_fracpos)
        # Boxcar width
        if 'BOX_RADIUS' in specobj.boxcar.keys():
            slit_pix = 2.0*specobj.boxcar['BOX_RADIUS']
            # Convert to arcsec
            binspatial, binspectral = parse.parse_binning(binning)
            boxsize.append(slit_pix*binspatial*spectrograph.detector[specobj.det-1]['platescale'])
        else:
            boxsize.append(0.)

        # Optimal profile (FWHM)
        opt_fwhm.append(np.median(specobj.fwhmfit)* binspatial*spectrograph.detector[specobj.det-1]['platescale'])
        # S2N -- default to boxcar
        #sext = (specobj.boxcar if (len(specobj.boxcar) > 0) else specobj.optimal)
        ivar = specobj.optimal['COUNTS_IVAR']
        is2n = np.median(specobj.optimal['COUNTS']*np.sqrt(ivar))
        s2n.append(is2n)

    # Generate the table, if we have at least one source
    if len(names) > 0:
        obj_tbl = Table()
        if spectrograph.pypeline == 'MultiSlit':
            obj_tbl['slit'] = slits
            obj_tbl['slit'].format = 'd'
        elif spectrograph.pypeline == 'Echelle':
            obj_tbl['order'] = slits
            obj_tbl['order'].format = 'd'
        obj_tbl['name'] = names
        obj_tbl['spat_pixpos'] = spat_pixpos
        obj_tbl['spat_pixpos'].format = '.1f'
        obj_tbl['spat_fracpos'] = spat_fracpos
        obj_tbl['spat_fracpos'].format = '.3f'
        obj_tbl['box_width'] = boxsize
        obj_tbl['box_width'].format = '.2f'
        obj_tbl['box_width'].unit = units.arcsec
        obj_tbl['opt_fwhm'] = opt_fwhm
        obj_tbl['opt_fwhm'].format = '.3f'
        obj_tbl['opt_fwhm'].unit = units.arcsec
        obj_tbl['s2n'] = s2n
        obj_tbl['s2n'].format = '.2f'
        # Write
        obj_tbl.write(science_dir+'/objinfo_{:s}.txt'.format(basename),
                      format='ascii.fixed_width', overwrite=True)


def save_2d_images(sci_output, rawfile, ext0, master_key_dict, mfdir, outfile, clobber=True, update_det=None):
    """ Write 2D images to the hard drive

    Args:
        sci_output: OrderedDict
        fitstbl: Table
        scidx: int
        ext0: int
        master_key_dict: str
        mfdir: str
        outfile: str
        clobber: bool, optional

    Returns:

    """
    hdus, prihdu = init_hdus(update_det, outfile)
    if hdus is None:
        # Original header
        head0 = fits.getheader(rawfile, ext=ext0)

        # Primary HDU for output
        prihdu = fits.PrimaryHDU()
        # Update with original header, skipping a few keywords
        hdus = [prihdu]
        hdukeys = ['BUNIT', 'COMMENT', '', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'HISTORY', 'EXTEND', 'DATASEC']
        for key in head0.keys():
            # Use new ones
            if key in hdukeys:
                continue
            # Update unused ones
            prihdu.header[key] = head0[key]
        # History
        if 'HISTORY' in head0.keys():
            # Strip \n
            tmp = str(head0['HISTORY']).replace('\n', ' ')
            prihdu.header.add_history(str(tmp))

        # PYPEIT
        # TODO Should the spectrograph be written to the header?
        prihdu.header['PIPELINE'] = str('PYPEIT')
        prihdu.header['DATE-RDX'] = str(datetime.date.today().strftime('%Y-%b-%d'))
        prihdu.header['FRAMMKEY'] = master_key_dict['frame']
        prihdu.header['BPMMKEY'] = master_key_dict['bpm']
        prihdu.header['BIASMKEY']  = master_key_dict['bias']
        prihdu.header['ARCMKEY']  = master_key_dict['arc']
        prihdu.header['TRACMKEY']  = master_key_dict['trace']
        prihdu.header['FLATMKEY']  = master_key_dict['flat']
        prihdu.header['PYPMFDIR'] = str(mfdir)


    # Fill in the images
    ext = len(hdus) - 1
    for key in sci_output.keys():
        if key in ['meta']:
            continue
        else:
            det = key
        sdet = parse.get_dnum(det, caps=True)  # e.g. DET02
        if 'sciimg' not in sci_output[det]:
            continue
        # Specified detector number?
        #if settings.argflag['reduce']['detnum'] is not None:
        #    if det not in map(int, settings.argflag['reduce']['detnum']):
        #        continue
        #    else:
        #        msgs.warn("Restricting the reduction to detector {:d}".format(det))

        # Processed frame
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = '{:s}-Processed'.format(sdet)
        hdu = fits.ImageHDU(sci_output[det]['sciimg']) #slf._sciframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

        # Raw Inverse Variance
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = '{:s}-IVARRAW'.format(sdet)
        hdu = fits.ImageHDU(sci_output[det]['sciivar']) #slf._modelvarframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

        # Background model
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = '{:s}-SKY'.format(sdet)
        hdu = fits.ImageHDU(sci_output[det]['skymodel']) #slf._modelvarframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

        # Object model
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = '{:s}-OBJ'.format(sdet)
        hdu = fits.ImageHDU(sci_output[det]['objmodel']) #slf._modelvarframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

        # Inverse Variance model
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = '{:s}-IVARMODEL'.format(sdet)
        hdu = fits.ImageHDU(sci_output[det]['ivarmodel'])  # slf._modelvarframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)

        # Final mask
        ext += 1
        keywd = 'EXT{:04d}'.format(ext)
        prihdu.header[keywd] = '{:s}-MASK'.format(sdet)
        hdu = fits.ImageHDU(sci_output[det]['outmask'])  # slf._modelvarframe[det-1])
        hdu.name = prihdu.header[keywd]
        hdus.append(hdu)



    # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfile, overwrite=clobber)
    msgs.info("Wrote: {:s}".format(outfile))


def init_hdus(update_det, outfile):
    hdus, prihdu = None, None
    if (update_det is not None) and os.path.isfile(outfile):
        hdus = fits.open(outfile)
        msgs.info("Using existing spec1d file, including the Header")
        msgs.info("Will only update the data extension for {} detector(s)".format(update_det))
        prihdu = hdus[0]
        # Names
        hdu_names = [hdu.name for hdu in hdus]
        # Remove the detector(s) being updated
        if not isinstance(update_det, list):
            update_det = [update_det]
        popme = []
        # Find em
        for ss,hdu_name in enumerate(hdu_names):
            for det in update_det:
                sdet = parse.get_dnum(det, prefix=False)
                idx = '{:s}{:s}'.format(specobjs.naming_model['det'], sdet)
                if idx in hdu_name:
                    popme.append(ss)
        # Remove em (and the bit in the Header too)
        for popthis in reversed(popme):
            hdus.pop(popthis)
            keywd = 'EXT{:04d}'.format(popthis)
            prihdu.header.remove(keywd)
    # Return
    return hdus, prihdu
