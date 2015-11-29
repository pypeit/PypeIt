import os
import astropy.io.fits as pyfits
import numpy as np

def save_arcids(slf, fname, pixels):
    # Setup the HDU
    hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([hdu]) # Insert the primary HDU (input model)
    for o in xrange(len(pixels)):
        hdulist.append(pyfits.ImageHDU(pixels[o])) # Add a new Image HDU
    ans = 'y'
    if os.path.exists(fname):
        if slf._argflag['out']['overwrite']:
            os.remove(fname)
        else:
            ans = ''
            while ans != 'y' and ans != 'n' and ans != 'r':
                msgs.warn("File %s exists!" % (fname), verbose=slf._argflag['out']['verbose'])
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
            msgs.error("Could not save extraction"+msgs.newline()+
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
    hdulist[0].header["PIXSIZE"] = (slf._argflag['reduce']['pixelsize'], 'ARMED: The size of each sampled pixel (km/s)')
    # Loop through all orders and write the wavelength into the header
    if wave is not None:
        for i in xrange(sciext.shape[1]):
            hdrname = "CDELT{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (np.log10(1.0 + slf._argflag['reduce']['pixelsize']/299792.458), 'ARMED: log10(1+pixsize/c)'.format(frametype))
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
        for j in xrange(len(kys)):
            hkey = kys[j][:5].upper()
            if np.ndim(extprops[kys[j]])==1 and np.size(extprops[kys[j]]==sciext.shape[1]):
                for i in xrange(sciext.shape[1]):
                    hdrname = "{0:s}{1:03d}".format(hkey,i+1)
                    hdulist[0].header[hdrname] = (extprops[kys[j]][i], 'ARMED: {0:s} for order {1:d}'.format(kys[j],i+1))
    # Write the file to disk
    if os.path.exists(filename):
        if slf._argflag['out']['overwrite'] == True:
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
                if rmfil == 'a': slf._argflag['run']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("{0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    else:
        hdulist.writeto(filename)
        msgs.info("{0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    return


def save_master(slf, data, filename="temp.fits", frametype="<None>", ind=[]):
    msgs.info("Saving master {0:s} frame as:".format(frametype)+msgs.newline()+filename)
    hdu = pyfits.PrimaryHDU(data)
    hdulist = pyfits.HDUList([hdu])
    msgs.info("Writing header information")
    for i in xrange(len(ind)):
        hdrname = "FRAME{0:03d}".format(i+1)
        hdulist[0].header[hdrname] = (slf._fitsdict['filename'][ind[i]], 'PYPIT: File used to generate Master {0:s}'.format(frametype))
    hdulist[0].header["FRAMETYP"] = (frametype, 'PYPIT: Master calibration frame type')
    # Write the file to disk
    if os.path.exists(filename):
        if slf._argflag['out']['overwrite'] == True:
            msgs.warn("Overwriting file:"+msgs.newline()+filename)
            os.remove(filename)
            hdulist.writeto(filename)
            msgs.info("Master {0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
        else:
            msgs.warn("This file already exists")
            rmfil=''
            while rmfil != 'n' and rmfil != 'y' and rmfil != 'a':
                rmfil=raw_input(msgs.input()+"Remove this file? ([y]es, [n]o, or [a]lways) - ")
            if rmfil == 'n':
                msgs.warn("Not saving master {0:s} frame:".format(frametype)+msgs.newline()+filename)
            else:
                os.remove(filename)
                if rmfil == 'a': slf._argflag['run']['overwrite'] = True
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
        if slf._argflag['out']['overwrite'] == True:
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
                if rmfil == 'a': slf._argflag['run']['overwrite'] = True
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
        if slf._argflag['out']['overwrite'] == True:
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
                if rmfil == 'a': slf._argflag['run']['overwrite'] = True
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
        if slf._argflag['out']['overwrite'] == True:
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
                if rmfil == 'a': slf._argflag['run']['overwrite'] = True
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
        if slf._argflag['out']['overwrite'] == True:
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
                if rmfil == 'a': slf._argflag['run']['overwrite'] = True
                hdulist.writeto(filename)
                msgs.info("Saved saturation mask for frame:"+msgs.newline()+fname)
    else:
        hdulist.writeto(filename)
        msgs.info("Saved saturation mask for frame:"+msgs.newline()+fname)

    return