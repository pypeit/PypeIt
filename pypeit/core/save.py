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


def save_all(sci_dict, master_key_dict, master_dir, spectrograph, head1d, head2d, scipath, basename,
                  only_1d=False, refframe='heliocentric', update_det=None, binning='None'):
    """
    Routine to save PypeIt 1d and 2d outputs
    Args:
        sci_dict: dict
            Dictionary containing extraction outputs
        master_key_dict: dict
            Dictionary with master key information for this reduction
        master_dir: str
            Directory where the master files live
        spectrograph: object, spectrograph
            Spectrograph object for the spectorgraph that was used
        head1d: dict
            fitstbl meta data dictionary that will become the header for the spec1d files
        head2d: dict
            rawfile header that will become the header for the spec2d files
        scipath: str
            path to which the outputs should be written
        basename: str
            the object basename
        only_1d: bool, default = False
            Only write out the
        refframe: str, default = 'heliocentric'
            Reference frame for the wavelengths
        update_det : int or list, default=None
            If provided, do not clobber the existing file but only update
            the indicated detectors.  Useful for re-running on a subset of detectors
        binning: str, default = None
          String indicating the binning of the data

    Returns:

    """

    # Filenames to write out
    objinfofile = scipath + '/objinfo_{:s}.txt'.format(basename)
    outfile1d = os.path.join(scipath, 'spec1d_{:s}.fits'.format(basename))
    outfile2d = scipath + '/spec2d_{:s}.fits'.format(basename)

    # TODO: Need some checks here that the exposure has been reduced

    # Build the final list of specobjs and vel_corr
    all_specobjs = specobjs.SpecObjs()

    vel_corr = 0.  # This will not be set for Standard stars, which is fine
    for key in sci_dict:
        if key in ['meta']:
            vel_corr = sci_dict['meta']['vel_corr']
            continue
        #
        try:
            all_specobjs.add_sobj(sci_dict[key]['specobjs'])
        except KeyError:  # No object extracted
            continue

    if len(all_specobjs) == 0:
        msgs.warn('No objects to save!')
        return

    # Create the helio_dict
    helio_dict = dict(refframe=refframe, vel_correction=vel_corr)

    save_1d_spectra_fits(all_specobjs, head1d, spectrograph, outfile1d,helio_dict=helio_dict, update_det=update_det)

    # 1D only?
    if only_1d:
        return
    # Obj info
    save_obj_info(all_specobjs, spectrograph, objinfofile, binning=binning)
    # Write 2D images for the Science Frame

    # TODO: Make sure self.det is correct!
    # master_key = self.fitstbl.master_key(frame, det=self.det)
    # TODO: Why is the raw file header needed?  Can the data be gotten from fitstbl?
    #  If not, is it worth adding the relevant information to fitstbl?
    # Original header
    save_2d_images(sci_dict, head2d, spectrograph.spectrograph, master_key_dict, master_dir, outfile2d, update_det=update_det)

    return


def save_1d_spectra_fits(specObjs, header, spectrograph, outfile, helio_dict=None, overwrite=True, update_det=None):
    """ Write 1D spectra to a multi-extension FITS file

    Args:
        specobjs : SpecObjs object
        header (dict or Row; dict-like):  Typically a Row from the fitstbl
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
          Name of PypeIt pipeline (e.g. 'MultiSlit')
        outfile (str):
        helio_dict (dict, optional):
        overwrite : bool, optional
        update_det : int or list, optional
          If provided, do not clobber the existing file but only update
          the indicated detectors.  Useful for re-running on a subset of detectors

    Returns:
        str: outfile

    """

    pypeline = spectrograph.pypeline
    instrume = spectrograph.spectrograph
    telescope = spectrograph.telescope
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



# TODO: (KBW) I don't think core algorithms should take class
# arguments...
def save_obj_info(all_specobjs, spectrograph, outfile, binning='None'):
    """

    Parameters
    ----------
    all_specobjs : list
    fitstbl : Table

    Returns
    -------

    """
    slits, names, spat_pixpos, spat_fracpos, boxsize, opt_fwhm, s2n = [], [], [], [], [], [], []  # Lists for a Table
    binspectral, binspatial = parse.parse_binning(binning)
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
            binspectral, binspatial = parse.parse_binning(binning)
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
        obj_tbl.write(outfile,format='ascii.fixed_width', overwrite=True)


#TODO 2d data model should be expanded to include:
# waveimage  --  flexure and heliocentric corrections should be applied to the final waveimage and since this is unique to
#                every exposure (i.e. it depneds on obstime, RA, DEC and the flexure incurred) it should be written out for
#                each science frame.
# tslits_dict -- flexure compensation implies that each frame will have a unique set of slit boundaries, so we probably need to
#                 write these for each file as well. Alternatively we could just write the offsets to the header.
def save_2d_images(sci_output, raw_header, spectrograph, master_key_dict, mfdir, outfile, clobber=True, update_det=None):
    """ Write 2D images to the hard drive

    Args:
        sci_output (OrderedDict):
        raw_header (astropy.fits.Header or dict):
        master_key_dict (str):
        mfdir (str):
        outfile (str):
        clobber: bool, optional

    Returns:

    """
    hdus, prihdu = init_hdus(update_det, outfile)
    if hdus is None:
        # Primary HDU for output
        prihdu = fits.PrimaryHDU()
        # Update with original header, skipping a few keywords
        hdus = [prihdu]
        hdukeys = ['BUNIT', 'COMMENT', '', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'HISTORY', 'EXTEND', 'DATASEC']
        for key in raw_header.keys():
            # Use new ones
            if key in hdukeys:
                continue
            # Update unused ones
            prihdu.header[key] = raw_header[key]
        # History
        if 'HISTORY' in raw_header.keys():
            # Strip \n
            tmp = str(raw_header['HISTORY']).replace('\n', ' ')
            prihdu.header.add_history(str(tmp))

        # PYPEIT
        # TODO Should the spectrograph be written to the header?
        prihdu.header['PIPELINE'] = str('PYPEIT')
        prihdu.header['SPECTROG'] = spectrograph
        prihdu.header['DATE-RDX'] = str(datetime.date.today().strftime('%Y-%b-%d'))
        prihdu.header['FRAMMKEY'] = master_key_dict['frame'][:-3]
        prihdu.header['BPMMKEY'] = master_key_dict['bpm'][:-3]
        prihdu.header['BIASMKEY']  = master_key_dict['bias'][:-3]
        prihdu.header['ARCMKEY']  = master_key_dict['arc'][:-3]
        prihdu.header['TRACMKEY']  = master_key_dict['trace'][:-3]
        prihdu.header['FLATMKEY']  = master_key_dict['flat'][:-3]
        prihdu.header['PYPMFDIR'] = str(mfdir)
        if sci_output['meta']['ir_redux']:
            prihdu.header['SKYSUB'] ='DIFF'
        else:
            prihdu.header['SKYSUB'] ='MODEL'



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
