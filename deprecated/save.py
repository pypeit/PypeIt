""" Output for PYPEIT
"""
import os
import datetime
import warnings

import numpy as np

from astropy import units
from astropy.io import fits
from astropy.table import Table
import copy

from IPython import embed


import linetools.utils

from pypeit import msgs
from pypeit.core import parse
from pypeit import io


def save_all(sci_dict, master_key_dict, master_dir, spectrograph, head1d, head2d, scipath, basename,
             update_det=None, binning='None'):
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
        refframe: str, default = 'heliocentric'
            Reference frame for the wavelengths
        update_det : int or list, default=None
            If provided, do not clobber the existing file but only update
            the indicated detectors.  Useful for re-running on a subset of detectors

    Returns:

    """
    from pypeit import specobjs # Hiding this to avoid circular import.  Will remove
    # Check for the directory
    if not os.path.isdir(scipath):
        os.makedirs(scipath)

    # Filenames to write out
    # TODO: These should be centrally defined so that they don't become
    # out of sync with what's in pypeit.PypeIt
    outfile1d = os.path.join(scipath, 'spec1d_{:s}.fits'.format(basename))
    outfile2d = os.path.join(scipath, 'spec2d_{:s}.fits'.format(basename))
    outfiletxt = os.path.join(scipath, 'spec1d_{:s}.txt'.format(basename))

    # TODO: Need some checks here that the exposure has been reduced

    # Build the final list of specobjs and vel_corr
    all_specobjs = specobjs.SpecObjs()

    for key in sci_dict:
        try:
            all_specobjs.add_sobj(sci_dict[key]['specobjs'])
        except KeyError:  # No object extracted
            continue

    if len(all_specobjs) == 0:
        msgs.warn('No objects to save. Only writing spec2d files!')
    else:
        # Build the spec1d output header.
        header = all_specobjs.build_header(head1d, head2d, spectrograph)
        all_specobjs.write_to_fits(header, outfile1d, update_det=update_det)
        # Txt file
        # TODO JFH: Make this a method in the specobjs class.
        save_obj_info(all_specobjs, spectrograph, outfiletxt, sci_dict, binning=binning)

    # Write 2D images for the Science Frame
    save_2d_images(sci_dict, head2d, spectrograph, master_key_dict, master_dir, outfile2d, update_det=update_det)

    return



# TODO: (KBW) I don't think core algorithms should take class
# arguments...
# TODO JFH: we make exceptions for core objects like specobjs
# TODO JXP: This will be elimiated as core and put back into specobjs once
#  we have proper datamodels for all these.
def save_obj_info(all_specobjs, spectrograph, outfile, sci_dict, binning='None'):
    """
    Write info to an ASCII file

    Args:
        all_specobjs (specobjs.SpecObjs):
        spectrograph (spectrograph.Spectrograph):
        outfile (str):
        binning (str, optional):

    Returns:

    """
    slits, names, spat_pixpos, spat_fracpos, boxsize, opt_fwhm, s2n = [], [], [], [], [], [], []  # Lists for a Table
    #binspectral, binspatial = parse.parse_binning(binning)
    for specobj in all_specobjs.specobjs:
        det = specobj.DET
        if specobj is None:
            continue
        # Detector items
        binspectral, binspatial = parse.parse_binning(sci_dict[det]['detector'].binning)
        platescale = sci_dict[det]['detector'].platescale
        # Append
        spat_pixpos.append(specobj.SPAT_PIXPOS)
        if spectrograph.pypeline == 'MultiSlit':
            spat_fracpos.append(specobj.SPAT_FRACPOS)
            slits.append(specobj.SLITID)
            names.append(specobj.NAME)
        elif spectrograph.pypeline == 'Echelle':
            spat_fracpos.append(specobj.ECH_FRACPOS)
            slits.append(specobj.ECH_ORDER)
            names.append(specobj.ECH_NAME)
        # Boxcar width
        if 'BOX_RADIUS' in specobj.keys():
            slit_pix = 2.0*specobj.BOX_RADIUS
            # Convert to arcsec
            binspectral, binspatial = parse.parse_binning(binning)
            # JFH TODO This should be using the order_platescale for each order. Furthermore, not all detectors
            # have the same platescale, i.e. with GNIRS it is the same detector but a different camera hence a
            # different attribute. platescale should be a spectrograph attribute determined on the fly.
            #boxsize.append(slit_pix*binspatial*spectrograph.detector[specobj.DET-1]['platescale'])
            boxsize.append(slit_pix*binspatial*platescale)
        else:
            boxsize.append(0.)

        # Optimal profile (FWHM)
        # S2N -- default to boxcar
        if hasattr(specobj, 'FWHMFIT'):
            #opt_fwhm.append(np.median(specobj.FWHMFIT)* binspatial*spectrograph.detector[specobj.DET-1]['platescale'])
            opt_fwhm.append(np.median(specobj.FWHMFIT)* binspatial*platescale)
            # S2N -- optimal
            ivar = specobj.OPT_COUNTS_IVAR
            is2n = np.median(specobj.OPT_COUNTS*np.sqrt(ivar))
            s2n.append(is2n)
        else: # Optimal is not required to occur
            opt_fwhm.append(0.)
            # S2N -- use boxcar
            ivar = specobj.BOX_COUNTS_IVAR
            is2n = np.median(specobj.BOX_COUNTS*np.sqrt(ivar))
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
    if os.path.isfile(outfile) and update_det is not None:
        hdus, prihdu = io.init_hdus(update_det, outfile)
    else:
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
        prihdu.header['PYPELINE'] = spectrograph.pypeline
        prihdu.header['SPECTROG'] = spectrograph.spectrograph
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



def save_sens_dict(sens_dict, outfile, overwrite=True):
    """
    Over-load the save_master() method in MasterFrame to write a FITS file

    Parameters
    ----------
    outfile : str, optional
      Use this input instead of the 'proper' (or unattainable) MasterFrame name

    Returns
    -------

    """

    if os.path.exists(outfile)and (not overwrite):
        msgs.warn("This file already exists.  Use overwrite=True to overwrite it")
        return
    #

    # jsonify has the annoying property that it modifies the objects when it jsonifies them so make a copy,
    # which converts lists to arrays, so we make a copy
    data_for_json = copy.deepcopy(sens_dict)
    gddict = linetools.utils.jsonify(data_for_json)
    linetools.utils.savejson(outfile, gddict, easy_to_read=True, overwrite=True)
    # Finish
    msgs.info("Sucessfuly save sensitivity function to file {:s}".format(outfile))


def write_fits(hdr, data, outfile, extnames=None, checksum=True):
    """
    Convenience method to write a set of data to a multi-extension FITS file

    Args:
        hdr (`astropy.io.fits.Header`):
            Header to be written to the primary image
        data (np.ndarray or list):
            One or more images to be written
        outfile (str):
        extnames (list, optional):
            Extension names to be used for each data item
        checksum (bool, optional):
            Add CHECKSUM to the header

    Returns:

    """
    warnings.warn("To be depracated!!")
    # Format the output
    ext = extnames if isinstance(extnames, list) else [extnames]
    if len(ext) > 1 and not isinstance(data, list):
        msgs.error('Input data type should be list, one numpy.ndarray per extension.')
    _data = data if isinstance(data, list) else [data]

    # Write the fits file
    fits.HDUList([fits.PrimaryHDU(header=hdr)]
                 + [fits.ImageHDU(data=d, name=n) for d, n in zip(_data, ext)]
                 ).writeto(outfile, overwrite=True, checksum=checksum)
