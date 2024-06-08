

# JFH This routine is deprecated. The preferred way to write out 1d coadds is now in the coadd1d class.
def save_coadd1d_to_fits(outfile, waves, fluxes, ivars, masks, telluric=None, obj_model=None,
                         header=None, ex_value='OPT', overwrite=True):
    '''
    Args:
        outfile (str): name of fitsfile you want to save to
        waves (ndarray): 1-D or 2-D (nspec by nexp/norder) wavelength array
        fluxes (ndarray): flux array
        ivars (ndarray): ivar array
        masks (ndarray): mask array
        header (dict): primary fits header
        ext_value (str): 'OPT' for optimal, and 'BOX' for boxcar
        overwrite (bool): if True, overwrite the old one, otherwise append it to the exist fits file.
    Returns:
        None
    '''

    # Estimate sigma from ivar
    sigs = np.sqrt(utils.inverse(ivars))

    if (os.path.exists(outfile)) and (np.invert(overwrite)):
        hdulist = fits.open(outfile)
        msgs.info("Reading primary HDU from existing file: {:s}".format(outfile))
    else:
        msgs.info("Creating an new primary HDU.")
        prihdu = fits.PrimaryHDU()
        if header is None:
            msgs.warn('The primary header is none')
        else:
            prihdu.header = header
        hdulist = fits.HDUList([prihdu])

    if waves.ndim == 1:
        wave_mask = waves > 1.0
        # Add Spectrum Table
        cols = []
        cols += [fits.Column(array=waves[wave_mask], name='{:}_WAVE'.format(ex_value), format='D')]
        cols += [fits.Column(array=fluxes[wave_mask], name='{:}_FLAM'.format(ex_value), format='D')]
        cols += [fits.Column(array=ivars[wave_mask], name='{:}_FLAM_IVAR'.format(ex_value), format='D')]
        cols += [fits.Column(array=sigs[wave_mask], name='{:}_FLAM_SIG'.format(ex_value), format='D')]
        cols += [fits.Column(array=masks[wave_mask].astype(float), name='{:}_MASK'.format(ex_value), format='D')]
        if telluric is not None:
            cols += [fits.Column(array=telluric[wave_mask], name='TELLURIC', format='D')]
        if obj_model is not None:
            cols += [fits.Column(array=obj_model[wave_mask], name='OBJ_MODEL', format='D')]

        coldefs = fits.ColDefs(cols)
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = 'OBJ0001-SPEC0001-{:}'.format(ex_value.capitalize())
        hdulist.append(tbhdu)
    else:
        nspec = waves.shape[1]

        for ispec in range(nspec):
            wave_mask = waves[:,ispec] > 1.0
            # Add Spectrum Table
            cols = []
            cols += [fits.Column(array=waves[:,ispec][wave_mask], name='{:}_WAVE'.format(ex_value), format='D')]
            cols += [fits.Column(array=fluxes[:,ispec][wave_mask], name='{:}_FLAM'.format(ex_value), format='D')]
            cols += [fits.Column(array=ivars[:,ispec][wave_mask], name='{:}_FLAM_IVAR'.format(ex_value), format='D')]
            cols += [fits.Column(array=sigs[:,ispec][wave_mask], name='{:}_FLAM_SIG'.format(ex_value), format='D')]
            cols += [fits.Column(array=masks[:,ispec][wave_mask].astype(float), name='{:}_MASK'.format(ex_value), format='D')]

            coldefs = fits.ColDefs(cols)
            tbhdu = fits.BinTableHDU.from_columns(coldefs)
            tbhdu.name = 'OBJ0001-SPEC{:04d}-{:}'.format(ispec+1, ex_value.capitalize())
            hdulist.append(tbhdu)

    if (os.path.exists(outfile)) and (np.invert(overwrite)):
        hdulist.writeto(outfile, overwrite=True)
        msgs.info("Appending 1D spectra to existing file {:s}".format(outfile))
    else:
        hdulist.writeto(outfile, overwrite=overwrite)
        msgs.info("Wrote 1D spectra to {:s}".format(outfile))

    return None

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


#def write_sensitivity():
    #sensfunc_name = "{0:s}/{1:s}/{2:s}_{3:03d}_{4:s}.yaml".format(os.getcwd(), settings.argflag['run']['directory']['master'], slf._fitsdict['target'][scidx[0]], 0, "sensfunc")
    #msgs.info("Writing sensfunc: {:s}".format(sensfunc_name))
    #with open(sensfunc_name, 'w') as yamlf:
    #    yamlf.write( yaml.dump(slf._sensfunc))
    #with io.open(sensfunc_name, 'w', encoding='utf-8') as f:
    #    f.write(unicode(json.dumps(slf._sensfunc, sort_keys=True, indent=4, separators=(',', ': '))))