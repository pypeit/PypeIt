# Shapely is needed if using the resample algorithm
try:
    import shapely
except ImportError:
    shapely = None


def generate_cube_resample(outfile, frame_wcs, slits, fluximg, ivarimg, raimg, decimg, waveimg, slitimg, gpm,
                           grid_nspat=5, grid_specsep=20,
                           overwrite=False, output_wcs=None, blaze_wave=None, blaze_spec=None, fluxcal=False,
                           sensfunc=None, specname=None, debug=False):
    """
    Save a datacube using the resample algorithm.

    This function takes the fully calibrated input data, and resamples
    all slices onto a regular 3D grid, while conserving flux. Note that
    the final datacube has correlations between voxels, and this covariance
    information is not saved.

    Args:
        outfile (`str`):
            Filename to be used to save the datacube
        frame_wcs (`astropy.wcs.WCS`_):
            World coordinate system for this frame.
        slits (:class:`pypeit.slittrace.SlitTraceSet`_)
            Information stored about the slits
        fluximg (`numpy.ndarray`_):
            Surface brightness of each pixel in the frame (units = erg/s/cm^2/A/arcsec^2)
        ivarimg (`numpy.ndarray`_):
            Inverse variance of each pixel in the frame
        raimg (`numpy.ndarray`_):
            Right ascension of each pixel in the frame (units = decimal degrees)
        decimg (`numpy.ndarray`_):
            Declination of each pixel in the frame (units = decimal degrees)
        waveimg (`numpy.ndarray`_):
            Wavelength of each pixel in the frame (units = Angstroms)
        slitimg (`numpy.ndarray`_):
            Slit image. -1 is not on a slit, and all other
            pixels are labelled with their spatial IDs.
        gpm (`numpy.ndarray`_):
            Good pixel mask (bool). True = good pixel
        grid_nspat (int, optional):
            Number of grid points in the spatial direction when evaluating the
            voxel geometry in detector coordinates. This should be an odd number
        grid_specsep (int, optional):
            Number of pixels between each grid point in the spectral direction
            when evaluating the voxel geometry in detector coordinates
        overwrite (bool, optional):
            If the output file exists, it will be overwritten if this parameter is True.
        output_wcs (`astropy.wcs.WCS`_, optional):
            World coordinate system for the output datacube. If None, frame_wcs will be used.
        blaze_wave (`numpy.ndarray`_, optional):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_, optional):
            Spectral blaze function
        fluxcal (bool, optional):
            Are the data flux calibrated? If True, the units are: erg/s/cm^2/Angstrom/arcsec^2
            multiplied by the PYPEIT_FLUX_SCALE. Otherwise, the units are: counts/s/Angstrom/arcsec^2")
        sensfunc (`numpy.ndarray`_, None, optional):
            Sensitivity function that has been applied to the datacube
        specname (str, None, optional):
            Name of the spectrograph
        debug (bool, optional):
            Debug the code by writing out a residuals cube?
    """
    # Set the output_wcs if it's not already set
    if output_wcs is None:
        output_wcs = frame_wcs
    # Check that grid_nspat is an odd number
    if grid_nspat % 2 == 0:
        msgs.warn(f"grid_nspat must be an odd number. Using grid_nspat={grid_nspat+1} instead")
        grid_nspat += 1
    debug = False
    # Get the grid spacing along the spatial direction
    frm_cd_spat = np.sqrt(frame_wcs.wcs.cd[1, 1] ** 2 + frame_wcs.wcs.cd[0, 1] ** 2)
    out_cd_spat = np.sqrt(output_wcs.wcs.cd[1, 1] ** 2 + output_wcs.wcs.cd[0, 1] ** 2)
    slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
    nvox_spat = int(np.ceil(slitlength*frm_cd_spat/out_cd_spat))
    crd_vox_spat = out_cd_spat * (np.arange(nvox_spat+1) - (nvox_spat+1)// 2)  # +1 to get bin edges
    # Get the grid spacing along the spectral direction
    out_cr_wave = output_wcs.wcs.crval[2]
    out_cd_wave = output_wcs.wcs.cd[2, 2]
    nvox_wave = int(np.ceil((np.max(waveimg)-out_cr_wave)/out_cd_wave))
    crd_vox_spec = out_cr_wave + out_cd_wave * np.arange(nvox_wave+1)  # +1 to get bin edges
    vox_shape = (nvox_wave+1, nvox_spat+1)

    # Detector spectal/spatial pixels and number of slices
    nspec, nspat, nslice = slits.nspec, slits.nspat, slits.spat_id.size

    # Generate the output datacube
    datcube = np.zeros((nslice, nvox_spat, nvox_wave), dtype=float)
    varcube = np.zeros((nslice, nvox_spat, nvox_wave), dtype=float)

    # Transform the voxel geometry to detector pixels
    grid_nspec = 1 + nspec // grid_specsep
    xgrid = np.zeros((grid_nspec, grid_nspat), dtype=int)
    ygridt = np.zeros(grid_nspec, dtype=int)
    ygridt[-1] = nspec - 1
    ygridt[1:-1] = (nspec % grid_specsep + 2 * grid_specsep) // 2 + np.arange(grid_nspec - 2) * grid_specsep
    ygrid = ygridt[:, np.newaxis].repeat(grid_nspat, axis=1)
    ra0, dec0 = np.zeros(nslice), np.zeros(nslice)
    offsimg = np.zeros_like(waveimg)
    varimgsq = utils.inverse(ivarimg ** 2)
    for sl, spat_id in enumerate(slits.spat_id):
        msgs.info(f"Calculating voxel geometry for slit {spat_id}")
        # Calculate RA and Dec of central traces
        wsl = np.where(slitimg == spat_id)
        this_ra, this_dec, this_wave = raimg[wsl], decimg[wsl], waveimg[wsl]
        _, spat_posn, _ = frame_wcs.wcs_world2pix(this_ra, this_dec, this_wave*1.0E-10, 0)
        asrt = np.argsort(spat_posn)
        ra0[sl] = np.interp(0.0, spat_posn[asrt], this_ra[asrt])
        dec0[sl] = np.interp(0.0, spat_posn[asrt], this_dec[asrt])
        # Generate the offsets
        cosdec = np.cos(dec0[sl] * np.pi / 180.0)
        diff_ra, diff_dec = (this_ra - ra0[sl]) * cosdec, this_dec - dec0[sl]
        msgs.bug("There is sometimes a sign error that needs to be resolved here...")
        msgs.error("Use another algorithm for the time being...")
        if np.max(diff_ra)-np.min(diff_ra) > np.max(diff_dec)-np.min(diff_dec):
            sgn = np.sign(diff_ra)
        else:
            sgn = np.sign(diff_dec)
        offsimg[wsl] = -sgn * np.sqrt(diff_ra**2 + diff_dec**2)
        # Update the xgrid values for this slice
        for yy in range(grid_nspec):
            wsl = np.where(slitimg == spat_id)
            allind = wsl[1][np.where(wsl[0] == ygridt[yy])]
            xgrid[yy, 0] = np.min(allind)
            xgrid[yy, -1] = np.max(allind)
            numpix = xgrid[yy, -1] - xgrid[yy, 0]
            sep = numpix // (grid_nspat - 1)
            xgrid[yy, 1:-1] = xgrid[yy, 0] + (numpix % sep + 2 * sep) // 2 + np.arange(grid_nspat - 2) * sep
        # Extract offset + wavelength information and estimate transform
        grid_coord = (ygrid.flatten(), xgrid.flatten())
        grid_offs = offsimg[grid_coord]
        grid_wave = waveimg[grid_coord]
        src = np.column_stack((grid_wave, grid_offs))
        dst = np.column_stack(grid_coord).astype(float)
        # Transform the voxel coordinates to detector coordinates
        evalpos = np.column_stack((crd_vox_spec[:,np.newaxis].repeat(crd_vox_spat.size, axis=1).flatten(),
                                   crd_vox_spat[np.newaxis,:].repeat(crd_vox_spec.size, axis=0).flatten()))
        # tform = LinearNDInterpolator(src, dst, rescale=True)
        # crd_det_tmp = tform(evalpos)

        src_off = np.min(src, axis=0)
        src_scl = np.max(src-src_off, axis=0)
        dst_off = np.min(dst, axis=0)
        dst_scl = np.max(dst-dst_off, axis=0)
        tform = RBFInterpolator((src-src_off)/src_scl, (dst-dst_off)/dst_scl, smoothing=0.01)
        crd_det = dst_off + dst_scl * tform((evalpos-src_off)/src_scl)
        if debug:
            plt.plot(crd_det[:, 0], crd_det[:, 1], 'rx')
            #plt.plot(crd_det_tmp[:, 0], crd_det_tmp[:, 1], 'bx')
            plt.plot(np.arange(slits.left_init.shape[0]), slits.left_init[:, 0], 'k-')
            plt.plot(np.arange(slits.right_init.shape[0]), slits.right_init[:, 0], 'k-')
            plt.show()

    # Calculate an "offsets" image, which indicates the offset in arcsec from (RA_0, DEC_0)
    # Create two splines of the offsets image: (1) offset predicts RA; (2) offset predicts Dec.
    # Use these splines to calculate the RA and DEC of the voxels, combine this with the output wavelength grid.
    # Generate all RA, DEC, WAVELENGTH triples (i.e. find the RA,DEC pairs along constant wavelength, for all wavelengths)
    # Use the WCS (which contains the astrometric transform) to go from world to pix
    #    i.e. need to invert this:
    #    world_ra, world_dec, _ = wcs.wcs_pix2world(slitID, evalpos, tilts[onslit_init]*(nspec-1), 0)
    # This gives us the x,y detector positions of the voxel geometry
        from shapely.geometry import Polygon, box as shapelyBox
        from shapely.strtree import STRtree

        crd_det_spec, crd_det_spat = crd_det[:, 0].reshape(vox_shape), crd_det[:, 1].reshape(vox_shape)
        # Generate a list of all detector pixels in this slice
        detpix_polys = []
        pix_spec, pix_spat = np.where(slitimg == spat_id)
        for ss in range(pix_spat.size):
            detpix_polys.append(shapely.geometry.box(pix_spat[ss], pix_spec[ss], pix_spat[ss]+1, pix_spec[ss]+1))
        # Create a Sort-Tile-Recursive tree of the detector pixels to quickly query overlapping voxels
        detgeom = shapely.strtree.STRtree(detpix_polys)
        # Loop through all voxels for this slice and calculate the overlapping area
        for wv in range(nvox_wave):
            for sp in range(nvox_spat):
                # Generate the voxel coordinates in detector pixel space (points must be counter-clockwise)
                voxel_geom = shapely.geometry.Polygon([(crd_det_spat[wv, sp],   crd_det_spec[wv,   sp]),
                                                       (crd_det_spat[wv, sp+1], crd_det_spec[wv,   sp]),
                                                       (crd_det_spat[wv, sp+1], crd_det_spec[wv+1, sp]),
                                                       (crd_det_spat[wv, sp],   crd_det_spec[wv+1, sp]),
                                                       (crd_det_spat[wv, sp],   crd_det_spec[wv,   sp])])
                # Find overlapping detector pixels
                result = detgeom.query(voxel_geom)
                # Sum all overlapping flux-weighted areas
                this_flx = 0
                this_var = 0
                this_area = 0
                for pp in range(len(result)):
                    area = voxel_geom.intersection(result[pp]).area
                    pix_spat = int(min(result[pp].exterior.coords[0][0], result[pp].exterior.coords[2][0]))
                    pix_spec = int(min(result[pp].exterior.coords[0][1], result[pp].exterior.coords[2][1]))
                    if ivarimg[pix_spec, pix_spat] != 0.0:
                        this_flx += area * fluximg[pix_spec, pix_spat]
                        this_var += area**2 * varimgsq[pix_spec, pix_spat]
                        this_area += area
                # Fill in the datacube
                this_area = 1 if this_area == 0 else this_area
                datcube[sl, sp, wv] = this_flx / this_area
                varcube[sl, sp, wv] = this_var / this_area**2

    # Generate a header
    hdr = output_wcs.to_header()

    # Add the unit of flux to the header
    if fluxcal:
        hdr['FLUXUNIT'] = (PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
    else:
        hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")

    # Save the final datacube
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    final_cube = DataCube(datcube.T, varcube.T, specname, blaze_wave, blaze_spec, sensfunc=sensfunc, fluxed=fluxcal)
    final_cube.to_file(outfile, hdr=hdr, overwrite=overwrite)


def generate_cube_ngp(outfile, hdr, all_sci, all_ivar, all_wghts, vox_coord, bins,
                      overwrite=False, blaze_wave=None, blaze_spec=None, fluxcal=False,
                      sensfunc=None, specname="PYP_SPEC", debug=False):
    """
    TODO :: Deprecate this routine once subpixellate works for combining cubes

    Save a datacube using the Nearest Grid Point (NGP) algorithm.

    Args:
        outfile (`str`):
            Filename to be used to save the datacube
        hdr (`astropy.io.fits.header_`):
            Header of the output datacube (must contain WCS)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used in the combination
        vox_coord (`numpy.ndarray`_):
            The voxel coordinates of each pixel in the spec2d frames. vox_coord is returned by the
            function `astropy.wcs.WCS.wcs_world2pix_` once a WCS is setup and every spec2d detector
            pixel has an RA, DEC, and WAVELENGTH.
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial and z wavelength coordinates
        overwrite (`bool`):
            If True, the output cube will be overwritten.
        blaze_wave (`numpy.ndarray`_):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_):
            Spectral blaze function
        fluxcal (bool):
            Are the data flux calibrated?
        sensfunc (`numpy.ndarray`_, None):
            Sensitivity function that has been applied to the datacube
        specname (str):
            Name of the spectrograph
        debug (bool):
            Debug the code by writing out a residuals cube?
    """
    # Add the unit of flux to the header
    if fluxcal:
        hdr['FLUXUNIT'] = (flux_calib.PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
    else:
        hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")

    # Use NGP to generate the cube - this ensures errors between neighbouring voxels are not correlated
    datacube, edges = np.histogramdd(vox_coord, bins=bins, weights=all_sci * all_wghts)
    normcube, edges = np.histogramdd(vox_coord, bins=bins, weights=all_wghts)
    nc_inverse = utils.inverse(normcube)
    datacube *= nc_inverse
    # Create the variance cube, including weights
    msgs.info("Generating variance cube")
    all_var = utils.inverse(all_ivar)
    var_cube, edges = np.histogramdd(vox_coord, bins=bins, weights=all_var * all_wghts ** 2)
    var_cube *= nc_inverse**2
    bpmcube = (normcube == 0).astype(np.uint8)

    # Save the datacube
    if debug:
        datacube_resid, edges = np.histogramdd(vox_coord, bins=bins, weights=all_sci * np.sqrt(all_ivar))
        normcube, edges = np.histogramdd(vox_coord, bins=bins)
        nc_inverse = utils.inverse(normcube)
        outfile_resid = "datacube_resid.fits"
        msgs.info("Saving datacube as: {0:s}".format(outfile_resid))
        hdu = fits.PrimaryHDU((datacube_resid*nc_inverse).T, header=hdr)
        hdu.writeto(outfile_resid, overwrite=overwrite)

    msgs.info("Saving datacube as: {0:s}".format(outfile))
    final_cube = DataCube(datacube.T, np.sqrt(var_cube.T), bpmcube.T, specname, blaze_wave, blaze_spec,
                          sensfunc=sensfunc, fluxed=fluxcal)
    final_cube.to_file(outfile, hdr=hdr, overwrite=overwrite)


def gaussian2D_cube(tup, intflux, xo, yo, dxdz, dydz, sigma_x, sigma_y, theta, offset):
    """
    Fit a 2D Gaussian function to a datacube. This function assumes that each
    wavelength slice of the datacube is well-fit by a 2D Gaussian. The centre of
    the Gaussian is allowed to vary linearly as a function of wavelength.

    .. note::

        The integrated flux does not vary with wavelength.

    Args:
        tup (:obj:`tuple`):
            A three element tuple containing the x, y, and z locations of each
            pixel in the cube
        intflux (float):
            The Integrated flux of the Gaussian
        xo (float):
            The centre of the Gaussian along the x-coordinate when z=0
        yo (float):
            The centre of the Gaussian along the y-coordinate when z=0
        dxdz (float):
            The change of xo with increasing z
        dydz (float):
            The change of yo with increasing z
        sigma_x (float):
            The standard deviation in the x-direction
        sigma_y (float):
            The standard deviation in the y-direction
        theta (float):
            The orientation angle of the 2D Gaussian
        offset (float):
            Constant offset

    Returns:
        `numpy.ndarray`_: The 2D Gaussian evaluated at the coordinate (x, y, z)
    """
    # Extract the (x, y, z) coordinates of each pixel from the tuple
    (x, y, z) = tup
    # Calculate the centre of the Gaussian for each z coordinate
    xo = float(xo) + z*dxdz
    yo = float(yo) + z*dydz
    # Account for a rotated 2D Gaussian
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    # Normalise so that the integrated flux is a parameter, instead of the amplitude
    norm = 1/(2*np.pi*np.sqrt(a*c-b*b))
    gtwod = offset + norm*intflux*np.exp(-(a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return gtwod.ravel()


def make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                               all_ivar=None, whitelightWCS=None, numra=None, numdec=None, trim=1):
    """
    Generate a whitelight image using the individual pixels of every input frame

    Args:
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the RA values of each pixel from all
            spec2d files
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the DEC values of each pixel from all
            spec2d files
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength values of each pixel
            from all spec2d files
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all
            spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights attributed to each pixel
            from all spec2d files
        all_idx (`numpy.ndarray`_):
            1D flattened array containing an integer identifier indicating which
            spec2d file each pixel originates from. For example, a 0 would
            indicate that a pixel originates from the first spec2d frame listed
            in the input file. a 1 would indicate that this pixel originates
            from the second spec2d file, and so forth.
        dspat (float):
            The size of each spaxel on the sky (in degrees)
        all_ivar (`numpy.ndarray`_, optional):
            1D flattened array containing of the inverse variance of each pixel
            from all spec2d files.  If provided, inverse variance images will be
            calculated and returned for each white light image.
        whitelightWCS (`astropy.wcs.WCS`_, optional):
            The WCS of a reference white light image. If supplied, you must also
            supply numra and numdec.
        numra (int, optional):
            Number of RA spaxels in the reference white light image
        numdec (int, optional):
            Number of DEC spaxels in the reference white light image
        trim (int, optional):
            Number of pixels to grow around a masked region

    Returns:
        tuple: two 3D arrays will be returned, each of shape [N, M, numfiles],
        where N and M are the spatial dimensions of the combined white light
        images.  The first array is a white light image, and the second array is
        the corresponding inverse variance image. If all_ivar is None, this will
        be an empty array.
    """
    # Determine number of files
    numfiles = np.unique(all_idx).size

    if whitelightWCS is None:
        # Generate a 2D WCS to register all frames
        coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
        coord_dlt = [dspat, dspat, np.max(all_wave) - np.min(all_wave)]
        whitelightWCS = generate_WCS(coord_min, coord_dlt)

        # Generate coordinates
        cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)
        numra = 1+int((np.max(all_ra) - np.min(all_ra)) * cosdec / dspat)
        numdec = 1+int((np.max(all_dec) - np.min(all_dec)) / dspat)
    else:
        # If a WCS is supplied, the numra and numdec must be specified
        if (numra is None) or (numdec is None):
            msgs.error("A WCS has been supplied to make_whitelight." + msgs.newline() +
                       "numra and numdec must also be specified")
    xbins = np.arange(1 + numra) - 1
    ybins = np.arange(1 + numdec) - 1
    spec_bins = np.arange(2) - 1
    bins = (xbins, ybins, spec_bins)

    whitelight_Imgs = np.zeros((numra, numdec, numfiles))
    whitelight_ivar = np.zeros((numra, numdec, numfiles))
    for ff in range(numfiles):
        msgs.info("Generating white light image of frame {0:d}/{1:d}".format(ff + 1, numfiles))
        ww = (all_idx == ff)
        # Make the cube
        pix_coord = whitelightWCS.wcs_world2pix(np.vstack((all_ra[ww], all_dec[ww], all_wave[ww] * 1.0E-10)).T, 0)
        wlcube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci[ww] * all_wghts[ww])
        norm, edges = np.histogramdd(pix_coord, bins=bins, weights=all_wghts[ww])
        nrmCube = (norm > 0) / (norm + (norm == 0))
        whtlght = (wlcube * nrmCube)[:, :, 0]
        # Create a mask of good pixels (trim the edges)
        gpm = grow_mask(whtlght == 0, trim) == 0  # A good pixel = 1
        whtlght *= gpm
        # Set the masked regions to the minimum value
        minval = np.min(whtlght[gpm == 1])
        whtlght[gpm == 0] = minval
        # Store the white light image
        whitelight_Imgs[:, :, ff] = whtlght.copy()
        # Now operate on the inverse variance image
        if all_ivar is not None:
            ivar_img, _ = np.histogramdd(pix_coord, bins=bins, weights=all_ivar[ww])
            ivar_img = ivar_img[:, :, 0]
            ivar_img *= gpm
            minval = np.min(ivar_img[gpm == 1])
            ivar_img[gpm == 0] = minval
            whitelight_ivar[:, :, ff] = ivar_img.copy()
    return whitelight_Imgs, whitelight_ivar, whitelightWCS


def make_sensfunc(ss_file, senspar, blaze_wave=None, blaze_spline=None, grating_corr=False):
    """
    Generate the sensitivity function from a standard star DataCube.

    Args:
        ss_file (:obj:`str`):
            The relative path and filename of the standard star datacube. It
            should be fits format, and for full functionality, should ideally of
            the form :class:`~pypeit.coadd3d.DataCube`.
        senspar (:class:`~pypeit.par.pypeitpar.SensFuncPar`):
            The parameters required for the sensitivity function computation.
        blaze_wave (`numpy.ndarray`_, optional):
            Wavelength array used to construct blaze_spline
        blaze_spline (`scipy.interpolate.interp1d`_, optional):
            Spline representation of the reference blaze function (based on the illumflat).
        grating_corr (:obj:`bool`, optional):
            If a grating correction should be performed, set this variable to True.

    Returns:
        `numpy.ndarray`_: A mask of the good sky pixels (True = good)
    """
    # TODO :: This routine has not been updated to the new spec1d plan of passing in a sensfunc object
    #      :: Probably, this routine should be removed and the functionality moved to the sensfunc object
    msgs.error("coding error - make_sensfunc is not currently supported.  Please contact the developers")
    # Check if the standard star datacube exists
    if not os.path.exists(ss_file):
        msgs.error("Standard cube does not exist:" + msgs.newline() + ss_file)
    msgs.info(f"Loading standard star cube: {ss_file:s}")
    # Load the standard star cube and retrieve its RA + DEC
    stdcube = fits.open(ss_file)
    star_ra, star_dec = stdcube[1].header['CRVAL1'], stdcube[1].header['CRVAL2']

    # Extract a spectrum of the standard star
    wave, Nlam_star, Nlam_ivar_star, gpm_star = extract_standard_spec(stdcube)

    # Extract the information about the blaze
    if grating_corr:
        blaze_wave_curr, blaze_spec_curr = stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data
        blaze_spline_curr = interp1d(blaze_wave_curr, blaze_spec_curr,
                                     kind='linear', bounds_error=False, fill_value="extrapolate")
        # Perform a grating correction
        grat_corr = correct_grating_shift(wave, blaze_wave_curr, blaze_spline_curr, blaze_wave, blaze_spline)
        # Apply the grating correction to the standard star spectrum
        Nlam_star /= grat_corr
        Nlam_ivar_star *= grat_corr ** 2

    # Read in some information above the standard star
    std_dict = flux_calib.get_standard_spectrum(star_type=senspar['star_type'],
                                                star_mag=senspar['star_mag'],
                                                ra=star_ra, dec=star_dec)
    # Calculate the sensitivity curve
    # TODO :: This needs to be addressed... unify flux calibration into the main PypeIt routines.
    msgs.warn("Datacubes are currently flux-calibrated using the UVIS algorithm... this will be deprecated soon")
    zeropoint_data, zeropoint_data_gpm, zeropoint_fit, zeropoint_fit_gpm = \
        flux_calib.fit_zeropoint(wave, Nlam_star, Nlam_ivar_star, gpm_star, std_dict,
                                 mask_hydrogen_lines=senspar['mask_hydrogen_lines'],
                                 mask_helium_lines=senspar['mask_helium_lines'],
                                 hydrogen_mask_wid=senspar['hydrogen_mask_wid'],
                                 nresln=senspar['UVIS']['nresln'],
                                 resolution=senspar['UVIS']['resolution'],
                                 trans_thresh=senspar['UVIS']['trans_thresh'],
                                 polyorder=senspar['polyorder'],
                                 polycorrect=senspar['UVIS']['polycorrect'],
                                 polyfunc=senspar['UVIS']['polyfunc'])
    wgd = np.where(zeropoint_fit_gpm)
    sens = np.power(10.0, -0.4 * (zeropoint_fit[wgd] - flux_calib.ZP_UNIT_CONST)) / np.square(wave[wgd])
    return interp1d(wave[wgd], sens, kind='linear', bounds_error=False, fill_value="extrapolate")
