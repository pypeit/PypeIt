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
        frame_wcs (`astropy.wcs.wcs.WCS`_):
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
        output_wcs (`astropy.wcs.wcs.WCS`_, optional):
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
