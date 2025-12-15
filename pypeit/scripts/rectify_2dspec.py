"""
Create an FITS file with rectified 2D spectra for all slits/orders.

.. include common links, assuming the primary doc root is up one directory
.. include:: ../include/links.rst
"""

import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt

from pypeit.scripts import scriptbase
from pypeit import spec2dobj, specobjs, msgs
from pypeit.core import coadd
from pypeit.core.wavecal import wvutils
from pypeit.core.moment import moment1d
from astropy.io import fits

from IPython import embed


class Rectify2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Create an FITS file with rectified '
                                                '2D sky-subtracted spectra for all slits/orders.', width=width)

        parser.add_argument('files', type = str, nargs='*', help = 'PypeIt spec2d file(s)')
        parser.add_argument('--no_rot', default=False, action='store_true',
                            help='Do not rotate the rectified image to have wavelength '
                                 'on the x-axis.')
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Embed in IPython shell in each detector loop, i.e., before saving to disk.')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        chk_version = not args.try_old

        for spec2file in args.files:
            msgs.info(f'Processing file: {spec2file}')
            # Get list of detectors
            hdr = fits.getheader(spec2file)
            detnames = hdr['HIERARCH ALLSPEC2D_DETS'].split(',')

            # Empty primary HDU
            hdu_list = [fits.PrimaryHDU()]

            for detname in detnames:
                msgs.info(f'DETECTOR: {detname}')
                spec2d = spec2dobj.Spec2DObj.from_file(spec2file, detname, chk_version=chk_version)
                pad = 10  # pixels to pad on each side
                slitmask = spec2d.slits.slit_img(pad=pad, flexure=spec2d.sci_spat_flexure)
                slit_ids = spec2d.slits.spat_id
                # this is just to print useful info in the terminal
                slitord_ids = spec2d.slits.slitord_id

                # get the wave grid
                # Do we have the spec1d file?
                spec1dfile = spec2file.replace('spec2d', 'spec1d')
                if Path(spec1dfile).is_file():
                    sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)
                    on_det = sobjs.DET == detname
                    waves, gpms = [], []
                    for sobj in sobjs[on_det]:
                        # we use the BOX extraction here for simplicity
                        if sobj.has_box_ext() and np.any(sobj.BOX_MASK):
                            waves.append(sobj.BOX_WAVE)
                            gpms.append(sobj.BOX_MASK)
                else:
                    # extract the wave and gpm from the spec2d object directly (taken from coadd2d)
                    waves, gpms = [], []
                    slits_left, slits_right, _ = spec2d.slits.select_edges()
                    row = np.arange(slits_left.shape[0])
                    box_radius = 3 # pixels
                    # Loop on the slits
                    for kk, spat_id in enumerate(slit_ids):
                        mask = slitmask == spat_id
                        # Create apertures at 5%, 50%, and 95% of the slit width to cover full range of wavelengths
                        # on this slit
                        trace_spat = slits_left[:, kk][:, np.newaxis] + np.outer((slits_right[:, kk] - slits_left[:, kk]),
                                                                                 [0.05, 0.5, 0.95])
                        box_denom = moment1d(spec2d.waveimg * mask > 0.0, trace_spat, 2 * box_radius, row=row)[0]
                        wave_box = moment1d(spec2d.waveimg * mask, trace_spat, 2 * box_radius,
                                            row=row)[0] / (box_denom + (box_denom == 0.0))
                        gpm_box = box_denom > 0.
                        waves += [wave for (wave, gpm) in zip(wave_box.T, gpm_box.T) if np.any(gpm)]
                        gpms += [(wave > 0.) & gpm for (wave, gpm) in zip(wave_box.T, gpm_box.T) if np.any(gpm)]


                if len(waves) == 0:
                    msgs.warn(f'There is a problem with the wavelengths on det {detname}. '
                              f'The RECTIFIED 2D spectral image will not be created.')
                    continue
                wmax = np.ceil(spec2d.waveimg[spec2d.waveimg>0].max())
                wmin = np.floor(spec2d.waveimg[spec2d.waveimg>0].min())
                wave_grid, wave_grid_mid, dsamp = wvutils.get_wave_grid(waves=waves, gpms=gpms,
                                                                        wave_grid_min=wmin, wave_grid_max=wmax)

                # Loop over slits and rectify each one
                imgrect_list = []
                for slitidx, slit_id in enumerate(slit_ids):
                    this_mask = (slitmask == slit_id)
                    # check if this slit was masked. If so, skip it.
                    slitord_id = slitord_ids[slitidx]
                    if not np.any(this_mask):
                        msgs.warn(f'Slit/order {slitord_id} on {detname} is fully masked. Skipping it.')
                        continue
                    msgs.info(f'Rectifying slit/order {slitord_id}')

                    slit_cen = spec2d.slits.center[:,slitidx]
                    mask = spec2d.bpmmask.mask == 0

                    imgrect_dict = coadd.compute_coadd2d([slit_cen], [spec2d.sciimg],
                                                         [spec2d.ivarmodel],[spec2d.skymodel],
                                                         [mask],[this_mask],
                                                         [spec2d.waveimg], wave_grid)
                    imgrect_list.append(imgrect_dict)


                # Get dimensions for each slit
                nspat_vec = np.array([cdict['nspat'] for cdict in imgrect_list])
                nspec_rect = len(wave_grid_mid)
                nspat_rect = int(np.sum(nspat_vec) + (spec2d.slits.nslits + 1) * pad)

                # Initialize output array
                # this is the rectified image on a common wavelength grid for all slits
                image_rect = np.zeros((nspec_rect, nspat_rect))
                # the outputs below are not really used, but are useful for double-checking
                # this is the wave image created by placing the rectified wavelengths for each slit in
                # their proper location in the rectified image
                waveimg_rect = np.zeros((nspec_rect, nspat_rect))
                # this is the common wave image for the whole rectified image, created by repeating
                # the wavelength grid for each spatial pixel
                wavegrid_image = np.repeat(wave_grid_mid[:, np.newaxis], nspat_rect, axis=1)

                # Loop over the rectified slits and place them in the output arrays
                spat_left = pad
                for islit, imgrect_dict in enumerate(imgrect_list):
                    spat_righ = spat_left + nspat_vec[islit]

                    # Get wavelength array for this slit
                    wave_slit = imgrect_dict['wave_mid']
                    nspec_slit = len(wave_slit)
                    # get the spectral pixel where the slit should start/stop
                    _wave_grid_mid = np.round(wave_grid_mid, 4)
                    _wave_slit = np.round(imgrect_dict['wave_mid'], 4)
                    wstart = np.where(_wave_grid_mid == _wave_slit[0])[0][0]
                    wend = np.where(_wave_grid_mid == _wave_slit[-1])[0][0] + 1

                    # Place the rectified data for this slit in the output arrays
                    for ispat_slit in range(nspat_vec[islit]):
                        # Get data for this spatial pixel
                        flux = imgrect_dict['imgminsky'][:nspec_slit, ispat_slit]

                        # Assign to output arrays
                        image_rect[wstart:wend, spat_left + ispat_slit] = flux
                        waveimg_rect[wstart:wend, spat_left + ispat_slit] = wave_slit

                    spat_left = spat_righ + pad

                if args.embed:
                    embed(header=f">>> Rectify2DSpec: useful variables are image_rect, waveimg_rect, wavegrid_image")

                # Rotate if desired
                if not args.no_rot:
                    image_rect = np.rot90(image_rect)
                    # Create HDU with WCS header
                    hdu = fits.ImageHDU(image_rect, name=detname)
                    hdu.header['CTYPE1'] = 'LAMBDA   '
                    hdu.header['CUNIT1'] = 'Angstrom'
                    hdu.header['CDELT1'] = dsamp
                    hdu.header['CRPIX1'] = 0
                    hdu.header['CRVAL1'] = wave_grid_mid[0]
                    hdu.header['CTYPE2'] = 'LINEAR  '
                    hdu.header['CUNIT2'] = 'pixel'
                    hdu.header['CDELT2'] = 1.
                    hdu.header['CRPIX2'] = 0
                    hdu.header['CRVAL2'] = 0.
                else:
                    # Create HDU with WCS header
                    hdu = fits.ImageHDU(image_rect, name=detname)
                    hdu.header['CTYPE1'] = 'LINEAR  '
                    hdu.header['CUNIT1'] = 'pixel'
                    hdu.header['CDELT1'] = 1.
                    hdu.header['CRPIX1'] = 0
                    hdu.header['CRVAL1'] = 0.
                    hdu.header['CTYPE2'] = 'LAMBDA   '
                    hdu.header['CUNIT2'] = 'Angstrom'
                    hdu.header['CDELT2'] = dsamp
                    hdu.header['CRPIX2'] = 0
                    hdu.header['CRVAL2'] = wave_grid_mid[0]
                # Add other keywords
                hdu.header['NSPEC'] = nspec_rect
                hdu.header['NSPAT'] = nspat_rect
                hdu.header['WAVEMIN'] = wave_grid_mid[0]
                hdu.header['WAVEMAX'] = wave_grid_mid[-1]
                hdu.header['DETNAME'] = detname
                hdu.header['PIPELINE'] = hdr['PIPELINE']
                hdu.header['PYPELINE'] = hdr['PYPELINE']
                hdu.header['PYP_SPEC'] = hdr['PYP_SPEC']

                hdu_list.append(hdu)

            # Write to file
            out_file = spec2file.replace('spec2d', 'rectified_spec2d')
            hdulist = fits.HDUList(hdu_list)
            hdulist.writeto(out_file, overwrite=True)
            msgs.info(f'Rectified images saved to {out_file}')