"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

import numpy as np

from IPython import embed

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import msgs
from pypeit import slittrace
from pypeit import specobjs
from pypeit import io

from pypeit.display import display
from pypeit.core.parse import get_dnum
from pypeit.images.imagebitmask import ImageBitMask
from pypeit import masterframe
from pypeit import spec2dobj
from pypeit.scripts import scriptbase


# TODO: Document this function...
def show_trace(specobjs, det, viewer, ch):

    if specobjs is None:
        return
    in_det = np.where(specobjs.DET == det)[0]
    for kk in in_det:
        trace = specobjs[kk]['TRACE_SPAT']
        obj_id = specobjs[kk].NAME
        maskdef_objname = specobjs[kk].MASKDEF_OBJNAME
        maskdef_extr_flag = specobjs[kk].MASKDEF_EXTRACT
        manual_extr_flag = specobjs[kk].hand_extract_flag
        if maskdef_objname is not None:
            trc_name = '{}     OBJNAME:{}'.format(obj_id, maskdef_objname)
        else:
            trc_name = obj_id
        if maskdef_extr_flag is not None and maskdef_extr_flag is True:
            display.show_trace(viewer, ch, trace, trc_name, color='gold') #hdu.name)
        elif manual_extr_flag is True:
            display.show_trace(viewer, ch, trace, trc_name, color='#33ccff') #hdu.name)
        else:
            display.show_trace(viewer, ch, trace, trc_name, color='orange') #hdu.name)


class Show2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display sky subtracted, spec2d image in a '
                                                 'Ginga viewer.  Run above the Science/ folder',
                                    width=width)

        parser.add_argument('file', type = str, default = None, help = 'PYPIT spec2d file')
        parser.add_argument('--list', default=False, help='List the extensions only?',
                            action='store_true')
        parser.add_argument('--det', default=1, type=int, help='Detector number')
        parser.add_argument('--showmask', default=False, help='Overplot masked pixels',
                            action='store_true')
        parser.add_argument('--removetrace', default=False, action="store_true",
                            help='Do not overplot traces in the skysub, sky_resid, and resid '
                                 'channels')
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        parser.add_argument('--ignore_extract_mask', default=False, action='store_true',
                            help='Ignore the extraction mask')
        parser.add_argument("--sensfunc", type=str, default=None,
                            help='Pass in a sensfunc to display the sky-subtracted image with a '
                                 'flux calibration')
        parser.add_argument('--channels', type=str,
                            help='Only show a subset of the channels (0-indexed), e.g. 1,3')
        return parser

    @staticmethod
    def main(args):

        # List only?
        if args.list:
            hdu = io.fits_open(args.file)
            hdu.info()
            return

        # Load it up -- NOTE WE ALLOW *OLD* VERSIONS TO GO FORTH
        spec2DObj = spec2dobj.Spec2DObj.from_file(args.file, args.det, chk_version=False)

        # Setup for PypeIt imports
        msgs.reset(verbosity=2)

        # Init
        # TODO: get_dnum needs to be deprecated...
        sdet = get_dnum(args.det, prefix=False)

        # Find the set of channels to show
        if args.channels is not None:
            show_channels = [int(item) for item in args.channels.split(',')]
        else:
            show_channels = [0,1,2,3]

        # Grab the slit edges
        slits = spec2DObj.slits
        if spec2DObj.sci_spat_flexure is not None:
            msgs.info("Offseting slits by {}".format(spec2DObj.sci_spat_flexure))
        all_left, all_right, mask = slits.select_edges(flexure=spec2DObj.sci_spat_flexure)
        # TODO -- This may be too restrictive, i.e. ignore BADFLTCALIB??
        gpm = mask == 0
        left = all_left[:, gpm]
        right = all_right[:, gpm]
        slid_IDs = spec2DObj.slits.slitord_id[gpm]
        maskdef_id = None if spec2DObj.slits.maskdef_id is None \
                            else spec2DObj.slits.maskdef_id[gpm]
        bitMask = ImageBitMask()

        # Object traces from spec1d file
        spec1d_file = args.file.replace('spec2d', 'spec1d')
        if args.file[-2:] == 'gz':
            spec1d_file = spec1d_file[:-3]
        if os.path.isfile(spec1d_file):
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=False)
        else:
            sobjs = None
            msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                      '                          No objects were extracted.')

        display.connect_to_ginga(raise_err=True, allow_new=True)

        # Now show each image to a separate channel

        # Show the bitmask?
        mask_in = None
        if args.showmask:
            viewer, ch_mask = display.show_image(spec2DObj.bpmmask, chname="BPM",
                                                 waveimg=spec2DObj.waveimg, clear=True)

        channel_names = []
        # SCIIMG
        if 0 in show_channels:
            image = spec2DObj.sciimg  # Processed science image
            mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
                                                   sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_sci = 'sciimg-det{:s}'.format(sdet)
            # Clear all channels at the beginning
            viewer, ch_sci = display.show_image(image, chname=chname_sci,
                                                waveimg=spec2DObj.waveimg, clear=True,
                                                cuts=(cut_min, cut_max))

            if sobjs is not None:
                show_trace(sobjs, args.det, viewer, ch_sci)
            display.show_slits(viewer, ch_sci, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_sci)

        # SKYSUB
        if 1 in show_channels:
            if args.ignore_extract_mask:
                # TODO -- Is there a cleaner way to do this?
                gpm = (spec2DObj.bpmmask == 0) | (spec2DObj.bpmmask == 2**bitMask.bits['EXTRACT'])
            else:
                gpm = spec2DObj.bpmmask == 0

            image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm
            mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
                                                   sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_skysub = 'skysub-det{:s}'.format(sdet)
            viewer, ch_skysub = display.show_image(image, chname=chname_skysub,
                                                   waveimg=spec2DObj.waveimg, bitmask=bitMask,
                                                   mask=mask_in, cuts=(cut_min, cut_max),
                                                   wcs_match=True)
            if not args.removetrace and sobjs is not None:
                    show_trace(sobjs, args.det, viewer, ch_skysub)
            display.show_slits(viewer, ch_skysub, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_skysub)

        # TODO Place holder for putting in sensfunc
        #if args.sensfunc:
        #    # Load the sensitivity function
        #    wave_sens, sfunc, _, _, _ = sensfunc.SensFunc.load(sensfunc_masterframe_name)
        #    # Interpolate the sensitivity function onto the wavelength grid of the data. Since the image is rectified
        #    # this is trivial and we don't need to do a 2d interpolation
        #    sens_factor = flux_calib.get_sensfunc_factor(
        #        pseudo_dict['wave_mid'][:,islit], wave_sens, sfunc, fits.getheader(files[0])['TRUITIME'],
        #        extrap_sens=parset['fluxcalib']['extrap_sens'])
        #    # Compute the median sensitivity and set the sensitivity to zero at locations 100 times the median. This
        #    # prevents the 2d image from blowing up where the sens_factor explodes because there is no throughput
        #    sens_gpm = sens_factor < 100.0*np.median(sens_factor)
        #    sens_factor_masked = sens_factor*sens_gpm
        #    sens_factor_img = np.repeat(sens_factor_masked[:, np.newaxis], pseudo_dict['nspat'], axis=1)
        #    imgminsky = sens_factor_img*pseudo_dict['imgminsky']
        #    imgminsky_gpm = sens_gpm[:, np.newaxis] & pseudo_dict['inmask']
        #else:
        #    imgminsky= pseudo_dict['imgminsky']

        # SKRESIDS
        if 2 in show_channels:
            # the block below is repeated because if showing this channel but
            # not channel 1 it will crash
            if args.ignore_extract_mask:
                # TODO -- Is there a cleaner way to do this?
                gpm = (spec2DObj.bpmmask == 0) | (spec2DObj.bpmmask == 2**bitMask.bits['EXTRACT'])
            else:
                gpm = spec2DObj.bpmmask == 0
            chname_skyresids = 'sky_resid-det{:s}'.format(sdet)
            image = (spec2DObj.sciimg - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) * gpm
            viewer, ch_sky_resids = display.show_image(image, chname_skyresids,
                                                       waveimg=spec2DObj.waveimg, cuts=(-5.0, 5.0),
                                                       bitmask=bitMask, mask=mask_in)
            if not args.removetrace and sobjs is not None:
                    show_trace(sobjs, args.det, viewer, ch_sky_resids)
            display.show_slits(viewer, ch_sky_resids, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_skyresids)

        # RESIDS
        if 3 in show_channels:
            chname_resids = 'resid-det{:s}'.format(sdet)
            # full model residual map
            image = (spec2DObj.sciimg - spec2DObj.skymodel - spec2DObj.objmodel) \
                        * np.sqrt(spec2DObj.ivarmodel) * (spec2DObj.bpmmask == 0)
            viewer, ch_resids = display.show_image(image, chname=chname_resids,
                                                   waveimg=spec2DObj.waveimg, cuts=(-5.0, 5.0),
                                                   bitmask=bitMask, mask=mask_in, wcs_match=True)
            if not args.removetrace and sobjs is not None:
                    show_trace(sobjs, args.det, viewer, ch_resids)
            display.show_slits(viewer, ch_resids, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_resids)


        # After displaying all the images sync up the images with WCS_MATCH
        shell = viewer.shell()
        shell.start_global_plugin('WCSMatch')
        shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [channel_names[-1]],
                                        {})

        if args.embed:
            embed()

        # Playing with some mask stuff
        #out = shell.start_operation('TVMask')
        #maskfile = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Shane_Kast_blue/600_4310_d55/shane_kast_blue_setup_A/crmask.fits'
        #out = shell.call_local_plugin_method(chname_resids, 'TVMask', 'load_file', [maskfile], {})


