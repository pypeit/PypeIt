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
from pypeit import utils
from pypeit import __version__
from pypeit.pypmsgs import PypeItError, PypeItDataModelError

from pypeit.display import display
from pypeit.images.imagebitmask import ImageBitMask
from pypeit.images.detector_container import DetectorContainer
from pypeit import spec2dobj
from pypeit.scripts import scriptbase


def show_trace(sobjs, det, viewer, ch):
    """
    Overplot the extracted object traces for this detector in the provided ginga
    channel.

    Args:
        sobjs (:class:`~pypeit.specobjs.SpecObjs`):
            Object holding the 1D spectral extractions.  If None, the function
            doesn't do anything.
        det (:obj:`str`):
            The string identifier for the detector or mosaic used to select the
            extractions to show.
        viewer (?):
        ch (?):
    """
    if sobjs is None:
        return
    in_det = np.where(sobjs.DET == det)[0]
    trace_list = []
    trc_name_list = []
    maskdef_extr_list = []
    manual_extr_list = []
    for kk in in_det:
        trace = sobjs[kk]['TRACE_SPAT']
        obj_id = sobjs[kk].NAME
        maskdef_objname = sobjs[kk].MASKDEF_OBJNAME
        maskdef_extr_flag = sobjs[kk].MASKDEF_EXTRACT
        manual_extr_flag = sobjs[kk].hand_extract_flag
        trc_name = '{}     OBJNAME:{}'.format(obj_id, maskdef_objname) if maskdef_objname is not None else obj_id

        trace_list.append(trace)
        trc_name_list.append(trc_name)
        maskdef_extr_list.append(maskdef_extr_flag is True)
        manual_extr_list.append(manual_extr_flag is True)

    display.show_trace(viewer, ch, np.swapaxes(trace_list, 1,0), np.array(trc_name_list),
                       maskdef_extr=np.array(maskdef_extr_list), manual_extr=np.array(manual_extr_list))


class Show2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display sky subtracted, spec2d image in a '
                                                'ginga viewer.',
                                    width=width)

        parser.add_argument('file', type=str, default=None, help='Path to a PypeIt spec2d file')
        parser.add_argument('--list', default=False, action='store_true',
                            help='List the extensions only?')
        parser.add_argument('--det', default='1', type=str,
                            help='Detector name or number.  If a number, the name is constructed '
                                 'assuming the reduction is for a single detector.  If a string, '
                                 'it must match the name of the detector object (e.g., DET01 for '
                                 'a detector, MSC01 for a mosaic).')
        parser.add_argument('--spat_id', type=int, default=None,
                            help='Restrict plotting to this slit (PypeIt ID notation)')
        parser.add_argument('--maskID', type=int, default=None,
                            help='Restrict plotting to this maskID')
        parser.add_argument('--showmask', default=False, help='Overplot masked pixels',
                            action='store_true')
        parser.add_argument('--removetrace', default=False, action='store_true',
                            help='Do not overplot traces in the skysub, sky_resid, and resid '
                                 'channels')
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        parser.add_argument('--ignore_extract_mask', default=False, action='store_true',
                            help='Ignore the extraction mask')
#        parser.add_argument('--sensfunc', type=str, default=None,
#                            help='Pass in a sensfunc to display the sky-subtracted image with a '
#                                 'flux calibration')
        parser.add_argument('--channels', type=str,
                            help='Only show a subset of the channels (0-indexed), e.g. 1,3')
        parser.add_argument('--prefix', type=str, default='',
                            help='Channel name prefix [lets you display more than one set]')
        parser.add_argument('--no_clear', dest='clear', default=True, 
                            action='store_false',
                            help='Do *not* clear all existing tabs')
        parser.add_argument('-v', '--verbosity', type=int, default=2,
                            help='Verbosity level between 0 [none] and 2 [all]')
        return parser

    @staticmethod
    def main(args):

        # List only?
        if args.list:
            io.fits_open(args.file).info()
            return

        # Setup for PypeIt imports
        msgs.reset(verbosity=args.verbosity)

        # Parse the detector name
        try:
            det = int(args.det)
        except:
            detname = args.det
        else:
            detname = DetectorContainer.get_name(det)

        # Find the set of channels to show
        show_channels = [0,1,2,3] if args.channels is None \
                            else [int(item) for item in args.channels.split(',')]

        # Try to read the Spec2DObj using the current datamodel, but allowing
        # for the datamodel version to be different
        try:
            spec2DObj = spec2dobj.Spec2DObj.from_file(args.file, detname, chk_version=False)
        except PypeItDataModelError:
            try:
                # Try to get the pypeit version used to write this file
                file_pypeit_version = fits.getval(args.file, 'VERSPYP', 0)
            except KeyError:
                file_pypeit_version = '*unknown*'
            msgs.warn(f'Your installed version of PypeIt ({__version__}) cannot be used to parse '
                      f'{args.file}, which was reduced using version {file_pypeit_version}.  You '
                      'are strongly encouraged to re-reduce your data using this (or, better yet, '
                      'the most recent) version of PypeIt.  Script will try to parse only the '
                      'relevant bits from the spec2d file and continue (possibly with more '
                      'limited functionality).')
            spec2DObj = None

        if spec2DObj is None:
            # Try to get the relevant elements directly from the fits file
            with io.fits_open(args.file) as hdu:
                names = [h.name for h in hdu]
                has_det = any([detname in n for n in names])
                if not has_det:
                    msgs.error(f'Provided file has no extensions including {detname}.')
                for ext in ['SCIIMG', 'SKYMODEL', 'OBJMODEL', 'IVARMODEL']:
                    _ext = f'{detname}-{ext}'
                    if _ext not in names:
                        msgs.error(f'{args.file} missing extension {_ext}.')

                sciimg = hdu[f'{detname}-SCIIMG'].data
                skymodel = hdu[f'{detname}-SKYMODEL'].data
                objmodel = hdu[f'{detname}-OBJMODEL'].data
                ivarmodel = hdu[f'{detname}-IVARMODEL'].data
                _ext = f'{detname}-BPMMASK'
                bpmmask = hdu[_ext].data if _ext in names else None
                _ext = f'{detname}-WAVEIMG'
                waveimg = hdu[_ext].data if _ext in names else None

                img_gpm = np.ones(sciimg.shape, dtype=bool) if bpmmask is None else bpmmask == 0
                model_gpm = img_gpm.copy()

                _ext = f'{detname}-SLITS'
                if _ext not in names:
                    msgs.warn(f'{args.file} missing extension {_ext}; cannot show slit edges.')
                else:
                    slit_columns = hdu[_ext].columns.names
                    slit_spat_id = hdu[_ext].data['spat_id'] if 'spat_id' in slit_columns else None
                    slit_left = hdu[_ext].data['left_tweak'].T if 'left_tweak' in slit_columns \
                                else (hdu[_ext].data['left_init'].T if 'left_init' in slit_columns
                                      else None)
                    slit_right = hdu[_ext].data['right_tweak'].T if 'right_tweak' in slit_columns \
                                else (hdu[_ext].data['right_init'].T if 'right_init' in slit_columns
                                      else None)
                    slit_mask = hdu[_ext].data['mask'] if 'mask' in slit_columns else None
                    slit_mask_id = hdu[_ext].data['maskdef_id'] \
                                        if 'maskdef_id' in slit_columns else None
                    if 'SCI_SPAT_FLEXURE' in hdu[f'{detname}-SCIIMG'].header \
                            and slit_left is not None and slit_right is not None:
                        sci_spat_flexure \
                                = float(hdu[f'{detname}-SCIIMG'].header['SCI_SPAT_FLEXURE'])
                        slit_left += sci_spat_flexure
                        slit_right += sci_spat_flexure
                        msgs.info(f'Offseting slits by {sci_spat_flexure} pixels.')
                    pypeline = hdu[f'{detname}-SCIIMG'].header['PYPELINE'] \
                                    if 'PYPELINE' in hdu[f'{detname}-SCIIMG'].header else None
                    if pypeline in ['MultiSlit', 'IFU']:
                        slit_slid_IDs = slit_spat_id
                    elif pypeline == 'Echelle':
                        slit_slid_IDs = hdu[_ext].data['ech_order'] \
                                            if 'ech_order' in slit_columns else None
                    else:
                        slit_slid_IDs = None

                # Ignore any spec1d file
                sobjs = None

        else:
            # Use the parsed SpecObjs object
            sciimg = spec2DObj.sciimg
            skymodel = spec2DObj.skymodel
            objmodel = spec2DObj.objmodel
            ivarmodel = spec2DObj.ivarmodel
            bpmmask = spec2DObj.bpmmask
            waveimg = spec2DObj.waveimg

            img_gpm = spec2DObj.select_flag(invert=True)
            model_gpm = img_gpm.copy()
            if args.ignore_extract_mask:
                model_gpm |= spec2DObj.select_flag(flag='EXTRACT')

            if spec2DObj.sci_spat_flexure is not None:
                msgs.info(f'Offseting slits by {spec2DObj.sci_spat_flexure}')
            slit_left, slit_right, slit_mask \
                    = spec2DObj.slits.select_edges(flexure=spec2DObj.sci_spat_flexure)
            slit_spat_id = spec2DObj.slits.spat_id
            slit_mask_id = spec2DObj.slits.maskdef_id
            slit_slid_IDs = spec2DObj.slits.slitord_id

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

        # TODO: This may be too restrictive, i.e. ignore BADFLTCALIB??
        slit_gpm = slit_mask == 0

        # Restrict on spat_id or maskid? 
        if slit_gpm is not None:
            if slit_spat_id is not None and args.spat_id is not None:
                slit_gpm &= slit_spat_id == args.spat_id
            # TODO: Should this be 'elif'?  Why not just 'if'?  I.e., are
            # spat_id and maskdef_id mutually exclusive?
            elif slit_mask_id is not None and args.maskID is not None:
                slit_gpm &= slit_mask_id == args.maskID

        left = None if slit_left is None else slit_left[:, slit_gpm]
        right = None if slit_right is None else slit_right[:, slit_gpm]
        slid_IDs = None if slit_slid_IDs is None else slit_slid_IDs[slit_gpm]
        maskdef_id = None if slit_mask_id is None else slit_mask_id[slit_gpm]

        # Connect to an open ginga window, or open a new one
        display.connect_to_ginga(raise_err=True, allow_new=True)

        # Show the bitmask?
        if args.showmask and bpmmask is not None:
            viewer, ch_mask = display.show_image(bpmmask, chname='MASK', waveimg=waveimg,
                                                 clear=args.clear)

        channel_names = []
        # SCIIMG
        if 0 in show_channels:
            mean, med, sigma = sigma_clipped_stats(sciimg[img_gpm], sigma_lower=5.0, sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_sci = args.prefix+f'sciimg-{detname}'
            # Clear all channels at the beginning
            viewer, ch_sci = display.show_image(sciimg, chname=chname_sci, waveimg=waveimg, 
                                                clear=args.clear, cuts=(cut_min, cut_max))
            if sobjs is not None:
                show_trace(sobjs, detname, viewer, ch_sci)
            display.show_slits(viewer, ch_sci, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_sci)

        # SKYSUB
        if 1 in show_channels:
            image = (sciimg - skymodel) * model_gpm.astype(float)
            mean, med, sigma = sigma_clipped_stats(image[model_gpm], sigma_lower=5.0,
                                                   sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_skysub = args.prefix+f'skysub-{detname}'
            viewer, ch_skysub = display.show_image(image, chname=chname_skysub,
                                                   waveimg=waveimg, cuts=(cut_min, cut_max),
                                                   wcs_match=True)
            if not args.removetrace and sobjs is not None:
                show_trace(sobjs, detname, viewer, ch_skysub)
            display.show_slits(viewer, ch_skysub, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_skysub)

        # TODO Place holder for putting in sensfunc
        #if args.sensfunc:
        #    # Load the sensitivity function
        #    wave_sens, sfunc, _, _, _ = sensfunc.SensFunc.load(sensfunc_name)
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
            chname_skyresids = args.prefix+f'sky_resid-{detname}'
            image = (sciimg - skymodel) * np.sqrt(ivarmodel) * model_gpm.astype(float)
            viewer, ch_sky_resids = display.show_image(image, chname_skyresids, waveimg=waveimg,
                                                       cuts=(-5.0, 5.0))
            if not args.removetrace and sobjs is not None:
                show_trace(sobjs, detname, viewer, ch_sky_resids)
            display.show_slits(viewer, ch_sky_resids, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_skyresids)

        # RESIDS
        if 3 in show_channels:
            chname_resids = args.prefix+f'resid-{detname}'
            # full model residual map
            image = (sciimg - skymodel - objmodel) * np.sqrt(ivarmodel) * img_gpm.astype(float)
            viewer, ch_resids = display.show_image(image, chname=chname_resids, waveimg=waveimg,
                                                   cuts=(-5.0, 5.0), wcs_match=True)
            if not args.removetrace and sobjs is not None:
                show_trace(sobjs, detname, viewer, ch_resids)
            display.show_slits(viewer, ch_resids, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_resids)


        # After displaying all the images sync up the images with WCS_MATCH
        shell = viewer.shell()
        shell.start_global_plugin('WCSMatch')
        shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [channel_names[-1]],
                                        {})

        if args.embed:
            embed(header=utils.embed_header())


        # Playing with some mask stuff
        #out = shell.start_operation('TVMask')
        #maskfile = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Shane_Kast_blue/600_4310_d55/shane_kast_blue_setup_A/crmask.fits'
        #out = shell.call_local_plugin_method(chname_resids, 'TVMask', 'load_file', [maskfile], {})


