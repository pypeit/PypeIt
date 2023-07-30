"""
Script for quick-look reductions for Multislit observations.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import copy
import os
import pathlib
import time

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

from pypeit import utils
from pypeit import data
from pypeit import par, msgs
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit import spec2dobj
from pypeit import coadd2d
from pypeit import specobjs
from pypeit import slittrace
from pypeit import extraction
from pypeit import find_objects
from pypeit import calibrations
from pypeit.display import display
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.parse import get_dnum, parse_binning
from pypeit import sensfunc
from pypeit.core import flux_calib
from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs

def parse_det(det, spectrograph):
    #TODO Need a utility routine to deal with this
    if 'mosaic' in det:
        mosaic = True
        _det = spectrograph.default_mosaic
        if _det is None:
            msgs.error(f'{spectrograph.name} does not have a known mosaic')
    else:
        try:
            _det = tuple(int(d) for d in det)
        except:
            msgs.error(f'Could not convert detector input to integer.')
        mosaic = len(_det) > 1
        if not mosaic:
            _det = _det[0]
    return _det

# TODO make the quicklook default_pypeit_par part of the spectrograph class
def config_lines(args):
    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format(args.spectrograph)]
    cfg_lines += ['    redux_path = {0}'.format(args.redux_path)]
    cfg_lines += ['    quicklook = True']
    cfg_lines += ['    scidir = Science_QL']
    # Calibrations
    cfg_lines += ['[baseprocess]']
    cfg_lines += ['    use_pixelflat = False']
    cfg_lines += ['    use_illumflat = False']
    cfg_lines += ['    use_biasimage = False'] # TODO determine how to make this instrument specific?

    cfg_lines += ['[calibrations]']
    cfg_lines += ['    [[wavelengths]]']
    cfg_lines += ['        refframe = observed']
    cfg_lines += ['[scienceframe]']
    cfg_lines += ['    [[process]]']
    cfg_lines += ['        spat_flexure_correct = False']
    if not args.mask_cr:
        cfg_lines += ['        mask_cr = False']
    cfg_lines += ['[reduce]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['        skip_optimal = True']
    if args.box_radius is not None:  # Boxcar radius
        cfg_lines += ['        boxcar_radius = {0}'.format(args.box_radius)]
    cfg_lines += ['    [[findobj]]']
    cfg_lines += ['        skip_second_find = True']
    cfg_lines += ['[flexure]']
    cfg_lines += ['        spec_method = skip']
    return cfg_lines

def print_offset_report(files, dither_pattern, dither_id, offset_arcsec, target, platescale):

    if len(np.unique(dither_pattern)) > 1:
        msgs.error('Script only supported for a single type of dither pattern.')
    A_files = files[dither_id == 'A']
    B_files = files[dither_id == 'B']
    nA = len(A_files)
    nB = len(B_files)

    # Print out a report on the offsets
    msg_string = msgs.newline() + '*******************************************************'
    msg_string += msgs.newline() + ' Summary of offsets for target {:s} with dither pattern:   {:s}'.format(target,
                                                                                                            dither_pattern[
                                                                                                                0])
    msg_string += msgs.newline() + '*******************************************************'
    msg_string += msgs.newline() + 'filename     Position         arcsec    pixels    '
    msg_string += msgs.newline() + '----------------------------------------------------'
    for iexp, file in enumerate(files):
        msg_string += msgs.newline() + '    {:s}    {:s}   {:6.2f}    {:6.2f}'.format(
            file.name, dither_id[iexp], offset_arcsec[iexp], offset_arcsec[iexp] / platescale)
    msg_string += msgs.newline() + '********************************************************'
    msgs.info(msg_string)

def build_calibrate(det, files, spectrograph, parset, bias_masterframe_name, slit_masterframe_name, wvcalib_masterframe_name,
                       tilts_masterframe_name):
    # Read in the bias
    msbias = buildimage.BiasImage.from_file(bias_masterframe_name) if bias_masterframe_name is not None else None
    # Read in the msbpm
    msbpm = spectrograph.bpm(files[0], det)
    # Read in the slits
    slits = slittrace.SlitTraceSet.from_file(slit_masterframe_name)
    # Reset the bitmask
    slits.mask = slits.mask_init.copy()
    # Read in the wv_calib
    wv_calib = wavecalib.WaveCalib.from_file(wvcalib_masterframe_name)
    # wv_calib.is_synced(slits)
    slits.mask_wvcalib(wv_calib)
    # Read in the tilts
    tilts_obj = wavetilts.WaveTilts.from_file(tilts_masterframe_name)
    tilts_obj.is_synced(slits)
    slits.mask_wavetilts(tilts_obj)

    # Build the Calibrate object
    caliBrate = calibrations.Calibrations(None, parset['calibrations'], spectrograph, None)
    caliBrate.msbias = msbias
    caliBrate.msbpm = msbpm
    caliBrate.slits = slits
    caliBrate.wavetilts = tilts_obj
    caliBrate.wv_calib = wv_calib
    caliBrate.binning = f'{slits.binspec},{slits.binspat}'
    caliBrate.det = det

    return caliBrate

def run(files, dither_id, offset_arcsec, caliBrate, spectrograph, platescale, parset, std_trace, show, bkg_redux=False):

    msg_string = msgs.newline() + '*******************************************************'
    spec2d_list = []
    offset_ref = offset_arcsec[0]
    offsets_dith_pix = []
    # Generalize to a multiple slits, doing one slit at a time?

    uniq_throws, uni_indx = np.unique(np.abs(offset_arcsec), return_inverse=True)
    # uniq_throws = uniq values of the dither throw
    # uni_indx = indices into the uniq_throws array needed to reconstruct the original array
    nuniq = uniq_throws.size
    for iuniq in range(nuniq):
        A_ind = (uni_indx == iuniq) & (dither_id == 'A')
        B_ind = (uni_indx == iuniq) & (dither_id == 'B')
        A_files_uni = files[A_ind]
        A_dither_id_uni = dither_id[A_ind]
        B_dither_id_uni = dither_id[B_ind]
        B_files_uni = files[B_ind]
        A_offset = offset_arcsec[A_ind]
        B_offset = offset_arcsec[B_ind]
        A_offset_pix = (np.mean(A_offset) - offset_ref) / platescale
        B_offset_pix = (np.mean(B_offset) - offset_ref) / platescale
        throw = np.abs(A_offset[0])
        if bkg_redux:
            msgs.info('Reducing A-B pairs for throw = {:}'.format(throw))
            if (len(A_files_uni) > 0) & (len(B_files_uni) > 0):
                spec2DObj_A, spec2DObj_B = reduce(A_files_uni, caliBrate, spectrograph, parset,
                                                  bkg_files=B_files_uni, show=show, std_trace=std_trace)
                spec2d_list += [spec2DObj_A, spec2DObj_B]
                offsets_dith_pix += [A_offset_pix, B_offset_pix]
            else:
                msgs.warn('Skpping files that do not have an A-B match with the same throw:')
                for iexp in range(len(A_files_uni)):
                    msg_string += msgs.newline() + '    {:s}    {:s}   {:6.2f}    {:6.2f}'.format(
                        A_files_uni[iexp].name, A_dither_id_uni[iexp],
                        A_offset[iexp], A_offset[iexp] / platescale)
                for iexp in range(len(B_files_uni)):
                    msg_string += msgs.newline() + '    {:s}    {:s}   {:6.2f}    {:6.2f}'.format(
                        B_files_uni[iexp].name, B_dither_id_uni[iexp],
                        B_offset[iexp], B_offset[iexp] / platescale)
        else:
            msgs.info('Reducing images for offset = {:}'.format(A_offset[0]))
            spec2DObj = reduce(A_files_uni, caliBrate, spectrograph, parset, show=show, std_trace=std_trace)
            spec2d_list += [spec2DObj]
            offsets_dith_pix += [A_offset_pix]

    offsets_dith_pix = np.array(offsets_dith_pix)

    return spec2d_list, offsets_dith_pix


def reduce(files, caliBrate, spectrograph, parset, bkg_files=None, show=False, std_trace=None):
    """
    Peform 2d extraction for a set of files at the same unique A-B offset location.

    Parameters
    ----------
    A_files (list of strings):
       Files at A position for this offset
    B_files (list of strings)
       Files at B position for this offeset
    caliBrate (object):
       CaliBrate object
    spectrograph (object):
       spectrograph object
    det (int):
       Detector number
    parset (parsect object)
       Parset
    show (bool, optional):
       Show 2d reduction outputs. Default=False
    std_trace (string, optional)
       Trace for standard star. Default=None

    Returns
    -------
    spec2DObj_A, spec2DObj_B

    spec2DObj_A (object, Spec2D):
       Spec2d Object for extraction at A position
    spec2DObj_B (object, Spec2D)
       Spec2d Object for extraction at B position

    """

    bkg_redux = bkg_files is not None
    # Build Science image
    sciImg = buildimage.buildimage_fromlist(
        spectrograph, caliBrate.det, parset['scienceframe'], list(files), bpm=caliBrate.msbpm, slits=caliBrate.slits,
        ignore_saturation=False)

    if bkg_files is not None:
        # Background Image?
        bgimg = buildimage.buildimage_fromlist(spectrograph, caliBrate.det, parset['scienceframe'],
                                               list(bkg_files), bpm=caliBrate.msbpm,
                                               slits=caliBrate.slits, ignore_saturation=False)
        sciImg = sciImg.sub(bgimg)

    # DP: Should find_negative be True here? JFH: For quicklook yes!
    objFind = find_objects.FindObjects.get_instance(sciImg, caliBrate.slits, spectrograph, parset, 'science',
                                                    waveTilts=caliBrate.wavetilts,
                                                    bkg_redux=bkg_redux, find_negative=bkg_redux, show=show)

    global_sky, sobjs_obj = objFind.run(std_trace=std_trace, show_peaks=show)

    # Instantiate Extract object
    extract = extraction.Extract.get_instance(sciImg, caliBrate.slits, sobjs_obj, spectrograph, parset, 'science',
                                              global_sky=global_sky, waveTilts=caliBrate.wavetilts,
                                              wv_calib=caliBrate.wv_calib,
                                              bkg_redux=bkg_redux, return_negative=bkg_redux, show=show)
    skymodel, objmodel, ivarmodel, outmask, sobjs, waveimg, tilts = extract.run()

    # TODO -- Do this upstream
    # Tack on detector
    for sobj in sobjs:
        sobj.DETECTOR = sciImg.detector

    # Construct table of spectral flexure
    spec_flex_table = Table()
    spec_flex_table['spat_id'] = caliBrate.slits.spat_id
    spec_flex_table['sci_spec_flexure'] = extract.slitshift

    # Construct the Spec2DObj with the positive image
    spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                    ivarraw=sciImg.ivar,
                                    skymodel=skymodel,
                                    objmodel=objmodel,
                                    ivarmodel=ivarmodel,
                                    scaleimg=None,
                                    waveimg=waveimg,
                                    bpmmask=outmask,
                                    detector=sciImg.detector,
                                    sci_spat_flexure=sciImg.spat_flexure,
                                    sci_spec_flexure=spec_flex_table,
                                    vel_corr=None,
                                    vel_type=parset['calibrations']['wavelengths']['refframe'],
                                    tilts=tilts,
                                    slits=copy.deepcopy(caliBrate.slits),
                                    wavesol=caliBrate.wv_calib.wave_diagnostics(print_diag=False),
                                    maskdef_designtab=None)
    spec2DObj.process_steps = sciImg.process_steps

    if not bkg_redux:
        return spec2DObj
    else:
        # Construct the Spec2DObj with the negative image
        spec2DObj_bkg = spec2dobj.Spec2DObj(sciimg=-sciImg.image,
                                           ivarraw=sciImg.ivar,
                                           skymodel=-skymodel,
                                           objmodel=-objmodel,
                                           ivarmodel=ivarmodel,
                                           scaleimg=None,
                                           waveimg=waveimg,
                                           bpmmask=outmask,
                                           detector=sciImg.detector,
                                           sci_spat_flexure=sciImg.spat_flexure,
                                           sci_spec_flexure=spec_flex_table,
                                           vel_corr=None,
                                           vel_type=parset['calibrations']['wavelengths']['refframe'],
                                           tilts=tilts,
                                           slits=copy.deepcopy(caliBrate.slits),
                                           wavesol=caliBrate.wv_calib.wave_diagnostics(print_diag=False),
                                           maskdef_designtab=None)
        return spec2DObj, spec2DObj_bkg




class QL_Multislit(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Script to produce quick-look multislit PypeIt reductions', width=width)
        parser.add_argument('spectrograph', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))
        parser.add_argument('full_rawpath', type=str, help='Full path to the raw files')
        parser.add_argument('files', type=str, nargs='+',
                            help='list of frames i.e. img1.fits img2.fits')
        parser.add_argument('--spec_samp_fact', default=1.0, type=float,
                            help='Make the wavelength grid finer (spec_samp_fact < 1.0) or '
                                 'coarser (spec_samp_fact > 1.0) by this sampling factor, i.e. '
                                 'units of spec_samp_fact are pixels.')
        parser.add_argument('--spat_samp_fact', default=1.0, type=float,
                            help='Make the spatial grid finer (spat_samp_fact < 1.0) or coarser '
                                 '(spat_samp_fact > 1.0) by this sampling factor, i.e. units of '
                                 'spat_samp_fact are pixels.')
        parser.add_argument("--bkg_redux", default=False, action='store_true',
                            help='If set the script will perform difference imaging quicklook. Namely it will identify '
                                 'sequences of AB pairs based on the dither pattern and perform difference imaging sky '
                                 'subtraction and fit for residuals')
        parser.add_argument("--flux", default=False, action='store_true',
                            help='This option will multiply in sensitivity function to obtain a '
                                 'flux calibrated 2d spectrum')
        parser.add_argument("--mask_cr", default=False, action='store_true',
                            help='This option turns on cosmic ray rejection. This improves the '
                                 'reduction but doubles runtime.')
        parser.add_argument("--writefits", default=False, action='store_true',
                            help="Write the ouputs to a fits file")
        parser.add_argument('--no_gui', default=False, action='store_true',
                            help="Do not display the results in a GUI")
        parser.add_argument('--box_radius', type=float,
                            help='Set the radius for the boxcar extraction')
        parser.add_argument('--offset', type=float, default=None,
                            help='Override the automatic offsets determined from the headers. '
                                 'Offset is in pixels.  This option is useful if a standard '
                                 'dither pattern was not executed.  The offset convention is '
                                 'such that a negative offset will move the (negative) B image '
                                 'to the left.')
        parser.add_argument("--redux_path", type=str, default='current working directory',
                            help="Location where reduction outputs should be stored.")
        parser.add_argument("--master_dir", type=str, default=os.getenv('QL_MASTERS'),
                            help="Location of PypeIt Master files used for the reduction.")
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        parser.add_argument("--show", default=False, action="store_true",
                            help='Show the reduction steps. Equivalent to the -s option when '
                                 'running pypeit.')
        parser.add_argument('--det', type=str, default='1', nargs='*',
                            help='Detector(s) to show.  If more than one, the list of detectors '
                                 'must be one of the allowed mosaics hard-coded for the selected '
                                 'spectrograph.')
        return parser


    @staticmethod
    def main(args):

        # Parse the detector this is taken from view_fits but this should be made into a utility function

        tstart = time.time()
        # Parse the files sort by MJD
        files = np.array([pathlib.Path(args.full_rawpath) / file for file in args.files])
        nfiles = len(files)


        # Read in the spectrograph, config the parset
        spectrograph = load_spectrograph(args.spectrograph)
        spectrograph_cfg_lines = spectrograph.config_specific_par(files[0]).to_config()
        parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines,
                                              merge_with=(config_lines(args),))
        _det = parse_det(args.det, spectrograph)

        target = spectrograph.get_meta_value(files[0], 'target')
        mjds = np.zeros(nfiles)
        for ifile, file in enumerate(files):
            mjds[ifile] = spectrograph.get_meta_value(file, 'mjd', ignore_bad_header=True,
                                                      no_fussing=True)
        files = files[np.argsort(mjds)]

        # Get the master path

        # Calibration Master directory
        master_dir = data.Paths.data / 'QL_MASTERS' if args.master_dir is None else pathlib.Path(args.master_dir)
        master_subdir = spectrograph.get_ql_master_dir(files[0])
        master_path = master_dir / master_subdir
        if not master_path.is_dir():
            msgs.error(f'{master_path} does not exist!  You must install the QL_MASTERS '
                       'directory; download the data from the PypeIt dev-suite Google Drive and '
                       'either define a QL_MASTERS environmental variable or use the '
                       'pypeit_install_ql_masters script.')

        bias_masterframe_name = \
            utils.find_single_file(master_path / "MasterBias*")
        slit_masterframe_name = \
            utils.find_single_file(master_path / "MasterSlits*")
        tilts_masterframe_name = \
            utils.find_single_file(master_path / "MasterTilts*")
        wvcalib_masterframe_name = \
            utils.find_single_file(master_path / "MasterWaveCalib*")
        std_spec1d_file = \
            utils.find_single_file(master_path / "spec1d_*")
        sensfunc_masterframe_name = \
            utils.find_single_file(master_path / "sens_*")


        # TODO Implement some kind of checking for minimal masters. If --flux is set check for sensfunc etc.
        #if (bias_masterframe_name is None or not os.path.isfile(bias_masterframe_name)) or \
        if (slit_masterframe_name is None or not slit_masterframe_name.is_file()) or \
                (tilts_masterframe_name is None or not tilts_masterframe_name.is_file()) or \
                (std_spec1d_file is None or not std_spec1d_file.is_file()):
            # or (sensfunc_masterframe_name is None or not os.path.isfile(sensfunc_masterframe_name)):
            msgs.error('Master frames not found.  Check that environment variable QL_MASTERS '
                       'points at the Master Calibs')

        det_container = spectrograph.get_detector_par(_det, hdu=fits.open(files[0]))
        binspectral, binspatial = parse_binning(det_container['binning'])
        platescale = det_container['platescale']*binspatial
        detname = det_container.name

        if std_spec1d_file is not None:
            std_trace = specobjs.get_std_trace(detname, std_spec1d_file, chk_version=False)
        else:
            std_trace = None

        # Parse the offset information out of the headers. TODO in the future
        # get this out of fitstable
        dither_pattern, dither_id, offset_arcsec = spectrograph.parse_dither_pattern(files)

        print_offset_report(files, dither_pattern, dither_id, offset_arcsec, target, platescale)
        caliBrate = build_calibrate(_det, files, spectrograph, parset, bias_masterframe_name,
                                        slit_masterframe_name, wvcalib_masterframe_name, tilts_masterframe_name)

        spec2d_list, offsets_dith_pix = run(files, dither_id, offset_arcsec, caliBrate, spectrograph,
                                            platescale, parset, std_trace, args.show, bkg_redux=args.bkg_redux)

        # Override offsets if they were passed in?
        if args.offset is not None:
            offsets_pixels = np.array([0.0, args.offset])
            msgs.info('Using user specified offsets instead: {:5.2f}'.format(args.offset))
        else:
            offsets_pixels = offsets_dith_pix


        # Instantiate Coadd2d
        coadd = coadd2d.CoAdd2D.get_instance(spec2d_list, spectrograph, parset, det=_det,
                                             offsets=offsets_pixels, weights='uniform',
                                             spec_samp_fact=args.spec_samp_fact,
                                             spat_samp_fact=args.spat_samp_fact,
                                             bkg_redux=True, debug=args.show)
        # Coadd the slits
        # TODO implement only_slits later
        islit = 0
        coadd_dict_list = coadd.coadd(only_slits=None, interp_dspat=False)
        # Create the pseudo images
        pseudo_dict = coadd.create_pseudo_image(coadd_dict_list)

        # Multiply in a sensitivity function to flux the 2d image
        if args.flux:
            # Load the sensitivity function
            #            wave_sens, sfunc, _, _, _ = sensfunc.SensFunc.load(sensfunc_masterframe_name)
            sens = sensfunc.SensFunc.from_file(sensfunc_masterframe_name)
            # Interpolate the sensitivity function onto the wavelength grid of
            # the data. Since the image is rectified this is trivial and we
            # don't need to do a 2d interpolation
            exptime = spectrograph.get_meta_value(files[0], 'exptime')
            sens_factor = flux_calib.get_sensfunc_factor(pseudo_dict['wave_mid'][:, islit],
                                                         sens.wave.flatten(), sens.zeropoint.flatten(), exptime,
                                                         extrap_sens=True) #parset['fluxcalib']['extrap_sens'])

            # Compute the median sensitivity and set the sensitivity to zero at
            # locations 100 times the median. This prevents the 2d image from
            # blowing up where the sens_factor explodes because there is no
            # throughput
            sens_gpm = sens_factor < 100.0 * np.median(sens_factor)
            sens_factor_masked = sens_factor * sens_gpm
            sens_factor_img = np.repeat(sens_factor_masked[:, np.newaxis], pseudo_dict['nspat'],
                                        axis=1)
            imgminsky = sens_factor_img * pseudo_dict['imgminsky']
            imgminsky_gpm = sens_gpm[:, np.newaxis] & pseudo_dict['inmask']
        else:
            imgminsky = pseudo_dict['imgminsky']
            imgminsky_gpm = pseudo_dict['inmask']

        ##########################
        # Now display the images #
        ##########################
        if not args.no_gui:
            sdet = get_dnum(_det, prefix=False)
            display.connect_to_ginga(raise_err=True, allow_new=True)

            # TODO: Bug in ginga prevents me from using cuts here for some
            # reason
            mean, med, sigma = sigma_clipped_stats(imgminsky[imgminsky_gpm], sigma_lower=3.0,
                                                   sigma_upper=3.0)
            chname_skysub = 'fluxed-skysub-det{:s}'.format(sdet) \
                if args.flux else 'skysub-det{:s}'.format(sdet)
            cuts_skysub = (med - 3.0 * sigma, med + 3.0 * sigma)
            cuts_resid = (-5.0, 5.0)
            # fits.writeto('/Users/joe/ginga_test.fits',imgminsky, overwrite=True)
            # fits.writeto('/Users/joe/ginga_mask.fits',imgminsky_gpm.astype(float), overwrite=True)
            # embed()

            # Clear all channels at the beginning
            # TODO: JFH For some reason Ginga crashes when I try to put cuts in here.
            viewer, ch_skysub = display.show_image(imgminsky, chname=chname_skysub,
                                                   waveimg=pseudo_dict['waveimg'], clear=True,
                                                   cuts=cuts_skysub)
            slit_left, slit_righ, _ = pseudo_dict['slits'].select_edges()
            slit_id = caliBrate.slits.slitord_id[0]
            display.show_slits(viewer, ch_skysub, slit_left, slit_righ, slit_ids=slit_id)

            # SKRESIDS
            chname_skyresids = 'sky_resid-det{:s}'.format(sdet)
            # sky residual map
            image = pseudo_dict['imgminsky'] * np.sqrt(pseudo_dict['sciivar']) * pseudo_dict['inmask']
            viewer, ch_skyresids = display.show_image(image, chname_skyresids,
                                                      waveimg=pseudo_dict['waveimg'],
                                                      cuts=cuts_resid)

            display.show_slits(viewer, ch_skyresids, slit_left, slit_righ, slit_ids=caliBrate.slits.slitord_id[0])
            shell = viewer.shell()
            out = shell.start_global_plugin('WCSMatch')
            out = shell.call_global_plugin_method('WCSMatch', 'set_reference_channel',
                                                  [chname_skysub], {})

        # TODO extract along a spatial position
        if args.writefits:
            head0 = fits.getheader(files[0])
            # TODO use meta tools for the object name in the future.
            outfile = target + '_specXspat_{:3.2f}X{:3.2f}.fits'.format(args.spec_samp_fact,
                                                                        args.spat_samp_fact)
            hdu = fits.PrimaryHDU(imgminsky, header=head0)
            hdu_resid = fits.ImageHDU(pseudo_dict['imgminsky'] \
                                      * np.sqrt(pseudo_dict['sciivar']) * pseudo_dict['inmask'])
            hdu_wave = fits.ImageHDU(pseudo_dict['waveimg'])
            hdul = fits.HDUList([hdu, hdu_resid, hdu_wave])
            msgs.info('Writing sky subtracted image to {:s}'.format(outfile))
            hdul.writeto(outfile, overwrite=True)

        msgs.info(utils.get_time_string(time.time()-tstart))


        if args.embed:
            embed()

        return 0
