"""
Script for quick-look reductions of Keck MOSFIRE observations.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import copy
from glob import glob
import time

from pkg_resources import resource_filename

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

from pypeit import utils
from pypeit import pypeit
from pypeit import par, msgs
from pypeit import pypeitsetup
from pypeit import wavecalib
from pypeit import wavetilts
from pypeit import spec2dobj
from pypeit import coadd2d
from pypeit import specobjs
from pypeit import slittrace
from pypeit import reduce
from pypeit import calibrations
from pypeit.display import display
from pypeit.images import buildimage
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.parse import get_dnum, parse_binning
from pypeit.core.wavecal import wvutils
from pypeit import sensfunc
from pypeit.core import flux_calib
from pypeit.scripts import scriptbase


def config_lines(args):
    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format('vlt_fors2')]
    cfg_lines += ['    redux_path = {0}'.format(args.redux_path)]
    cfg_lines += ['    scidir = Science_QL']
    # Calibrations
    cfg_lines += ['[baseprocess]']
    cfg_lines += ['    use_pixelflat = False']
    cfg_lines += ['    use_illumflat = False']
    cfg_lines += ['[calibrations]']
    cfg_lines += ['    [[wavelengths]]']
    cfg_lines += ['        refframe = observed']
    if not args.mask_cr:
        cfg_lines += ['[scienceframe]']
        cfg_lines += ['    [[process]]']
        cfg_lines += ['        mask_cr = False']
    cfg_lines += ['[reduce]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['        skip_optimal = True']
    if args.box_radius is not None:  # Boxcar radius
        cfg_lines += ['        boxcar_radius = {0}'.format(args.box_radius)]
    cfg_lines += ['    [[findobj]]']
    cfg_lines += ['        skip_second_find = True']

    return cfg_lines


def run(files, caliBrate, spectrograph, det, parset, show=False, std_trace=None):
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

    # Build Science image
    sciImg = buildimage.buildimage_fromlist(
        spectrograph, det, parset['scienceframe'], list(files),  bias=caliBrate.msbias,
        bpm=caliBrate.msbpm, slits=caliBrate.slits, ignore_saturation=False)

    # Instantiate Reduce object
    # Required for pypeline specific object
    # At instantiaton, the fullmask in self.sciImg is modified
    redux = reduce.Reduce.get_instance(sciImg, spectrograph, parset, caliBrate, 'science', ir_redux=True, show=show, det=det)

    # skymodel, objmodel, ivarmodel, outmask, sobjs, scaleimg, waveimg, tilts = redux.run(
    #     std_trace=std_trace, return_negative=True, show_peaks=show)

    global_sky, sobjs_obj, skymask = redux.run_objfind(std_trace=std_trace, show_peaks=show)
    skymodel, objmodel, ivarmodel, outmask, sobjs, scaleimg, waveimg, tilts = redux.run_extraction(
        global_sky, sobjs_obj, skymask)

    # TODO -- Do this upstream
    # Tack on detector
    for sobj in sobjs:
        sobj.DETECTOR = sciImg.detector

    # Construct table of spectral flexure
    spec_flex_table = Table()
    spec_flex_table['spat_id'] = caliBrate.slits.spat_id
    spec_flex_table['sci_spec_flexure'] = redux.slitshift

    # Construct the Spec2DObj with the positive image
    spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                    ivarraw=sciImg.ivar,
                                    skymodel=skymodel,
                                    objmodel=objmodel,
                                    ivarmodel=ivarmodel,
                                    scaleimg=scaleimg,
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
    return spec2DObj


class QLVLTFORS2(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Script to produce quick-look PypeIt reductions '
                                                'on a pair of MOSFIRE files (A-B)', width=width)
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
        parser.add_argument("--redux_path", type=str, default=os.getcwd(),
                            help="Location where reduction outputs should be stored.")
        parser.add_argument("--master_dir", type=str, default=os.getenv('QL_MASTERS'),
                            help="Location of PypeIt Master files used for the reduction.")
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        parser.add_argument("--show", default=False, action="store_true",
                            help='Show the reduction steps. Equivalent to the -s option when '
                                 'running pypeit.')
        return parser

    @staticmethod
    def main(args):

        tstart = time.time()
        # Parse the files sort by MJD
        files = np.array([os.path.join(args.full_rawpath, file) for file in args.files])
        nfiles = len(files)

        # Read in the spectrograph, config the parset
        spectrograph = load_spectrograph('vlt_fors2')
        #spectrograph_def_par = spectrograph.default_pypeit_par()
        spectrograph_cfg_lines = spectrograph.config_specific_par(files[0]).to_config()
        parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines,
                                              merge_with=config_lines(args))
        science_path = os.path.join(parset['rdx']['redux_path'], parset['rdx']['scidir'])

        target = spectrograph.get_meta_value(files[0], 'target')
        mjds = np.zeros(nfiles)
        for ifile, file in enumerate(files):
            mjds[ifile] = spectrograph.get_meta_value(file, 'mjd', ignore_bad_header=True,
                                                      no_fussing=True)
        files = files[np.argsort(mjds)]

        # Calibration Master directory
        #TODO hardwired for now
        master_dir ='./'
        #master_dir = resource_filename('pypeit', 'data/QL_MASTERS') \
        #    if args.master_dir is None else args.master_dir
        if not os.path.isdir(master_dir):
            msgs.error(f'{master_dir} does not exist!  You must install the QL_MASTERS '
                       'directory; download the data from the PypeIt dev-suite Google Drive and '
                       'either define a QL_MASTERS environmental variable or use the '
                       'pypeit_install_ql_masters script.')

        # Define some hard wired master files here to be later parsed out of the directory
        fors2_grism = spectrograph.get_meta_value(files[0], 'dispname')
        fors2_masters = os.path.join(master_dir, 'FORS2_MASTERS', fors2_grism)


        bias_masterframe_name = \
            utils.find_single_file(os.path.join(fors2_masters, "MasterBias*"))
        slit_masterframe_name \
            = utils.find_single_file(os.path.join(fors2_masters, "MasterSlits*"))
        tilts_masterframe_name \
            = utils.find_single_file(os.path.join(fors2_masters, "MasterTilts*"))
        wvcalib_masterframe_name \
            = utils.find_single_file(os.path.join(fors2_masters, 'MasterWaveCalib*'))
        std_spec1d_file = utils.find_single_file(os.path.join(fors2_masters, 'spec1d_*'))
        sensfunc_masterframe_name = utils.find_single_file(os.path.join(fors2_masters, 'sens_*'))

        # TODO make and impelement sensfunc
        if (bias_masterframe_name is None or not os.path.isfile(bias_masterframe_name)) or \
                (slit_masterframe_name is None or not os.path.isfile(slit_masterframe_name)) or \
                (tilts_masterframe_name is None or not os.path.isfile(tilts_masterframe_name)) or \
                (std_spec1d_file is None or not os.path.isfile(std_spec1d_file)):
            # or (sensfunc_masterframe_name is None or not os.path.isfile(sensfunc_masterframe_name)):
            msgs.error('Master frames not found.  Check that environment variable QL_MASTERS '
                       'points at the Master Calibs')

        # We need the platescale

        # Get detector (there's only one)
        #det = 1 # MOSFIRE has a single detector
        #detector = spectrograph.get_detector_par(det)
        #detname = detector.name

        # We need the platescale
        det_container = spectrograph.get_detector_par(1, hdu=fits.open(files[0]))
        binspectral, binspatial = parse_binning(det_container['binning'])
        platescale = det_container['platescale']*binspatial
        # Parse the offset information out of the headers.
        _, _, offset_arcsec = spectrograph.parse_dither_pattern(files)

        # Print out a report on the offsets
        msg_string = msgs.newline()  + '*******************************************************'
        msg_string += msgs.newline() + ' Summary of offsets for target {:s}:                   '
        msg_string += msgs.newline() + '*******************************************************'
        msg_string += msgs.newline() + '           filename                arcsec   pixels    '
        msg_string += msgs.newline() + '----------------------------------------------------'
        for iexp, file in enumerate(files):
            msg_string += msgs.newline() + '    {:s}    {:6.2f}    {:6.2f}'.format(
                os.path.basename(file), offset_arcsec[iexp], offset_arcsec[iexp] / platescale)
        msg_string += msgs.newline() + '********************************************************'
        msgs.info(msg_string)

        ## Read in the master frames that we need
        ##
        det = 1  # Currently CHIP1 is supported
        if std_spec1d_file is not None:
            # Get the standard trace if need be
            sobjs = specobjs.SpecObjs.from_fitsfile(std_spec1d_file)
            this_det = sobjs.DET == det
            if np.any(this_det):
                sobjs_det = sobjs[this_det]
                sobjs_std = sobjs_det.get_std()
                std_trace = None if sobjs_std is None else sobjs_std.TRACE_SPAT.flatten()
            else:
                std_trace = None
        else:
            std_trace = None

        # Read in the bias
        msbias = buildimage.BiasImage.from_file(bias_masterframe_name)
        # Read in the msbpm
        sdet = get_dnum(det, prefix=False)
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

        # Find the unique offsets. This is a bit of a kludge, i.e. we are considering offsets within
        # 0.1 arcsec of each other to be the same throw, but I should like to be able to specify a tolerance here,
        # but then I need a version of unique that accepts a tolerance
        uniq_offsets, uni_indx = np.unique(np.around(offset_arcsec), return_inverse=True)
        nuniq = uniq_offsets.size
        spec2d_list = []
        offset_ref = offset_arcsec[0]
        offsets_dith_pix = []
        # Generalize to a multiple slits, doing one slit at a time?
        islit = 0

        # Loop over the unique throws and create a spec2d_A and spec2D_B for
        # each, which are then fed into coadd2d with the correct offsets

        # TODO Rework the logic here so that we can print out a unified report
        # on what was actually reduced.

        for iuniq in range(nuniq):
            indx = uni_indx == iuniq
            files_uni = files[indx]
            offsets = offset_arcsec[indx]
            msgs.info('Reducing images for offset = {:}'.format(offsets[0]))
            spec2DObj = run(files_uni, caliBrate, spectrograph, det, parset, show=args.show, std_trace=std_trace)
            spec2d_list += [spec2DObj]
            offsets_dith_pix += [np.mean(offsets)/platescale]

        offsets_dith_pix = np.array(offsets_dith_pix)

        if args.offset is not None:
            offsets_pixels = np.array([0.0, args.offset])
            msgs.info('Using user specified offsets instead: {:5.2f}'.format(args.offset))
        else:
            offsets_pixels = offsets_dith_pix


        # Instantiate Coadd2d
        coadd = coadd2d.CoAdd2D.get_instance(spec2d_list, spectrograph, parset, det=det,
                                             offsets=offsets_pixels, weights='uniform',
                                             spec_samp_fact=args.spec_samp_fact,
                                             spat_samp_fact=args.spat_samp_fact,
                                             ir_redux=True, debug=args.show)
        # Coadd the slits
        # TODO implement only_slits later
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
                                                         sens.wave, sens.zeropoint, exptime,
                                                         extrap_sens=parset['fluxcalib']['extrap_sens'])

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
            slit_id = slits.slitord_id[0]
            display.show_slits(viewer, ch_skysub, slit_left, slit_righ, slit_ids=slit_id)

            # SKRESIDS
            chname_skyresids = 'sky_resid-det{:s}'.format(sdet)
            # sky residual map
            image = pseudo_dict['imgminsky'] * np.sqrt(pseudo_dict['sciivar']) * pseudo_dict['inmask']
            viewer, ch_skyresids = display.show_image(image, chname_skyresids,
                                                      waveimg=pseudo_dict['waveimg'],
                                                      cuts=cuts_resid)

            display.show_slits(viewer, ch_skyresids, slit_left, slit_righ,
                               slit_ids=slits.slitord_id[0])
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


