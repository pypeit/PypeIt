"""  Base module for organizing calib frame generation
Allows arms, armed, etc. to call recipes from here instead
of duplicating them within their codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np

from pypit import msgs

from pypit.core import arsort

from pypit import armasters

from pypit import arcimage
from pypit import biasframe
from pypit import bpmimage
from pypit import flatfield
from pypit import traceslits
from pypit import traceimage
from pypit import wavecalib
from pypit import wavetilts
from pypit import waveimage

from pypit import ardebug as debugger


def get_msarc(det, setup, sci_ID, spectrograph, fitstbl, tsettings, msbias):
    # Instantiate with everything needed to generate the image (in case we do)
    arcImage = arcimage.ArcImage([], spectrograph=spectrograph,
                                 settings=tsettings, det=det, setup=setup,
                                 sci_ID=sci_ID, msbias=msbias, fitstbl=fitstbl)
    # Load the MasterFrame (if it exists and is desired)?
    msarc = arcImage.master()
    if msarc is None:  # Otherwise build it
        msgs.info("Preparing a master {0:s} frame".format(arcImage.frametype))
        msarc = arcImage.build_image()
        # Save to Masters
        arcImage.save_master(msarc, raw_files=arcImage.file_list, steps=arcImage.steps)
    # Return
    return msarc, arcImage


def get_msbias(det, setup, sci_ID, fitstbl, tsettings):
    # Init
    biasFrame = biasframe.BiasFrame(settings=tsettings, setup=setup,
                                    det=det, fitstbl=fitstbl, sci_ID=sci_ID)
    # Load the MasterFrame (if it exists and is desired) or the command (e.g. 'overscan')
    msbias = biasFrame.master()
    if msbias is None:  # Build it and save it
        msbias = biasFrame.build_image()
        biasFrame.save_master(msbias, raw_files=biasFrame.file_list, steps=biasFrame.steps)
    # Return
    return msbias, biasFrame


def get_mspbm(det, spectrograph, tsettings, shape, binning=None, reduce_badpix=None, msbias=None):
    bpmImage = bpmimage.BPMImage(spectrograph=spectrograph,
                                 settings=tsettings, det=det,
                                 shape=shape,
                                 binning=binning,
                                 reduce_badpix=reduce_badpix,
                                 msbias=msbias)
    msbpm = bpmImage.build()
    # Return
    return msbpm, bpmImage


def get_msflat(det, setup, sci_ID, fitstbl, tslits_dict, datasec_img,
               flat_settings, msbias, mstilts):
    # Instantiate
    pixflat_image_files = arsort.list_of_files(fitstbl, 'pixelflat', sci_ID)
    flatField = flatfield.FlatField(file_list=pixflat_image_files, msbias=msbias,
                                    settings=flat_settings,
                                    tslits_dict=tslits_dict,
                                    tilts=mstilts, det=det, setup=setup)

    # Load from disk (MasterFrame)?
    mspixflatnrm = flatField.master()
    if mspixflatnrm is None:
        # Use mstrace if the indices are identical
        #if np.all(arsort.ftype_indices(fitstbl,'trace',1) ==
        #                  arsort.ftype_indices(fitstbl, 'pixelflat', 1)) and (traceSlits.mstrace is not None):
        #    flatField.mspixelflat = traceSlits.mstrace.copy()
        # Run
        mspixflatnrm, slitprof = flatField.run(datasec_img, armed=False)
        # Save to Masters
        flatField.save_master(mspixflatnrm, raw_files=pixflat_image_files, steps=flatField.steps)
        flatField.save_master(slitprof, raw_files=pixflat_image_files, steps=flatField.steps,
                              outfile=armasters.core_master_name('slitprof', setup, flat_settings['masters']['directory']))
    else:
        slitprof, _, _ = flatField.load_master_slitprofile()
    # Return
    return mspixflatnrm, slitprof, flatField

def get_mswave(setup, tslits_dict, wvimg_settings, mstilts, wv_calib, maskslits):
    # Instantiate
    waveImage = waveimage.WaveImage(mstilts, wv_calib, settings=wvimg_settings,
                                    setup=setup, maskslits=maskslits,
                                    slitpix=tslits_dict['slitpix'])
    # Attempt to load master
    mswave = waveImage.master()
    if mswave is None:
        mswave = waveImage._build_wave()
    # Save to hard-drive
    waveImage.save_master(mswave, steps=waveImage.steps)
    # Return
    return mswave, waveImage


def get_tslits_dict(det, setup, spectrograph, sci_ID, ts_settings, tsettings,
               fitstbl, pixlocn, msbias, msbpm, trim=True):
    # Instantiate (without mstrace)
    traceSlits = traceslits.TraceSlits(None, pixlocn, settings=ts_settings,
                                       det=det, setup=setup, binbpx=msbpm)

    # Load via masters, as desired
    if not traceSlits.master():
        # Build the trace image first
        trace_image_files = arsort.list_of_files(fitstbl, 'trace', sci_ID)
        Timage = traceimage.TraceImage(trace_image_files,
                                       spectrograph=spectrograph,
                                       settings=tsettings, det=det)
        mstrace = Timage.process(bias_subtract=msbias, trim=trim)

        # Load up and get ready
        traceSlits.mstrace = mstrace
        _ = traceSlits.make_binarr()
        # Now we go forth
        traceSlits.run(arms=True)
        # QA
        traceSlits._qa()
        # Save to disk
        traceSlits.save_master()
    # Return
    return traceSlits.tslits_dict, traceSlits


def get_wv_calib(det, setup, spectrograph, sci_ID, wvc_settings, fitstbl,
                 tslits_dict, pixlocn, msarc, nonlinear=None):

    # Instantiate
    waveCalib = wavecalib.WaveCalib(msarc, spectrograph=spectrograph,
                                    settings=wvc_settings, det=det,
                                    setup=setup, fitstbl=fitstbl, sci_ID=sci_ID)
    # Load from disk (MasterFrame)?
    wv_calib = waveCalib.master()
    # Build?
    if wv_calib is None:
        wv_calib, _ = waveCalib.run(tslits_dict['lcen'], tslits_dict['rcen'],
                                    pixlocn, nonlinear=nonlinear)
        # Save to Masters
        waveCalib.save_master(waveCalib.wv_calib)
    else:
        waveCalib.wv_calib = wv_calib
    # Mask
    wv_maskslits = waveCalib._make_maskslits(tslits_dict['lcen'].shape[1])

    # Return
    return wv_calib, wv_maskslits, waveCalib

def get_wv_tilts(det, setup, tilt_settings, settings_det, tslits_dict,
                 pixlocn, msarc, wv_calib, maskslits):
    # Instantiate
    waveTilts = wavetilts.WaveTilts(msarc, settings=tilt_settings,
                                    det=det, setup=setup,
                                    tslits_dict=tslits_dict, settings_det=settings_det,
                                    pixlocn=pixlocn)
    # Master
    mstilts = waveTilts.master()
    if mstilts is None:
        mstilts, wt_maskslits = waveTilts.run(maskslits=maskslits,
                                              wv_calib=wv_calib)
        waveTilts.save_master()
    else:
        wt_maskslits = np.zeros_like(maskslits, dtype=bool)
    # Return
    return mstilts, wt_maskslits, waveTilts
