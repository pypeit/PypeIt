"""  Base module for organizing calib frame generation
Allows arms, armed, etc. to call recipes from here instead
of duplicating them within their codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect

import numpy as np
import os

from pypit import msgs

from pypit.core import arsort

from pypit import arcimage
from pypit import biasframe
from pypit import bpmimage
from pypit import traceslits
from pypit import traceimage

from pypit import ardebug as debugger


def msarc(det, setup, sci_ID, spectrograph, fitstbl, tsettings, msbias):
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


def msbias(det, setup, sci_ID, fitstbl, tsettings):
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


def mspbm(det, spectrograph, tsettings, shape, binning=None, reduce_badpix=None, msbias=None):
    bpmImage = bpmimage.BPMImage(spectrograph=spectrograph,
                                 settings=tsettings, det=det,
                                 shape=shape,
                                 binning=binning,
                                 reduce_badpix=reduce_badpix,
                                 msbias=msbias)
    msbpm = bpmImage.build()
    # Return
    return msbpm, bpmImage


def tslits_dict(det, setup, spectrograph, sci_ID, ts_settings, tsettings,
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

