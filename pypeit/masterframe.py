""" Abstract class for Master Image frames
      This could go inside masters.py
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import os
import warnings

from pypeit import msgs
from astropy.io import fits
from pypeit.par import pypeitpar
import os

from abc import ABCMeta


class MasterFrame(object):
    """
    This class is designed to gather a set of like methods
    for Master calibration frames

    Args:
        frametype (str):  Frametype of the parent class
        master_key (str or None):
          Root of the MasterFrame names, e.g. 'A_1_01'.
          If None, set to 'master'
        master_dir (str or None):
          Name of the MasterFrame folder, e.g. MF_keck_deimos
          If None, set to './'
        reuse_masters (bool, optional):
          Reuse already created master files from disk.

    """
    __metaclass__ = ABCMeta

    def __init__(self, frametype, master_key, master_dir, reuse_masters=False):

        # Output names
        self.frametype = frametype
        if master_dir is None:
            master_dir = './'
        self.master_dir = master_dir

        if master_key is None:
            master_key = 'master'
        self.master_key = master_key

        # Other parameters
        self.reuse_masters=reuse_masters
        self.msframe = None

    @property
    def ms_name(self):
        """ Default filenames for MasterFrames

        Returns:
            str:  Master filename
        """
        return master_name(self.frametype, self.master_key, self.master_dir)

    @property
    def mdir(self):
        """

        Returns:
            str: Master frames folder

        """
        return self.master_dir

    def master(self, prev_build=False):
        """
        Load the master frame from disk, as settings allows. This routine checks the the mode of master usage
        then calls the load_master method. This method should not be overloaded by children of this class. Instead
        one should overload the load_master method below.

        Args:
            prev_build (bool, optional):
                If True, try to load master from disk

        Returns:
            ndarray or None:  Master image

        """
        # Are we loading master files from disk?
        if self.reuse_masters or prev_build:
            self.msframe, head = self.load_master(self.ms_name)
            return self.msframe
        else:
            return None

    def load_master(self, filename, exten=0):
        """
        Generic master file reader. This function should mostly be replaced by specific load_master methods in the children
        of this class.

        Args:
            filename (str):  Master filename
            exten (int, optional):  Extension of the file to load from

        Returns:
            ndarray or None:  master frame image

        """

        # Does the master file exist?
        if not os.path.isfile(filename):
            msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, filename))
            return None, None
        else:
            msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(self.frametype) + " from filename: {:}".format(filename))
            hdu = fits.open(filename)
            # msgs.info("Master {0:s} frame loaded successfully:".format(hdu[0].header['FRAMETYP'])+msgs.newline()+name)
            head0 = hdu[0].header
            data = hdu[exten].data.astype(np.float)
            # List of files used to generate the Master frame (e.g. raw file frames)
            file_list = []
            for key in head0:
                if 'FRAME' in key:
                    file_list.append(head0[key])
            #return data , head0, file_list
            return data, head0

    def save_master(self, data, outfile=None, raw_files=None, steps=None, overwrite=True, extensions=None, names=None):
        """
        Core function to write a MasterFrame image.  Intended for simple images only; more complex masters
        need their own method and may be written by their own Class (e.g. TraceSlits) and should thus overwrite this method.

        Args:
            data (ndarray or dict):  Object to save
            outfile (str, optional): Outfile name
            raw_files (list, optional): List of raw filenames
            steps (list, optional): list of steps executed in the master file creation
            extensions (list, optional): Additional data images to write
            names (list, optional): Names of the extensions

        """
        _outfile = self.ms_name if outfile is None else outfile
        # Additional keywords for the Header
        keywds = None if steps is None else dict(steps=','.join(steps))
        # Check for existing
        if os.path.exists(_outfile) and (not overwrite):
            msgs.warn("This file already exists.  Use overwrite=True to overwrite it")
            return
        #
        msgs.info("Saving master {0:s} frame as:".format(self.frametype) + msgs.newline() + _outfile)

        hdu = fits.PrimaryHDU(data)
        hlist = [hdu]
        # Extensions
        if extensions is not None:
            for kk,exten in enumerate(extensions):
                hdu = fits.ImageHDU(exten)
                if names is not None:
                    hdu.name = names[kk]
                hlist.append(hdu)
        # HDU list
        hdulist = fits.HDUList(hlist)
        # Header
        msgs.info("Writing header information")
        if raw_files is not None:
            for i in range(len(raw_files)):
                hdrname = "FRAME{0:03d}".format(i+1)
                hdulist[0].header[hdrname] = (raw_files[i], 'PyepIt: File used to generate Master {0:s}'.format(self.frametype))
        hdulist[0].header["FRAMETYP"] = (self.frametype, 'PyepIt: Master calibration frame type')
        if keywds is not None:
            for key in keywds.keys():
                hdulist[0].header[key] = keywds[key]
        # Write the file to disk
        if os.path.exists(_outfile):
            msgs.warn("Overwriting file:" + msgs.newline() + _outfile)

        hdulist.writeto(_outfile, overwrite=True)

        # Finish
        msgs.info("Master {0:s} frame saved successfully:".format(self.frametype) + msgs.newline() + _outfile)
        return


def master_name(ftype, master_key, mdir):
    """ Default filenames for MasterFrames

    Args:
        ftype (str):
          Frame type
        master_key (str):
          Setup name, e.b. A_1_01
        mdir (str):
          Master directory

    Returns:
        str: masterframe filename
    """
    name_dict = dict(bias='{:s}/MasterBias_{:s}.fits'.format(mdir, master_key),
                     badpix='{:s}/MasterBadPix_{:s}.fits'.format(mdir, master_key),
                     trace='{:s}/MasterTrace_{:s}.fits'.format(mdir, master_key),
                     pinhole='{:s}/MasterPinhole_{:s}.fits'.format(mdir, master_key),
                     pixelflat='{:s}/MasterPixelFlat_{:s}.fits'.format(mdir, master_key),
                     illumflat='{:s}/MasterIllumFlat_{:s}.fits'.format(mdir, master_key),
                     arc='{:s}/MasterArc_{:s}.fits'.format(mdir, master_key),
                     wave='{:s}/MasterWave_{:s}.fits'.format(mdir, master_key),
                     wv_calib='{:s}/MasterWaveCalib_{:s}.json'.format(mdir, master_key),
                     tilts='{:s}/MasterTilts_{:s}.fits'.format(mdir, master_key),
                     # sensfunc='{:s}/MasterSensFunc_{:s}_{:s}.yaml'.format(mdir, master_key[0], master_key[-2:]),
                     sensfunc='{:s}/MasterSensFunc_{:s}_{:s}.fits'.format(mdir, master_key[0], master_key[-2:]),
                     )
    return name_dict[ftype]


def set_master_dir(redux_path, spectrograph, par):
    """
    Set the master directory auto-magically

    Args:
        redux_path: str or None
        spectrograph: Spectrograph or None
        par: ParSet or None

    Returns:
        master_dir : str
          Path of the MasterFrame directory

    """
    # Parameters
    if par is None:
        tmppar = pypeitpar.CalibrationsPar()
    else:
        if 'caldir' not in par.keys():
            tmppar = pypeitpar.CalibrationsPar()
        else:
            tmppar = par
    # Redux path
    if redux_path is None:
        redux_path = os.getcwd()
    master_dir = os.path.join(redux_path, tmppar['caldir'])
    # Spectrograph
    if spectrograph is not None:
        master_dir += '_'+spectrograph.spectrograph
    # Return
    return master_dir

