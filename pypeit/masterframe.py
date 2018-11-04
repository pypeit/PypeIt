""" Abstract class for Master Image frames
      This could go inside masters.py
"""
from __future__ import absolute_import, division, print_function

import os
import warnings

from pypeit import msgs
from astropy.io import fits
from pypeit.par import pypeitpar
import os

from pypeit import debugger

from abc import ABCMeta


class MasterFrame(object):
    """
    This class is designed to gather a set of like methods
    for Master calibration frames

    Parameters
    ----------
    frametype : str
    setup : str
    master_dir : str, optional
    redux_path : str, optional
      Path for reduction
    spectrograph : Spectrograph, optional
      Only used for directory_path;  should be Deprecated

    Attributes
    ----------
    """
    __metaclass__ = ABCMeta

    def __init__(self, frametype, setup, spectrograph=None,
                 master_dir=None, mode=None, par=None, redux_path=None):

        # Output path
        if master_dir is None:
            self.master_dir = self._set_master_dir(redux_path, spectrograph, par)
        else:
            self.master_dir = master_dir

        # Other parameters
        self.frametype = frametype
        self.setup = setup
        self.mode = mode
        self.msframe = None

    @property
    def ms_name(self):
        """ Default filenames for MasterFrames

        Parameters
        ----------

        Returns
        -------
        msname : str
        """
        name_dict = dict(bias='{:s}/MasterBias_{:s}.fits'.format(self.master_dir, self.setup),
                         badpix='{:s}/MasterBadPix_{:s}.fits'.format(self.master_dir, self.setup),
                         trace='{:s}/MasterTrace_{:s}'.format(self.master_dir, self.setup),  # Just a root as FITS+JSON are generated
                         pinhole='{:s}/MasterPinhole_{:s}.fits'.format(self.master_dir, self.setup),
                         pixelflat='{:s}/MasterPixelFlat_{:s}.fits'.format(self.master_dir, self.setup),
                         illumflat='{:s}/MasterIllumFlat_{:s}.fits'.format(self.master_dir, self.setup),
                         arc='{:s}/MasterArc_{:s}.fits'.format(self.master_dir, self.setup),
                         wave='{:s}/MasterWave_{:s}.fits'.format(self.master_dir, self.setup),
                         wv_calib='{:s}/MasterWaveCalib_{:s}.json'.format(self.master_dir, self.setup),
                         tilts='{:s}/MasterTilts_{:s}.fits'.format(self.master_dir, self.setup),
                         # sensfunc='{:s}/MasterSensFunc_{:s}_{:s}.yaml'.format(self.master_dir, self.setup[0], self.setup[-2:]),
                         sensfunc='{:s}/MasterSensFunc_{:s}_{:s}.fits'.format(self.master_dir, self.setup[0], self.setup[-2:]),
                         )

        return name_dict[self.frametype]

    @property
    def mdir(self):
        return self.master_dir

    def _masters_load_chk(self):
        # Logic on whether to load the masters frame
        return self.mode == 'reuse' or self.mode == 'force'

    def _set_master_dir(self, redux_path, spectrograph, par):
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
            master_dir += '_' + spectrograph.spectrograph
        # Return
        return master_dir


    def master(self, force = False):
        """
        Load the master frame from disk, as settings allows. This routine checks the the mode of master usage
        then calls the load_master method. This method should not be overloaded by children of this class. Instead
        one should overload the load_master method below.

        Returns
        -------
        msframe : ndarray or None
         master image

        """
        # Are we loading master files from disk?
        if self._masters_load_chk() or force:
            # Does the master file exist?
            if not os.path.isfile(self.ms_name):
                msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, self.ms_name))
                if force:
                    msgs.error("Crashing out because reduce-masters-force=True:" + msgs.newline() + self.ms_name)
                return None
            else:
                # If master file exists then load and return it
                self.msframe = self.load_master(self.ms_name)
                return self.msframe.copy()
        else:
            return None


    def load_master(self, filename, exten = 0):
        """
        Generic master file reader. This function should mostly be replaced by specific load_master methods in the children
        of this class.

        Returns
        -------
        msframe : ndarray or None
          master frame image

        """

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
        return data

    def save_master(self, data, outfile=None, raw_files=None, steps=None, overwrite=True, extensions=None, names=None):
        """
        Core function to write a MasterFrame image.  Intended for simple images only; more complex masters
        need their own method and may be written by their own Class (e.g. TraceSlits) and should thus overwrite this method.

        Parameters
        ----------
        data : ndarray or dict
        outfile : str (optional)
        raw_files : list (optional)
        steps : list (optional)
          steps executed in the master file creation
        extensions : list, optional
           Additional data images to write
        names : list, optional
           Names of the extensions


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
                hdulist[0].header[hdrname] = (raw_files[i], 'PYPIT: File used to generate Master {0:s}'.format(self.frametype))
        hdulist[0].header["FRAMETYP"] = (self.frametype, 'PYPIT: Master calibration frame type')
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