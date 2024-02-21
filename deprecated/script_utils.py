"""
This file contains a series of utility functions
that can be used with the PypeIt scripts.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import numpy as np
from configobj import ConfigObj

from astropy.table import Table

from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import PypeItPar
from pypeit.metadata import PypeItMetaData
from pypeit.core import procimg
from pypeit import msgs
from pypeit import slittrace
from pypeit.pypeitsetup import PypeItSetup
from pypeit.inputfiles import PypeItFile


class Utilities:
    """
    Class with convenience functions for pypeit scripts.
    
    None of the functions in this script should be used by the reduction code.

    Args:
        pypeit_file (:obj:`str`, optional):
            Name of the pypeit file to read
        spectrograph (:obj:`str`, optional):
            Name of the spectrograph to load.  If None, no spectrograph class is
            instantiated.
        iFile (:obj:`int`, optional):
            0-indexed number of a selected frame in the pypeit file
        det (:obj:`int`, :obj:`tuple`, optional):
            The detector/mosaic to process.

    Attributes:
        pypeit_file (:obj:`str`):
            The name of the input .pypeit file used for reduction.
        pypeitFile (:obj:`~pypeit.inputfiles.PypeItFile`):
            Object containing the data read from the pypeit file
        det (:obj:`int`, :obj:`tuple`):
            Detector/Mosoic identifier
        iFile (:obj:`int`):
            0-indexed number of a selected frame in the pypeit file
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            PypeIt parameters used to set the code behavior.
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            Object with the metadata pulled from all headers of all the files in
            the pypeit file.
        calib_dir (`Path`_):
            Path where the processed calibration files are stored.
    """
    def __init__(self, pypeit_file=None, spectrograph=None, iFile=None, det=0):
        self.load_pypeit(pypeit_file)
        self.det = det
        if spectrograph is None:
            return
        self.load_spectrograph(spectrograph)
        if iFile is None:
            return
        self.load_par(iFile=iFile)
        self.load_metadata()

    def _reinit(self):
        """Reset all the attributes."""
        self.pypeit_file = None
        self.pypeitFile = None
        self.det = None
        self.iFile = None
        self.spectrograph = None
        self.par = None
        self.fitstbl = None
        self.calib_dir = None

    def load_pypeit(self, pypeit_file):
        """
        Read the pypeit file.

        Args:
            pypeit_file (:obj:`str`):
                Name of the file to read.  If None, the object is reinitialized!
        """
        if pypeit_file is None:
            self._reinit()
            return
        self.pypeit_file = pypeit_file
        msgs.info(f'Loading the PypeIt file: {self.pypeit_file}')
        self.pypeItFile = PypeItFile.from_file(self.pypeit_file) 

    def load_spectrograph(self, spectrograph):
        """
        Instantiate the spectrograph object.

        Args:
            spectrograph (:obj:`str`):
                Name of the spectrograph to load.  If None, will try to pull the
                name of the spectrograph from the ingested pypeit file
                configuration parameters.  If that also fails, an error is
                raised.
        """
        if spectrograph is None:
            try:
                spectrograph = ConfigObj(self.pypeitFile.cfg_lines)['rdx']['spectrograph']
            except:
                spectrograph = None
        if spectrograph is None:
            msgs.error('Spectrograph not provided and could not be loaded from pypeit file.')
        self.spectrograph = load_spectrograph(spectrograph)

    def load_par(self, iFile=None):
        """
        Load the pypeit parameters.

        Args:
            iFile (:obj:`int`, optional):
                0-indexed number of a selected frame in the pypeit file.  This
                frame is used to set the configuration-specific parameters.  If
                None, the first file in the pypeit file is used.
        """
        if self.pypeitFile is None or self.spectrograph is None:
            msgs.error('Must first load pypeit file and spectrograph to get parameters')
        self.iFile = iFile
        if self.iFile is None:
            msgs.warn('No example file selected.  Configuration specific parameters will use '
                      'first file listed in the pypeit file.')
            self.iFile = 0
        spectrograph_cfg_lines = self.spectrograph.config_specific_par(
                                    self.pypeitFile.filenames[self.iFile]).to_config()
        # Load the parset
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, 
                                            merge_with=(self.pypeitFile.cfg_lines,))        

    def load_metadata(self):
        """
        Use the ingested pypeit file data to build the metadata
        (:attr:`fitstbl`).
        """
        if self.spectrograph is None or self.par is None or self.pypeitFile is None:
            msgs.error('To get metadata, must set spectrograph, par, and read pypeit file.')
        msgs.info('Compiling metadata')
        self.fitstbl = PypeItMetaData(self.spectrograph, self.par, files=self.pypeItFile.filenames,
                                      usrdata=self.pypeItFile.data, strict=True)
        self.fitstbl.finalize_usr_build(self.pypeItFile.frametypes, self.pypeItFile.setup_name)

    def load_calib_dir(self):
        """
        Load the directory for/with the processed calibration frames
        (:attr:`calib_dir`).
        """
        if self.par is None:
            msgs.error('Must first get the pypeit parameters.')
        self.calib_dir = Path(self.par['rdx']['redux_path']).resolve() \
                            / self.par['calibrations']['calib_dir']

        if not self.calib_dir.exists():
            calib_dir_base = Path().resolve() / self.calib_dir.name
            msgs.warn(f'Calibration directory set by parameters does not exist: {self.calib_dir}.  '
                      f'Trying {calib_dir_base}')
            if not calib_dir_base.exists():
                msgs.error(f'{calib_dir_base} also does not exist.  Cannot set valid calibration '
                           'directory.')
            self.calib_dir = calib_dir_base

    def set_index(self, iFile):
        """
        Set the index of the example file from the pypeit file (:attr:`iFile`).

        Args:
            iFile (:obj:`int`):
                The index of the file in the pypeit file.  If None, use
                :attr:`iFile`.  If that is also None, an error is raised.  An
                error is also raised if the file index is less than 0 or more
                than the number of files in the pypeit file.
        """
        if iFile is None:
            iFile = self.iFile
        if iFile is None:
            msgs.error('Must select file index')
        nfiles = len(self.pypeitFile.filenames)
        if iFile < 0 or iFile >= nfiles:
            msgs.error('Selected file index not within range of ingested pypeit file: '
                       f'{self.pypeit_file} contains {nfiles} files.')
        self.iFile = iFile

    def get_basename(self, iFile=None):
        """
        Get the name of the file at the specified index in the pypeit file.

        Only the file name is returned, not its full path.

        Args:
            iFile (:obj:`int`):
                The index of the file in the pypeit file; see :func:`set_index`.

        Returns:
            :obj:`str`: The name of the file.
        """
        if self.pypeitFile is None:
            msgs.error('Must first load pypeit file')
        self.set_index(iFile)
        return Path(self.pypeitFile.filenames[self.iFile]).revolve().name

    def load_frame(self, iFile=None):
        """
        Load, trim, and re-orient the specified file.

        Args:
            iFile (:obj:`int`):
                The index of the file in the pypeit file; see :func:`set_index`.

        Returns:
            `numpy.ndarray`_: The image data.
        """
        if self.pypeitFile is None or self.spectrograph is None:
            msgs.error('Must first load pypeit file and spectrograph')
        if self.det is None:
            msgs.error('Must provide detector number(s).')
        self.set_index(iFile)
        detpar, rawimage, _, _, datasec, _ \
                = self.spectrograph.get_rawimage(self.pypeitFile.filenames[self.iFile], self.det)
        rawimage = procimg.trim_frame(rawimage, datasec < 1)
        return self.spectrograph.orient_image(detpar, rawimage)

    def select_science_frame(self, use_first=False, standard=False):
        """
        Find all of the indices that correspond to science frames.
        """
        if self.pypeitFile is None:
            msgs.error('Must first load pypeit file')
        sciidx = np.array([], dtype=int)
        cntr = 0
        print("\nList of science frames:")
        for tt in range(len(self.pypeitFile.data)):
            ftype = self.pypeitFile.data['frametype'][tt].split(",")
            for ff in range(len(ftype)):
                if ftype[ff] == "science" or (standard and ftype[ff] == "standard"):
                    sciidx = np.append(sciidx, tt)
                    print("  ({0:d}) {1:s}    ({2:s} frame)".format(cntr+1,
                           self.pypeitFile.data['filename'][tt], ftype[ff]))
                    cntr += 1
                    break
        # Determine which science frame the user wants
        if cntr == 1:
            msgs.info("Only one frame listed in .pypeit file - using that frame")
            return sciidx[0]
        if use_first:
            return sciidx[0]

        # Interactively select the frame to use
        ans = ''
        while True:
            ans = input(f' Which frame would you like to select [1-{cntr}]: ')
            try:
                ans = int(ans)
                if 1 <= ans <= cntr:
                    return sciidx[ans-1]
            except ValueError:
                msgs.info("That is not a valid option!")

    def get_calib_key(self, iFile=None):
        """
        Get the calibration key for a given file.

        Args:
            iFile (:obj:`int`, optional):
                Index of the file in the ingested pypeit file.  If None, use
                :attr:`iFile`.

        Returns:
            :obj:`str`: Calibration key.
        """
        if self.fitstbl is None:
            msgs.error('Must first be able to create PypeItMetaData object.  Make sure you have '
                       'provided the pypeit file and spectrograph.')
        if self.det is None:
            msgs.error('Must provide detector number(s) to find slit calibration frame file.')
        self.set_index(iFile)
        setup = self.fitstbl['setup'][self.iFile]
        calib_id = CalibFrame.ingest_calib_id(self.fitstbl['calib'][self.iFile])
        detname = self.spectrograph.get_det_name(self.det)
        return CalibFrame.construct_calib_key(setup, calib_id, detname)

    def get_slits(self, iFile=None):
        """
        Find the slits calibration frame associated with a given file.

        Args:
            iFile (:obj:`int`, optional):
                Index of the file in the ingested pypeit file.  If None, use
                :attr:`iFile`.

        Returns:
            :class:`~pypeit.slittrace.SlitTraceSet`: The object containing the
            slit information.
        """
        if self.fitstbl is None:
            msgs.error('Must first be able to create PypeItMetaData object.  Make sure you have '
                       'provided the pypeit file and spectrograph.')
        if self.det is None:
            msgs.error('Must provide detector number(s) to find slit calibration frame file.')
        if self.calib_dir is None:
            self.load_calib_dir()
        calib_key = self.get_calib_key(iFile=iFile)

        # Construct the expected calibration frame file name
        slits_file = Path(slittrace.SlitTraceSet.construct_file_name(calib_key,
                            calib_dir=self.calib_dir)).resolve()
        if not slits_file.exists():
            msgs.error(f'Slit calibration frame does not exist: {slits_file}')
        return slittrace.SlitTraceSet.from_file(slits_file)

