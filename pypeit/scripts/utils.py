#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This file contains a series of utility functions
that can be used with the PypeIt scripts.
"""

import os
import numpy as np
from configobj import ConfigObj

from astropy.table import Table

from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import PypeItPar
from pypeit.par.util import parse_pypeit_file
from pypeit.metadata import PypeItMetaData
from pypeit.core import procimg
from pypeit import msgs
from pypeit import slittrace
from pypeit import masterframe


class Utilities:
    """This class contains a series of convenience functions
    that are used by many of the pypeit scripts. None of the
    functions in this script should be used by the reduction
    code.

    Attributes:
        pypeit_file (:obj:`str`):
            The name of the input .pypeit file used for reduction.
        det (:obj:`int`):
            Detector index
        iFile (:obj:`int`):
            Index of the frame being displayed/loaded.
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.
        par (:class:`pypeit.par.pypeitpar.PypeItPar`):
            PypeIt parameters used to set the code behavior.  If not
            provided, the default parameters specific to the provided
            spectrograph are used.
        cfg_lines (:obj:`list`):
            List of configuration lines
        data_files (:obj:`list`):
            List of data files to read
        frametype (:obj:`list`):
            List of frametypes for each file
        usrdata (:obj:`astropy.table.Table`):
            Table of user supplied info on data files
        setups (:obj:`list`):
            List of setup lines
    """
    def __init__(self, pypeit_file, det=0):
        self.pypeit_file = pypeit_file
        self.det = det
        self.iFile = None  # if a single frame is being used, set the index, otherwise None
        self.spectrograph = None
        self.par = None

        # Load the pypeit file
        self.cfg_lines, self.data_files, self.frametype, self.usrdata, self.setups =\
            parse_pypeit_file(pypeit_file, runtime=False)

    def check_index(self, iFile):
        if self.iFile is None and iFile is None:
            msgs.error("file index must be selected to get base name")
        elif iFile is None:
            iFile = self.iFile
        return iFile

    def get_basename(self, iFile=None):
        idx = self.check_index(iFile)
        # Grab the filename
        fname = self.data_files[idx]
        return os.path.basename(fname)

    def get_master_dirkey(self, iFile=None):
        # Do some checks
        if self.spectrograph is None or self.par is None:
            msgs.error()
        idx = self.check_index(iFile)
        # Get the master dir
        if self.par['calibrations']['master_dir'] == 'default':
            mdir = os.path.join(self.par['rdx']['redux_path'], 'Masters')
        else:
            mdir = self.par['calibrations']['master_dir']

        if not os.path.exists(mdir):
            mdir_base = os.path.join(os.getcwd(), os.path.basename(mdir))
            msgs.warn('Master file dir: {0} does not exist. Using {1}'.format(mdir, mdir_base))
            mdir = mdir_base
        # Get the base name
        file_base = self.get_basename(iFile=iFile)
        ftdict = dict({file_base: 'science'})
        fitstbl = PypeItMetaData(self.spectrograph, self.par,
                                 files=[self.data_files[idx]], usrdata=Table(self.usrdata[idx]), strict=True)
        fitstbl.finalize_usr_build(ftdict, self.setups[0])
        mkey = fitstbl.master_key(0, det=self.det)
        # Return the result
        return mdir, mkey

    def load_frame(self, iFile=None):
        idx = self.check_index(iFile)
        # Grab the filename
        fname = self.data_files[idx]
        detpar, rawimage, _, _, datasec, _ = self.spectrograph.get_rawimage(fname, self.det)
        rawimage = procimg.trim_frame(rawimage, datasec < 1)
        return self.spectrograph.orient_image(detpar, rawimage)

    def load_spectrograph_parset(self, iFile=None):
        # Do some checks
        idx = self.check_index(iFile)
        # Grab the filename
        fname = self.data_files[idx]
        # Setup the configuration
        cfg = ConfigObj(self.cfg_lines)
        # Load the spectrograph
        spectrograph_name = cfg['rdx']['spectrograph']
        self.spectrograph = load_spectrograph(spectrograph_name)
        msgs.info('Loaded spectrograph {0}'.format(self.spectrograph.name))
        spectrograph_cfg_lines = self.spectrograph.config_specific_par(fname).to_config()
        # Load the parset
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=self.cfg_lines)
        return

    def select_science_frame(self, use_first=False, standard=False):
        """Find all of the indices that correspond to science frames
        """
        sciidx = np.array([], dtype=np.int)
        cntr = 0
        print("\nList of science frames:")
        for tt in range(len(self.usrdata)):
            ftype = self.usrdata['frametype'][tt].split(",")
            for ff in range(len(ftype)):
                if ftype[ff] == "science" or (standard and ftype[ff] == "standard"):
                    sciidx = np.append(sciidx, tt)
                    print("  ({0:d}) {1:s}    ({2:s} frame)".format(cntr+1, self.usrdata['filename'][tt], ftype[ff]))
                    cntr += 1
                    break
        # Determine which science frame the user wants
        if cntr == 1:
            msgs.info("Only one frame listed in .pypeit file - using that frame")
            idx = sciidx[0]
        elif use_first:
            idx = sciidx[0]
        else:
            ans = ''
            while True:
                ans = input(" Which frame would you like to select [1-{0:d}]: ".format(cntr))
                try:
                    ans = int(ans)
                    if 1 <= ans <= cntr:
                        idx = sciidx[ans-1]
                        break
                except ValueError:
                    pass
                msgs.info("That is not a valid option!")
        # Set the selected filename
        self.iFile = idx
        return idx


def get_slits(mkey, mdir):
    slit_masterframe_name = masterframe.construct_file_name(slittrace.SlitTraceSet,
                                                            master_key=mkey,
                                                            master_dir=mdir)
    return slittrace.SlitTraceSet.from_file(slit_masterframe_name)
