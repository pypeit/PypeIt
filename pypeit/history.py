"""
Module for managing the history of PypeIt output files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os.path
import numpy as np
from IPython import embed

from astropy.time import Time
from astropy.io import fits
from pypeit.core.framematch import FrameTypeBitMask

class History:
    """
    Holds and creates history entries for FITS files.

    Args:

        header (`astropy.io.fits.Header`_, optional):
            Header from a fits file. The
            history keyword entries in this header will be used to populate this 
            History object. Defaults to None.

    Attributes:
    
        history (:obj:`list` of `str`): List of history entries.
    """

    def __init__(self, header=None):
        """
        Initializes history.
        """
        self.history = []
        if header is not None and 'HISTORY' in header:
            for card in header['HISTORY']:
                self.history.append(str(card))

    def add_reduce(self, calib_id, metadata, frames, bg_frames):
        """
        Add history entries for reducing a frame. For example::

            HISTORY 2021-03-05T23:56 PypeIt Reducing target HIP15339                        
            HISTORY Combining frames:                                                       
            HISTORY "S20161108S0087.fits.gz"                                                
            HISTORY "S20161108S0090.fits.gz"                                                
            HISTORY Subtracted background from frames:                                      
            HISTORY "S20161108S0088.fits.gz"                                                
            HISTORY "S20161108S0089.fits.gz"                                                
            HISTORY Calibration frames:                                                    
            HISTORY arc,science,tilt "S20161108S0069.fits.gz"                               
            HISTORY arc,science,tilt "S20161108S0070.fits.gz"                               
            HISTORY arc,science,tilt "S20161108S0071.fits.gz"                               
            HISTORY arc,science,tilt "S20161108S0072.fits.gz"                               
            HISTORY pixelflat,trace "S20161108S0078.fits.gz"                                

        Args:
            calib_id (int): The calibration id being reduced.

            metadata (:class:`pypeit.metadata.PypeItMetaData`): The metadata for all
                the fits files PypeIt knows of.

            frames (`numpy.ndarray`_): Array of indexes into metadata of the 
                frames being combined in the reduction.

            bg_frames (`numpy.ndarray`_): Array of indexes into metadata of the 
                frames being subtracted in the reduction.
        """
        self.append(f'PypeIt Reducing target {metadata["target"][frames[0]]}')
        self.append('Combining frames:', add_date=False)
        for frame in frames:
            self.append(f'"{metadata["filename"][frame]}"', add_date=False)
        if len(bg_frames) > 0:
            self.append('Subtracted background from frames:', add_date=False)
            for frame in bg_frames:
                self.append(f'"{metadata["filename"][frame]}"', add_date=False)

        frametype_bitmask = FrameTypeBitMask()
        calibration_types = [x for x in frametype_bitmask.keys() if x not in ['science', 'standard']]

        calib_frames = metadata[metadata.find_frames(calibration_types, calib_id)]
        if len(calib_frames) > 0:
            self.append('Calibration frames:', add_date=False)

            for frame in calib_frames:
                self.append(f'{frame["frametype"]} "{frame["filename"]}"', add_date=False)

    def add_coadd1d(self, spec1d_files, objids, gpm_exp=None):
        """
        Add history entries for 1D coadding.
        
        The history shows what files and objects were used for coadding.
        For example::
            
            HISTORY 2021-01-23T02:12 PypeIt Coadded 4 objects from 3 spec1d files           
            HISTORY File 0 "spec1d_DE.20170425.53065-dra11_DEIMOS_2017Apr25T144418.240.fits"  
            HISTORY File 1 "spec1d_DE.20170425.51771-dra11_DEIMOS_2017Apr25T142245.350.fits"  
            HISTORY File 2 "spec1d_DE.20170425.50487-dra11_DEIMOS_2017Apr25T140121.014.fits"  
            HISTORY Object ID SPAT0692-SLIT0704-DET08 from file 0                           
            HISTORY Object ID SPAT0695-SLIT0706-DET04 from file 2                           
            HISTORY Object ID SPAT0691-SLIT0704-DET08 from file 2                           
            HISTORY Object ID SPAT0695-SLIT0706-DET04 from file 1

        Args:
            spec1d_files (:obj:`list`): List of the spec1d files used for coadding.
            objids (:obj:`list`): List of the PypeIt object ids used in coadding.
            gpm_exp (:obj:`list`, optional): List of boolean indicating which exposures were coadded.
        """

        if gpm_exp is not None:
            # Not coadded files and objids
            notcoadded_spec1d_files = [spec1d_file for (spec1d_file, gpm_exp) in zip(spec1d_files, gpm_exp) if not gpm_exp]
            notcoadded_objids = [objid for (objid, gpm_exp) in zip(objids, gpm_exp) if not gpm_exp]
            combined_notcoadd_files_objids = list(zip(notcoadded_spec1d_files, notcoadded_objids))

            # Coadded files and objids
            coadded_spec1d_files = [spec1d_file for (spec1d_file, gpm_exp) in zip(spec1d_files, gpm_exp) if gpm_exp]
            coadded_objids = [objid for (objid, gpm_exp) in zip(objids, gpm_exp) if gpm_exp]
            combined_files_objids = list(zip(coadded_spec1d_files, coadded_objids))
        else:
            combined_files_objids = list(zip(spec1d_files, objids))
            combined_notcoadd_files_objids = None

        files_objids = [combined_files_objids, combined_notcoadd_files_objids]
        # add history
        for file_objid in files_objids:
            if file_objid is None:
                continue
            elif file_objid == combined_files_objids:
                self.append(f'PypeIt Coadded {len(file_objid)} objects '
                            f'from {np.unique([f[0] for f in file_objid]).size} spec1d files')
            elif file_objid == combined_notcoadd_files_objids and len(file_objid) > 0:
                self.append(f'PypeIt DID NOT COADD {len(file_objid)} objects '
                            f'from {np.unique([f[0] for f in file_objid]).size} spec1d files', add_date=False)

            current_spec1d = ""
            for (spec1d, objid) in file_objid:
                if spec1d != current_spec1d:
                    current_spec1d = spec1d

                    self.append(f'From "{os.path.basename(spec1d)}"', add_date=False)
                    header = fits.getheader(spec1d)
                    additional_info = None
                    if 'SEMESTER' in header:
                        additional_info = f"Semester: {header['SEMESTER']}"
                    if 'PROGID' in header:
                        additional_info += f" Program ID: {header['PROGID']}"
                    if additional_info is not None:
                        self.append(additional_info, add_date=False)
                obj_info = objid
                # get extension names
                hnames = [h.name for h in fits.open(spec1d)]
                # find the extension name that include objid
                ind_ext = np.where([objid in h for h in hnames])[0]
                if ind_ext.size > 0:
                    # get the header for this extension
                    this_ext_header = fits.getheader(spec1d, ext=ind_ext[0])
                    if 'MASKDEF_ID' in this_ext_header:
                        obj_info += f" {this_ext_header['MASKDEF_ID']}"
                    if 'MASKDEF_OBJNAME' in this_ext_header:
                        obj_info += f" {this_ext_header['MASKDEF_OBJNAME']}"
                self.append(obj_info, add_date=False)

    def append(self, history, add_date=True):
        """Append a new history entry.

        Args: 
            history (str): The history text to add.

            add_date (bool): If true a isot formatted date willbe prepended
                to the history entry. Defaults to True.
        """

        if add_date:
            self.history.append(f'{Time.now().to_value("isot", subfmt="date_hm")} {history}')
        else:
            self.history.append(history)
        
    def write_to_header(self, header):
        """Write history entries to a FITS header.

        Args:
            header (`astropy.io.fits.Header`_): The header to write to.
        """
        
        for line in self.history:
            header['history'] = line

