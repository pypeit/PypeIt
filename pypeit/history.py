"""
Module for performing two-dimensional coaddition of spectra.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from astropy.time import Time
from pypeit.core.framematch import FrameTypeBitMask

class History:
    """
    Holds and creates history entries for FITS files.

    Attributes:
    history (list of str): List of history entries.
    """
    def __init__(self, header=None):
        """
        Initializes history.

        Args:
        header (:obj:`astropy.io.fits.Header`):
            Header from a fits file. The history tags in ths header will be
            used to populate this History object. Defaults to None.
        """
        self.history = []
        if header is not None and 'HISTORY' in header:
            for card in header['HISTORY']:
                self.history.append(str(card))

    def add_reduce(self, calib_id, metadata, frames, bg_frames):
        """
        Add history entries for reducing a frame.

        Args:
        calib_id (int):  The calibration id being reduced.
        metadata (:obj:`pypeit.metadata.PypeItMetadata1):
            The metadata for all the fits files PypeIt knows of.
        frames (`numpy.ndarray'_): Array of indexes into metadata of the 
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
            self.append('Callibration frames:', add_date=False)

            for frame in calib_frames:
                self.append(f'{frame["frametype"]} "{frame["filename"]}"', add_date=False)

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
        """Write history entries to a FITS header
        Args:
        header (`astropy.io.fits.Header`): The header to write to.
        """
        
        for line in self.history:
            header['history'] = line

