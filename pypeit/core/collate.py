#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script collates multiple 1d spectra in multiple files by object, 
runs flux calibration/coadding on them, and produces files suitable
for KOA archiving.
"""

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle

from pypeit import specobjs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spec2dobj import AllSpec2DObj
from pypeit import msgs
from pypeit.slittrace import SlitTraceBitMask

class SourceObject:
    """ A group of reduced spectra from the same source object. This contains
    the information needed to coadd the spectra and archive the metadata.

    An instance is initiated with the first spectra of the group. Additional
    spectra can be compared with this object to see if it matches using the
    match method, and are added to it if they do.

    Args:
    spec1d_obj (:obj:`pypeit.specobj.SpecObj`):
        The initial spectra of the group as a SpecObj.
    spec1d_header (:obj:`astropy.io.fits.Header`): 
        The header for the first spec1d file in the group.
    spec1d_file (str): Filename of the first spec1d file in the group.
    spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`): 
        The spectrograph that was used to take the data.
    match_type (str): How spectra should be compared. 'ra/dec' means the
        spectra should be compared using the sky coordinates in RA and DEC.
        'pixel' means the spectra should be compared by the spatial pixel
        coordinates in the image.

    Attributes:
    spec_obj_list (list of :obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
        The list of spectra in the group as SpecObj objects.
    spec1d_file_list (list of str): 
        The pathnames of the spec1d files in the group.
    spec1d_header_list: (list of :obj:`astropy.io.fits.Header`):
        The headers of the spec1d files in the group
    """
    def __init__(self, spec1d_obj, spec1d_header, spec1d_file, spectrograph, match_type):
        self.spec_obj_list = [spec1d_obj]
        self.spec1d_file_list = [spec1d_file]
        self.spec1d_header_list = [spec1d_header]
        self._spectrograph = spectrograph
        self.match_type = match_type

        if (match_type == 'ra/dec'):
            self.coord = SkyCoord(spec1d_obj.RA, spec1d_obj.DEC, unit='deg')
        else:
            self.coord = spec1d_obj['SPAT_PIXPOS']

        self.coaddfile = self.build_coadd_file_name()

    def build_coadd_file_name(self):
        """Build the output file name for coadding.
        The filename convention is J<hmsdms+dms>_<instrument name>_<YYYYMMDD>.fits
        when matching by RA/DEC and SPAT_<spatial position>_<instrument name>_<YYYYMMDD>.fits
        when matching by pixel position..

        Currently instrument_name is taken from spectrograph.camera

        Return: str  The name of the coadd output file.
        """
        time_portion = Time(self.spec1d_header_list[0]['MJD'], format="mjd").strftime('%Y%m%d')
        if self.match_type == 'ra/dec':
            coord_portion = 'J' + self.coord.to_string('hmsdms', sep='', precision=2).replace(' ', '')
        else:
            coord_portion = self.spec_obj_list[0]['NAME'].split('_')[0]
        instrument_name = self._spectrograph.camera

        return f'{coord_portion}_{instrument_name}_{time_portion}.fits'

    def _config_key_match(self, header):
        """
        Check to see if the configuration keys from a spec1d file match the
        ones for this SourceObject.

        Args:
        header (:obj:`astropy.io.fits.Header`):
            Header from a spec1d file.

        Returns (bool): True if the configuration keys match, 
            false if they do not.
        """
        # Make sure the spectrograph matches
        if 'PYP_SPEC' not in header or header['PYP_SPEC'] != self._spectrograph.name:
            return False

        first_header = self.spec1d_header_list[0]

        for key in self._spectrograph.configuration_keys():
            # Ignore "decker" because it's valid to coadd spectra with different slit masks
            if key != "decker":
                if key not in first_header and key not in header:
                    # Both are missing the key, this is ok
                    continue
                elif key not in first_header or key not in header:
                    # One has a value and the other doesn't so they don't match
                    return False

                if first_header[key] != header[key]:
                    return False

        # No mismatches were found
        return True

    def match(self, spec_obj, spec1d_header, tolerance, unit = u.arcsec):
        """Determine if a SpecObj matches this group within the given tolerance.
        This will also compare the configuration keys to make sure the SpecObj
        is compatible with the ones in this SourceObject.

        Args:
        spec_obj (:obj:`pypeit.specobj.SpecObj`): 
            The SpecObj to compare with this SourceObject.

        spec1d_header (:obj:`astropy.io.fits.Header`):
            The header from the spec1d that dontains the SpecObj.
            
        tolerance (float): 
            Maximum distance that two spectra can be from each other to be 
            considered to be from the same source. Measured in floating
            point pixels or as an angular distance (see ``unit1`` argument).

        unit (:obj:`astropy.units.Unit`):
            Units of ``tolerance`` argument if match_type is 'ra/dec'. 
            Defaults to arcseconds. Igored if match_type is 'pixel'.

        Returns (bool): True if the SpecObj matches this group,
            False otherwise.
        """

        if not self._config_key_match(spec1d_header):
            return False

        if self.match_type == 'ra/dec':
            coord2 = SkyCoord(ra=spec_obj.RA, dec=spec_obj.DEC, unit='deg')
            return self.coord.separation(coord2) <= Angle(tolerance, unit=unit)
        else:
            coord2 =spec_obj['SPAT_PIXPOS'] 
            return np.fabs(coord2 - self.coord) <= tolerance

def find_slits_to_exclude(spec2d_files, par):
    """Find slits that should be excluded according to the input parameters.
    The slit mask ids are returned in a map alongside the text labels for the
    flags that caused the slit to be excluded.

    Args:
        spec2d_files (:obj:`list` of str): 
            List of spec2d files to build the map from.
        par (:obj:`Collate1DPar):
            Parameters from .collate1d file

    Return:
        :obj:`dict` Mapping of slit mask ids to the flags that caused the
        slit to be excluded.
    """

    # Get the types of slits to exclude from our parameters
    exclude_flags = par['collate1d']['slit_exclude_flags']
    if isinstance(exclude_flags, str):
        exclude_flags = [exclude_flags]

    # Go through the slit_info of all spec2d files and find
    # which slits should be excluded based on their flags
    bit_mask = SlitTraceBitMask()
    exclude_map = dict()
    for spec2d_file in spec2d_files:

        allspec2d = AllSpec2DObj.from_fits(spec2d_file)
        for sobj2d in [allspec2d[det] for det in allspec2d.detectors]:
            for (slit_id, mask, slit_mask_id) in sobj2d['slits'].slit_info:
                for flag in exclude_flags:
                    if bit_mask.flagged(mask, flag):
                        if slit_mask_id not in exclude_map:
                            exclude_map[slit_mask_id] = {flag}
                        else:
                            exclude_map[slit_mask_id].add(flag)

    return exclude_map

def group_spectra_by_source(spec1d_files, exclude_map, match_type, tolerance, unit=u.arcsec):
    """Given a list of spec1d files from PypeIt, group the spectra within the
    files by their source object. The grouping is done by comparing the 
    position of each spectra (using either pixel or RA/DEC) using a given tolerance.

    Args:
        spec1d_files (list of str): A list of spec1d files created by PypeIt.
        exclude_map (dict): Mapping of excluded slit mask ids to the reason
            they were excluded.  Any spectra for these slits will be ignored.
        match_type (str): How to match the positions of spectra. 'pixel' for
            matching by spatial pixel distance or 'ra/dec' for matching by
            angular distance.
        tolerance (float): 
            Maximum distance that two spectra can be from each other to be 
            considered to be from the same source. Measured in floating
            point pixels or as an angular distance (see ``unit1`` argument).
        unit (:obj:`astropy.units.Unit`):
            Units of ``tolerance`` argument if match_type is 'ra/dec'. 
            Defaults to arcseconds. Igored if match_type is 'pixel'.

    Returns (list of `obj`:SourceObject): The grouped spectra as SourceObjects.

    """
    source_list = []

    for spec1d_file in spec1d_files:
        sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)

        for sobj in sobjs.specobjs:

            if sobj.MASKDEF_ID in exclude_map:
                msgs.info(f'Excluding {sobj.MASKDEF_ID} in {spec1d_file} because of flags {exclude_map[sobj.MASKDEF_ID]}')
                continue
            if sobj.OPT_COUNTS is None and sobj.BOX_COUNTS is None:
                msgs.info(f'Excluding {sobj.NAME} in {spec1d_file} because of missing both OPT_COUNTS and BOX_COUNTS')
                continue

            # Search for a SourceObject that matches this SpecObj.
            # If one can't be found, trat this as a new SourceObject.
            found = False
            for source in source_list:
                if  source.match(sobj, sobjs.header, tolerance, unit):
                    source.spec_obj_list.append(sobj)
                    source.spec1d_file_list.append(spec1d_file)
                    source.spec1d_header_list.append(sobjs.header)
                    found = True

            if not found:
                spectrograph = load_spectrograph(sobjs.header['PYP_SPEC'])
                source_list.append(SourceObject(sobj, sobjs.header, spec1d_file, spectrograph, match_type))

    return source_list
