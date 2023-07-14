#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This module contains code for collating multiple 1d spectra source object.
"""

#TODO -- Consider moving this into the top-level, i.e. out of core

import copy
import os.path

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle

from pypeit import specobjs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import msgs

class SourceObject:

    """A group of reduced spectra from the same source object. This contains
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
            try:
                self.coord = SkyCoord(spec1d_obj.RA, spec1d_obj.DEC, unit='deg')
            except Exception as e:
                msgs.error(f"Cannot do ra/dec matching on {spec1d_obj.NAME}, could not read RA/DEC.")
        else:
            self.coord = spec1d_obj['SPAT_PIXPOS']

    @classmethod
    def build_source_objects(cls, spec1d_files, match_type):
        """Build a list of SourceObjects from a list of spec1d files. There will be one SourceObject per
        SpecObj in the resulting list (i.e. no combining or collating is done by this method).

        Args:
            spec1d_files (list of str): List of spec1d filenames
            match_type (str):           What type of matching the SourceObjects will be configured for.
                                        Must be either 'ra/dec' or 'pixel'
        Returns: 
            list of :obj:`SourceObject`: A list of uncollated SourceObjects with one SpecObj per SourceObject.
        """
        result = []
        for spec1d_file in spec1d_files:            
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
            spectrograph = load_spectrograph(sobjs.header['PYP_SPEC'])
            for sobj in sobjs:
                result.append(SourceObject(sobj, sobjs.header, spec1d_file, spectrograph, match_type))
    
        return result

    def _config_key_match(self, header):
        """
        Check to see if the configuration keys from a spec1d file match the
        ones for this SourceObject.

        Args:
            header (:obj:`astropy.io.fits.Header`):
                Header from a spec1d file.

        Returns:
            bool: True if the configuration keys match, 
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

        Returns:
            bool: True if the SpecObj matches this group, False otherwise.
        """

        if not self._config_key_match(spec1d_header):
            return False

        if self.match_type == 'ra/dec':
            coord2 = SkyCoord(ra=spec_obj.RA, dec=spec_obj.DEC, unit='deg')
            return self.coord.separation(coord2) <= Angle(tolerance, unit=unit)
        else:
            coord2 =spec_obj['SPAT_PIXPOS'] 
            return np.fabs(coord2 - self.coord) <= tolerance

    def combine(self, other_source_object):
        """Combine this SourceObject with another. The two objects must be from the
        same spectrograph and use the same match type.

        Args:
            other_source_object (:obj:`SourceObject`): The other object to combine with.

        Returns:
            (:obj:`SourceObject`): This SourceObject, now combined with other_source_object.
        """
        
        if other_source_object._spectrograph.name != self._spectrograph.name or \
           other_source_object.match_type != self.match_type:
           msgs.error(f"Can't append incompatible source objects. {self.spectrograph.name}/{self.match_type} does not match {other_source_object.spectrograph.name}/{other_source_object.match_type}")

        self.spec_obj_list += other_source_object.spec_obj_list
        self.spec1d_file_list += other_source_object.spec1d_file_list
        self.spec1d_header_list += other_source_object.spec1d_header_list

        return self

def collate_spectra_by_source(source_list, tolerance, unit=u.arcsec):
    """Given a list of spec1d files from PypeIt, group the spectra within the
    files by their source object. The grouping is done by comparing the 
    position of each spectra (using either pixel or RA/DEC) using a given tolerance.

    Args:
        source_list (list of :obj:`SourceObject`): A list of source objects, one
            SpecObj per object, ready for collation.
        tolerance (float): 
            Maximum distance that two spectra can be from each other to be 
            considered to be from the same source. Measured in floating
            point pixels or as an angular distance (see ``unit`` argument).
        unit (:obj:`astropy.units.Unit`):
            Units of ``tolerance`` argument if match_type is 'ra/dec'. 
            Defaults to arcseconds. Ignored if match_type is 'pixel'.

    Returns:
        list: The collated spectra as SourceObjects.

    """

    collated_list = []
    for source in source_list:

        # Search for a collated SourceObject that matches this one.
        # If one can't be found, treat this as a new collated SourceObject.
        found = False
        for collated_source in collated_list:
            if  collated_source.match(source.spec_obj_list[0], 
                                      source.spec1d_header_list[0],
                                      tolerance, unit):

                collated_source.combine(source)
                found = True

        if not found:
            collated_list.append(copy.deepcopy(source))

    return collated_list
