"""
Implements an object to handle manual object extraction.

.. include:: ../include/links.rst
"""
import inspect

from IPython import embed

import numpy as np

from pypeit import datamodel 

class ManualExtractionObj(datamodel.DataContainer):
    """
    A data container holding the arguments for how to perform the
    manual extraction of a spectrum.

    A list of these objects is generated in pypeit.py
    to perform a set of user-defined extractions.

    For an example of how to define a series of manual extractions in
    the pypeit input file, see :ref:`pypeit_file`.

    Args:
        frame (:obj:`str`):
            The name of the fits file for a manual extraction
        spat (np.ndarray): Array of spatial positions to hand extract
        spec (np.ndarray): Array of spectral positions to hand extract
        det (np.ndarray): Array of detectors for hand extraction. 
            This must be a aligned with spec and spat .
            The values can be negative (for negative images)
        fwhm (np.ndarray): Array of FWHM for hand extraction. 
            This must be aligned with spec and spat.

    """
    version = '1.0.0'

    datamodel = {
        'frame': dict(otype=str,
                    descr='The name of the fits file for a manual extraction'),
        'det': dict(otype=np.ndarray, atype=np.integer, 
                    descr='detectors for hand extraction. These can be negative for the negative image'),
        'spec': dict(otype=np.ndarray, atype=np.floating, 
                    descr='spectral positions to hand extract'),
        'spat': dict(otype=np.ndarray, atype=np.floating, 
                    descr='spatial positions to hand extract'),
        'fwhm': dict(otype=np.ndarray, atype=np.floating, 
                    descr='FWHMs for hand extractions'),
    }

    @classmethod
    def by_fitstbl_input(cls, frame: str, inp: str):
        """Generate the object from an entry in the fitstbl

        Args:
            frame (str): filename
            inp (str): det:spat:spec:fwhm
                1:1181.8:3820.6:3.

        Returns:
            ManualExtractionObj:
        """
        # Generate a dict
        idict = dict(spat=[], spec=[], det=[], fwhm=[])
        m_es = inp.split(',')
        for m_e in m_es:
            parse = m_e.split(':')
            idict['det'] += [int(parse[0])]
            idict['spat'] += [float(parse[1])]
            idict['spec'] += [float(parse[2])]
            idict['fwhm'] += [float(parse[3])]

        # Build me
        return cls(frame=frame, spat=np.array(idict['spat']), 
                   spec=np.array(idict['spec']),
                   fwhm=np.array(idict['fwhm']),
                   det=np.array(idict['det']))


    def __init__(self, frame=None, spat=None, spec=None, det=None, fwhm=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _validate(self):
        """Validate

        A couple of quick checks..

        Raises:
            ValueError: Raised if one of the arrays is not set or if they don't have the same length
        """
        if len(self.spec) != len(self.spat):
            raise ValueError("spec and spat not of the same length")
        if len(self.fwhm) != len(self.det):
            raise ValueError("FWHM and not det not of the same length")

    def dict_for_objfind(self, det, neg=False):
        """
        Repackage into a dict for the extraction code

        Args:
            det (int):
                Detector under consideration
            neg (bool, optional):
                If True, return the negative image requests

        Returns:
            dict or None: To be passed into reduce.find_objects()

        """
        # Find the ones we want
        if neg:
            gd_det = self.det == -1*det
        else:
            gd_det = self.det == det
        # None?
        if not np.any(gd_det):
            return None
        # Fill 
        manual_extract_dict = {}
        for key in ['spec', 'spat', 'det', 'fwhm']:
            manual_extract_dict[key] = self[key][gd_det]
        if neg:
            manual_extract_dict['det'] = -1 * manual_extract_dict['det']
        # Return
        return manual_extract_dict

