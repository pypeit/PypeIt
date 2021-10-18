import numpy as np
import inspect

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
        spat (list): List of spatial positions to hand extract
            e.g. "1243.3:1200," or "1243.3:1200,1345:1200'
        spec (list): List of spectral positions to hand extract
            e.g. "1243.3:1200," or "1243.3:1200,1345:1200'
        det (list): List of detectors for hand extraction. 
            This must be a list aligned with spec and spat lists.
        fwhm (list): List of FWHM for hand extraction. 
            This must be a list aligned with spec and spat lists.

    """
    version = '1.0.0'

    datamodel = {
        'frame': dict(otype=str,
                    descr='The name of the fits file for a manual extraction'),
        'det': dict(otype=np.ndarray, atype=np.integer, 
                    descr='detectors for hand extraction.'),
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
            ValueError: [description]
            ValueError: [description]
        """
        if len(self.spec) != len(self.spat):
            raise ValueError("spec and spat not of the same length")
        if len(self.fwhm) != len(self.det):
            raise ValueError("FWHM and not det not of the same length")

    def dict_for_objfind(self, neg=False):
        """
        Repackage into a dict for the extraction code

        Args:
            neg (bool, optional):
                Negative image

        Returns:
            dict: To be passed into reduce.find_objects()

        """
        manual_dict =  dict(hand_extract_spec=self.spec, 
                    hand_extract_spat=self.spat,
                    hand_extract_det=self.det, 
                    hand_extract_fwhm=self.fwhm)
        #
        dets = np.atleast_1d(manual_dict['hand_extract_det'])
        # Grab the ones we want
        gd_det = dets > 0
        if not neg:
            gd_det = np.invert(gd_det)
        # Any?
        if not np.any(gd_det):
            return manual_dict
        # Fill
        manual_extract_dict = {}
        for key in manual_dict.keys():
            sgn = 1
            if key == 'hand_extract_det':
                sgn = -1
            manual_extract_dict[key] = sgn*np.atleast_1d(manual_dict[key])[gd_det]
        # Return
        return manual_extract_dict