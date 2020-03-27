""" Simple object to hold + process a single image.
"""
import numpy as np
import inspect

from pypeit import datamodel

from IPython import embed


class DetectorContainer(datamodel.DataContainer):
    """
    Class to hold a Detector

    Args:

    Attributes:

    """
    # Set the version of this class
    version = '1.0.0'
    # Be careful.  None of these can match default FITS header cards
    datamodel_v100 = {
        'dataext': dict(otype=int, desc='Index of fits extension containing data'),
        'specaxis': dict(otype=int, desc='Spectra are dispersed along this axis. Allowed values are 0 ' \
                        '(first dimension for a numpy array shape) or 1 (second dimension ' \
                                         'for numpy array shape)'),
        'specflip': dict(otype=bool, desc='If this is True then the dispersion dimension (specificed by ' \
                        'the specaxis) will be flipped.  PypeIt expects wavelengths to ' \
                        'increase with increasing pixel number.  If this is not the case ' \
                        'for this instrument, set specflip to True.'),
        'spatflip': dict(otype=bool, desc='If this is True then the spatial dimension will be flipped.  ' \
                        'PypeIt expects echelle orders to increase with increasing pixel ' \
                        'number.  I.e., setting spatflip=True can reorder images so that ' \
                        'blue orders appear on the left and red orders on the right.'),
        'xgap': dict(otype=(int, float), desc='Gap between the square detector pixels (expressed as a fraction of the ' \
                                              'x pixel size -- x is predominantly the spatial axis)'),
        'ygap': dict(otype=(int, float), desc='Gap between the square detector pixels (expressed as a fraction of the ' \
                                              'y pixel size -- y is predominantly the spectral axis)'),
        'ysize': dict(otype=(int, float), desc='The size of a pixel in the y-direction as a multiple of the x pixel ' \
                     'size (i.e. xsize = 1.0 -- x is predominantly the dispersion axis)'),
        'platescale': dict(otype=(int, float), desc='arcsec per pixel in the spatial dimension for an unbinned pixel'),
        'darkcurr': dict(otype=(int, float), desc='Dark current (e-/hour)'),
        'saturation': dict(otype=(int, float), desc='The detector saturation level'),
        'mincounts': dict(otype=(int, float), desc='Counts in a pixel below this value will be ignored as being unphysical'),
        'nonlinear': dict(otype=(int, float), desc='Percentage of detector range which is linear (i.e. everything ' \
                         'above nonlinear*saturation will be flagged as saturated)'),
        'numamplifiers': dict(otype=int, desc='Number of amplifiers'),
        'gain': dict(otype=np.ndarray, atype=np.floating, desc='Inverse gain (e-/ADU). A list should be provided if a detector ' \
                    'contains more than one amplifier.'),
        'ronoise': dict(otype=np.ndarray, atype=np.floating, desc='Read-out noise (e-). A list should be provided if a detector ' \
                       'contains more than one amplifier.'),
        'datasec': dict(otype=np.ndarray, atype=str, desc='Either the data sections or the header keyword where the valid ' \
                       'data sections can be obtained, one per amplifier. If defined ' \
                       'explicitly should be in FITS format (e.g., [1:2048,10:4096]).'),
        'oscansec': dict(otype=np.ndarray, atype=str, desc='Either the overscan section or the header keyword where the valid ' \
                        'data sections can be obtained, one per amplifier. If defined ' \
                        'explicitly should be in FITS format (e.g., [1:2048,10:4096]).'),
        'det': dict(otype=int, desc='PypeIt designation for detector number.  1-based indexing'),
        'binning': dict(otype=str, desc='Binning in PypeIt orientation (not the original)'),
    }

    datamodel = datamodel_v100.copy()

    def __init__(self, dataext, specaxis, specflip, spatflip, platescale, saturation,
                 mincounts, nonlinear, numamplifiers, gain, ronoise, det,
                 binning,  # Up to here are required
                 xgap=None, ygap=None, ysize=None, darkcurr=None,
                 datasec=None, oscansec=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])

        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _bundle(self):
        """
        Overload for the HDU name

        Returns:
            list:

        """
        return super(DetectorContainer, self)._bundle(ext='DETECTOR')


