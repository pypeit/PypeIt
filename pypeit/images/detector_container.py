""" Simple object to hold + process a single image.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect

from IPython import embed

import numpy as np

from astropy.table import Table

from pypeit import datamodel
from pypeit import msgs


class DetectorContainer(datamodel.DataContainer):
    """
    Class to hold a Detector

    Args:

    Attributes:

    """
    # Set the version of this class
    version = '1.0.1'
    # Force the full datamodel into a single row of an astropy Table
    one_row_table = True
    # Be careful.  None of these can match default FITS header cards
    datamodel = {'dataext': dict(otype=int, descr='Index of fits extension containing data'),
                 'specaxis': dict(otype=int,
                                  descr='Spectra are dispersed along this axis. Allowed '
                                        'values are 0 (first dimension for a numpy array '
                                        'shape) or 1 (second dimension for numpy array '
                                        'shape).'),
                 'specflip': dict(otype=bool,
                                  descr='If this is True then the dispersion dimension '
                                        '(specified by the specaxis) will be flipped.  '
                                        'PypeIt expects wavelengths to increase with '
                                        'increasing pixel number.  If this is not the case '
                                        'for this instrument, set specflip to True.'),
                 'spatflip': dict(otype=bool,
                                  descr='If this is True then the spatial dimension will be '
                                        'flipped.  PypeIt expects echelle orders to increase '
                                        'with increasing pixel number.  I.e., setting '
                                        'spatflip=True can reorder images so that blue '
                                        'orders appear on the left and red orders on the '
                                        'right.'),
                 'xgap': dict(otype=(int, float),
                              descr='Gap between the square detector pixels (expressed as a '
                                    'fraction of the x pixel size -- x is predominantly the '
                                    'spatial axis)'),
                 'ygap': dict(otype=(int, float),
                              descr='Gap between the square detector pixels (expressed as a '
                                    'fraction of the y pixel size -- y is predominantly the '
                                    'spectral axis)'),
                 'ysize': dict(otype=(int, float),
                               descr='The size of a pixel in the y-direction as a multiple '
                                     'of the x pixel size (i.e. xsize = 1.0 -- x is '
                                     'predominantly the dispersion axis)'),
                 'platescale': dict(otype=(int, float),
                                    descr='arcsec per pixel in the spatial dimension for an '
                                          'unbinned pixel'),
                 'darkcurr': dict(otype=(int, float), descr='Dark current (e-/pixel/hour)'),
                 'saturation': dict(otype=(int, float), descr='The detector saturation level'),
                 'mincounts': dict(otype=(int, float),
                                   descr='Counts in a pixel below this value will be ignored '
                                         'as being unphysical'),
                 'nonlinear': dict(otype=(int, float),
                                   descr='Percentage of detector range which is linear '
                                         '(i.e. everything above ``nonlinear*saturation`` will '
                                         'be flagged as saturated)'),
                 'numamplifiers': dict(otype=int, descr='Number of amplifiers'),
                 'gain': dict(otype=np.ndarray, atype=np.floating,
                              descr='Inverse gain (e-/ADU). A list should be provided if a '
                                    'detector contains more than one amplifier.'),
                 'ronoise': dict(otype=np.ndarray, atype=np.floating,
                                 descr='Read-out noise (e-). A list should be provided if a '
                                       'detector contains more than one amplifier. If any '
                                       'element of this list is <=0, the readout noise will '
                                       'be determined from the overscan regions defined by '
                                       'oscansec.'),
                 'datasec': dict(otype=np.ndarray, atype=str,
                                 descr='Either the data sections or the header keyword '
                                       'where the valid data sections can be obtained, one '
                                       'per amplifier. If defined explicitly should be in '
                                       'FITS format (e.g., [1:2048,10:4096]).'),
                 'oscansec': dict(otype=np.ndarray, atype=str,
                                  descr='Either the overscan section or the header keyword '
                                        'where the valid data sections can be obtained, one '
                                        'per amplifier. If defined explicitly should be in '
                                        'FITS format (e.g., [1:2048,10:4096]).'),
                 'det': dict(otype=int,
                             descr='PypeIt designation for detector number (1-based).'),
                 'binning': dict(otype=str,
                                 descr='Binning in PypeIt orientation (not the original)')}

    def __init__(self, dataext, specaxis, specflip, spatflip, platescale, saturation,
                 mincounts, nonlinear, numamplifiers, gain, ronoise, det,
                 binning,  # Up to here are required
                 xgap=None, ygap=None, ysize=None, darkcurr=None,
                 datasec=None, oscansec=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])

        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)
        if self.darkcurr is None:
            # Use of darkcurr in RawImage means that it cannot be None.
            self.darkcurr = 0.

    def _bundle(self):
        """
        Overload base class bundling so that detector data is all forced into a
        single BinTableHDU row.

        Returns:
            :obj:`list`: List of dictionaries to write to HDUs.
        """
        return super(DetectorContainer, self)._bundle(ext='DETECTOR') #, one_row_table=True)

#    @classmethod    
#    def from_hdu(cls, hdu, **kwargs):
#        """
#        Overload the base class instantiation from an HDU to force the reader to
#        extract the data from the single-row table.
#
#        .. warning::
#            ``kwargs`` are passed directly to
#            :func:`~pypeit.datamodel.DataContainer.from_hdu`, which will raise
#            an exception if any of the keywords are not recognized.
#        """
#        if 'one_row_table' in kwargs:
#            if not kwargs['one_row_table']:
#                msgs.warn(f'{cls.__name__} must always use one_row_table=True!')
#            del kwargs['one_row_table']
#        return super(DetectorContainer, cls).from_hdu(hdu, one_row_table=True, **kwargs)

    @property
    def name(self):
        """
        Return a string identifier for the detector.
        """
        return f'DET{self.det:02}'

    @classmethod
    def parse_name(cls, name):
        """
        Parse the string identifier of the detector into its integer index.
        """
        return int(name[3:])



