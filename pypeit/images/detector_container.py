"""
:class:`~pypeit.datamodel.DataContainer` object to hold detector properties.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect

from IPython import embed

import numpy as np

from pypeit import datamodel
from pypeit import msgs
from pypeit.core import procimg


class DetectorContainer(datamodel.DataContainer):
    """
    Class to hold a detector properties.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_detectorcontainer.rst

    """
    # Set the version of this class
    version = '1.0.1'
    # Force the full datamodel into a single row of an astropy Table
    one_row_table = True
    # Be careful.  None of these can match default FITS header cards
    datamodel = {'dataext': dict(otype=int,
                                 descr='Index of fits extension containing data'),
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
                # TODO: There are actually two types of "saturation": (1) the
                # point at which the amplifier A/D converter reaches the upper
                # limit of the bit representation (e.g., 65535 = 2**16-1) or (2)
                # the point at which the electron well depth is filled.  The
                # reason this matters is that detection of the former should be
                # done using the ADU/DN value in the *raw* frame --- before
                # subtracting the bias, applying the gain, etc. --- and
                # detection of the latter should be done using counts in the
                # *bias-subtracted* frame.  Looking across all our instruments,
                # it looks like we're mixing how we define this number...
                 'saturation': dict(otype=(int, float),
                                    descr='The detector saturation level in ADU/DN'),
                 'mincounts': dict(otype=(int, float),
                                   descr='Counts (e-) in a pixel below this value will be ignored '
                                         'as being unphysical.'),
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

    name_prefix = 'DET'
    """
    Prefix for the name of the detector.
    """

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
        Overload base class bundling to select appropriate extension name.

        Returns:
            :obj:`list`: List of dictionaries to write to HDUs.
        """
        # Run the base-level method, which will force all the components of the
        # datamodel into a single-row table.
        d = super()._bundle(ext='DETECTOR')
        # Add the name to the table metadata so that it gets added to the header
        # of the detector extension
        d[0]['DETECTOR'].meta['name'] = self.name
        return d

    @property
    def name(self):
        """
        Return a string identifier for the detector.  This is a simple wrapper
        for :func:`get_name` using :attr:`det`.
        """
        return self.get_name(self.det)

    @staticmethod
    def get_det_str(det):
        """
        Return a string identifier for the detector.  Currently a zero-padded
        two character string with the detector number.

        Args:
            det (:obj:`int`):
                1- indexed detector number.

        Returns:
            :obj:`str`: String representation of the detector number used in
            constructing the detector name.
        """
        return f'{det:02}'

    @staticmethod
    def get_name(det):
        """
        Return a string identifier for the detector.  Currently, e.g., DET01 for
        det=1.

        Args:
            det (:obj:`int`):
                1- indexed detector number.

        Returns:
            :obj:`str`: Detector name.
        """
        return f'{DetectorContainer.name_prefix}{DetectorContainer.get_det_str(det)}'

    @staticmethod
    def parse_name(name):
        """
        Parse the string identifier of the detector into its integer index.

        Args:
            name (:obj:`str`):
                Detector name.  Assumed to have been created by
                :func:`get_name`.

        Returns:
            :obj:`int`: The parsed detector number.  For example, returns 2 when
            the name is ``DET02``.
        """
        return int(name[len(DetectorContainer.name_prefix):])

    def nonlinear_counts(self, datasec_img=None, apply_gain=True):
        """
        Return the ADU/DN or counts at which the detector response becomes
        non-linear.

        Args:
            datasec_img (`numpy.ndarray`_, optional):
                An image identifying the amplifier used to read each pixel in
                the detector data section.  If provided, the returned object is
                an image giving the non-linear counts for each pixel.
            apply_gain (:obj:`bool`, optional):
                Apply gain in the calculation. I.e., convert the value to
                counts. If only a float is returned, (i.e. ``datasec_img`` is
                not provided), the mean of the gains for all amplifiers is
                used.

        Returns:
            :obj:`float`, `numpy.ndarray`_: Counts at which the detector
            response becomes nonlinear. If ``datasec_img`` is provided, an
            image of the same shape is returned with the pixel-specific
            nonlinear-count threshold.
        """
        if apply_gain:
            gain = np.mean(self.gain) if datasec_img is None \
                    else procimg.gain_frame(datasec_img, self.gain)
        else:
            gain = 1.
        return self.saturation * self.nonlinear * gain



