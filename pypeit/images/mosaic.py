"""
:class:`~pypeit.datamodel.DataContainer` object to hold mosaic properties.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
from IPython import embed

import numpy as np

from astropy import table
from astropy.io import fits

from pypeit import datamodel
from pypeit import io
from pypeit.images.detector_container import DetectorContainer
from pypeit import msgs


class Mosaic(datamodel.DataContainer):
    """
    Class to hold mosaic parameters and the details of the detectors used to
    construct the mosaic.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_mosaic.rst

    """

    # Set the version of this class
    version = '1.0.0'

    # WARNING: `binning` and `platescale` have the same names as datamodel
    # components in pypeit.images.detector_container.DetectorContainer.  This is
    # to ease use of either DetectorContainer or Mosaic objects thoughout the
    # code.  But this requires special treatment for I/O (see _parse()), so
    # beware of adding new datamodel components with the same names as those in
    # DetectorContainer!
    datamodel = {'id': dict(otype=int, descr='Mosaic ID number'),
                 'detectors': dict(otype=np.ndarray, atype=DetectorContainer,
                                   descr='List of objects with detector parameters.'),
                 'binning': dict(otype=str, descr='On-chip binning'),
                 'platescale': dict(otype=float, descr='Detector platescale in arcsec/pixel'),
                 'shape': dict(otype=tuple,
                               descr='Shape of each processed detector image'),
                 'shift': dict(otype=np.ndarray, atype=float,
                               descr='Raw, hard-coded pixel shifts for each unbinned detector'),
                 'rot': dict(otype=np.ndarray, atype=float,
                             descr='Raw, hard-coded rotations (counter-clockwise in degrees) for '
                                   'each unbinned detector'),
                 'tform': dict(otype=np.ndarray, atype=float,
                               descr='The full transformation matrix for each detector used to '
                                     'construct the mosaic.'),
                 'msc_order': dict(otype=int, descr='Order of the interpolation used to construct the mosaic.')}

    name_prefix = 'MSC'
    """
    Prefix for the name of the mosaic.
    """

    def __init__(self, id, detectors, shape, shift, rot, tform, msc_order):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])

        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _validate(self):
        """
        Validate the mosaic.
        """
        if self.detectors is not None:
            self.platescale = self.detectors[0].platescale
            self.binning = self.detectors[0].binning
            for i in range(1,self.ndet):
                if self.detectors[i].platescale != self.platescale:
                    msgs.error('Platescale difference between detectors in mosaic.')
                if self.detectors[i].binning != self.binning:
                    msgs.error('Binning difference between detectors in mosaic.')

    def _bundle(self):
        """
        Overload base class bundling so that the mosaic data is collected into
        an astropy Table with one row for each detector.

        Returns:
            :obj:`list`: List of dictionaries to write to HDUs.
        """
        # Collate all the detector parameters into a single astropy table

        # NOTE: This *requires* that the DetectorContainer _bundle method forces
        # all the detector parameters into a single row of an astropy Table.  It
        # also requires that datamodel components *cannot* be None for some
        # detectors in the mosaic and not others.
        try:
            tbl = table.vstack([d._bundle()[0]['DETECTOR'] for d in self.detectors],
                                join_type='exact')
        except:
            msgs.error('CODING ERROR: Could not stack detector parameter tables when writing '
                       'mosaic metadata.')
        if self.shift is not None:
            tbl['shift'] = self.shift
        if self.rot is not None:
            tbl['rot'] = self.rot
        if self.tform is not None:
            tbl['tform'] = self.tform
        if self.msc_order is not None:
            tbl.meta['msc_order'] = self.msc_order
        if self.id is not None:
            tbl.meta['id'] = self.id
        if self.shape is not None:
            tbl.meta['shape'] = str(self.shape)
        # Set the name
        tbl.meta['name'] = self.name
        # Keep the version of the DetectorContainer data model, needed for
        # version checking when loading files.
        tbl.meta['DETMODV'] = DetectorContainer.version
        return [{'MOSAIC':tbl}]

    @classmethod
    def _parse(cls, hdu, hdu_prefix=None, **kwargs):
        """
        Parse the data from the provided HDU.

        See :func:`pypeit.datamodel.DataContainer._parse` for the
        argument descriptions.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            hdu_prefix (:obj:`str`, optional):
                Only parse HDUs with extension names matched to this prefix. If
                None, :attr:`ext` is used. If the latter is also None, all HDUs
                are parsed. See :func:`~pypeit.io.hdu_iter_by_ext`.
            kwargs (:obj:`dict`):
                Used to collect extra keywords, but **has no affect on code
                flow**.
        """
        # Running the default parser should collect everything but the
        # DetectorContainer data
        d, version_passed, type_passed, parsed_hdus = super()._parse(hdu, hdu_prefix=hdu_prefix)

        # This should only ever read one hdu!
        if len(parsed_hdus) > 1:
            msgs.error('CODING ERROR: Parsing saved Mosaic instances should only parse 1 HDU.')

        # These are the same as the attributes for the detectors, so we need to
        # get rid of them.  We'll get them back via the _validate function.
        del d['platescale']
        del d['binning']

        # Now fake out the DetectorContainer method to get the detector metadata
        _hdu = hdu[parsed_hdus[0]] if hasattr(hdu, '__len__') else hdu
        ndet = _hdu.data.shape[0]
        tbl = table.Table(_hdu.data)
        hdr = fits.Header()
        hdr['DMODCLS'] = DetectorContainer.__name__
        hdr['DMODVER'] = _hdu.header['DETMODV']
        d['detectors'] = np.array([DetectorContainer.from_hdu(
                                        fits.BinTableHDU(data=table.Table(tbl[i]),
                                                         name='DETECTOR', header=hdr))
                                    for i in range(ndet)])

        return d, version_passed, type_passed, parsed_hdus


    @property
    def det(self):
        """
        Return a tuple with the set of detector ids in the mosaic.
        """
        return tuple() if self.detectors is None else tuple(d.det for d in self.detectors)

    @property
    def ndet(self):
        """
        Return the number of detectors in the mosaic.
        """
        return self.detectors.size

    @property
    def name(self):
        """
        Return a name for the mosaic.
        """
        return self.get_name(self.id)

    @staticmethod
    def get_id_str(mosaic_id):
        """
        Return a string identifier for the detector.
        """
        return f'{mosaic_id:02}'

    @staticmethod
    def get_name(mosaic_id):
        """
        Return a name for the mosaic.
        """
        return f'{Mosaic.name_prefix}{Mosaic.get_id_str(mosaic_id)}'

    @staticmethod
    def parse_name(name):
        """
        Parse the numerical identifier of the mosaic from its name.
        """
        return int(name[len(Mosaic.name_prefix):])



