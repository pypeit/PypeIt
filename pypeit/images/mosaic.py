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
    """

    # Set the version of this class
    version = '1.0.0'

    # WARNING: None of the keywords in this datamodel should be the same as in
    # pypeit.images.detector_container.DetectorContainer.
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
                                     'construct the mosaic.')}

    def __init__(self, id, detectors, shape, shift, rot, tform):

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
    def from_hdu(cls, hdu, chk_version=True, **kwargs):
        """
        Instantiate the object from an HDU extension.

        This is primarily a wrapper for :func:`_parse`.

        Args:
            hdu (`astropy.io.fits.HDUList`_, `astropy.io.fits.ImageHDU`_, `astropy.io.fits.BinTableHDU`_):
                The HDU(s) with the data to use for instantiation.
            chk_version (:obj:`bool`, optional):
                If True, raise an error if the datamodel version or
                type check failed. If False, throw a warning only.
            kwargs (:obj:`dict`):
                Used to collect extra keywords, but **has no affect on code
                flow**.
        """
        # TODO: Unsure about the generality of this for additional code
        # development...  It may cause issues if the hdu has many extensions.
        # Specifically, it won't work if the hdu has multiple mosaic extensions.
        # Find the mosaic extension
        _ext, _hdu = io.hdu_iter_by_ext(hdu)
        _ext = np.atleast_1d(np.array(_ext, dtype=object))
        ext_to_use = None
        for e in _ext:
            if 'DMODCLS' in hdu[e].header and hdu[e].header['DMODCLS'] == cls.__name__:
                ext_to_use = e
                break
        if ext_to_use is None:
            msgs.error('The HDU(s) cannot be parsed by a {0} object!'.format(cls.__name__))

        # Parsing this extension alone should collect everything but the
        # DetectorContainer data
        d, dm_version_passed, dm_type_passed, parsed_hdus = cls._parse(hdu[ext_to_use])

        # Now fake out the DetectorContainer method to get the detector metadata
        ndet = hdu[ext_to_use].data.shape[0]
        tbl = table.Table(hdu[ext_to_use].data)
        hdr = fits.Header()
        hdr['DMODCLS'] = DetectorContainer.__name__
        hdr['DMODVER'] = hdu[ext_to_use].header['DETMODV']
        d['detectors'] = np.array([DetectorContainer.from_hdu(
                                        fits.BinTableHDU(data=table.Table(tbl[i]),
                                                         name='DETECTOR', header=hdr))
                                    for i in range(ndet)])

        # Check version and type?
        if not dm_type_passed:
            msgs.error('The HDU(s) cannot be parsed by a {0} object!'.format(cls.__name__))
        if not dm_version_passed:
            _f = msgs.error if chk_version else msgs.warn
            _f('Current version of {0} object in code (v{1})'.format(cls.__name__, cls.version)
               + ' does not match version used to write your HDU(s)!')

        # Instantiate
        return cls.from_dict(d=d)

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
        return f'MSC{Mosaic.get_id_str(mosaic_id)}'

    @staticmethod
    def parse_name(name):
        """
        Parse the numerical identifier of the mosaic from its name.
        """
        return int(name[3:])



