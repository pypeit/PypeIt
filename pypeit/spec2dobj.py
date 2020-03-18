"""
Module for the Spec2DObj class

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
import datetime
from IPython import embed

import numpy as np

from scipy import interpolate

from astropy import units
from astropy.io import fits

from pypeit import msgs
from pypeit import io
from pypeit import datamodel
from pypeit.images import detector_container


def spec2d_hdu_prefix(det):
    return 'DET{:02d}-'.format(det)


class Spec2DObj(datamodel.DataContainer):
    """Class to handle 2D spectral image outputs of PypeIt

    One generates one of these Objects for each detector in the exposure.

    Args:

    Attributes:
        See datamodel
        head0 (`astropy.fits.Header`):
            Primary header if instantiated from a FITS file

    """
    version = '1.0.0'

    # TODO 2d data model should be expanded to include:
    # waveimage  --  flexure and heliocentric corrections should be applied to the final waveimage and since this is unique to
    #                every exposure (i.e. it depneds on obstime, RA, DEC and the flexure incurred) it should be written out for
    #                each science frame.
    # tslits_dict -- flexure compensation implies that each frame will have a unique set of slit boundaries, so we probably need to
    #                 write these for each file as well. Alternatively we could just write the offsets to the header.

    datamodel = {
        'sciimg': dict(otype=np.ndarray, atype=np.floating, desc='2D processed science image'),
        'ivarraw': dict(otype=np.ndarray, atype=np.floating, desc='2D processed inverse variance image'),
        'skymodel': dict(otype=np.ndarray, atype=np.floating, desc='2D sky model image'),
        'objmodel': dict(otype=np.ndarray, atype=np.floating, desc='2D object model image'),
        'ivarmodel': dict(otype=np.ndarray, atype=np.floating, desc='2D ivar model image'),
        'mask': dict(otype=np.ndarray, atype=np.integer, desc='2D mask image'),
        'detector': dict(otype=detector_container.DetectorContainer, desc='Detector DataContainer'),
        'det': dict(otype=int, desc='Detector index'),
    }

    @classmethod
    def from_file(cls, file, det):
        hdul = fits.open(file)
        slf = super(Spec2DObj, cls).from_hdu(hdul, hdu_prefix=spec2d_hdu_prefix(det))
        slf.head0 = hdul[0].header
        return slf

    def __init__(self, det, sciimg, ivarraw, skymodel, objmodel, ivarmodel, mask, detector):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)


    def _init_internals(self):
        """
        All modifiable internals go here
        """
        self.head0 = None

    def _vaildate(self):
        """
        Assert that the detector has been set

        Returns:

        """
        assert self.det is not None, 'Must set det at instantiation!'

    def _bundle(self):
        """
        Over-write default _bundle() method to separate the DetectorContainer
        into its own HDU

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """
        d = []
        # Rest of the datamodel
        for key in self.keys():
            # Skip Nones
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray:
                tmp = {}
                tmp[key] = self[key]
                d.append(tmp)
            # Deal with the detector
            elif key == 'detector':
                d.append(dict(detector=self.detector))
            else: # Add to header of the primary image
                d[0][key] = self[key]
        # Return
        return d

    @property
    def hdu_prefix(self):
        """
        Provides for a dynamic hdu_prefix based on our naming model

        see :func:`spec2d_hdu_prefix`

        Returns:
            str

        """
        return spec2d_hdu_prefix(self.det)


class AllSpec2DObj(object):
    """
    Simple object to hold Spec2DObj objects
    and perform I/O

    Anything that goes into self['meta'] must be parseable into a FITS Header

    Restrict keys to be type int or 'meta'
    and items to be :class:`Spec2DObj`

    """
    hdr_prefix = 'ALLSPEC2D_'
    @classmethod
    def from_fits(cls, filename):
        """

        Args:
            filename (str):

        Returns:
            :class:`AllSpec2DObj`:

        """
        # Instantiate
        slf = cls()
        # Open
        hdul = fits.open(filename)
        # Meta
        hkeys = list(hdul[0].header.keys())
        for key in hkeys:
            if slf.hdr_prefix in key:
                slf['meta'][key.split(slf.hdr_prefix)[-1]] = hdul[0].header[key]
        # Detectors included
        detectors = hdul[0].header[slf.hdr_prefix+'DETS']
        for det in [int(item) for item in detectors.split(',')]:
            obj = Spec2DObj.from_hdu(hdul, hdu_prefix=spec2d_hdu_prefix(det))
            slf[det] = obj
        # Header
        slf['meta']['head0'] = hdul[0].header
        return slf

    def __init__(self):
        # Init meta
        self['meta'] = {}

    # TODO -- Turn off attribute setting

    def keys(self):
        return self.__dict__.keys()

    def __setitem__(self, item, value):
        # Check item
        if not isinstance(item, int) and item != 'meta':
            raise KeyError('Key must be an integer, i.e. detector number or "meta"')
        # Check value
        if isinstance(item, int):
            assert isinstance(value, Spec2DObj), 'Item must be a Spec2DObj'
        self.__dict__[item] = value

    def __getitem__(self, item):
        """Get an item directly from the internal dict."""
        return self.__dict__[item]

    def build_primary_hdr(self, raw_header, spectrograph, master_key_dict=None, master_dir=None):
        hdr = io.initialize_header(primary=True)

        hdukeys = ['BUNIT', 'COMMENT', '', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'HISTORY', 'EXTEND', 'DATASEC']
        for key in raw_header.keys():
            # Use new ones
            if key in hdukeys:
                continue
            # Update unused ones
            hdr[key] = raw_header[key]
        # History
        if 'HISTORY' in raw_header.keys():
            # Strip \n
            tmp = str(raw_header['HISTORY']).replace('\n', ' ')
            hdr.add_history(str(tmp))

        # PYPEIT
        # TODO Should the spectrograph be written to the header?
        hdr['PIPELINE'] = str('PYPEIT')
        hdr['PYPELINE'] = spectrograph.pypeline
        hdr['SPECTROG'] = spectrograph.spectrograph
        hdr['DATE-RDX'] = str(datetime.date.today().strftime('%Y-%b-%d'))
        # MasterFrame info
        # TODO -- Should this be in the header of the individual HDUs ?
        if master_key_dict is not None:
            hdr['FRAMMKEY'] = master_key_dict['frame'][:-3]
            hdr['BPMMKEY'] = master_key_dict['bpm'][:-3]
            hdr['BIASMKEY'] = master_key_dict['bias'][:-3]
            hdr['ARCMKEY'] = master_key_dict['arc'][:-3]
            hdr['TRACMKEY'] = master_key_dict['trace'][:-3]
            hdr['FLATMKEY'] = master_key_dict['flat'][:-3]
        if master_dir is not None:
            hdr['PYPMFDIR'] = str(master_dir)
        # Sky sub mode
        if self['meta']['ir_redux']:
            hdr['SKYSUB'] = 'DIFF'
        else:
            hdr['SKYSUB'] = 'MODEL'
        #
        return hdr

    def write_to_fits(self, outfile, pri_hdr=None, update_det=None, overwrite=True):

        # Start the HDU list (or read from disk) and Primary header
        if os.path.isfile(outfile) and update_det is not None:
            hdus, prihdu = io.init_hdus(update_det, outfile)
        else:
            # Primary HDU for output
            prihdu = fits.PrimaryHDU()
            # Header
            if pri_hdr is not None:
                prihdu.header = pri_hdr
            # Update with original header, skipping a few keywords
            hdus = [prihdu]

        # Add meta to Primary Header
        for key in self['meta']:
            prihdu.header[self.hdr_prefix+key.upper()] = self['meta'][key]

        # Loop on em (in order of detector)
        keys = list(self.keys())
        keys.remove('meta')
        try:
            keys.sort()
        except TypeError:
            embed(header='244 of spec2dobj')
        extnum = 1
        for key in keys:
            #
            hdul = self[key].to_hdu()
            # TODO -- Make adding EXT000X a default of DataContainer
            for hdu in hdul:
                keywd = 'EXT{:04d}'.format(extnum)
                prihdu.header[keywd] = hdu.name
                extnum += 1
            # Add em in
            hdus += hdul

        # Detectors included
        detectors = str(keys)[1:-1]
        prihdu.header[self.hdr_prefix+'DETS'] = detectors

        # Finish
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(outfile, overwrite=overwrite)
        msgs.info("Wrote: {:s}".format(outfile))


