""" Simple object to hold + process a single image.
"""
import numpy as np
import os

from astropy.io import fits

from pypeit import msgs
from pypeit import ginga

from pypeit.images import maskimage
from pypeit.images import detector_container
from pypeit import datamodel

from IPython import embed


class PypeItImage(datamodel.DataContainer):
    """
    Class to hold a single image from a single detector in PypeIt
    and its related images (e.g. ivar, mask).

    Oriented in its spec,spat format

    The intent is to keep this object as light-weight as possible.

    PypeItImage has an optional, internal attribute `mask` which is intended to hold
    a maskimage.ImageMask DataContainer

    Args:
        image (np.ndarray):
        ivar (np.ndarray, optional):
        rn2img (np.ndarray, optional):
        binning (tuple, optional):
        bpm (np.ndarray, optional):
            Passed to self.mask
        crmask (np.ndarray, optional):
            Passed to self.mask
        fullmask (np.ndarray, optional):
            Passed to self.mask

    Attributes:
        mask (class:`pypeit.images.maskimage.ImageMask`):
            Image mask(s)
        hdu_prefix (str, optional):
            Appended to the HDU name, if provided.
            Mainly used to enable output of multiple PypeItImage objects
            in more complex DataContainers
        head0 (astropy.io.fits.Header):
        detector_par (:class:`pypeit.par.pypeitpar.DetectorPar`):

    """
    # Set the version of this class
    version = '1.0.0'
    #
    datamodel_v100 = {
        'image': dict(otype=np.ndarray, atype=np.floating, desc='Main data image'),
        'ivar': dict(otype=np.ndarray, atype=np.floating, desc='Main data inverse variance image'),
        'rn2img': dict(otype=np.ndarray, atype=np.floating, desc='Read noise squared image'),
        'BIN_SPEC': dict(otype=(int, np.integer), desc='Binning in spectral dimension'),
        'BIN_SPAT': dict(otype=(int, np.integer), desc='Binning in spatial dimension'),
        'HEAD0': dict(otype=fits.header.Header, desc='Image header of primary HDU'),
        'mask': dict(otype=maskimage.ImageMask, desc='Mask DataContainer'),
        'detector': dict(otype=detector_container.Detector, desc='Detector DataContainer'),
    }

    datamodel = datamodel_v100.copy()

    @classmethod
    def from_file(cls, file, hdu_prefix=None):
        """
        Instantiate from a file on disk (FITS file)

        Overloading DataContainer method to deal with mask

        Args:
            file (str):
            hdu_prefix (str, optional):

        Returns:
            :class:`pypeit.images.pypeitimage.PypeItImage`:
                Loaded up PypeItImage with the primary Header attached

        """
        # Open
        hdul = fits.open(file)

        slf = super(PypeItImage, cls).from_hdu(hdul, hdu_prefix=hdu_prefix)

        # Header
        slf.HEAD0 = hdul[0].header
        # Mask
        slf.mask = maskimage.ImageMask.from_hdu(hdul, hdu_prefix=hdu_prefix)

        # Return
        return slf

    @classmethod
    def from_pypeitimage(cls, pypeitImage):
        slf = cls(pypeitImage.image, ivar=pypeitImage.ivar, rn2img=pypeitImage.rn2img,
                  binning=pypeitImage.binning)
        if pypeitImage.mask is not None:
            slf.mask = maskimage.ImageMask(pypeitImage.mask.bpm, crmask=pypeitImage.mask.crmask,
                                   fullmask=pypeitImage.mask.fullmask)
        # Return
        return slf

    def __init__(self, image, ivar=None, rn2img=None, bpm=None,
                 binning=None, crmask=None, fullmask=None, prefix=None,
                 detector=None):

        self.prefix = prefix

        # Setup the DataContainer
        super(PypeItImage, self).__init__({'image': image, 'ivar': ivar, 'rn2img': rn2img,
                                          'mask': maskimage.ImageMask(bpm, crmask=crmask, fullmask=fullmask),
                                          'detector': detector}
                                          )

        # Internals need to come after
        self.binning = binning

    def _init_internals(self):

        #self.mask = None
        self.binning = None


    def _bundle(self):
        """
        Over-write default _bundle() method to write one
        HDU per image.  Any extras are in the HDU header of
        the primary image.

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """
        d = []
        # Primary image
        d.append(dict(image=self.image))

        # Rest of the datamodel
        for key in self.keys():
            if key in ['image', 'HEAD0']:
                continue
            # Skip None
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray:
                tmp = {}
                tmp[key] = self[key]
                d.append(tmp)
            # Deal with the mask
            elif key == 'mask':
                d.append(dict(mask=self.mask))
            else: # Add to header of the primary image
                d[0][key] = self[key]
        # Return
        return d

    @property
    def shape(self):
        return () if self.image is None else self.image.shape

    '''
    def to_file(self, ofile, overwrite=True, checksum=True, primary_hdr=None, hdr=None, hdu_prefix=None):
        """ Overload to_file to add in the mask
        """
        # Get PypeitImage hdul
        hdul = self.to_hdu(primary_hdr=primary_hdr, add_primary=True, hdu_prefix=hdu_prefix)

        ## Mask HDU
        #if self.mask is not None:
        #    mask_hdul = self.mask.to_hdu(hdu_prefix=hdu_prefix)
        #else:
        #    mask_hdul = []

        ## Combine
        #for ihdu in mask_hdul:
        #    hdul.append(ihdu)

        # Write
        write_to_fits(hdul, ofile, overwrite=overwrite, checksum=checksum, hdr=hdr)
    '''


    def show(self):
        """
        Show the image in a ginga viewer.
        """
        if self.image is None:
            # TODO: This should fault.
            msgs.warn("No image to show!")
            return
        ginga.show_image(self.image, chname='image')

    def __repr__(self):
        repr = '<{:s}: '.format(self.__class__.__name__)
        # Image
        rdict = {}
        for attr in ['image', 'ivar', 'rn2img']:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                rdict[attr] = True
            else:
                rdict[attr] = False
        repr += ' images={}'.format(rdict)
        repr = repr + '>'
        return repr

