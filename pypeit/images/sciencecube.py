""" Object to hold + process a single image.
This module also includes the build_from_list() method
which is how the ScienceImage is most frequently generated. """

import numpy as np

# NOTE: This is currently the only use of shapely in pypeit! See
# ScienceCube.calculate_area().
try:
    from shapely.geometry import Polygon
except:
    Polygon = None

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.images import pypeitimage

from IPython import embed


class ScienceCube(pypeitimage.PypeItImage):
    """
    Class to generate and hold a science cube

    Child of PypeItImage

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.
        det (:obj:`int`):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.ProcessImagesPar`):
            Parameters that dictate the processing of the images.  See
            :class:`pypeit.par.pypeitpar.ProcessImagesPar` for the
            defaults.
        image (np.ndarray):
        ivar (np.ndarray):
        bpm (np.ndarray):
            Bad pixel mask.  Held in ImageMask
        rn2img (np.ndarray, optional):
        crmask (np.ndarray, optional):
        mask (np.ndarray, optional):
        files (list, optional):
            List of filenames that went into the loaded image

    """
    frametype = 'cube'

    def __init__(self, spectrograph, det, par, cube, ivar, bpm, rn2img=None,
                 crmask=None, mask=None, files=None):

        # Init me
        pypeitimage.PypeItImage.__init__(self, cube, ivar=ivar, rn2img=rn2img,
                                         bpm=bpm, crmask=crmask, fullmask=mask)

        if files is None:
            files = []

        # Required attribs
        self.spectrograph = spectrograph
        if not isinstance(par, pypeitpar.CubePar):
            msgs.error('Provided ParSet for must be type CubePar.')
        self.par = par
        self.det = det

        # Not required
        self.files = files

    def calculate_area(self):
        """Calculate the overlapping area of a pixel and voxel"""
        if Polygon is None:
            msgs.error('Shapley python package is required to calculate spaxel area.')
        voxel = Polygon([(0.1, 0.2), (1.2, 1.23), (2.5, 0.8), (0.5, 0.12)])
        pixel = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
        area = voxel.intersection(pixel).area


