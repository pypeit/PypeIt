""" Object to hold + process a single image"""

from pypeit import msgs

from pypeit.core import procimg


class SingleImage(object):

    def __init__(self, filename, spectrograph, par, det, frametype=None):

        # Required parameters
        self.spectrograph = spectrograph
        self.par = par  # ProcessImagesPar
        self.det = det
        self.filename = filename

        # Optional parameters
        self.frametype = frametype

        # Attributes
        self.disk_image = None       # Image loaded from disk, typically "raw"
        self.header = None           # Image header
        self.processed_image = None  # Image processed in one or more ways

        # All possible process steps
        self.steps = dict(bias_subtracted=False,
                          overscan_subtracted=False,
                          dark_subtracted=False,
                          trimmed=False,
                          flattened=False,
                          )

    @property
    def bpm(self):
        if self.processed_image is None:
            bpm = None
        else:
            bpm = self.spectrograph.bpm(shape=self.processed_image.shape,
                                        filename=self.filename,
                                        det=self.det)
        return bpm

    @property
    def datasec_img(self):
        dimg = self.spectrograph.get_datasec_img(self.filename, self.det)
        return dimg

    @property
    def oscansec_img(self):
        oimg = self.spectrograph.get_oscansec_img(self.filename, self.det)
        return oimg

    def load(self):
        """
        Load the image from disk

        Returns:

        """
        self.disk_image, self.header \
            = self.spectrograph.load_raw_frame(self.filename, det=self.det)
        #
        return self.disk_image, self.header

    def reset_steps(self):
        for key in self.steps.keys():
            self.steps[key] = False

    def subtract_overscan(self, image, force=False):
        # Checks
        assert image.shape == self.disk_image.shape

        numamplifiers = self.spectrograph.detector[self.det-1]['numamplifiers']
        temp = procimg.subtract_overscan(image, numamplifiers, self.datasec[kk],
                                         self.oscansec[kk],
                                         method=self.par['overscan'],
                                         params=self.par['overscan_par'])

        # Return
        return self.processed_image

    def trim(self, image, force=False):
        # Checks
        assert image.shape == self.disk_image.shape
        if self.steps['trimmed'] and (not force):
            msgs.warn("Image was already trimmed")
        # Do it
        trim_image = procimg.trim_frame(image, self.datasec_img < 1)
        # Save
        self.processed_image = trim_image
        self.steps['trimmed'] = True
        # Return
        return self.processed_image



