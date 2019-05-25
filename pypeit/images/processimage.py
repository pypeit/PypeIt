""" Object to hold + process a single image"""

from pypeit import msgs

from pypeit.core import procimg
from pypeit.images import pypeitimage


class ProcessImage(pypeitimage.PypeItImage):

    def __init__(self, filename, spectrograph, par, det, frametype=None):

        # Init me
        pypeitimage.PypeItImage.__init__(filename, spectrograph, det)
        # Required parameters
        self.par = par  # ProcessImagesPar

        # Optional parameters
        self.frametype = frametype

        # All possible processing steps
        self.steps = dict(bias_subtracted=False,
                          overscan_subtracted=False,
                          dark_subtracted=False,
                          trimmed=False,
                          flattened=False,
                          )

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



