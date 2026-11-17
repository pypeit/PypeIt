"""
Ginga image type that allows the value of a wavelength map to be displayed
for any image.

.. note::

    This module used to also define the ``SlitWavelength`` ginga plugin, which
    received images sent to a separate ginga process through ginga's remote
    control interface.  PypeIt now displays images with a ginga viewer that
    runs in the same process (see :mod:`pypeit.display.viewer`), so the plugin
    is no longer needed.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from ginga.AstroImage import AstroImage


class SlitImage(AstroImage):
    """
    Custom image type for a 2D slit image.

    This is a child of ginga's `AstroImage`_ and primarily enables
    the display and cross-image registration of the slit wavelengths.

    TODO: Also include the slit/order spatial position?

    Parameters
    ----------
    wav_np : `numpy.ndarray`_, optional
        Wavelength map image.
    rms_np : `numpy.ndarray`_, optional
        Per-pixel image of the wavelength-fit RMS.  Each slit has a single
        RMS, so every pixel within a given slit holds the same value.  Must
        have the same shape as *wav_np* when provided.  Pixels outside all
        slits should be 0 or ``nan``; those values are suppressed in the
        readout.
    **kwargs:
        Keyword arguments for `AstroImage`_
    """
    def __init__(self, wav_np=None, rms_np=None, **kwargs):
        super(SlitImage, self).__init__(**kwargs)
        self.wav_np = wav_np
        self.rms_np = rms_np

    def info_xy(self, data_x, data_y, settings):
        """
        Function to show info for data.

        Parameters
        ----------
        data_x: float
            image x coordinate
        data_y: float
            image y coordinate
        settings: ``ginga.misc.Settings.SettingGroup``
            ginga settings group

        Returns
        -------
        info : ``ginga.misc.Bunch.Bunch``
            Metadata for this coordinate
        """
        info = super(SlitImage, self).info_xy(data_x, data_y, settings)
        if self.wav_np is None:
            # Allow the wavelengths to be undefined and just return the
            # same result as would be provided by AstroImage.
            return info

        # We report the value across the pixel, even though the coords
        # change halfway across the pixel
        _d_x = int(np.floor(data_x + 0.5))
        _d_y = int(np.floor(data_y + 0.5))

        # Convert wavelength value to a string. Return an empty
        # string if the provided coordinate is off the image
        _ht, _wd = self.wav_np.shape
        in_bounds = 0 <= _d_y < _ht and 0 <= _d_x < _wd
        wav_s = '{:<14.6g}'.format(self.wav_np[_d_y, _d_x]) if in_bounds else ''

        # Build the RMS readout when an RMS image has been supplied and the
        # cursor is inside a slit (rms value > 0 and finite).
        rms_lbl = ''
        rms_s = ''
        if self.rms_np is not None and in_bounds:
            rms_val = self.rms_np[_d_y, _d_x]
            if np.isfinite(rms_val) and rms_val > 0:
                rms_lbl = 'RMS'
                rms_s = '{:<10.4g}'.format(rms_val)

        # TODO: Not sure we need this try/except block
        try:
            info.update(dict(ra_lbl="\u03bb", ra_txt=wav_s, dec_lbl=rms_lbl, dec_txt=rms_s))
        except Exception as e:
            self.logger.error('Error including wavelength value: {0}'.format(e), exc_info=True)
        return info
