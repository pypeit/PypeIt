"""
Pluging for ginga that allows the value of a wavelength map to be
displayed for any image.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from ginga import GingaPlugin
from ginga.AstroImage import AstroImage

class SlitImage(AstroImage):
    """
    Custom image type for a 2D slit image.

    This is a child of ginga's `AstroImage`_ and primarily enables
    the display and cross-image registration of the slit wavelengths.

    TODO: Also include the slit/order spatial position?

    Parameters
    ----------
    wav_np : `numpy.ndarray`_
        Wavelength map image
    **kwargs:
        Keyword arguments for `AstroImage`_
    """
    def __init__(self, wav_np=None, **kwargs):
        super(SlitImage, self).__init__(**kwargs)
        self.wav_np = wav_np

    def info_xy(self, data_x, data_y, settings):
        """
        Function to show info for data.

        Parameters
        ----------
        data_x: float
            image x coordinate
        data_y: float
            image y coordinate
        settings: :class:`ginga.misc.Settings.SettingGroup`
            ginga settings group

        Returns
        -------
        info : :class:`ginga.misc.Bunch.Bunch`
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
        wav_s = '{:<14.6g}'.format(self.wav_np[_d_y, _d_x]) \
                    if 0 <= _d_y < _ht and 0 <= _d_x < _wd else ''

        # TODO: Not sure we need this try/except block
        try:
            info.update(dict(ra_lbl="\u03bb", ra_txt=wav_s, dec_lbl='', dec_txt=''))
        except Exception as e:
            self.logger.error('Error including wavelength value: {0}'.format(e), exc_info=True)
        return info


class SlitWavelength(GingaPlugin.GlobalPlugin):
    """
    ginga plugin that enables display and registration of slit
    wavelength coordinates for a 2D slit image.

    See the documentation for `ginga GlobalPlugin`_.
    """
    def __init__(self, fv):
        super(SlitWavelength, self).__init__(fv)

    def load_buffer(self, imname, chname, img_buf, dims, dtype, header, wav_buf, wav_dtype,
                    metadata):
        """
        Load and display the 2D slit image.

        Parameters
        ----------
        imname : :obj:`str`
            Name to use for the image in Ginga
        chname : :obj:`str`
            Channel in which to load the image
        img_buf : bytes
            Image data, as a buffer
        dims : :obj:`tuple`
            Image dimensions in pixels, typically (height, width)
        dtype : :obj:`str`
            numpy data type of encoding (e.g. 'float32')
        header : :obj:`dict`
            Fits file header as a dictionary
        wav_buf : bytes
            Wavelength image data, as a buffer. Ultimate 2D shape
            once parsed must match the input image data.
        wav_dtype : :obj:`str`
            numpy data type of wav_buf array encoding (e.g. 'float32')
        metadata : :obj:`dict`
            other metadata about image to attach to image

        Returns
        -------
        status : :obj:`int`
            Load status number.  Currently always returns 0.

        Notes
        -----

        * Get array dims: data.shape
        * Get array dtype: str(data.dtype)
        * Make a string from a numpy array: buf = grc.Blob(data.tobytes())

        """
        self.logger.info("received image data len=%d" % (len(img_buf)))

        # Unpack the data
        # TODO: Should trim down what's in this try/except block
        try:
            # dtype string works for most instances
            if dtype == '':
                dtype = float

            byteswap = metadata.get('byteswap', False)

            # unpack the auxillary wavelength file
            data = np.fromstring(wav_buf, dtype=wav_dtype)
            if byteswap:
                data.byteswap(True)
            wav_np = data.reshape(dims)

            # Create image container
            image = SlitImage(wav_np=wav_np, logger=self.logger)
            image.load_buffer(img_buf, dims, dtype, byteswap=byteswap, metadata=metadata)
            image.update_keywords(header)
            image.set(name=imname, path=None)

        except Exception as e:
            # Some kind of error unpacking the data
            errmsg = 'Error creating image data for {0}: {1}'.format(imname, e)
            self.logger.error(errmsg)
            raise GingaPlugin.PluginError(errmsg)

        # Display the image
        channel = self.fv.gui_call(self.fv.get_channel_on_demand, chname)

        # Note: this little hack needed to let window resize in time for
        # file to auto-size properly
        self.fv.gui_do(self.fv.change_channel, channel.name)
        self.fv.gui_do(self.fv.add_image, imname, image, chname=channel.name)
        return 0

    def __str__(self):
        return 'slitwavelength'
