import numpy as np

from ginga import GingaPlugin
from ginga.AstroImage import AstroImage


class PypeItImage(AstroImage):
    """
    Custom image type for PypeIt
    """
    def __init__(self, wav_np=None, **kwargs):

        AstroImage.__init__(self, **kwargs)

        self.wav_np = wav_np

    def info_xy(self, data_x, data_y, settings):
        info = super(PypeItImage, self).info_xy(data_x, data_y, settings)

        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            _d_x, _d_y = (int(np.floor(data_x + 0.5)),
                          int(np.floor(data_y + 0.5)))

            _ht, _wd = self.wav_np.shape
            if 0 <= _d_y < _ht and 0 <= _d_x < _wd:
                # spectral wavelength is stored in auxillary array
                wavelength = self.wav_np[_d_y, _d_x]
                # choose your best formatting here...
                wav_s = "{:<14.6g}".format(wavelength)
            else:
                wav_s = ''
                
            info.update(dict(ra_lbl="\u03bb", ra_txt=wav_s,
                             dec_lbl='', dec_txt=''))

        except Exception as e:
            self.logger.error("Error getting wavelength value: {}".format(e),
                              exc_info=True)

        return info


class PypeIt(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        super(PypeIt, self).__init__(fv)

    def load_buffer(self, imname, chname, img_buf, dims, dtype,
                    header, wav_buf, wav_dtype, metadata):
        """Display a FITS image buffer.

        Parameters
        ----------
        imname : string
            a name to use for the image in Ginga
        chname : string
            channel in which to load the image
        img_buf : bytes
            the image data, as a buffer
        dims : tuple
            image dimensions in pixels (usually (height, width))
        dtype : string
            numpy data type of encoding (e.g. 'float32')
        header : dict
            fits file header as a dictionary
        wav_buf : bytes
            the wavelength data, as a buffer
        wav_dtype : string
            numpy data type of wav_buf array encoding (e.g. 'float32')
        metadata : dict
            other metadata about image to attach to image

        Returns
        -------
        0

        Notes
        -----

        * Get array dims: data.shape
        * Get array dtype: str(data.dtype)
        * Make a string from a numpy array: buf = grc.Blob(data.tobytes())

        """
        self.logger.info("received image data len=%d" % (len(img_buf)))

        # Unpack the data
        try:
            # dtype string works for most instances
            if dtype == '':
                dtype = np.float

            byteswap = metadata.get('byteswap', False)

            # unpack the auxillary wavelength file
            data = np.fromstring(wav_buf, dtype=wav_dtype)
            if byteswap:
                data.byteswap(True)
            wav_np = data.reshape(dims)

            # Create image container
            image = PypeItImage(logger=self.logger, wav_np=wav_np)
            image.load_buffer(img_buf, dims, dtype, byteswap=byteswap,
                              metadata=metadata)
            image.update_keywords(header)
            image.set(name=imname, path=None)

        except Exception as e:
            # Some kind of error unpacking the data
            errmsg = "Error creating image data for '%s': %s" % (
                imname, str(e))
            self.logger.error(errmsg)
            raise GingaPlugin.PluginError(errmsg)

        # Display the image
        channel = self.fv.gui_call(self.fv.get_channel_on_demand, chname)

        # Note: this little hack needed to let window resize in time for
        # file to auto-size properly
        self.fv.gui_do(self.fv.change_channel, channel.name)

        self.fv.gui_do(self.fv.add_image, imname, image,
                       chname=channel.name)
        return 0

    def __str__(self):
        return 'pypeit'
