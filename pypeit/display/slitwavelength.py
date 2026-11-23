"""
Plug-in for ginga that allows the value of a wavelength map to be
displayed for any image.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from ginga import GingaPlugin

from pypeit.display.ginga_draw import SlitImage


class SlitWavelength(GingaPlugin.GlobalPlugin):
    """
    ginga plugin that enables display and registration of slit
    wavelength coordinates for a 2D slit image.

    See the documentation for `ginga GlobalPlugin`_.
    """
    def __init__(self, fv):
        super(SlitWavelength, self).__init__(fv)

    def load_buffer(self, imname, chname, img_buf, dims, dtype, header, wav_buf, wav_dtype,
                    metadata, rms_buf=None, rms_dtype=None):
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
        rms_buf : bytes, optional
            Per-pixel wavelength-fit RMS data, as a buffer.  Each slit has a
            single RMS, so every pixel within a given slit holds the same
            value.  Must have the same length as *wav_buf* once interpreted
            as *rms_dtype*.  When ``None`` (default) the DEC readout field is
            left blank, preserving the original behaviour.
        rms_dtype : :obj:`str`, optional
            numpy data type of *rms_buf* (e.g. ``'float32'``).  Required when
            *rms_buf* is provided; ignored otherwise.

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
            data = np.frombuffer(wav_buf, dtype=wav_dtype)
            if byteswap:
                data.byteswap(True)
            wav_np = data.reshape(dims)

            # Unpack the optional per-pixel RMS image
            rms_np = None
            if rms_buf is not None and rms_dtype is not None:
                rms_data = np.frombuffer(rms_buf, dtype=rms_dtype)
                if byteswap:
                    rms_data.byteswap(True)
                rms_np = rms_data.reshape(dims)

            # Create image container
            image = SlitImage(wav_np=wav_np, rms_np=rms_np, logger=self.logger)
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
