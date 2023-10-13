"""
Implements the objects used to hold slit edge data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import datamodel
from pypeit import calibframe
from pypeit.display import display


class ScatteredLight(calibframe.CalibFrame):
    """
    Defines a generic class for generating a model of the scattered light.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_scatteredlight.rst

    """
    calib_type = 'ScatteredLight'
    """Name for type of calibration frame."""

    calib_file_format = 'fits.gz'
    """File format for the calibration frame file."""

    version = '1.0.0'
    """Scattered Light data model version."""

    internals = calibframe.CalibFrame.internals
    """
    Attributes kept separate from the datamodel.
    """

    # Define the data model
    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'pypeline': dict(otype=str, descr='PypeIt pypeline name'),
                 'detname': dict(otype=str, descr='Identifier for detector or mosaic'),
                 'nspec': dict(otype=int, descr='Number of pixels in the image spectral direction.'),
                 'nspat': dict(otype=int, descr='Number of pixels in the image spatial direction.'),
                 'binning': dict(otype=str, descr='Binning in PypeIt orientation (not the original)'),
                 'pad': dict(otype=int, descr='Integer number of pixels to mask beyond the slit edges.'),
                 'scattlight_raw': dict(otype=np.ndarray, atype=np.floating, descr='Processed, combined scattered light image'),
                 'scattlight_model': dict(otype=np.ndarray, atype=np.floating, descr='Model of the scattered light in scattlight_raw'),
                 'scattlight_param': dict(otype=np.ndarray, atype=np.floating, descr='Model parameters that define the scattered light model')}
    """Provides the class data model."""

    # TODO: Allow tweaked edges to be arguments?
    # TODO: May want nspat to be a required argument.
    # The INIT must contain every datamodel item or risk fail on I/O when it is a nested container
    def __init__(self, pypeline=None, detname=None, nspec=None, nspat=None, PYP_SPEC=None, binning=None, pad=0,
                 scattlight_raw=None, scattlight_model=None, scattlight_param=None):

        # Instantiate the DataContainer
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # The dictionary passed to DataContainer.__init__ does not
        # contain self.
        # TODO: Does it matter if the calling function passes the
        # keyword arguments in a different order? No.
        datamodel.DataContainer.__init__(self, d=_d)

    def _validate(self):
        """
        Validate the slit traces.
        """
        # Allow the object to be empty
        if self.scattlight_model is None:
            return

        self.nspec, self.nslits = self.left_init.shape

        if self.PYP_SPEC is None:
            self.PYP_SPEC = 'unknown'

    def get_model(self, image):
        """Generate a model of the scattered light, based on an input image. This routine
        requires that the scattered light has already been predefined.

        Parameters
        ----------
        image : `numpy.ndarray`_
            A 2D image that you want to determine the amount of scattered light

        Returns
        -------
        model : `numpy.ndarray`_
            A model of the expected scattered light in the input image. Shape is (nspec, nspat).
        """
        msgs.info("Generating a scattered light image")
        if self.scattlight_param is None:
            msgs.warn("No scattered light parameters are available")
            return np.zeros_like(image)
        # Return the model of the scattered light
        return self.spec.scattered_light_model(self.scattlight_param, image)

    def show(self, image=None, slits=None, wcs_match=True):
        """ Display the master scattered light frame, the model, and data-model.

        Parameters
        ----------
        image : `numpy.ndarray`_, optional
            A 2D image that you want to display (including a scattered light model). If None,
            the master Scattered Light frame wil be displayed by default
        slits : :class:`~pypeit.slittrace.SlitTraceSet`, optional
            The current slit traces
        wcs_match : :obj:`bool`, optional
            Use a reference image for the WCS and match all image in other channels to it.
        """
        # Prepare the frames
        _data = self.scattlight_raw if image is None else image
        _model = self.scattlight_model if image is None else self.get_model(image)
        _resid = _data - _model
        modmax = np.max(_model)
        image_list = zip([_data, _model, _resid],
                             ['data', 'model', 'data-model'],
                             [(0.0, modmax), (0.0, modmax), (-modmax/2,modmax/2)])

        # Display frames
        show_scattered_light(image_list, slits=slits, wcs_match=wcs_match)


def show_scattered_light(image_list, slits=None, wcs_match=True):
    """
    Interface to ginga to show the quality of the Scattered Light subtraction

    Parameters
    ----------
    image_list : zip
        A zip of the images to show, their names, and the scales
    slits : :class:`~pypeit.slittrace.SlitTraceSet`, optional
        The current slit traces
    wcs_match : :obj:`bool`, optional
        Use a reference image for the WCS and match all image in other channels to it.
    """
    display.connect_to_ginga(raise_err=True, allow_new=True)
    if slits is not None:
        left, right, mask = slits.select_edges()
        gpm = mask == 0
    # Loop me
    clear = True
    for img, name, cut in image_list:
        if img is None:
            continue
        viewer, ch = display.show_image(img, chname=name, cuts=cut, wcs_match=wcs_match, clear=clear)
        if slits is not None:
            display.show_slits(viewer, ch, left[:, gpm], right[:, gpm], slit_ids=slits.spat_id[gpm])
        # Turn off clear
        if clear:
            clear = False
