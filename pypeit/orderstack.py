"""
Module for the SpecObj classes

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect

from IPython import embed

import numpy as np

from pypeit import utils
from pypeit import datamodel
from pypeit import io
from pypeit.spectrographs.util import load_spectrograph


class OrderStack(datamodel.DataContainer):
    """
    Class to handle the coadded array of orders from a single setup of an Echelle spectrum.

    One generates one of these objects for each echelle setup in the exposure. They
    are instantiated by the echelle coadding routine, and then all coadding and setup information
    for the object are added as attributes

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_orderstack.rst

    Args:
        wave
        wave_grid_mid
        flux
        PYP_SPEC
        ivar
        mask
        telluric
        obj_model
        ext_mode
        fluxed

    Attributes:
        head0 (`astropy.io.fits.Header`):
            Primary header
        spect_meta (:obj:`dict`):
            Parsed meta from the header
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Build from PYP_SPEC
    """

    version = '1.0.0'
    """
    Current datamodel version number.
    """

    datamodel = {'wave_stack': dict(otype=np.ndarray, atype=np.floating,
                              descr='Wavelength array from individual, coadded orders'),
                 'flux_stack': dict(otype=np.ndarray, atype=np.floating,
                              descr='Flux array from coadded orders, in units of counts/s or 10^-17 erg/s/cm^2/Ang; '
                                    'see ``fluxed``'),
                 'ivar_stack': dict(otype=np.ndarray, atype=np.floating,
                              descr='Inverse variance array of coadded orders (matches units of flux)'),
                 'sigma_stack': dict(otype=np.ndarray, atype=np.floating,
                              descr='One sigma noise array of coadded orders, equivalent to 1/sqrt(ivar) (matches units of flux)'),
                 'mask_stack': dict(otype=np.ndarray, atype=np.integer,
                              descr='Mask array of coadded orders (1=Good,0=Bad)'),
                 'PYP_SPEC': dict(otype=str, descr='``PypeIt`` spectrograph designation'),
                 'ext_mode': dict(otype=str, descr='Extraction mode (options: BOX, OPT)'),
                 'fluxed': dict(otype=bool, descr='Boolean indicating if the spectrum is fluxed.'),
                 'spect_meta': dict(otype=dict, descr='header dict'), 
                 'setup_name': dict(otype=str, descr='Echelle spectrograph setup'),
                 }

    internals = ['head0',
                 'filename',
                 'spectrograph',
                 'spect_meta',
                 'history']

    """
    Defines the current datmodel.
    """

    def __init__(self, wave_stack, flux_stack, ivar_stack=None, mask_stack=None, PYP_SPEC=None, sigma_stack=None, 
                 ext_mode=None, fluxed=None, setup_name = None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

    def _bundle(self):
        """
        Override the base class method simply to set the HDU extension name.
        """
        return super()._bundle(ext='SPECTRUM')

    def to_file(self, ofile, primary_hdr=None, history=None, **kwargs):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.to_file`
        to deal with the header

        Args:
            ofile (:obj:`str`): Filename
            primary_hdr (`astropy.io.fits.Header`_, optional):
            **kwargs:  Passed to super.to_file()

        """
        if primary_hdr is None:
            primary_hdr = io.initialize_header()
        # Build the header
        if self.head0 is not None and self.PYP_SPEC is not None:
            spectrograph = load_spectrograph(self.PYP_SPEC)
            subheader = spectrograph.subheader_for_spec(self.head0, self.head0,
                                                        extra_header_cards = ['RA_OBJ', 'DEC_OBJ'])
        else:
            subheader = {}
        # Add em in
        for key in subheader:
            primary_hdr[key] = subheader[key]

        # Add history
        if history is not None:
            history.write_to_header(primary_hdr)

        # Do it
        super().to_file(ofile, primary_hdr=primary_hdr, **kwargs)
