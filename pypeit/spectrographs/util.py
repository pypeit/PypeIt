"""
Spectrograph utility methods.
"""
from IPython import embed

import pathlib

from astropy.io import fits

from pypeit import spectrographs
from pypeit import msgs


def load_spectrograph(
    spec:str|spectrographs.spectrograph.Spectrograph, pypeit_fits:bool=False
) -> spectrographs.spectrograph.Spectrograph:
    """
    Instantiate a spectrograph from the available subclasses of
    :class:`~pypeit.spectrographs.spectrograph.Spectrograph`.

    Args:
        spec (:obj:`str`, :class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph to instantiate. If the input object is ``None``
            or has :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
            as a base class, the instance is simply returned. If it is a
            string, the string is used to instantiate the relevant
            spectrograph instance.
        pypeit_fits (:obj:`bool`, optional):
            The spectrograph loader is being called from a post-processing
            script where the expected input files are PypeIt-written FITS files
            only.  This has the effect of overriding the :attr:`allowed_extensions`
            attribute to be ``[".fits"]``.

    Returns:
        :class:`~pypeit.spectrographs.spectrograph.Spectrograph`: The
        spectrograph used to obtain the data to be reduced.

    Raises:
        PypeItError:
            Raised if the input is a string that does not select a recognized
            PypeIt spectrograph.
    """
    # If given None, return None
    if spec is None:
        return None

    # If provided a full spectrograph class, update if necessary and return
    if isinstance(spec, spectrographs.spectrograph.Spectrograph):
        if pypeit_fits:
            spec.allowed_extensions = [".fits"]
        return spec

    # The function was provided the name of a spectrograph; return the class
    classes = spectrographs.spectrograph_classes()
    if spec in classes.keys():
        spectrograph = classes[spec]()
        if pypeit_fits:
            spectrograph.allowed_extensions = [".fits"]
        return spectrograph

    # Check if we were given a file, and if so try to read the spectrograph type from its header
    if pathlib.Path(spec).is_file():
        header = fits.getheader(spec)
        if 'PYP_SPEC' in header:
            pyp_spec = header['PYP_SPEC']
            if pyp_spec in classes.keys():
                spectrograph = classes[pyp_spec]()
                if 'DISPNAME' in header:
                    spectrograph.dispname = header['DISPNAME']
                if pypeit_fits:
                    spectrograph.allowed_extensions = [".fits"]
                return spectrograph
            else:
                msgs.error(f'Unknown PYP_SPEC {pyp_spec} found in {spec}')
        else:
            msgs.error(f'{spec} did not contain PYP_SPEC in its header')
            
    msgs.error('{0} is not a supported spectrograph.'.format(spec))


