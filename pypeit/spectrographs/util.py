"""
Spectrograph utility methods.
"""
from IPython import embed

import os.path

import numpy as np

from astropy.io import fits

from pypeit import spectrographs
from pypeit import msgs


def load_spectrograph(spec):
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

    Returns:
        :class:`~pypeit.spectrographs.spectrograph.Spectrograph`: The
        spectrograph used to obtain the data to be reduced.

    Raises:
        PypeItError:
            Raised if the input is a string that does not select a recognized
            PypeIt spectrograph.
    """
    if spec is None or isinstance(spec, spectrographs.spectrograph.Spectrograph):
        return spec

    classes = spectrographs.spectrograph_classes()
    if spec in classes.keys():
        return classes[spec]()

    # Check if we were given a file, and if so try to read the spectrograph type from its header
    if os.path.isfile(spec):
        header = fits.getheader(spec)
        if 'PYP_SPEC' in header:
            pyp_spec = header['PYP_SPEC']
            if pyp_spec in classes.keys():
                spectrograph = classes[pyp_spec]()
                if 'DISPNAME' in header:
                    spectrograph.dispname = header['DISPNAME']
                return spectrograph
            else:
                msgs.error(f'Unknown PYP_SPEC {pyp_spec} found in {spec}')
        else:
            msgs.error(f'{spec} did not contain PYP_SPEC in its header')
            
    msgs.error('{0} is not a supported spectrograph.'.format(spec))


