"""
Spectrograph utility methods.
"""
from IPython import embed

import numpy as np

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
            ``PypeIt`` spectrograph.
    """
    if spec is None or isinstance(spec, spectrographs.spectrograph.Spectrograph):
        return spec

    classes = spectrographs.spectrograph_classes()
    if spec in classes.keys():
        return classes[spec]()

    msgs.error('{0} is not a supported spectrograph.'.format(spec))


