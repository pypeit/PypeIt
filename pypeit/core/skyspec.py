from pypeit import dataPaths
from pypeit import onespec


def load_sky_spectrum(sky_file: str) -> onespec.OneSpec:
    """
    Load a sky spectrum from the PypeIt data directory into an XSpectrum1D
    object.

    .. todo::

        Try to eliminate the XSpectrum1D dependancy

    Args:
        sky_file (:obj:`str`):
            The filename (NO PATH) of the sky file to use.

    Returns:
        :class:`~pypeit.onespec.OneSpec`: Sky spectrum
    """
    path = dataPaths.sky_spec.get_file_path(sky_file)
    return onespec.OneSpec.from_xspec_file(path)
