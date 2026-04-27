
try:
    from specutils import Spectrum, SpectrumList
except ModuleNotFoundError:
    Spectrum = None
    SpectrumList = None
else:
    from pypeit.specutils import pypeit_loaders


