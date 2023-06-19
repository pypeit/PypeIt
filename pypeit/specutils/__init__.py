
try:
    from specutils import Spectrum1D, SpectrumList
except ModuleNotFoundError:
    Spectrum1D = None
    SpectrumList = None
else:
    from pypeit.specutils import pypeit_loaders


