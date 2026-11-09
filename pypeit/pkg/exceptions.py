"""
Provides pypeit specific exceptions.
"""

__all__ = [
    'PypeItError',
    'PypeItWavelengthCalibrationError',
    'PypeItBitMaskError',
    'PypeItDataModelError',
    'PypeItPathError'
]

class PypeItError(Exception):
    pass

class PypeItWavelengthCalibrationError(PypeItError):
    """A wavelength fit failed because the data cannot support a solution."""

class PypeItBitMaskError(PypeItError):
    pass

class PypeItDataModelError(PypeItError):
    pass

class PypeItPathError(PypeItError):
    pass
