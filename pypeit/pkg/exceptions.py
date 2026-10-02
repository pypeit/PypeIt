"""
Provides pypeit specific exceptions.
"""

__all__ = [
    'PypeItError',
    'PypeItBitMaskError',
    'PypeItDataModelError',
    'PypeItPathError'
]

class PypeItError(Exception):
    pass

class PypeItBitMaskError(PypeItError):
    pass

class PypeItDataModelError(PypeItError):
    pass

class PypeItPathError(PypeItError):
    pass
