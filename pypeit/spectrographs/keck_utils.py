"""
Utilities shared by the Keck spectrograph classes.
"""


def koa_qlview_header_fields(hdr):
    """
    Read the header fields common to Keck/KOA raw files for display in the
    quicklook viewer.

    Used by the ``qlview_raw_info`` methods of the Keck spectrographs, which
    override individual entries where an instrument uses different keywords.

    Args:
        hdr (`astropy.io.fits.Header`_):
            Primary header of a raw file.

    Returns:
        :obj:`dict`: Values for the keys ``OBJECT``, ``FRAMENO``, ``IMTYPE``
        (from ``KOAIMTYP``), ``MASKNAME``, ``OBSMODE``, and ``EXPTIME``.
        ``EXPTIME`` falls back through ``TTIME``, ``ITIME``, ``ETIME``, and
        ``ELAPTIME``, in that order.  Missing keywords are ``'N/A'``.
    """
    exptime = None
    for key in ['EXPTIME', 'TTIME', 'ITIME', 'ETIME', 'ELAPTIME']:
        exptime = hdr.get(key, None)
        if exptime is not None:
            break
    return {
        'OBJECT': hdr.get('OBJECT', 'N/A'),
        'FRAMENO': hdr.get('FRAMENO', 'N/A'),
        'IMTYPE': hdr.get('KOAIMTYP', 'N/A'),
        'MASKNAME': hdr.get('MASKNAME', 'N/A'),
        'OBSMODE': hdr.get('OBSMODE', 'N/A'),
        'EXPTIME': 'N/A' if exptime is None else exptime,
    }
