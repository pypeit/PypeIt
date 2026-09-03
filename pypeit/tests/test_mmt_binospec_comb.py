"""
Unit tests for automatic combination-group assignment on MMT/Binospec IFU
data.

The Binospec IFU has no nod partner (sky comes from dedicated sky fibers via
a joint fit), so :meth:`get_comb_group` never sets ``bkg_id``; it only groups
science/standard frames that share a pointing onto a common ``comb_id`` so
same-pointing frames coadd while distinct dither positions each reduce to
their own spec1d/spec2d.
"""
import numpy as np
from astropy.table import Table

from pypeit.spectrographs.util import load_spectrograph


def _comb_table(rows, setups=None):
    """
    Build a minimal metadata table shaped like the one
    :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups` hands to
    ``get_comb_group``: every science/standard frame already carries a unique
    ``comb_id`` (1..k in table order) and ``bkg_id`` is -1 everywhere.  The
    ``dithoff``/``dithpos`` columns are seeded with their card defaults, as
    ``PypeItMetaData._build`` does from the ``init_meta`` definitions.

    Args:
        rows (list): sequence of ``(frametype, ra_deg, dec_deg, mjd)`` tuples.
        setups (list, optional): per-row setup label; defaults to all ``'A'``.
    """
    n = len(rows)
    ftype = [r[0] for r in rows]
    comb = np.full(n, -1, dtype=int)
    sci = [i for i, f in enumerate(ftype) if 'science' in f or 'standard' in f]
    for k, i in enumerate(sci):
        comb[i] = k + 1
    t = Table()
    t['filename'] = [f'f{i}.fits' for i in range(n)]
    t['frametype'] = ftype
    t['setup'] = setups if setups is not None else ['A'] * n
    t['ra'] = [r[1] for r in rows]
    t['dec'] = [r[2] for r in rows]
    t['mjd'] = [r[3] for r in rows]
    t['comb_id'] = comb
    t['bkg_id'] = np.full(n, -1, dtype=int)
    t['dithoff'] = np.zeros(n, dtype=float)
    t['dithpos'] = np.full(n, 'None')
    return t


def _offset(ra0, dec0, east_as, north_as):
    """Offset a pointing by (east, north) arcsec, returning (ra_deg, dec_deg)."""
    ra = ra0 + (east_as / 3600.0) / np.cos(np.radians(dec0))
    dec = dec0 + north_as / 3600.0
    return ra, dec


# JADES field pointing used throughout (12:36:52.6 +62:07:56.6).
RA0, DEC0 = 189.21916666666664, 62.13238888888889


def test_same_pointing_shares_comb_id():
    # Eight science frames at one identical pointing must coadd: a single
    # shared comb_id, no background pairing, ~zero dither offset.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('science', RA0, DEC0, 60819.14 + 0.01 * i) for i in range(8)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert np.all(comb == comb[0])
    assert comb[0] > 0
    assert np.all(np.asarray(t['bkg_id']) == -1)
    assert np.allclose(np.asarray(t['dithoff'], dtype=float), 0.0, atol=1e-6)
    assert len(set(str(x) for x in t['dithpos'])) == 1


def test_two_position_dither_splits_comb_id():
    # A two-position dither (0.5" apart, well above the 0.1" tolerance) taken
    # ABAB: frames at the same position share a comb_id, the two positions get
    # different comb_ids, and nothing is background-paired.
    spec = load_spectrograph('mmt_binospec_ifu')
    posA = (RA0, DEC0)
    posB = _offset(RA0, DEC0, 0.5, 0.0)
    seq = [posA, posB, posA, posB]
    rows = [('science', p[0], p[1], 60819.14 + 0.01 * i)
            for i, p in enumerate(seq)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == comb[2]            # both A frames
    assert comb[1] == comb[3]            # both B frames
    assert comb[0] != comb[1]            # A vs B distinct
    assert len(set(comb)) == 2
    assert np.all(np.asarray(t['bkg_id']) == -1)
    assert len(set(str(x) for x in t['dithpos'])) == 2


def test_subtolerance_jitter_stays_one_group():
    # Pointing wander below the tolerance (0.02") must not split the group.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = []
    for i in range(4):
        ra, dec = _offset(RA0, DEC0, 0.02 * (i % 2), 0.0)
        rows.append(('science', ra, dec, 60819.14 + 0.01 * i))
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert np.all(comb == comb[0])


def test_large_offset_gets_own_comb_id():
    # A large sky-acquisition offset (60") lands in its own cluster: the six
    # on-object frames share one comb_id, the offset frame gets its own. No
    # special sky handling -- it is just another cluster.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('science', RA0, DEC0, 60819.14 + 0.01 * i) for i in range(6)]
    sky = _offset(RA0, DEC0, 60.0, 0.0)
    rows.append(('science', sky[0], sky[1], 60819.14 + 0.07))
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert np.all(comb[:6] == comb[0])
    assert comb[6] != comb[0]
    assert len(set(comb)) == 2


def test_two_setups_not_merged():
    # Identical pointing but two different configurations must never share a
    # comb_id -- combination cannot cross a setup.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('science', RA0, DEC0, 60819.14 + 0.01 * i) for i in range(4)]
    setups = ['A', 'A', 'B', 'B']
    t = spec.get_comb_group(_comb_table(rows, setups=setups))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == comb[1]
    assert comb[2] == comb[3]
    assert comb[0] != comb[2]


def test_standard_at_own_pointing_separate():
    # A standard star observed at a different pointing than the science field
    # gets its own comb_id; the science frames still coadd together.
    spec = load_spectrograph('mmt_binospec_ifu')
    std = _offset(RA0, DEC0, 200.0, 120.0)
    rows = [('science', RA0, DEC0, 60819.14),
            ('science', RA0, DEC0, 60819.15),
            ('science', RA0, DEC0, 60819.16),
            ('standard', std[0], std[1], 60819.20)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == comb[1] == comb[2]
    assert comb[3] != comb[0]
    assert np.all(np.asarray(t['bkg_id']) == -1)


def test_calibration_frames_untouched():
    # Non science/standard frames keep comb_id/bkg_id == -1.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('arc,tilt', RA0, DEC0, 60819.10),
            ('pixelflat,illumflat,trace', RA0, DEC0, 60819.11),
            ('science', RA0, DEC0, 60819.14),
            ('science', RA0, DEC0, 60819.15)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == -1 and comb[1] == -1
    assert comb[2] == comb[3] and comb[2] > 0
