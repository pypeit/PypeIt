"""
Module to run tests on Alignment frames
"""
from pathlib import Path

from IPython import embed

import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator

from astropy.io import fits

from pypeit import alignframe
from pypeit.tests.tstutils import data_output_path


def test_alignments():
    nspat, nspec, nslits, nalign = 100, 1000, 10, 5
    tmp = np.ones((nspec, nspat)) * 10.
    traces = np.zeros((nspec, nalign, nslits))
    instant_dict = dict(alignframe=tmp,
                        nspec=nspec,
                        nalign=nalign,
                        nslits=nslits,
                        traces=traces,
                        PYP_SPEC='keck_kcwi',
                        spat_id=np.arange(nslits))

    alignments = alignframe.Alignments(**instant_dict)
    alignments.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(alignments.get_path()).absolute()

    # I/O
    alignments.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'
    with fits.open(ofile) as hdu:
        assert hdu[1].name == 'ALIGN', 'Extension name changed'
    _alignments = alignframe.Alignments.from_file(str(ofile))

    # Test
    for key in instant_dict.keys():
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(alignments[key], _alignments[key])
        else:
            assert alignments[key] == _alignments[key]

    # Clean-up
    ofile.unlink()



def test_interp1d_rows():
    rng = np.random.default_rng(1234)
    nrow, n, m = 50, 5, 30
    x = rng.uniform(0., 10., size=(nrow, n))
    y = rng.normal(size=(nrow, n))
    # Includes values outside the range of x, to test the extrapolation
    x_new = rng.uniform(-3., 13., size=(nrow, m))
    x_new[:, 0] = x[:, 2]   # A value equal to one of the abscissae
    result = alignframe.interp1d_rows(x, y, x_new)
    for i in range(nrow):
        expected = interp1d(x[i], y[i], kind='linear', bounds_error=False,
                            fill_value='extrapolate')(x_new[i])
        assert np.array_equal(result[i], expected)
    # Integer ordinates are cast to float, as in interp1d
    yint = np.arange(nrow*n).reshape(nrow, n)
    result = alignframe.interp1d_rows(x, yint, x_new)
    assert np.array_equal(result[3], interp1d(x[3], yint[3], bounds_error=False,
                                              fill_value='extrapolate')(x_new[3]))


def reference_alignment_transforms(traces, locations, tilts):
    """
    The per-pixel loop that :func:`pypeit.alignframe.AlignmentSplines.build_splines`
    used before it was vectorised.  Returns the grid and values of the
    transform of each slit.
    """
    nspec, ntrace, nslit = traces.shape
    locations = np.asarray(locations)
    spl_fulltilts = RegularGridInterpolator((np.arange(tilts.shape[0]), np.arange(tilts.shape[1])),
                                            tilts * (nspec - 1), method='linear')
    spldict = dict(kind='linear', bounds_error=False, fill_value="extrapolate")
    ycoord = np.arange(nspec)
    out = []
    for sl in range(nslit):
        xlr, tlr = np.zeros((nspec, 2)), np.zeros((nspec, 2))
        spl_loc = nspec * [None]
        eval_trim = 2
        for sp in range(nspec):
            xlr[sp, :] = interp1d(locations, traces[sp,:,sl], **spldict)([0.0, 1.0])
            tmptilt = spl_fulltilts([[sp, xlr[sp,0] + eval_trim],
                                     [sp, xlr[sp,0] + eval_trim+1],
                                     [sp, xlr[sp,1] - eval_trim],
                                     [sp, xlr[sp,1] - eval_trim-1]])
            tlr[sp, :] = [tmptilt[0]-eval_trim*(tmptilt[1]-tmptilt[0]),
                          tmptilt[2]-eval_trim*(tmptilt[3]-tmptilt[2])]
            spl_loc[sp] = interp1d(traces[sp,:,sl], locations, **spldict)
        tilt_ypos = interp1d(tlr[:,1], ycoord, **spldict)
        ypos_xpos = interp1d(ycoord, xlr[:,1], **spldict)
        ypos = tilt_ypos(tlr[:,0])
        xpos = ypos_xpos(ypos)
        slitlen = np.sqrt((ypos - np.arange(ypos.size)) ** 2 + (xpos - xlr[:, 0]) ** 2)
        spl_slen = interp1d(ycoord, slitlen, **spldict)
        xcoord = np.arange(np.floor(np.min(xlr)), np.ceil(np.max(xlr))+1, 1.0)
        out_transform = np.zeros((nspec, xcoord.size))
        for sp in range(nspec):
            out_transform[sp,:] = spl_loc[sp](xcoord) * spl_slen(sp)
        out.append((ycoord, xcoord, out_transform))
    return out


def test_alignment_splines():
    nspec, nspat = 300, 260
    spec = np.arange(nspec, dtype=float)
    spat = np.arange(nspat, dtype=float)
    # Curved, tilted arc lines
    tilts = (spec[:,None] + 1e-3*(spat[None,:]-130.)**2 + 0.05*(spat[None,:]-130.)) / (nspec-1)
    # Three curved slits
    cen = np.array([50., 120., 200.])[None,:] + 2e-4*(spec[:,None]-150.)**2 + 0.02*spec[:,None]
    width = np.array([30., 36., 25.]) + 0.01*spec[:,None]

    # Slit edges, as used by coadd3d when there are no alignment traces
    edges = np.stack([cen - width/2, cen + width/2], axis=1)
    # Alignment traces that do not reach the slit edges (so the edges are
    # extrapolated), given in a random order
    align_loc = np.array([0.5, 0.1, 0.9, 0.3, 0.7])
    align = cen[:,None,:] + (align_loc[None,:,None] - 0.5) * width[:,None,:] \
                + 0.3*np.sin(spec/40.)[:,None,None]

    pixspec, pixspat = np.meshgrid(np.linspace(0, nspec-1, 77), np.linspace(10, nspat-10, 91),
                                   indexing='ij')
    for traces, locations in [(edges, [0.0, 1.0]), (align, align_loc)]:
        splines = alignframe.AlignmentSplines(traces, locations, tilts)
        for sl, (ycoord, xcoord, values) in \
                enumerate(reference_alignment_transforms(traces, locations, tilts)):
            assert np.array_equal(splines.spl_transform[sl].grid[0], ycoord)
            assert np.array_equal(splines.spl_transform[sl].grid[1], xcoord)
            assert np.array_equal(splines.spl_transform[sl].values, values)
            expected = RegularGridInterpolator((ycoord, xcoord), values, method='linear',
                                               bounds_error=False,
                                               fill_value=None)((pixspec, pixspat))
            assert np.array_equal(splines.transform(sl, pixspat, pixspec), expected)
