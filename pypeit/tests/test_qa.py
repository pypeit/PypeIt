"""
Module to run tests on arqa
"""
import matplotlib
matplotlib.use('Agg', force=True)
from matplotlib import pyplot as plt

import numpy as np
from PIL import Image
import pytest

import pypeit
from pypeit import log
from pypeit import qa
from pypeit.pkg.qawriter import QAWriter


def _draw(seed=2718):
    """Build a deterministic test figure."""
    rng = np.random.default_rng(seed)
    fig, ax = plt.subplots(figsize=(3, 3))
    ax.imshow(rng.normal(size=(64, 64)))
    ax.set_title('QA')
    return fig


def test_package_writer_is_serial():
    """The package-level writer must default to the serial path."""
    assert isinstance(pypeit.qaWriter, QAWriter)
    assert not pypeit.qaWriter.parallel, \
        'the default QA writer must write figures serially'


def test_save_figure_serial(tmp_path):
    writer = QAWriter()
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    out = tmp_path / 'serial.png'
    writer.save_figure(fig, out, dpi=50)
    assert out.exists(), 'serial save_figure must write in-line'
    assert len(plt.get_fignums()) == 0, 'figure must be closed'


def test_save_figure_threaded_matches_serial(tmp_path):
    """Deferred writes must be pixel-identical to in-line writes."""
    serial = QAWriter()
    ref = tmp_path / 'ref.png'
    serial.save_figure(_draw(), ref, dpi=80)
    serial.flush()

    threaded = QAWriter(ncpu=4)
    assert threaded.parallel
    outs = []
    for i in range(8):
        o = tmp_path / f'par{i}.png'
        threaded.save_figure(_draw(), o, dpi=80)
        outs.append(o)
    threaded.flush()

    # Compare decoded pixel arrays, not raw bytes: matplotlib embeds a
    # 'Software' PNG chunk carrying its version string.
    ref_px = np.asarray(Image.open(ref))
    for o in outs:
        assert np.array_equal(np.asarray(Image.open(o)), ref_px), \
            'threaded QA write differs from the serial write'


def test_save_figure_pending_isolated_from_pyplot(tmp_path):
    """
    A queued figure must be deregistered from pyplot so that later
    pyplot-state plotting (e.g. a bare plt.plot) cannot draw into it.
    """
    ref = tmp_path / 'ref.png'
    QAWriter().save_figure(_draw(), ref, dpi=80)

    writer = QAWriter(ncpu=2)
    out = tmp_path / 'pending.png'
    writer.save_figure(_draw(), out, dpi=80)
    # This must land on a fresh figure, not the still-pending one above
    plt.plot([0, 1], [5, 5], color='red', linewidth=10)
    plt.close('all')
    writer.flush()

    assert np.array_equal(np.asarray(Image.open(out)), np.asarray(Image.open(ref))), \
        'pyplot-state plotting contaminated a queued QA figure'


def test_flush_reraises(tmp_path):
    """A failed background write must surface on the calling thread."""
    writer = QAWriter(ncpu=2)
    fig, _ = plt.subplots()
    writer.save_figure(fig, tmp_path / 'nonexistent_dir' / 'x.png', dpi=50)
    with pytest.raises(Exception):
        writer.flush()


def test_init_is_repeatable(tmp_path):
    """
    init() must be safe to call repeatedly; a forked worker uses it to drop the
    pool inherited from its parent.
    """
    writer = QAWriter(ncpu=4)
    assert writer.parallel
    writer.save_figure(_draw(), tmp_path / 'queued.png', dpi=80)
    # Re-initializing drops the inherited pool and any queued encodes
    writer.init(ncpu=1)
    assert not writer.parallel
    assert len(writer.pending) == 0
    writer.init(ncpu=2)
    assert writer.parallel
    out = tmp_path / 'after.png'
    writer.save_figure(_draw(), out, dpi=80)
    writer.flush()
    assert out.exists()


def test_get_dimen():
    """ Get the plotting dimensions
    Returns
    -------

    """
    npanels, maxp = 1, 25
    pages, npp = qa.get_dimen(npanels, maxp=maxp)
    assert len(pages) == 1 and pages[0][0]*pages[0][1] == 1 and len(npp) == 1 and npp[0] == npanels
    npanels, maxp = 5, 5
    pages, npp = qa.get_dimen(npanels, maxp=maxp)
    assert len(pages) == 1 and pages[0][0] * pages[0][1] == 6 and len(npp) == 1 and npp[0] == npanels
    npanels, maxp = 22, 8
    pages, npp = qa.get_dimen(npanels, maxp=maxp)
    assert (len(pages) == 3) and (pages[0][0] * pages[0][1] == maxp) and (pages[1][0] * pages[1][1] == maxp)
    assert (len(npp) == 3) and (npp[0] == maxp) and (npp[1] == maxp) and (npp[2] == 6)
    npanels, maxp = 22, 7
    pages, npp = qa.get_dimen(npanels, maxp=maxp)
    assert (len(pages) == 4) and (pages[0][0] * pages[0][1] == maxp+1) and (pages[1][0] * pages[1][1] == maxp+1) \
        and (pages[2][0] * pages[2][1] == maxp + 1) and (pages[3][0] * pages[3][1] == 1)
    assert (len(npp) == 4) and (npp[0] == maxp) and (npp[1] == maxp) and (npp[2] == maxp) and (npp[3] == 1)
