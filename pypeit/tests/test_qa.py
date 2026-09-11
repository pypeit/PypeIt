"""
Module to run tests on arqa
"""
import matplotlib
matplotlib.use('Agg', force=True)
from matplotlib import pyplot as plt

import numpy as np
from PIL import Image
import pytest

from pypeit import log
from pypeit import qa


def test_save_figure_serial(tmp_path):
    qa.init_qa_pool(1)
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    out = tmp_path / 'serial.png'
    qa.save_figure(fig, out, dpi=50)
    assert out.exists(), 'serial save_figure must write in-line'
    assert len(plt.get_fignums()) == 0, 'figure must be closed'


def test_save_figure_threaded_matches_serial(tmp_path):
    """Deferred writes must be pixel-identical to in-line writes."""
    rng = np.random.default_rng(2718)          # fixed seed: deterministic
    data = rng.normal(size=(64, 64))

    def _draw():
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.imshow(data)
        ax.set_title('QA')
        return fig

    qa.init_qa_pool(1)
    ref = tmp_path / 'ref.png'
    qa.save_figure(_draw(), ref, dpi=80)
    qa.flush_qa()

    qa.init_qa_pool(4)
    outs = []
    for i in range(8):
        o = tmp_path / f'par{i}.png'
        qa.save_figure(_draw(), o, dpi=80)
        outs.append(o)
    qa.flush_qa()
    qa.init_qa_pool(1)

    # Compare decoded pixel arrays, not raw bytes: matplotlib embeds a
    # 'Software' PNG chunk carrying its version string.
    ref_px = np.asarray(Image.open(ref))
    for o in outs:
        assert np.array_equal(np.asarray(Image.open(o)), ref_px), \
            'threaded QA write differs from the serial write'


def test_save_figure_pending_isolated_from_pyplot(tmp_path):
    """
    A figure queued by save_figure must be deregistered from pyplot so that
    later pyplot-state plotting (e.g. a bare plt.plot) cannot draw into it.
    """
    def _draw():
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.plot([0, 1], [1, 0])
        return fig

    qa.init_qa_pool(1)
    ref = tmp_path / 'ref.png'
    qa.save_figure(_draw(), ref, dpi=80)

    qa.init_qa_pool(2)
    out = tmp_path / 'pending.png'
    qa.save_figure(_draw(), out, dpi=80)
    # This must land on a fresh figure, not the still-pending one above
    plt.plot([0, 1], [5, 5], color='red', linewidth=10)
    plt.close('all')
    qa.flush_qa()
    qa.init_qa_pool(1)

    assert np.array_equal(np.asarray(Image.open(out)), np.asarray(Image.open(ref))), \
        'pyplot-state plotting contaminated a queued QA figure'


def test_flush_qa_reraises(tmp_path):
    """A failed background write must surface on the main thread."""
    qa.init_qa_pool(2)
    fig, _ = plt.subplots()
    qa.save_figure(fig, tmp_path / 'nonexistent_dir' / 'x.png', dpi=50)
    with pytest.raises(Exception):
        qa.flush_qa()
    qa.init_qa_pool(1)


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
