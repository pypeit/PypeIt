"""
Module to run tests on arqa
"""
import pathlib

from pypeit import log
from pypeit import qa

def test_set_qa_filename_bare_name_only():
    # set_qa_filename returns only a bare file name -- no directory
    # components at all (callers join it with outputPaths.qa_pngs
    # themselves). Every method still supported after the dead-branch
    # removal (§4.1) is checked here.
    methods = ['slit_trace_qa', 'slit_profile_qa', 'arc_fit_qa', 'arc_fwhm_qa',
               'arc_fit2d_global_qa', 'arc_fit2d_orders_qa', 'arc_tilts_spec_qa',
               'arc_tilts_spat_qa', 'arc_tilts_2d_qa', 'pca_arctilt',
               'plot_orderfits_Blaze', 'obj_trace_qa', 'obj_profile_qa',
               'spat_flexure_qa_corr', 'spec_flexure_qa_corr', 'spec_flexure_qa_sky',
               'spatillum_finecorr', 'detector_structure']
    for method in methods:
        name = qa.set_qa_filename('test', method, det='DET01', slit=1,
                                  prefix='pre', mode='global')
        assert len(pathlib.Path(name).parts) == 1, \
            f'{method} returned more than a bare file name: {name}'
        assert 'QA' not in name and 'PNGs' not in name, \
            f'{method} embedded a directory name in the file name: {name}'


def test_set_qa_filename_dead_methods_removed():
    # These two methods have no live call site anywhere in the codebase
    # (confirmed by a subsequent audit -- three of the five originally
    # flagged as dead in §4.1 turned out to be reachable indirectly, via
    # `html_mf_pngs`'s `fname=` dictionary values, and were restored).
    # They should raise, not silently succeed.
    for method in ['plot_orderfits_Arc', 'pca_plot']:
        try:
            qa.set_qa_filename('test', method, slit=1, prefix='pre')
        except IOError:
            pass
        else:
            raise AssertionError(f'{method} should no longer be a valid QA method')


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
