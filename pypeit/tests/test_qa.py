"""
Module to run tests on arqa
"""
from pypeit import msgs
from pypeit.core import qa

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
