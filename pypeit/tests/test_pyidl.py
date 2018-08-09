import numpy as np
from pypeit.core.pydl import bspline

def test_bsplinetodict():
    """ Test for writing a bspline onto a dict
    (and also reading it out).
    """

    x = np.random.rand(500)
    # Create bspline
    
    init_bspline = bspline(x, bkspace=0.01*(np.max(x)-np.min(x)))
    bspline_todict = init_bspline.to_dict()
    
    bspline_fromdict = bspline.from_dict(bspline_todict)
    
    print(' ')
    print('Returns 0. if the breakpoints are the same before and after')
    print('reading from the dictionary:')
    print(np.max(np.array(bspline_todict['breakpoints'])-bspline_fromdict.breakpoints))

