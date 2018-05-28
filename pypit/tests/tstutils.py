# Odds and ends in support of tests
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

from pypit import traceslits
from pypit import arcimage


def load_kast_blue_masters(aimg=False, tslits=False):
    settings = dict(masters={})
    settings['masters']['directory'] = os.getenv('PYPIT_DEV') + '/Cooked/MF_shane_kast_blue'
    settings['masters']['reuse'] = True
    setup = 'A_01_aa'
    # Load up the Masters
    ret = []
    if aimg:
        AImg = arcimage.ArcImage(setup=setup, settings=settings)
        msarc, header, _ = AImg.load_master_frame()
        ret.append(msarc)
    #
    if tslits:
        TSlits = traceslits.TraceSlits.from_master_files(settings['masters']['directory'] + '/MasterTrace_A_01_aa')
        TSlits._make_pixel_arrays()
        ret.append(TSlits)
    # Return
    return ret
