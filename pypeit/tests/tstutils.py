# Odds and ends in support of tests
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import pytest

from pypit import arcimage
from pypit import arpixels
from pypit import traceslits
from pypit import wavecalib
from pypit import wavetilts
from pypit import processimages
from pypit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def load_kast_blue_masters(get_spectrograph=False, aimg=False, tslits=False, tilts=False,
                           datasec=False, wvcalib=False):

    spectrograph = load_spectrograph(spectrograph='shane_kast_blue')
    spectrograph.naxis = (2112,350)     # Image shape with overscan

    root_path = data_path('MF') if os.getenv('PYPIT_DEV') is None \
                    else os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'MF')
    directory_path = root_path+'_'+spectrograph.spectrograph
                    
    mode = 'reuse'

    # Load up the Masters
    ret = []

    if get_spectrograph:
        ret.append(spectrograph)

    setup = 'A_01_aa'
    if aimg:
        AImg = arcimage.ArcImage(spectrograph, setup=setup, root_path=root_path, mode=mode)
        msarc, header, _ = AImg.load_master_frame()
        ret.append(msarc)

    if tslits:
        TSlits = traceslits.TraceSlits.from_master_files(os.path.join(directory_path,
                                                                      'MasterTrace_A_01_aa'))
        TSlits._make_pixel_arrays()
        _ = TSlits._fill_tslits_dict()
        ret.append(TSlits)

    if tilts:
        wvTilts = wavetilts.WaveTilts(None, spectrograph=spectrograph, setup=setup,
                                      root_path=root_path, mode=mode)
        tilts = wvTilts.master()
        ret.append(tilts)

    if datasec:
        datasec_img = spectrograph.get_datasec_img(data_path('b1.fits.gz'), 1)
        ret.append(datasec_img)

    if wvcalib:
        Wavecalib = wavecalib.WaveCalib(None, spectrograph=spectrograph, setup=setup,
                                        root_path=root_path, mode=mode)
        wv_calib = Wavecalib.master()
        ret.append(wv_calib)

    # Return
    return ret
