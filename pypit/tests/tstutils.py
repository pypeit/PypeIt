# Odds and ends in support of tests
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

from pypit import arcimage
from pypit import arpixels
from pypit import traceslits
from pypit import wavecalib
from pypit import wavetilts
from pypit import processimages
from pypit.spectrographs import io

pypdev_path = os.getenv('PYPIT_DEV')


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def load_kast_blue_masters(get_settings=False, aimg=False, tslits=False, tilts=False,
                           datasec=False, wvcalib=False):
    spectrograph = 'shane_kast_blue'
    settings = processimages.default_settings()
    settings['masters'] = {}
    if pypdev_path is not None:
        settings['masters']['directory'] = pypdev_path + '/Cooked/MF_shane_kast_blue'
    else:
        settings['masters']['directory'] = data_path('MF')+'_shane_kast_blue'
    settings['masters']['reuse'] = True
    settings['masters']['loaded'] = []

    settings['detector']['dataext'] = 0
    settings['detector']['datasec01'] = [[0, 1024], [0, 0]]
    settings['detector']['datasec02'] = [[1024, 2048], [0, 0]]
    settings['detector']['oscansec01'] = [[2049, 2080], [0, 0]]
    settings['detector']['oscansec02'] = [[2080, 2111], [0, 0]]
    settings['detector']['naxis0'] = 2112  # Raw frame, with overscan
    settings['detector']['naxis1'] = 350
    settings['detector']['numamplifiers'] = 2
    settings['detector']['gain'] = [1.2, 1.2]

    # Load up the Masters
    ret = []

    if get_settings:
        ret.append(settings)
    #
    setup = 'A_01_aa'
    if aimg:
        AImg = arcimage.ArcImage(setup=setup, settings=settings)
        msarc, header, _ = AImg.load_master_frame()
        ret.append(msarc)
    #
    if tslits:
        TSlits = traceslits.TraceSlits.from_master_files(settings['masters']['directory'] + '/MasterTrace_A_01_aa')
        TSlits._make_pixel_arrays()
        _ = TSlits._fill_tslits_dict()
        ret.append(TSlits)
    if tilts:
        wvTilts = wavetilts.WaveTilts(None, settings=settings, setup=setup)
        tilts = wvTilts.master()
        ret.append(tilts)
    if datasec:
        datasec, _ = io.get_datasec(spectrograph, filename=None, det_settings=settings['detector'],
                                    numamplifiers=settings['detector']['numamplifiers'], det=1)
        datasec_img = arpixels.pix_to_amp(settings['detector']['naxis0'],
                                          settings['detector']['naxis1'],
                                          datasec, settings['detector']['numamplifiers'])
        ret.append(datasec_img)
    if wvcalib:
        Wavecalib = wavecalib.WaveCalib(None, settings=settings, setup=setup)
        wv_calib = Wavecalib.master()
        ret.append(wv_calib)

    # Return
    return ret
