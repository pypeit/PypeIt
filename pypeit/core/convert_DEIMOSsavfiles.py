from scipy.io import readsav
import numpy as np
import os.path

from pypeit import io
from pkg_resources import resource_filename


def sav_to_fits(savfile):
    """
    Simple method to convert the DEIMOS .sav files, which contain the
    optical model maps, to .fits files.
    ToDO: This is specific for keck_deimos `static_calib` data, since the path
    is explicitly mentioned. If needed, this method could be generalized.
    TODO: move it to `pypeit.io`.
    Args:
        savfile: path to the .sav file

    Returns:
        Save the content of the .sav file into a .fits file in `data/static_calibs/keck_deimos/`

    """
    savfile_name = os.path.splitext(os.path.basename(savfile))[0]
    sav = readsav(savfile, python_dict=True)

    list_keys = list(sav.keys())
    for k in list_keys:
        if type(sav[k]) is not np.ndarray:
            sav[k] = np.asarray([sav[k]])

    to_path = resource_filename('pypeit', 'data/static_calibs/keck_deimos/')
    io.write_to_fits(sav, to_path + savfile_name + '.fits', overwrite=True)

    return
