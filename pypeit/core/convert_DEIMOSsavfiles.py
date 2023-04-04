from scipy.io import readsav
import numpy as np
import pathlib

from pypeit import io
from pypeit import data


def sav_to_fits(savfile):
    """
    Simple method to convert the DEIMOS .sav files, which contain the
    optical model maps, to .fits files.  The new file has the
    same name with a .fits extension.

    ToDO: This is specific for keck_deimos `static_calib` data, since the path
    is explicitly mentioned. If needed, this method could be generalized.

    TODO: move it to `pypeit.io`.

    Args:
        savfile (str): full path to the .sav file


    """
    savfile_name = pathlib.Path(savfile).stem
    sav = readsav(savfile, python_dict=True)

    list_keys = list(sav.keys())
    for k in list_keys:
        if type(sav[k]) is not np.ndarray:
            sav[k] = np.asarray([sav[k]])

    io.write_to_fits(
        sav,
        data.Paths.static_calibs / 'keck_deimos', f'{savfile_name}.fits',
        overwrite=True
    )

    return
