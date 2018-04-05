#!/usr/bin/env python3

import time
from pypit.par import ParSet

#-----------------------------------------------------------------------------

class PypitLoadPar(ParSet):
    def __init__(self, settings=None, spect=None):
        par_types = [ParSet, dict]

        pars = [ 'settings', 'spect' ]

        values = [ settings, spect ]

        # No defaults

        # No fixed options

        dtypes =   [ par_types, par_types ]

        # Can-call is all False

        super(PypitLoadPar, self).__init__(pars, values=values, dtypes=dtypes)


class PypitDirectoryPar(ParSet):
    def __init__(self, master=None, science=None, qa=None):
        pars = [ 'master', 'science', 'qa' ]

        values = [ master, science, qa ]

        defaults = [ './master', './science', 'qa' ]

        # No fixed options

        dtypes =   [ str, str, str ]

        # Can-call is all False

        super(PypitDirectoryPar, self).__init__(pars, values=values, defaults=defaults,
                                                dtypes=dtypes)


class PypitRunPar(ParSet):
    def __init__(self, spectrograph=None, pypitdir=None, setup=None, ncpus=None, calcheck=None,
                 preponly=False, redname=None, useIDname=None, load=None, directory=None):
       
        par_types = [ParSet, dict]
        spectrograph_options = PypitRunPar._available_spectrographs()

        pars = [ 'spectrograph', 'pypitdir', 'setup', 'ncpus', 'calcheck', 'preponly', 'redname',
                 'useIDname', 'load', 'directory' ]

        values = [ spectrograph, pypitdir, setup, ncpus, calcheck, preponly, redname, useIDname,
                   load, directory ]

        defaults = [ 'shane_kast_blue', None, False, 1, False, False, None, False,
                     PypitLoadPar(), PypitDirectoryPar() ]

        options =  [ spectrograph_options, None, None, None, None, None, None, None, None, None ]

        dtypes =   [ str, str, bool, int, bool, bool, str, bool, par_types, par_types ]

        # Can-call is all False

        super(PypitRunPar, self).__init__(pars, values=values, defaults=defaults, options=options,
                                          dtypes=dtypes)


    @staticmethod
    def _available_spectrographs():
        """
        Return the list of allowed spectrographs for pypit reductions.

        .. todo::
            - Listed explicitly for now.  We can do something like what
              is currently done in the code by trolling through a file
              directory.
        """
        return [ 'apf_levy', 'keck_deimos', 'keck_hires', 'keck_lris_blue', 'keck_lris_red',
                 'shane_kast_blue', 'shane_kast_red', 'tng_dolores', 'wht_isis_blue' ]

    

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    test = PypitRunPar()

    print(test)

    try:
        test['spectrograph'] = 'new_spectrograph'
    except ValueError as e:
        print(e)

    try:
        test['calcheck'] = 5
    except TypeError as e:
        print(e)

    try:
        test['ncpus'] = 3.2
    except TypeError as e:
        print(e)
    
    try:
        test['load']['settings'] = 'test'
    except TypeError as e:
        print(e)
    
    try:
        test['directory']['master'] = 5.3
    except TypeError as e:
        print(e)
    

    print('Elapsed time: {0} seconds'.format(time.clock() - t))



