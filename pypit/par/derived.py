
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import inspect
from .parset import ParSet
from configobj import ConfigObj

#-----------------------------------------------------------------------------
# Helper functions

def _recursive_dict_evaluate(d):
    """
    will not work for lists of dicts...
    """
    for k in d.keys():
        if isinstance(d[k], dict):
           d[k] = _recursive_dict_evaluate(d[k])
        elif isinstance(d[k], list):
            try:
                d[k] = [ eval(e) for e in d[k] ]
            except:
                pass
        else:
            try:
                d[k] = eval(d[k])
            except:
                pass

    return d


#-----------------------------------------------------------------------------
# Lower-level ParSets

class LoadPar(ParSet):
    def __init__(self, settings=None, spect=None):
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        dtypes['settings'] = [ParSet, dict]
        descr['settings'] = 'Load a reduction settings file (Note: this command overwrites ' \
                            'all default settings'

        dtypes['spect'] = [ParSet, dict]
        descr['spect'] = 'Load a spectrograph settings file (Note: this command overwrites ' \
                         'all default settings'

        # Instantiate the parameter set
        super(LoadPar, self).__init__(list(pars.keys()),
                                      values=list(pars.values()),
                                      dtypes=list(dtypes.values()),
                                      descr=list(descr.values()))

    @classmethod
    def from_dict(cls, cfg):
        return cls(settings=cfg['settings'], spect=cfg['spect'])


class DirectoryPar(ParSet):
    def __init__(self, master=None, science=None, qa=None):
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['master'] = 'master'
        dtypes['master'] = str
        descr['master'] = 'Directory relative to calling directory to write master files.'
        
        defaults['science'] = 'science'
        dtypes['science'] = str
        descr['science'] = 'Directory relative to calling directory to write science files.'
        
        defaults['qa'] = 'qa'
        dtypes['qa'] = str
        descr['qa'] = 'Directory relative to calling directory to write qa files.'
        
        # Instantiate the parameter set
        super(DirectoryPar, self).__init__(list(pars.keys()),
                                           values=list(pars.values()),
                                           defaults=list(defaults.values()),
                                           dtypes=list(dtypes.values()),
                                           descr=list(descr.values()))

    @classmethod
    def from_dict(cls, cfg):
        return cls(master=cfg['master'], science=cfg['science'], qa=cfg['qa'])


#-----------------------------------------------------------------------------
# Top-level ParSets

class RunPar(ParSet):
    def __init__(self, spectrograph=None, pypitdir=None, setup=None, ncpus=None, calcheck=None,
                 preponly=False, redname=None, useIDname=None, load=None, directory=None):
        # Grab the parameter names and values from the function
        # arguments
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        pars = dict([(k,values[k]) for k in args[1:]])

        # Initialize the other used specifications for this parameter
        # set
        defaults = dict.fromkeys(pars.keys())
        options = dict.fromkeys(pars.keys())
        dtypes = dict.fromkeys(pars.keys())
        descr = dict.fromkeys(pars.keys())

        # Fill out parameter specifications.  Only the values that are
        # *not* None (i.e., the ones that are defined) need to be set
        defaults['spectrograph'] = 'shane_kast_blue'
        options['spectrograph']  = RunPar._available_spectrographs()
        dtypes['spectrograph']   = str
        descr['spectrograph']    = 'Spectrograph used to obtain the data to be reduced.  ' \
                                   'Options are: {0}'.format(
                                       ', '.join(RunPar._available_spectrographs()))

        dtypes['pypitdir']   = str
        descr['pypitdir']    = 'Directory with pypit python code.  Should not be set by the user.'
        
        defaults['setup'] = False
        dtypes['setup']   = bool
        descr['setup']    = 'If True, run in setup mode.  Useful to parse files when starting ' \
                            'reduction on a large set of data'

        defaults['ncpus'] = 1
        dtypes['ncpus']   = int
        descr['ncpus']    = 'Number of CPUs to use (-1 means all bar one CPU, -2 means all bar ' \
                            'two CPUs)'

        defaults['calcheck'] = False
        dtypes['calcheck']   = bool
        descr['calcheck']    = 'Flag to skip the data reduction and just checks to make sure all ' \
                               'calibration data are present'

        defaults['preponly'] = False
        dtypes['preponly']   = bool
        descr['preponly']    = 'If True, pypit will prepare the calibration frames; if false, ' \
                               'will only reduce the science frames'

        dtypes['redname']   = str
        descr['redname']    = 'Name of the reduction configuration file provided by the user.  ' \
                              'Set by the code internally, not by the user.'

        defaults['useIDname'] = False
        dtypes['useIDname']   = bool
        descr['useIDname']    = 'If True, file sorting will ensure that the idname is made'

        defaults['load'] = LoadPar()
        dtypes['load']   = [ParSet, dict]
        descr['load']    = 'Additional settings files to load'

        defaults['directory'] = DirectoryPar()
        dtypes['directory']   = [ParSet, dict]
        descr['directory']    = 'Directories for output'

        # Instantiate the parameter set
        super(RunPar, self).__init__(list(pars.keys()),
                                     values=list(pars.values()),
                                     defaults=list(defaults.values()),
                                     options=list(options.values()),
                                     dtypes=list(dtypes.values()),
                                     descr=list(descr.values()),
                                     cfg_section='run', cfg_comment='Execution options')


    @classmethod
    def from_cfg_file(cls, cfg_file, merge_with=None, evaluate=True):
        # Read the configuration file
        cfg = ConfigObj(cfg_file)
        # Get the list of other files to merge it with
        _merge_with = [] if merge_with is None else \
                        ([merge_with] if isinstance(merge_with, str) else merge_with)
        for f in _merge_with:
            cfg.merge(ConfigObj(f))

        # Evaluate the strings if requested
        if evaluate:
            cfg['run'] = _recursive_dict_evaluate(cfg['run'])
        
        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(cfg['run'])


    @classmethod
    def from_dict(cls, cfg):
        return cls(spectrograph=cfg['spectrograph'], pypitdir=cfg['pypitdir'], setup=cfg['setup'],
                   ncpus=cfg['ncpus'], calcheck=cfg['calcheck'], preponly=cfg['preponly'],
                   redname=cfg['redname'], useIDname=cfg['useIDname'],
                   load=LoadPar.from_dict(cfg['load']),
                   directory=DirectoryPar.from_dict(cfg['directory']))


    @staticmethod
    def _available_spectrographs():
        """
        Return the list of allowed spectrographs for pypit reductions.

        .. todo::
            - Listed explicitly for now.  We can do something like what
              is currently done in the code by trolling through a file
              directory.
            - Should allow for a 'user' spectrograph
        """
        return [ 'apf_levy', 'keck_deimos', 'keck_hires', 'keck_lris_blue', 'keck_lris_red',
                 'shane_kast_blue', 'shane_kast_red', 'tng_dolores', 'wht_isis_blue' ]

    

