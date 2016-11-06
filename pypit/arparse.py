from __future__ import (print_function, absolute_import, division, unicode_literals)
from future.utils import iteritems

import collections
import inspect
from os.path import exists as pathexists
from multiprocessing import cpu_count
from os.path import dirname, basename, isfile
from astropy.time import Time
from textwrap import wrap as wraptext
from glob import glob

# Logging
from pypit import ardebug
from pypit import armsgs
debug = ardebug.init()
msgs = armsgs.get_logger()

# Initialize the settings variables
argflag, spect = None, None

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger


class NestedDict(dict):
    """
    A class to generate nested dicts
    """
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


class BaseFunctions(object):
    def __init__(self, defname, savname):
        """ Initialise a settings class. These base functions are used by both spect and argflag settings.

        Parameters
        ----------
        defname : str
          Name of the default settings file to load.
        savname : str
          Name of the file that should save the complete list of the loaded settings.
        """
        self._defname = defname
        self._afout = open(savname, 'w')

    def load_file(self, filename=None):
        """ Load a settings file

        Parameters
        ----------
        filename : str
          Name of the settings file to load. If None, the default settings will be loaded.

        Returns
        -------
        linesarr : list
          Each element of this list contains a line that is read in from a settings file.
        """
        if filename is None:
            msgs.info("Loading default settings")
        else:
            msgs.info("Loading settings")
        try:
            if filename is None:
                lines = open(self._defname, 'r').readlines()
            else:
                lines = open(filename, 'r').readlines()
        except IOError:
            if filename is None:
                msgs.error("Default settings file does not exist:" + msgs.newline() +
                           self._defname)
            else:
                msgs.error("Settings file does not exist:" + msgs.newline() +
                           self._defname)
        linesarr = self.load_lines(lines)
        return linesarr

    def load_lines(self, lines):
        """ Load a lines of a settings file.
        Ignore comment lines (those that start with a #)
        Replace special characters.

        Parameters
        ----------
        lines : list
          A list containing all settings lines to be parsed.

        Returns
        -------
        linesarr : list
          Each element of this list contains a line that is read in from a settings file.
        """
        linesarr = []
        for ll in lines:
            ll = ll.replace("\t", " ").replace("\n", " ")
            if len(ll.strip()) == 0:
                # Nothing on a line
                continue
            elif ll.strip()[0] == '#':
                # A comment line
                continue
            # Remove comments
            ll = ll.split("#")[0].strip()
            linesarr.append(ll.split())
        return linesarr


class BaseArgFlag(BaseFunctions):
    def __init__(self, defname, savname):
        """ Initialize the base functions for the arguments and flags.
        This class contains the functions to load the arguments and flags
        that are common to all reduction programs (i..e ARMED, ARMLSD, etc.)

        Parameters
        ----------
        defname : str
          Name of the default settings file to load.
        savname : str
          Name of the file that should save the complete list of the loaded settings.
        """
        super(BaseArgFlag, self).__init__(defname, savname)
        self._argflag = NestedDict()

    def save(self):
        """ Save the arguments and flags settings used for a given reduction
        """
        def savedict(dct, keylst, keys):
            for (key, value) in iteritems(dct):
                keys += [str(key)]
                if isinstance(value, dict):
                    savedict(value, keylst, keys)
                else:
                    keylst += [str(' ').join(keys) + str(" ") +
                               str("{0}\n".format(value).replace(" ", ""))]
                del keys[-1]

        keylst = []
        savedict(self._argflag.copy(), keylst, [])
        # Sort the list
        keylst = sorted(keylst, key=str.lower)
        # Write the list out in the set order
        for i in range(len(keylst)):
            lstsplt = keylst[i].split(" ")[0]
            # Include a newline after a new keyword
            if i == 0:
                prev = lstsplt
            if prev != lstsplt:
                self._afout.write(str("\n"))
            self._afout.write(keylst[i])
            prev = lstsplt
        self._afout.close()

    def set_param(self, lst, value=None):
        """ Save a single parameter to the argflag dictionary

        Parameters
        ----------
        lst : str or list
          Either a string containing the keyword argument (e.g. 'run redname ARMLSD')
          or a list containing the elements of the keyword argument (e.g. ['run', 'redname']).
          If lst is a list, value must be specified.
        value : any type
          The value of the keyword argument provided by lst (when lst is of type list).
        """
        members = [x for x, y in inspect.getmembers(self, predicate=inspect.ismethod)]
        if type(lst) is str:
            lst = lst.split()
            value = None  # Force the value to be None
        else:
            if type(lst) is list and value is None:
                msgs.error("Couldn't read the value from one of the keyword arguments:" + msgs.newline() +
                           " ".join(lst))
        if value is None:
            func = "_".join(lst[:-1])
            value = "{0:s}".format(str(lst[-1]))
        else:
            func = "_".join(lst)
        if func in members:
            func = "self." + func + "('{0:s}')".format(value)
            eval(func)
        else:
            msgs.error("There appears to be an error on the following parameter:" + msgs.newline() +
                       " ".join(lst) + " {0:s}".format(str(value)))

    def set_paramlist(self, lstall):
        """ Save a list of parameters to the argflag dictionary

        Parameters
        ----------
        lstall : list
          Each element of the lstall is a list containing a full line of a setting
          (e.g. a single element of lstall might look like ['run', 'redname', 'ARMLSD'])
        """
        for ll in range(len(lstall)):
            lst = lstall[ll]
            cnt = 1
            succeed = False
            members = [x for x, y in inspect.getmembers(self, predicate=inspect.ismethod)]
            while cnt < len(lst):
                func = "_".join(lst[:-cnt])
                # Determine if there are options that need to be passed to this function
                options = ""
                nmbr = [[],   # Suffix on 1st arg
                        [],    # Suffix on 2nd arg
                        ["manual"]]    # Suffix on 3rd arg
                ltr = "a"
                for nn in range(len(nmbr)):
                    if nn == 0:
                        ltr = "a"
                    elif nn == 1:
                        ltr = "b"
                    elif nn == 2:
                        ltr = "c"
                    anmbr = nmbr[nn]
                    for aa in anmbr:
                        fspl = func.split("_")
                        if len(fspl) <= nn:
                            continue
                        aatmp = func.split("_")[nn]
                        if aa in aatmp:
                            try:
                                aanum = int(aatmp.lstrip(aa))
                                options += ", {0:s}nmbr={1:d}".format(ltr, aanum)
                            except ValueError:
                                msgs.error("There must be an integer suffix on the {0:s} keyword argument:".format(aa) +
                                           msgs.newline() + " ".join(lst))
                            func = func.replace(aatmp, aa)
                if func in members:
                    func = "self." + func + "('{0:s}'".format(" ".join(lst[-cnt:]))
                    func += options
                    func += ")"
                    eval(func)
                    succeed = True
                    break
                else:
                    cnt += 1
            if not succeed:
                msgs.error("There appears to be an error on the following input line:" + msgs.newline() +
                           " ".join(lst))

    def update(self, v, ll=None):
        """ Update an element of the argflag dictionary

        Parameters
        ----------
        v : any type
          The value of a keyword argument
        ll : list (optional)
          A list containing the keyword (i.e. without the value).
          In general, ll is determined by traceback to the argument
          that called update.
        """
        def ingest(dct, upd):
            """
            Ingest the upd dictionary into dct
            """
            for (kk, vv) in iteritems(upd):
                if isinstance(vv, collections.Mapping):
                    r = ingest(dct.get(kk, {}), vv)
                    dct[kk] = r
                else:
                    dct[kk] = upd[kk]
            return dct
        # First derive a list of the arguments for the keyword to be updated
        if ll is None:
            # update() is called from within this class,
            # so grab the name of the parent function
            ll = inspect.currentframe().f_back.f_code.co_name.split('_')
        # Store a copy of the dictionary to be updated
        dstr = self._argflag.copy()
        # Formulate a dictionary that lists the argument to be updated
        udct = dict({ll[-1]: v})
        for ii in range(1, len(ll)):
            udct = dict({ll[-ii-1]: udct.copy()})
        # Update the master dictionary
        self._argflag = ingest(dstr, udct).copy()

    def arc_combine_match(self, v=-1.0):
        """ Match similar arc frames together? A successful match is found when the frames
        are similar to within N-sigma, where N is the argument of this expression. If v<0,
        arc frames will not be matched.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def arc_combine_method(self, v='weightmean'):
        """ What method should be used to combine the arc frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def arc_combine_reject_cosmics(self, v=-1.0):
        """ Specify the rejection threshold (in standard deviations) for
        cosmic rays when combining the arc frames. If v<0, cosmic rays
        will not be rejected.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def arc_combine_reject_lowhigh(self, v=[0, 0]):
        """ Specify the number of low/high pixels to be rejected when combining
        the arc frames, in the format: [low,high].

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0 or v[1] < 0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def arc_combine_reject_level(self, v=[3.0, 3.0]):
        """ Specify the significance threshold (in standard deviations)
        used to reject deviant pixels when combining the arc frames,
        in the format: [low,high].

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0.0 or v[1] < 0.0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def arc_combine_reject_replace(self, v='maxnonsat'):
        """ What should be done if all pixels are rejected when
        combining the arc frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def arc_combine_satpix(self, v='reject'):
        """ What should be done to saturated pixels when combining the arc frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    def arc_extract_binby(self, v=1.0):
        """ Binning factor to use when extracting 1D arc spectrum. A value of
        1 means that no binning will be performed. This argument does not need
        to be an integer, but it should be >=1.0.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        if v < 1.0:
            msgs.error("The argument of {0:s} must be >= 1.0".format(get_current_name()))
        self.update(v)

    def arc_load_extracted(self, v):
        v = key_bool(v)
        self.update(v)

    def arc_load_calibrated(self, v):
        v = key_bool(v)
        self.update(v)

    def arc_calibrate_detection(self, v):
        v = key_float(v)
        self.update(v)

    def arc_calibrate_IDpixels(self, v):
        v = key_list(v)
        self.update(v)

    def arc_calibrate_IDwaves(self, v):
        v = key_list(v)
        self.update(v)

    def arc_calibrate_lamps(self, v):
        allowed = ['ArI', 'CdI', 'HgI', 'HeI', 'KrI', 'NeI', 'XeI', 'ZnI', 'ThAr']
        v = key_list_allowed(v, allowed)
        self.update(v)

    def arc_calibrate_method(self, v):
        allowed = ['fit', 'simple']
        v = key_allowed(v, allowed)
        self.update(v)

    def arc_calibrate_nfitpix(self, v):
        v = key_int(v)
        if v % 2 == 0:
            msgs.warn("An odd integer is recommended for the argument of {0:s}".format(get_current_name()))
        self.update(v)

    def arc_calibrate_numsearch(self, v):
        v = key_int(v)
        self.update(v)

    def arc_useframe(self, v):
        allowed = ["arc"]
        v = key_none_allowed_filename(v, allowed)
        self.update(v)

    def bias_combine_method(self, v):
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def bias_combine_reject_cosmics(self, v):
        v = key_float(v)
        self.update(v)

    def bias_combine_reject_lowhigh(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0 or v[1] < 0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def bias_combine_reject_level(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0.0 or v[1] < 0.0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def bias_combine_reject_replace(self, v):
        """ What should be done if all pixels are rejected when
        combining the bias frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def bias_combine_satpix(self, v):
        """ What should be done to saturated pixels when combining the arc frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    def bias_useoverscan(self, v):
        v = key_bool(v)
        self.update(v)

    def bias_useframe(self, v):
        allowed = ['bias', 'overscan', 'dark', 'none']
        if "," in v:
            # Must be a list - multiple options
            v = load_list(v)
            if "none" in v:
                msgs.error("'none' cannot be a list element for the argument of {0:s}".format(get_current_name()))
            if ("bias" in v) and ("dark" in v):
                msgs.error("'bias' and 'dark' cannot both be a list elements of {0:s}".format(get_current_name()))
            for ll in v:
                if ll not in allowed:
                    msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                               ", ".join(allowed) + " or a list containing these options")
            if "overscan" in v:
                self.bias_useoverscan("true")
            if "bias" in v:
                v = "bias"
            elif "dark" in v:
                v = "dark"
            else:
                msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                           ", ".join(allowed) + " or a list containing these options")
        elif v.lower() == "none":
            v = None
        elif v.lower() not in allowed:
            msgs.info("Assuming the following is the name of a bias frame:" + msgs.newline() + v)
        self.update(v)

    def output_overwrite(self, v):
        v = key_bool(v)
        self.update(v)

    def output_sorted(self, v):
        """
        The argument of this keyword argument should have no extension
        """
        v = key_none(v)
        self.update(v)

    def output_verbosity(self, v):
        v = key_int(v)
        if (v < 0) or (v > 2):
            msgs.error("The verbosity can only take values between 0 (minimum) and 2 (maximum)" + msgs.newline() +
                       "Please change the argument of {0:s}".format(get_current_name()))
        self.update(v)

    def pixelflat_combine_match(self, v):
        v = key_float(v)
        self.update(v)

    def pixelflat_combine_method(self, v):
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def pixelflat_combine_reject_cosmics(self, v):
        v = key_float(v)
        self.update(v)

    def pixelflat_combine_reject_lowhigh(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0 or v[1] < 0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def pixelflat_combine_reject_level(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0.0 or v[1] < 0.0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def pixelflat_combine_reject_replace(self, v):
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def pixelflat_combine_satpix(self, v):
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    def pixelflat_norm_recnorm(self, v):
        v = key_bool(v)
        self.update(v)

    def pixelflat_useframe(self, v):
        allowed = ['pixelflat']
        v = key_none_allowed_filename(v, allowed)
        self.update(v)

    def reduce_badpix(self, v):
        v = key_bool(v)
        self.update(v)

    def reduce_calibrate_nonlinear(self, v):
        v = key_bool(v)
        self.update(v)

    def reduce_calibrate_refframe(self, v):
        allowed = ['geocentric', 'heliocentric', 'barycentric']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_calibrate_wavelength(self, v):
        allowed = ['air', 'vacuum', 'none']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_flatfield_method(self, v):
        allowed = ['polyscan']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_flatfield_params(self, v):
        v = key_list(v)
        self.update(v)

    def reduce_flatfield_perform(self, v):
        v = key_bool(v)
        self.update(v)

    def reduce_flatfield_useframe(self, v):
        allowed = ['pixelflat', 'slitflat']
        v = key_allowed_filename(v, allowed)
        self.update(v)

    def reduce_trace_useframe(self, v):
        allowed = ['trace', 'slitflat', 'science']
        v = key_allowed_filename(v, allowed)
        self.update(v)

    def reduce_flexure_maxshift(self, v):
        v = key_int(v)
        self.update(v)

    def reduce_flexure_method(self, v):
        allowed = ['none', 'boxcar', 'slitcen']
        v = key_none_allowed(v, allowed)
        self.update(v)

    def reduce_flexure_spectrum(self, v):
        if v.lower() == 'none':
            v = None
        else:
            if not pathexists(self._argflag['run']['pypitdir'] + 'data/sky_spec/' + v):
                files_sky = glob(self._argflag['run']['pypitdir'] + 'data/sky_spec/*.fits')
                skyfiles = ""
                for i in files_sky:
                    skyfiles += msgs.newline() + "  - " + str(i.split("/")[-1])
                msgs.error("The following archive sky spectrum file does not exist:" + msgs.newline() +
                           "  " + v + msgs.newline() + msgs.newline() + "Please use one of the files listed below:" +
                           skyfiles)
        self.update(v)

    def reduce_masters_file(self, v):
        if v.lower() == 'none':
            v = ''
        self.update(v)

    def reduce_masters_loaded(self, v):
        v = key_list(v)
        self.update(v)

    def reduce_masters_reuse(self, v):
        v = key_bool(v)
        self.update(v)

    def reduce_masters_setup(self, v):
        if v.lower() == 'none':
            v = ''
        self.update(v)

    def reduce_overscan_method(self, v):
        allowed = ['polynomial', 'savgol', 'median']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_overscan_params(self, v):
        v = key_list(v)
        self.update(v)

    def reduce_pixel_locations(self, v):
        if v.lower() == "none":
            v = None
        elif v.split(".")[-1] == "fits":
            pass
        elif v.split(".")[-2] == "fits" and v.split(".")[-1] == "gz":
            pass
        else:
            msgs.error("The argument of {0:s} must be 'None' or a fits file".format(get_current_name()))
        self.update(v)

    def reduce_pixel_size(self, v):
        v = key_float(v)
        self.update(v)

    def reduce_skysub_bspline_everyn(self, v):
        v = key_int(v)
        self.update(v)

    def reduce_skysub_method(self, v):
        allowed = ['bspline', 'polyscan']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_skysub_perform(self, v):
        v = key_bool(v)
        self.update(v)

    def reduce_trim(self, v):
        v = key_bool(v)
        self.update(v)

    def run_calcheck(self, v):
        v = key_bool(v)
        self.update(v)

    def run_directory_master(self, v):
        self.update(v)

    def run_directory_qa(self, v):
        self.update(v)

    def run_directory_science(self, v):
        self.update(v)

    def run_load_settings(self, v):
        if v.lower() == "none":
            v = None
        elif not isfile(v):
                msgs.error("The argument of {0:s} must be a PYPIT settings file".format(get_current_name()) +
                           msgs.newline() + "or 'None'. The following file does not exist:" + msgs.newline() + v)
        self.update(v)

    def run_load_spect(self, v):
        if v.lower() == "none":
            v = None
        elif isfile(v):
            pass
        elif isfile(v+".spect"):
            v += ".spect"
        else:
             msgs.error("The argument of {0:s} must be a PYPIT spectrograph settings".format(get_current_name()) +
                        msgs.newline() + "file or 'None'. The following file does not exist:" + msgs.newline() + v)
        self.update(v)

    def run_ncpus(self, v):
        if 'ncpus' in self._argflag['run'].keys():
            curcpu = self._argflag['run']['ncpus']
        else:
            curcpu = 0
        cpucnt = cpu_count()
        if v == 'all':
            # Use all available cpus
            v = cpucnt
            if v != curcpu:
                msgs.info("Setting {0:d} CPUs".format(v))
        elif v is None:
            # Use all but 1 available cpus
            v = cpucnt-1
            if v != curcpu:
                msgs.info("Setting {0:d} CPUs".format(v))
        else:
            try:
                v = int(v)
                if v > cpucnt:
                    msgs.warn("You don't have {0:d} CPUs!".format(v))
                    v = cpucnt
                elif v < 0:
                    v += cpucnt
                if v != curcpu:
                    msgs.info("Setting {0:d} CPUs".format(v))
            except ValueError:
                msgs.error("Incorrect argument given for number of CPUs" + msgs.newline() +
                           "Please choose from -" + msgs.newline() +
                           "all, 1..."+str(cpucnt))
                if cpucnt == 1:
                    if cpucnt != curcpu:
                        msgs.info("Setting 1 CPU")
                    v = 1
                else:
                    v = cpu_count()-1
                    if v != curcpu:
                        msgs.info("Setting {0:d} CPUs".format(v))
        self.update(v)

    def run_preponly(self, v):
        v = key_bool(v)
        self.update(v)

    def run_progname(self, v):
        self.update(v)

    def run_pypitdir(self, v):
        self.update(v)

    def run_qa(self, v):
        v = key_bool(v)
        self.update(v)

    def run_redname(self, v):
        self.update(v)

    def run_spectrograph(self, v):
        # Check that v is allowed
        stgs_arm = glob(dirname(__file__)+"/settings/settings.arm*")
        stgs_all = glob(dirname(__file__)+"/settings/settings.*")
        stgs_spc = list(set(stgs_arm) ^ set(stgs_all))
        spclist = [basename(stgs_spc[0]).split(".")[-1].lower()]
        for i in range(1, len(stgs_spc)):
            spclist += [basename(stgs_spc[i]).split(".")[-1].lower()]
        # Check there are no duplicate names
        if len(spclist) != len(set(spclist)):
            msgs.bug("Duplicate settings files found")
            msgs.error("Cannot continue with an ambiguous settings file")
        # Check the settings file exists
        if v.lower() not in spclist:
            msgs.error("Settings do not exist for the {0:s} spectrograph".format(v.lower()) + msgs.newline() +
                       "Please use one of the following spectrograph settings:" + msgs.newline() +
                       wraptext(", ".join(spclist), width=60))
        self.update(v)
        return

    def run_stopcheck(self, v):
        v = key_bool(v)
        self.update(v)

    def run_useIDname(self, v):
        v = key_bool(v)
        self.update(v)

    def science_extraction_manual(self, cnmbr=1, frame="none", params="[1,1000,500,[10,10]]"):
        # Send parameters away to individual arguments
        self.science_extraction_manual_frame(frame, cnmbr=cnmbr)
        self.science_extraction_manual_params(params, cnmbr=cnmbr)

    def science_extraction_manual_frame(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        if v.lower() == "none":
            v = None
        elif ".fits" not in v:
            msgs.error("The argument of {0:s} must be a fits file".format(cname))
        self.update(v, ll=cname.split('_'))

    def science_extraction_manual_params(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        if v.lower() == "none":
            v = None
        else:
            try:
                v = eval(v)
            except:
                msgs.error("The argument of {0:s} must be a 4 parameter list.".format(cname))
            if len(v) != 4:
                msgs.error("The argument of {0:s} must be a 4 parameter list.".format(cname))
        self.update(v, ll=cname.split('_'))

    def science_extraction_maxnumber(self, v):
        v = key_int(v)
        self.update(v)

    def science_extraction_profile(self, v):
        v = key_allowed(v)
        self.update(v)

    def science_extraction_reuse(self, v):
        v = key_bool(v)
        self.update(v)

    def slitflat_combine_match(self, v):
        v = key_float(v)
        self.update(v)

    def slitflat_combine_method(self, v):
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def slitflat_combine_reject_cosmics(self, v):
        v = key_float(v)
        self.update(v)

    def slitflat_combine_reject_lowhigh(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0 or v[1] < 0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def slitflat_combine_reject_level(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0.0 or v[1] < 0.0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def slitflat_combine_reject_replace(self, v):
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def slitflat_combine_satpix(self, v):
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    def slitflat_norm_recnorm(self, v):
        v = key_bool(v)
        self.update(v)

    def slitflat_useframe(self, v):
        allowed = ['slitflat']
        v = key_none_allowed_filename(v, allowed)
        self.update(v)

    def trace_combine_match(self, v):
        v = key_float(v)
        self.update(v)
        return

    def trace_combine_method(self, v):
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_combine_reject_cosmics(self, v):
        v = key_float(v)
        self.update(v)

    def trace_combine_reject_lowhigh(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0 or v[1] < 0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_combine_reject_level(self, v):
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] <= 0.0 or v[1] <= 0.0:
            msgs.error("The list values of argument {0:s} must be > 0.0".format(get_current_name()))
        self.update(v)

    def trace_combine_reject_replace(self, v):
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_combine_satpix(self, v):
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_slits_diffpolyorder(self, v):
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_dispersion_direction(self, v):
        v = key_int(v)
        if v != 0 and v != 1:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       "0 or 1 (if the dispersion axis is along a row or column respectively)")
        self.update(v)

    def trace_slits_fracignore(self, v):
        v = key_float(v)
        if v < 0.0 or v > 1.0:
            msgs.error("The argument of {0:s} must be between 0 and 1".format(get_current_name()))
        self.update(v)

    def trace_slits_function(self, v):
        allowed = ['polynomial', 'legendre', 'chebyshev']
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_slits_maxgap(self, v):
        if v.lower() == "none":
            v = None
        else:
            try:
                v = int(v)
            except ValueError:
                msgs.error("The argument of {0:s} must be of type int, or set to 'none'".format(get_current_name()))
            if v <= 1:
                msgs.error("The argument of {0:s} must be > 1 to set the maximum slit gap".format(get_current_name()))
        self.update(v)

    def trace_slits_number(self, v):
        if v.lower() == "auto":
            v = -1
        else:
            try:
                v = int(v)
            except ValueError:
                msgs.error("The argument of {0:s} must be of type int, or set to 'auto'".format(get_current_name()))
            if v == 0 or v == -1:
                v = -1
            elif v < -1:
                msgs.error("The argument of {0:s} must be >= 1 to manually set the number of slits,".format(get_current_name()) + msgs.newline() +
                           "or can be set to -1 (or 'auto') if you wish PYPIT to find slits automatically.")
        self.update(v)

    def trace_slits_polyorder(self, v):
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_pca_type(self, v):
        allowed = ['pixel', 'order']
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_slits_pca_params(self, v):
        v = key_list(v)
        self.update(v)

    def trace_slits_pca_extrapolate_pos(self, v):
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_pca_extrapolate_neg(self, v):
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_sigdetect(self, v):
        v = key_float(v)
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0".format(get_current_name()))
        self.update(v)

    def trace_slits_single(self, v):
        v = key_list(v)
        self.update(v)

    def trace_slits_tilts_idsonly(self, v):
        v = key_bool(v)
        self.update(v)

    def trace_slits_tilts_method(self, v):
        allowed = ['PCA', 'spline', 'spca', 'interp', 'perp', 'zero']
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_slits_tilts_params(self, v):
        v = key_list(v)
        self.update(v)

    def trace_slits_tilts_disporder(self, v):
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_tilts_order(self, v):
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_useframe(self, v):
        allowed = ['trace', 'science']
        v = key_none_allowed_filename(v, allowed)
        self.update(v)


class BaseSpect(BaseFunctions):
    def __init__(self, defname, savname):
        super(BaseSpect, self).__init__(defname, savname)
        self._spect = NestedDict()
        self._settings = []
        self.set_default()

    def save(self):
        """
        Save the settings used for this reduction
        """
        def savedict(dct, keylst, keys):
            for (key, value) in iteritems(dct):
                keys += [str(key)]
                if isinstance(value, dict):
                    savedict(value, keylst, keys)
                else:
                    keylst += [str(' ').join(keys) + str(" ") +
                               str("{0}\n".format(value).replace(" ", ""))]
                del keys[-1]

        keylst = []
        savedict(self._spect.copy(), keylst, [])
        # Sort the list
        keylst = sorted(keylst, key=str.lower)
        # Write the list out in the set order
        for i in range(len(keylst)):
            lstsplt = keylst[i].split(" ")[0]
            # Include a newline after a new keyword
            if i == 0:
                prev = lstsplt
            if prev != lstsplt:
                self._afout.write(str("\n"))
            self._afout.write(keylst[i])
            prev = lstsplt
        self._afout.close()
        return

    def set_default(self):
        """ Set some arguments that are not used in the settings file
        """
        self.update([], ll="set_arc".split("_"))
        self.update([], ll="set_bias".split("_"))
        self.update([], ll="set_pixelflat".split("_"))
        self.update([], ll="set_science".split("_"))
        self.update([], ll="set_slitflat".split("_"))
        self.update([], ll="set_standard".split("_"))
        self.update([], ll="set_trace".split("_"))
        self.update([], ll="arc_index".split("_"))
        self.update([], ll="bias_index".split("_"))
        self.update([], ll="pixelflat_index".split("_"))
        self.update([], ll="science_index".split("_"))
        self.update([], ll="slitflat_index".split("_"))
        self.update([], ll="standard_index".split("_"))
        self.update([], ll="trace_index".split("_"))
        return

    def set_param(self, lst, value=None):
        members = [x for x, y in inspect.getmembers(self, predicate=inspect.ismethod)]
        if type(lst) is str:
            lst = lst.split()
        if value is None:
            func = "_".join(lst[:-1])
            value = "{0:s}".format(str(lst[-1]))
        else:
            func = "_".join(lst)
        if func in members:
            func = "self." + func + "('{0:s}')".format(value)
            eval(func)
        else:
            msgs.error("There appears to be an error on the following parameter:" + msgs.newline() +
                       " ".join(lst) + " {0:s}".format(str(value)))
        return

    def set_paramlist(self, lstall):
        frmtyp = ["standard", "bias", "pixelflat", "trace", "slitflat", "arc"]
        for ll in range(len(lstall)):
            lst = lstall[ll]
            cnt = 1
            succeed = False
            members = [x for x, y in inspect.getmembers(self, predicate=inspect.ismethod)]
            while cnt < len(lst):
                func = "_".join(lst[:-cnt])
                # Determine if there are options that need to be passed to this function
                options = ""
                nmbr = [["det"],   # Suffix on 1st arg
                        ["datasec", "oscansec", "lampname", "lampstat", "headext"],    # Suffix on 2nd arg
                        ["condition"]]    # Suffix on 3rd arg
                ltr = "a"
                for nn in range(len(nmbr)):
                    if nn == 0:
                        ltr = "a"
                    elif nn == 1:
                        ltr = "b"
                    elif nn == 2:
                        ltr = "c"
                    anmbr = nmbr[nn]
                    for aa in anmbr:
                        fspl = func.split("_")
                        if len(fspl) <= nn:
                            continue
                        aatmp = func.split("_")[nn]
                        if aa in aatmp:
                            try:
                                aanum = int(aatmp.lstrip(aa))
                                options += ", {0:s}nmbr={1:d}".format(ltr, aanum)
                            except ValueError:
                                msgs.error("There must be an integer suffix on the {0:s} keyword argument:".format(aa) +
                                           msgs.newline() + " ".join(lst))
                            func = func.replace(aatmp, aa)
                # Now test if this is a function
                if func in members:
                    func = "self." + func + "('{0:s}'".format(" ".join(lst[-cnt:]))
                    func += options
                    func += ")"
                    eval(func)
                    succeed = True
                    break
                else:
                    cnt += 1
            if not succeed:
                # Try a few manual options
                if lst[0] == "check":
                    self.update(" ".join(lst[2:]), ll=lst[:2])
                elif lst[0] in frmtyp and lst[1] == "match":
                    self.update(lst[-1], ll=lst[:-1])
                else:
                    msgs.error("There appears to be an error on the following input line:" + msgs.newline() +
                               " ".join(lst))
        return

    def settings(self, v):
        self._settings.append(v.split())
        return

    def update(self, v, ll=None):
        """
        Update an element in spect
        """
        def ingest(dct, upd):
            """
            Ingest the upd dictionary into dct
            """
            for (kk, vv) in iteritems(upd):
                if isinstance(vv, collections.Mapping):
                    r = ingest(dct.get(kk, {}), vv)
                    dct[kk] = r
                else:
                    dct[kk] = upd[kk]
            return dct

        # First derive a list of the arguments for the keyword to be updated
        if ll is None:
            # update() is called from within this class,
            # so grab the name of the parent function
            ll = inspect.currentframe().f_back.f_code.co_name.split('_')
        # Store a copy of the dictionary to be updated
        dstr = self._spect.copy()
        # Formulate a dictionary that lists the argument to be updated
        udct = dict({ll[-1]: v})
        for ii in range(1, len(ll)):
            udct = dict({ll[-ii - 1]: udct.copy()})
        # Update the master dictionary
        self._spect = ingest(dstr, udct).copy()
        return

    def arc_canbe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if "," in v:
                v = v.strip("()[]").split(",")
            else:
                v = [v]
        # Update argument
        self.update(v)

    def arc_check_condition(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        # Check that v is allowed
        text = v.strip().replace('_', ' ')
        if ',' in text and text[0:2] != '%,':
            # There are multiple possibilities - split the text
            v = text.split(',')
        else:
            v = text
        # Update argument
        self.update(v, ll=cname.split('_'))

    def arc_idname(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def arc_number(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)

    def bias_canbe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if "," in v:
                v = v.strip("()[]").split(",")
            else:
                v = [v]
        # Update argument
        self.update(v)

    def bias_check_condition(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        # Check that v is allowed
        text = v.strip().replace('_', ' ')
        if ',' in text and text[0:2] != '%,':
            # There are multiple possibilities - split the text
            v = text.split(',')
        else:
            v = text
        # Update argument
        self.update(v, ll=cname.split('_'))

    def bias_idname(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def bias_number(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)

    def det_datasec(self, v, anmbr=1, bnmbr=1):
        cname = get_nmbr_name(anmbr=anmbr, bnmbr=bnmbr)
        # Check that v is allowed
        try:
            v = load_sections(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be detector section".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_oscansec(self, v, anmbr=1, bnmbr=1):
        cname = get_nmbr_name(anmbr=anmbr, bnmbr=bnmbr)
        # Check that v is allowed
        try:
            v = load_sections(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be detector section".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_darkcurr(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_gain(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = v.split(",")
            for i in range(len(v)):
                v[i] = float(v[i])
                if v[i] <= 0.0:
                    msgs.error("Each argument of {0:s} must be > 0.0".format(cname))
        except ValueError:
            msgs.error("Each argument of {0:s} must be of type float".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_ronoise(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = v.split(",")
            for i in range(len(v)):
                v[i] = float(v[i])
                if v[i] <= 0.0:
                    msgs.error("Each argument of {0:s} must be > 0.0".format(cname))
        except ValueError:
            msgs.error("Each argument of {0:s} must be of type float".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_nonlinear(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v <= 0.0 or v > 1.0:
            msgs.error("The argument of {0:s} must be > 0.0 and <= 1.0".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_numamplifiers(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(cname))
        if v <= 0:
            msgs.error("The argument of {0:s} must be >= 1".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_saturation(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0.0".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_suffix(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed

        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_xgap(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_ygap(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def det_ysize(self, v, anmbr=1):
        cname = get_nmbr_name(anmbr=anmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0.0".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def fits_calwin(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        if v <= 0.0:
            msgs.error("The calibration time window must be > 0.0")
        # Update argument
        self.update(v)

    def fits_dataext(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The fits data extension number must be >= 0")
        # Update argument
        self.update(v)

    def fits_headext(self, v, bnmbr=1):
        cname = get_nmbr_name(bnmbr=bnmbr)
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(cname))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(cname))
        # Update argument
        self.update(v, ll=cname.split('_'))

    def fits_numhead(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v <= 0:
            msgs.error("The number of fits headers must be >= 1")
        # Update argument
        self.update(v)

    def fits_numlamps(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The number of lamps must be >= 0")
        # Update argument
        self.update(v)

    def fits_timeunit(self, v):
        # Check that v is allowed
        allowed = ['s', 'm', 'h']
        astropy_allowed = Time.FORMATS.keys()
        if v not in allowed and v not in astropy_allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed) + "or one of the astropy Time formats:" + msgs.newline() +
                       ", ".join(astropy_allowed))
        # Update argument
        self.update(v)

    def keyword_ra(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_dec(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_target(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_airmass(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_binning(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_binningspatial(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_binningspectral(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_date(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_decker(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_detrot(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_dispname(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)


    def keyword_dispangle(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_equinox(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_exptime(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_filter1(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_filter2(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_hatch(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_idname(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_lamps(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_lampname(self, v, bnmbr=1):
        cname = get_nmbr_name(bnmbr=bnmbr)
        # Check that v is allowed
        if "." not in v:
            # User must have passed the name of the lamp
            pass
        else:
            # Get the name of the lamp from the header
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(cname) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v, ll=cname.split('_'))

    def keyword_lampstat(self, v, bnmbr=1):
        cname = get_nmbr_name(bnmbr=bnmbr)
        # Check that v is allowed
        if "." not in v:
            # User must have specified the status
            pass
        else:
            # Get the lamp status from the header
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(cname) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v, ll=cname.split('_'))

    def keyword_naxis0(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_naxis1(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_slitwid(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_slitlen(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def keyword_time(self, v):
        # Check that v is allowed
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        # Update argument
        self.update(v)

    def mosaic_camera(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def mosaic_elevation(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)

    def mosaic_latitude(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)

    def mosaic_longitude(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)

    def mosaic_ndet(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)

    def mosaic_minexp(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)

    def mosaic_reduction(self, v):
        # Check that v is allowed
        allowed = ['ARMLSD', 'ARMED']
        if v.upper() not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v.upper())

    def pixelflat_canbe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if "," in v:
                v = v.strip("()[]").split(",")
            else:
                v = [v]
        # Update argument
        self.update(v)

    def pixelflat_check_condition(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        # Check that v is allowed
        text = v.strip().replace('_', ' ')
        if ',' in text and text[0:2] != '%,':
            # There are multiple possibilities - split the text
            v = text.split(',')
        else:
            v = text
        # Update argument
        self.update(v, ll=cname.split('_'))

    def pixelflat_idname(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def pixelflat_lscomb(self, v):
        # Check that v is allowed
        v = v.lower()
        if v == "true":
            v = True
        elif v == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} must be True or False".format(get_current_name()))
        # Update argument
        self.update(v)

    def pixelflat_number(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)

    def science_canbe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if "," in v:
                v = v.strip("()[]").split(",")
            else:
                v = [v]
        # Update argument
        self.update(v)

    def science_check_condition(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        # Check that v is allowed
        text = v.strip().replace('_', ' ')
        if ',' in text and text[0:2] != '%,':
            # There are multiple possibilities - split the text
            v = text.split(',')
        else:
            v = text
        # Update argument
        self.update(v, ll=cname.split('_'))

    def science_idname(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def set_arc(self, v):
        # Check that v is allowed
        v = load_list(v)
        v = self._spect['set']['arc'] + v
        # Update argument
        self.update(v)
        return

    def set_bias(self, v):
        # Check that v is allowed
        v = load_list(v)
        v = self._spect['set']['bias'] + v
        # Update argument
        self.update(v)
        return

    def set_pixelflat(self, v):
        # Check that v is allowed
        v = load_list(v)
        v = self._spect['set']['pixelflat'] + v
        # Update argument
        self.update(v)
        return

    def set_science(self, v):
        # Check that v is allowed
        v = load_list(v)
        v = self._spect['set']['science'] + v
        # Update argument
        self.update(v)
        return

    def set_slitflat(self, v):
        # Check that v is allowed
        v = load_list(v)
        v = self._spect['set']['slitflat'] + v
        # Update argument
        self.update(v)
        return

    def set_standard(self, v):
        # Check that v is allowed
        v = load_list(v)
        v = self._spect['set']['standard'] + v
        # Update argument
        self.update(v)
        return

    def set_trace(self, v):
        # Check that v is allowed
        v = load_list(v)
        v = self._spect['set']['trace'] + v
        # Update argument
        self.update(v)
        return

    def slitflat_canbe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if "," in v:
                v = v.strip("()[]").split(",")
            else:
                v = [v]
        # Update argument
        self.update(v)

    def slitflat_check_condition(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        # Check that v is allowed
        text = v.strip().replace('_', ' ')
        if ',' in text and text[0:2] != '%,':
            # There are multiple possibilities - split the text
            v = text.split(',')
        else:
            v = text
        # Update argument
        self.update(v, ll=cname.split('_'))

    def slitflat_idname(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def slitflat_lscomb(self, v):
        # Check that v is allowed
        v = v.lower()
        if v == "true":
            v = True
        elif v == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} must be True or False".format(get_current_name()))
        # Update argument
        self.update(v)

    def slitflat_number(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)

    def standard_canbe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if "," in v:
                v = v.strip("()[]").split(",")
            else:
                v = [v]
        # Update argument
        self.update(v)

    def standard_check_condition(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        # Check that v is allowed
        text = v.strip().replace('_', ' ')
        if ',' in text and text[0:2] != '%,':
            # There are multiple possibilities - split the text
            v = text.split(',')
        else:
            v = text
        # Update argument
        self.update(v, ll=cname.split('_'))

    def standard_idname(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def standard_number(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)

    def trace_canbe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if "," in v:
                v = v.strip("()[]").split(",")
            else:
                v = [v]
        # Update argument
        self.update(v)

    def trace_check_condition(self, v, cnmbr=1):
        cname = get_nmbr_name(cnmbr=cnmbr)
        # Check that v is allowed
        text = v.strip().replace('_', ' ')
        if ',' in text and text[0:2] != '%,':
            # There are multiple possibilities - split the text
            v = text.split(',')
        else:
            v = text
        # Update argument
        self.update(v, ll=cname.split('_'))

    def trace_idname(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def trace_lscomb(self, v):
        # Check that v is allowed
        v = v.lower()
        if v == "true":
            v = True
        elif v == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} must be True or False".format(get_current_name()))
        # Update argument
        self.update(v)

    def trace_number(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)


class ARMLSD(BaseArgFlag):

    def reduce_calibrate_flux(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)

    def reduce_flexure_maxshift(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)

    def reduce_flexure_spec(self, v):
        # Check that v is allowed
        allowed = ['boxcar', 'slit_cen', 'none']
        v = v.lower()
        if v == "none":
            v = None
        elif v in allowed:
            pass
        else:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)


class ARMLSD_spect(BaseSpect):

    def keyword_dichroic(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)


class ARMED_spect(BaseSpect):

    def keyword_echangle(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                vspl = v.split(".")
                int(vspl[0])
            except ValueError:
                msgs.error("The argument of {0:s} must be of the form:".format(get_current_name()) + msgs.newline() +
                           "##.NAME" + msgs.newline() +
                           "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                           "and NAME is the header keyword name")
        # Update argument
        self.update(v)


def init(afclass, spclass):
    """
    Initialize the settings
    ----------
    afclass : class
      Class of arguments and flags
    spclass : class
      Class of spectrograph settings
    """
    global argflag
    global spect
    argflag = afclass.__dict__['_argflag']
    spect = spclass.__dict__['_spect']
    return


def get_argflag_class(init=None):
    """
    Get the Arguments and Flags
    ----------
    init : tuple
      For instantiation

    Returns
    -------
    argflag : Arguments and Flags
    """
    try:
        defname = glob(dirname(__file__))[0] + "/settings/settings." + init[0].lower()
        return eval(init[0]+"(defname='{0:s}', savname='{1:s}.settings')".format(defname, init[1]))
    except RuntimeError:
        msgs.error("Reduction type '{0:s}' is not allowed".format(init))


def get_spect_class(init):
    """
    Get the Spectrograph settings class
    ----------
    init : tuple
      For instantiation

    Returns
    -------
    spect_class : Class of spectrograph settings
    """
    try:
        defname = glob(dirname(__file__))[0] + "/settings/settings." + init[1].lower()
        return eval(init[0]+"_spect(defname='{0:s}', savname='{1:s}.spect')".format(defname, init[2]))
    except RuntimeError:
        msgs.error("{0:s} is not implemented yet".format(init[1]))


def get_current_name():
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    return "'" + " ".join(ll) + "'"


def get_nmbr_name(anmbr=None, bnmbr=None, cnmbr=None):
    cspl = inspect.currentframe().f_back.f_code.co_name.split('_')
    if anmbr is not None:
        cspl[0] += "{0:02d}".format(anmbr)
    if bnmbr is not None:
        cspl[1] += "{0:02d}".format(bnmbr)
    if cnmbr is not None:
        cspl[2] += "{0:02d}".format(cnmbr)
    return "_".join(cspl)


def load_sections(string):
    """
    From the input string, return the coordinate sections

    Parameters
    ----------
    string : str
      character string of the form [x1:x2,y1:y2]
      x1 = left pixel
      x2 = right pixel
      y1 = bottom pixel
      y2 = top pixel

    Returns
    -------
    sections : list (or None)
      the detector sections
    """
    xyrng = string.strip('[]()').split(',')
    if xyrng[0] == ":":
        xyarrx = [0, 0]
    else:
        xyarrx = xyrng[0].split(':')
    if xyrng[1] == ":":
        xyarry = [0, 0]
    else:
        xyarry = xyrng[1].split(':')
    return [[int(xyarrx[0]), int(xyarrx[1])], [int(xyarry[0]), int(xyarry[1])]]


def get_dnum(det):
    """Convert a detector index into a string used by the settings dictionary

    Parameters
    ----------
    det : int
      Detector index

    Returns
    -------
    dnum : str
      A string used by the settings dictionary
    """
    return 'det{0:02d}'.format(det)


def key_allowed(v, allowed):
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    v = v.lower()
    if v not in allowed:
        msgs.error("The argument of {0:s} must be one of".format(func_name) + msgs.newline() +
                   ", ".join(allowed))
    return v


def key_allowed_filename(v, allowed):
    vt = v.lower()
    if vt not in allowed:
        msgs.warn("Assuming the following is the name of a file:" + msgs.newline() + v)
    else:
        v = vt
    return v


def key_bool(v):
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if v.lower() == "true":
        v = True
    elif v.lower() == "false":
        v = False
    else:
        msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(func_name))
    return v


def key_float(v):
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    try:
        v = float(v)
    except ValueError:
        msgs.error("The argument of {0:s} must be of type float".format(func_name))
    return v


def key_int(v):
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    try:
        v = int(v)
    except ValueError:
        msgs.error("The argument of {0:s} must be of type int".format(func_name))
    return v


def key_list(strlist):
    # Check if the input array is a null list
    if strlist == "[]" or strlist == "()":
        return []
    # Remove outer brackets and split by commas
    temp = strlist.lstrip('([').rstrip(')]').split(',')
    addarr = []
    # Find the type of the array elements
    for i in temp:
        if i.lower() == 'none':
            # None type
            addarr += [None]
        elif i.lower() == 'true' or i.lower() == 'false':
            # bool type
            addarr += [i.lower() in ['true']]
        elif ',' in i:
            # a list
            addarr += i.lstrip('([').rstrip('])').split(',')
            msgs.bug("nested lists could cause trouble if elements are not strings!")
        elif '.' in i:
            try:
                # Might be a float
                addarr += [float(i)]
            except ValueError:
                # Must be a string
                addarr += [i]
        else:
            try:
                # Could be an integer
                addarr += [int(i)]
            except ValueError:
                # Must be a string
                addarr += [i]
    return addarr


def key_list_allowed(v, allowed):
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    v = key_list(v)
    for ll in v:
        if ll not in allowed:
            msgs.error("The allowed list does not include: {0:s}".format(ll) + msgs.newline() +
                       "Please choose one of the following:" + msgs.newline() +
                       ", ".join(allowed) + msgs.newline() +
                       "for the argument of {0:s}".format(func_name))
    return v


def key_none_allowed(v, allowed):
    """ First check if a keyword is None, then see if it is in the allowed list
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if v.lower() == "none":
        v = None
    elif v.lower() in allowed:
        for i in allowed:
            if v.lower() == i:
                v = i
                break
    else:
        msgs.error("The argument of {0:s} must be one of".format(func_name) + msgs.newline() +
                   ", ".join(allowed))
    return v


def key_none_allowed_filename(v, allowed):
    """ First check if a keyword is None, then see if it is in the allowed list,
    and finally assume that it is the name of a file if not in the list
    """
    if v.lower() == "none":
        v = None
    elif v.lower() in allowed:
        for i in allowed:
            if v.lower() == i:
                v = i
                break
    else:
        msgs.info("Assuming the following is the name of a file:" + msgs.newline() + v)
    return v


def key_none(v):
    if v.lower() == "none":
        v = None
    return v


def combine_methods():
    """ The methods that can be used to combine a set of frames into a master frame
    """
    methods = ['mean', 'median', 'weightmean']
    return methods


def combine_replaces():
    """ The options that can be used to replace rejected pixels when combining a set of frames
    """
    methods = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
    return methods


def combine_satpixs():
    methods = ['reject', 'force', 'nothing']
    return methods