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

    def load_file(self, filename=None, base=False):
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
                if base:
                    if isinstance(self, BaseArgFlag):
                        basefile = glob(dirname(__file__))[0] + "/data/settings/settings.baseargflag"
                    elif isinstance(self, BaseSpect):
                        basefile = glob(dirname(__file__))[0] + "/data/settings/settings.basespect"
                    else:
                        msgs.error("No base for this class")
                    msgs.info("Loading base settings from {:s}".format(basefile.split('/')[-1]))
                    lines = open(basefile, 'r').readlines()
                else:
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

    def init_param(self):
        """ Initialze the parameter list
        """
        # Base
        lines = self.load_file(base=True)
        self.set_paramlist(lines)
        # Pipeline specific
        lines = self.load_file()
        self.set_paramlist(lines)

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

    def arc_load_calibrated(self, v):
        """ If the extracted arc have previously been calibrated and saved, load the calibration files?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def arc_load_extracted(self, v):
        """ If the master arc has previously been extracted and saved, load the 1D extractions?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def arc_calibrate_detection(self, v):
        """ How significant should the arc line detections be (in units of a standard deviation)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def arc_calibrate_IDpixels(self, v):
        """ Manually set the pixels to be identified

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def arc_calibrate_IDwaves(self, v):
        """ Manually set the corresponding ID wavelengths

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def arc_calibrate_lamps(self, v):
        """ Name of the ions used for the wavelength calibration

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['ArI', 'CdI', 'HgI', 'HeI', 'KrI', 'NeI', 'XeI', 'ZnI', 'ThAr']
        v = key_list_allowed(v, allowed)
        self.update(v)

    def arc_calibrate_method(self, v):
        """ What method should be used to fit the individual arc lines.
        The 'fit' option is perhaps the most accurate; the 'simple' method
        uses a polynomial fit (to the log of a gaussian), is the fastest
        and is reliable

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['fit', 'simple', 'arclines']
        v = key_allowed(v, allowed)
        self.update(v)

    def arc_calibrate_nfitpix(self, v):
        """ Number of pixels to fit when deriving the centroid of the
        arc lines (an odd number is best)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v % 2 == 0:
            msgs.warn("An odd integer is recommended for the argument of {0:s}".format(get_current_name()))
        self.update(v)

    def arc_calibrate_numsearch(self, v):
        """ Number of brightest arc lines to search for preliminary identification

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        self.update(v)

    def arc_useframe(self, v):
        """ What filetype should be used for wavelength calibration (arc),
        you can also specify a master calibrations file if it exists.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ["arc"]
        v = key_none_allowed_filename(v, allowed)
        self.update(v)

    def bias_combine_method(self, v):
        """ What method should be used to combine the bias frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def bias_combine_reject_cosmics(self, v):
        """ Specify the rejection threshold (in standard deviations) for
        cosmic rays when combining the bias frames. If v<0, cosmic rays
        will not be rejected.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def bias_combine_reject_lowhigh(self, v):
        """ Specify the number of low/high pixels to be rejected when combining
        the bias frames, in the format: [low,high].

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

    def bias_combine_reject_level(self, v):
        """  Specify the significance threshold (in standard deviations)
        used to reject deviant pixels when combining the bias frames,
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
        """ What should be done to saturated pixels when combining the bias frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    '''
    def bias_useoverscan(self, v):
        """ Subtract the bias level using the overscan region?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)
    '''

    def bias_useframe(self, v):
        """ How to subtract the detector bias (bias, overscan, dark, none),
        you can also specify a master calibrations file if it exists.
        Alternatively, you can select more than one option (e.g. dark and
        overscan or bias and overscan), by providing a list of values
        (e.g. bias useframe [bias,overscan]).

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['bias', 'overscan', 'dark', 'none']
        if "," in v:
            # Must be a list - multiple options
            v = key_list(v)
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
        """ Overwrite any existing output files?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def output_sorted(self, v):
        """ A filename given to output the details of the sorted files.
        If no value is set, the default is the Settings File with .pypit
        removed.
        User can turn off this output by setting: output sorted off
        The value of this keyword argument should have no extension name

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none(v)
        if v is None:
            try:
                v = self._argflag['run']['redname'].replace('.pypit', '')
            except (AttributeError, KeyError):
                pass
        elif v == 'off':
            v = None
        self.update(v)

    def output_verbosity(self, v):
        """ Level of screen output (0 is No screen output, 1 is low level output, 2 is output everything)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if (v < 0) or (v > 2):
            msgs.error("The verbosity can only take values between 0 (minimum) and 2 (maximum)" + msgs.newline() +
                       "Please change the argument of {0:s}".format(get_current_name()))
        self.update(v)


    def pinhole_combine_match(self, v):
        """ Match similar pinhole frames together? A successful match is found when the frames
        are similar to within N-sigma, where N is the argument of this expression. If v<0,
        pinhole frames will not be matched.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def pinhole_combine_method(self, v):
        """ What method should be used to combine the pinhole frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def pinhole_combine_reject_cosmics(self, v):
        """ Specify the rejection threshold (in standard deviations) for
        cosmic rays when combining the pinhole frames. If v<0, cosmic rays
        will not be rejected.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def pinhole_combine_reject_lowhigh(self, v):
        """ Specify the number of low/high pixels to be rejected when combining
        the pinhole frames, in the format: [low,high].

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

    def pinhole_combine_reject_level(self, v):
        """ Specify the significance threshold (in standard deviations)
        used to reject deviant pixels when combining the pinhole frames,
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

    def pinhole_combine_reject_replace(self, v):
        """ What should be done if all pixels are rejected when
        combining the pinhole frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def pinhole_combine_satpix(self, v):
        """ What should be done to saturated pixels when combining the pinhole frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    def pinhole_useframe(self, v):
        """ What filetype should be used to identify the slit edges?
        you can also specify a master calibrations file if it exists.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['pinhole', 'science']
        v = key_none_allowed_filename(v, allowed)
        self.update(v)

    def pixelflat_combine_match(self, v):
        """ Match similar pixel flat frames together? A successful match is found when the frames
        are similar to within N-sigma, where N is the argument of this expression. If v<0,
        pixel flat frames will not be matched.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def pixelflat_combine_method(self, v):
        """ What method should be used to combine the pixel flat frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def pixelflat_combine_reject_cosmics(self, v):
        """ Specify the rejection threshold (in standard deviations) for
        cosmic rays when combining the pixel flat frames. If v<0, cosmic rays
        will not be rejected.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def pixelflat_combine_reject_lowhigh(self, v):
        """ Specify the number of low/high pixels to be rejected when combining
        the pixel flat frames, in the format: [low,high].

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

    def pixelflat_combine_reject_level(self, v):
        """ Specify the significance threshold (in standard deviations)
        used to reject deviant pixels when combining the pixel flat frames,
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

    def pixelflat_combine_reject_replace(self, v):
        """ What should be done if all pixels are rejected when
        combining the pixel flat frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def pixelflat_combine_satpix(self, v):
        """ What should be done to saturated pixels when combining the pixelflat frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)


    def pixelflat_useframe(self, v):
        """ What filetype should be used for pixel-to-pixel calibration (flat),
        you can also specify a master calibrations file if it exists.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['pixelflat']
        v = key_none_allowed_filename(v, allowed)
        self.update(v)

    def reduce_badpix(self, v):
        """ Make a bad pixel mask? (This step requires bias frames)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def reduce_calibrate_nonlinear(self, v):
        """ Perform a non-linear correction? This step requires a series of
        pixel flat frames of the same lamp and setup, and with a variety of
        exposure times and count rates in every pixel.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def reduce_calibrate_refframe(self, v):
        """ Which reference frame do you want the data in?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['geocentric', 'heliocentric', 'barycentric']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_calibrate_wavelength(self, v):
        """ Wavelength calibrate the data? The data will not be wavelength
        calibrated if the value of the keyword is set to None.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['air', 'vacuum', 'none']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_flatfield_method(self, v):
        """ Specify the method that should be used to normalize the flat field

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['polyscan', 'bspline']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_flatfield_params(self, v):
        """ Flat field method parameters, where the parameters relate to the method
        specified by the 'reduce flatfield method' keyword:

        polyscan:  [Polynomial order, Number of pixels, Number of repeats]
        bspline:   [Number of pixels in the dispersion direction between each knot]

        Note: if the bspline argument is 0 < number < 1, it will be assumed to be a fraction of the pixels in the dispersion direction

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def reduce_flatfield_perform(self, v):
        """ Flatfield the data?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def reduce_flatfield_useframe(self, v):
        """ What frame should be used to flat field the data? You can also
        specify a master calibrations file if it exists.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['pixelflat', 'trace']
        v = key_allowed_filename(v, allowed)
        self.update(v)

    def reduce_slitcen_useframe(self, v):
        """ What frame should be used to trace the slit centroid? You can also
        specify a master calibrations file if it exists.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['trace', 'pinhole', 'science']
        v = key_allowed_filename(v, allowed)
        self.update(v)

    def reduce_trace_useframe(self, v):
        """ What frame should be used to trace the slit edges? You can also
        specify a master calibrations file if it exists.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['trace']
        v = key_allowed_filename(v, allowed)
        self.update(v)

    def reduce_masters_file(self, v):
        """

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        if v.lower() == 'none':
            v = ''
        self.update(v)

    def reduce_masters_loaded(self, v):
        """

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def reduce_masters_reuse(self, v):
        """

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def reduce_masters_setup(self, v):
        """

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        if v.lower() == 'none':
            v = ''
        self.update(v)

    def reduce_overscan_method(self, v):
        """ Specify the method that should be used to fit the overscan

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['polynomial', 'savgol', 'median']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_overscan_params(self, v):
        """ Parameters of the overscan subtraction, where the parameters
        relate to the method specified by the 'reduce overscan method' keyword:

        polynomial:  [Polynomial order, Number of pixels, Number of repeats]
        savgol:  [Polynomial order, Window Size]   (note: window size should be odd)
        median:  The median method does not require any parameters.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def reduce_pixel_locations(self, v):
        """ If desired, a fits file can be specified (of the appropriate form) to
        specify the locations of the pixels on the detector (in physical space)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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
        """ The size of the extracted pixels (as an scaled number of Arc FWHM), -1 will not resample

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def reduce_skysub_bspline_everyn(self, v):
        """ bspline fitting parameters

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        self.update(v)

    def reduce_skysub_method(self, v):
        """ Method used for the sky subtraction

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['bspline', 'polyscan']
        v = key_allowed(v, allowed)
        self.update(v)

    def reduce_skysub_perform(self, v):
        """ Subtract the sky background from the data?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def reduce_slitprofile_perform(self, v):
        """ Determine the spatial slit profile?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def reduce_trim(self, v):
        """ Trim the frame to isolate the data?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def run_calcheck(self, v):
        """ If True, PYPIT will not reduce the data, it will just check to
        make sure all calibration data are present

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def run_directory_master(self, v):
        """ Child Directory name for master calibration frames

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def run_directory_qa(self, v):
        """ Child Directory name for quality assurance

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def run_directory_science(self, v):
        """ Child Directory name for extracted science frames

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def run_load_settings(self, v):
        """ Load a reduction settings file (Note: this command overwrites all default settings)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        if v.lower() == "none":
            v = None
        elif not isfile(v):
                msgs.error("The argument of {0:s} must be a PYPIT settings file".format(get_current_name()) +
                           msgs.newline() + "or 'None'. The following file does not exist:" + msgs.newline() + v)
        self.update(v)

    def run_load_spect(self, v):
        """ Load a spectrograph settings file (Note: this command overwrites all default settings)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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
        """ Number of CPUs to use (-1 means all bar one CPU available,
        -2 means all bar two CPUs available)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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
        """ If True, PYPIT will prepare the calibration frames and will
        only reduce the science frames when preponly is set to False

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def run_progname(self, v):
        """ A variable that is set by PYPIT during execution. This parameter
        is not available for user input.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)


    def run_pypitdir(self, v):
        """ A variable that is set by PYPIT during execution. This parameter
        is not available for user input.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def run_qa(self, v):
        """ Run quality control in real time? Setting this keyword to False will
        still produce the checks, but won't display the results during the
        reduction.

        Not currently implemented.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def run_redname(self, v):
        """ A variable that is set by PYPIT during execution. This parameter
        is not available for user input.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def run_setup(self, v):
        """ If True, run in setup mode.  Useful to parse files when starting
        reduction on a large set of data
        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def run_spectrograph(self, v):
        """ The name of the spectrograph data that should be reduced.
        A corresponding settings file must be available.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        # Check that v is allowed
        stgs_arm = glob(dirname(__file__)+"/data/settings/settings.arm*")
        stgs_all = glob(dirname(__file__)+"/data/settings/settings.*")
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
        """ If True, PYPIT will stop and require a user carriage
        return at every quality control check

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def run_useIDname(self, v):
        """ If True, file sorting will ensure that the idname is made

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def science_extraction_manual(self, cnmbr=1, frame="none", params="[1,1000,500,[10,10]]"):
        """ See documentation for the child parameters of this function

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        # Send parameters away to individual arguments
        self.science_extraction_manual_frame(frame, cnmbr=cnmbr)
        self.science_extraction_manual_params(params, cnmbr=cnmbr)

    def science_extraction_manual_frame(self, v, cnmbr=1):
        """ Specify the name of the fits file that a manual extraction will be performed on

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        if v.lower() == "none":
            v = None
        elif ".fits" not in v:
            msgs.error("The argument of {0:s} must be a fits file".format(cname))
        self.update(v, ll=cname.split('_'))

    def science_extraction_manual_params(self, v, cnmbr=1):
        """ Provide the parameters of the manual extraction in the format:
        [1,1000,500,[10,10]], where in this example '1' is the detector number,
        '1000' is the spatial location that the trace must go through, '500' is
        the spectral location that the trace must go through, '[10,10]' is the
        width around the stated (spatial,spectral) location specified above that
        should also be in the trace.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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
        """ Maximum number of objects to extract in a science frame

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        self.update(v)

    def science_extraction_profile(self, v):
        """ Fitting function used to extract science data, only if the extraction
        is 2D. Note, the available options of this argument that have a suffix 'func'
        will fit a function to the pixels whereas the options without this suffix take
        into account the integrated function within each pixel (and is closer to truth).

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['gaussian', 'gaussfunc', 'moffat', 'moffatfunc']
        v = key_allowed(v, allowed)
        self.update(v)

    def science_extraction_reuse(self, v):
        """ If the science frame has previously been extracted and saved, load the extractions

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def trace_combine_match(self, v):
        """ Match similar trace flat frames together? A successful match is found when the frames
        are similar to within N-sigma, where N is the argument of this expression. If v<0,
        trace flat frames will not be matched.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)
        return

    def trace_combine_method(self, v):
        """ What method should be used to combine the trace frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_methods()
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_combine_reject_cosmics(self, v):
        """ Specify the rejection threshold (in standard deviations) for
        cosmic rays when combining the trace frames. If v<0, cosmic rays
        will not be rejected.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def trace_combine_reject_lowhigh(self, v):
        """ Specify the number of low/high pixels to be rejected when combining
        the trace frames, in the format: [low,high].

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

    def trace_combine_reject_level(self, v):
        """ Specify the significance threshold (in standard deviations)
        used to reject deviant pixels when combining the trace frames,
        in the format: [low,high].

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] <= 0.0 or v[1] <= 0.0:
            msgs.error("The list values of argument {0:s} must be > 0.0".format(get_current_name()))
        self.update(v)

    def trace_combine_reject_replace(self, v):
        """ What should be done if all pixels are rejected when
        combining the trace frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_replaces()
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_combine_satpix(self, v):
        """ What should be done to saturated pixels when combining the trace frames?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = combine_satpixs()
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_dispersion_direction(self, v):
        """ Specify the primary dispersion direction of the raw data (0 for row, 1 for column)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v != 0 and v != 1:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       "0 or 1 (if the dispersion axis is along a row or column respectively)")
        self.update(v)

    def trace_slits_diffpolyorder(self, v):
        """ What is the order of the 2D function that should be used to fit
        the 2D solution for the spatial size of all slits?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_expand(self, v):
        """ If you are tracing the slit edges with a pinhole frame (i.e. a pinhole/science frame),
        you should expand the slit edges to the edges defined by the trace frame, which
        should be a flatfield exposure taken with the same slit length as the science frame.
        If the slits are traced with a trace frame, there is no need to expand the slits.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def trace_slits_fracignore(self, v):
        """ If a slit spans less than this fraction over the spectral size of the detector,
        it will be ignored (and reconstructed when/if an 'order' PCA analysis is performed).

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        if v < 0.0 or v > 1.0:
            msgs.error("The argument of {0:s} must be between 0 and 1".format(get_current_name()))
        self.update(v)

    def trace_slits_function(self, v):
        """ What function should be used to trace the slits?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['polynomial', 'legendre', 'chebyshev']
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_slits_maxgap(self, v):
        """ Maximum gap between slits. Use 'None' if the neighbouring
        slits are far apart, or of similar illumination.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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
        """ Manually set the number of slits to identify (>=1).
        'auto' or -1 will automatically identify the number of slits.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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
        """ What is the order of the function (specified by 'trace slits function')
        that should be used to trace the slits ?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_pad(self, v):
        """ How many pixels should be considered beyond the automatic slit
        edge trace. Note that this parameter does not change the location
        of the slit edges. This parameter allows for a smooth model to be
        fit past the automatically detected slit edges.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_pca_type(self, v):
        """ Should the PCA be performed using pixel position (pixel) or by spectral order (order).
        The latter is used for echelle spectroscopy, or for slits where the slit separation is a
        smooth function of the slit number. The former option can be used for multi-object spectroscopy
        where the gap between slits is irregular.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['pixel', 'order']
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_slits_pca_params(self, v):
        """ What order polynomials should be used to fit the principle components

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def trace_slits_pca_extrapolate_pos(self, v):
        """ How many extra echelle orders to predict in the positive direction

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_pca_extrapolate_neg(self, v):
        """ How many extra orders to predict in the negative direction

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_sigdetect(self, v):
        """ Sigma detection threshold for edge detection

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0".format(get_current_name()))
        self.update(v)

    def trace_slits_single(self, v):
        """ Add a user-defined slit?
        Syntax is a list of values, 2 per detector that define the slit
        according to column values.  The 2nd value (for the right edge)
        must be >0 to be applied.  Example for LRISr [-1, -1, 7, 295]
        which means the code skips user-definition for the first detector
        but adds one for the 2nd.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def trace_slits_tilts_idsonly(self, v):
        """ Use only the arc lines that have an identified wavelength
        to trace the spectral tilt

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def trace_slits_tilts_method(self, v):
        """ What method should be used to trace the spectral tilt of the
        slit along an order?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['pca', 'spline', 'spca', 'interp', 'perp', 'zero']
        v = key_allowed(v, allowed)
        self.update(v)

    def trace_slits_tilts_params(self, v):
        """ Parameters that should be used for the 'trace slits tilts method' arguement.
        Options include:

        pca, spca :  A list containing the order of the polynomials that should be used to fit the tilt principle components

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        self.update(v)

    def trace_slits_tilts_order(self, v):
        """ What is the order of the polynomial function to be used for the tilt of an individual arc line

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_useframe(self, v):
        """ What frame should be used to trace the slit edges, based on the
        average of the left/right edges.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['trace']
        v = key_none_allowed_filename(v, allowed)
        self.update(v)


class BaseSpect(BaseFunctions):
    def __init__(self, defname, savname):
        super(BaseSpect, self).__init__(defname, savname)
        self._spect = NestedDict()
        self._settings = []
        self.set_default()

    def load_ftype(self, ftype_dict):
        """ Parse the dict generated from a .pypit file on frametypes
        Parameters
        ----------
        ftype_dict : dict

        Returns
        -------
        linesarr : list
          Each element of this list is a line equivalent to that in a PYPIT file
          e.g.  arc number 1

        """
        # Save
        self.__dict__['_ftdict'] = ftype_dict.copy()
        # Dict to hold values
        fdict = {}
        for key,value in ftype_dict.items():
            ftypes = value.split(',')
            for ftype in ftypes:
                if ftype == 'science':
                    continue
                if ftype not in fdict.keys():
                    fdict[ftype] = 1
                else:
                    fdict[ftype] += 1
        # Generate the lines
        linesarr = []
        for key,value in fdict.items():
            if value > self.__dict__['_spect'][key]['number']:
                linesarr.append(' {:s} number {:d}\n'.format(key,value))
        # Return
        return linesarr

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
        self.update([], ll="set_pinhole".split("_"))
        self.update([], ll="set_pixelflat".split("_"))
        self.update([], ll="set_science".split("_"))
        self.update([], ll="set_standard".split("_"))
        self.update([], ll="set_trace".split("_"))
        self.update([], ll="arc_index".split("_"))
        self.update([], ll="bias_index".split("_"))
        self.update([], ll="pinhole_index".split("_"))
        self.update([], ll="pixelflat_index".split("_"))
        self.update([], ll="science_index".split("_"))
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
        frmtyp = ["standard", "bias", "pixelflat", "trace", "pinhole", "arc"]
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
        """ Adjust a PYPIT setting. This can be used to force certain default
        reduction options for a given spectrograph. For example, to use a set
        number of cpus when reducing a given instrument, you could specify in
        the settings.instrument_name file: 'settings run ncpus 3'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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
        """ If there are frames that will be an arc in addition to other frame types,
        include the other frame types here.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none_list(v)
        self.update(v)

    def arc_check_condition(self, v, cnmbr=1):
        """ Check that a frame satisfies a series of conditions before it is
        labelled as an arc frame. Multiple conditions can be specified,
        where each new condition has a different integer suffix appended to
        the condition variable.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        v = key_check(v)
        self.update(v, ll=cname.split('_'))

    def arc_idname(self, v):
        """ Header key value of arc frames for header keyword: 'keyword idname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def arc_number(self, v):
        """ Number of arc frames to use

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def bias_canbe(self, v):
        """ If there are frames that will be a bias in addition to other frame types,
        include the other frame types here.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none_list(v)
        self.update(v)

    def bias_check_condition(self, v, cnmbr=1):
        """ Check that a frame satisfies a series of conditions before it is
        labelled as a bias frame. Multiple conditions can be specified,
        where each new condition has a different integer suffix appended to
        the condition variable.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        v = key_check(v)
        self.update(v, ll=cname.split('_'))

    def bias_idname(self, v):
        """ Header key value of bias frames for header keyword: 'keyword idname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def bias_number(self, v):
        """ Number of bias frames to use

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        #key_min_val(v,-1)
        self.update(v)

    def det_datasec(self, v, anmbr=1, bnmbr=1):
        """ Either the data sections or the header keyword where the
        valid data sections can be obtained.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr, bnmbr=bnmbr)
        try:
            v = load_sections(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be a detector section".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_oscansec(self, v, anmbr=1, bnmbr=1):
        """ Either the overscan sections or the header keyword where the
        valid overscan sections can be obtained.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr, bnmbr=bnmbr)
        try:
            v = load_sections(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be detector section".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_darkcurr(self, v, anmbr=1):
        """ Dark current (e-/hour)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        v = key_float(v)
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_gain(self, v, anmbr=1):
        """ Inverse gain (e-/ADU). A list should be provided if a detector contains more than one amplifier.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        try:
            v = v.split(",")
            for i in range(len(v)):
                v[i] = float(v[i])
                if v[i] <= 0.0:
                    msgs.error("Each argument of {0:s} must be > 0.0".format(cname))
        except ValueError:
            msgs.error("Each argument of {0:s} must be of type float".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_ronoise(self, v, anmbr=1):
        """ Read-out noise (e-). A list should be provided if a detector contains more than one amplifier.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        try:
            v = v.split(",")
            for i in range(len(v)):
                v[i] = float(v[i])
                if v[i] <= 0.0:
                    msgs.error("Each argument of {0:s} must be > 0.0".format(cname))
        except ValueError:
            msgs.error("Each argument of {0:s} must be of type float".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_nonlinear(self, v, anmbr=1):
        """ Percentage of detector range which is linear (i.e. everything above
        nonlinear*saturation will be flagged as saturated)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        v = key_float(v)
        if v <= 0.0 or v > 1.0:
            msgs.error("The argument of {0:s} must be > 0.0 and <= 1.0".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_numamplifiers(self, v, anmbr=1):
        """ Number of amplifiers for each detector.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        v = key_int(v)
        if v <= 0:
            msgs.error("The argument of {0:s} must be >= 1".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_saturation(self, v, anmbr=1):
        """ The detector saturation level

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        v = key_float(v)
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0.0".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_suffix(self, v, anmbr=1):
        """ Suffix to be appended to all saved calibration and extraction frames

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        self.update(v, ll=cname.split('_'))

    def det_xgap(self, v, anmbr=1):
        """ Gap between the square detector pixels (expressed as a fraction
        of the pixel size along the dispersion axis)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        v = key_float(v)
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_ygap(self, v, anmbr=1):
        """ Gap between the square detector pixels (expressed as a fraction
        of the pixel size along the spatial axis)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        v = key_float(v)
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        self.update(v, ll=cname.split('_'))

    def det_ysize(self, v, anmbr=1):
        """ The size of a pixel in the spatial direction as a multiple of the
        pixel size along the spectral direction (i.e. assume xsize = 1.0)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(anmbr=anmbr)
        v = key_float(v)
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0.0".format(cname))
        self.update(v, ll=cname.split('_'))

    def fits_calwin(self, v):
        """ The window of time in hours to search for matching calibration
        frames for a science frame.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        #if v <= 0.0:
        #    msgs.error("The calibration time window must be > 0.0")
        self.update(v)

    def fits_dataext(self, v):
        """ Extension number of data

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The fits data extension number must be >= 0")
        self.update(v)

    def fits_headext(self, v, bnmbr=1):
        """ How many headers need to be read in for a given file

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(bnmbr=bnmbr)
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(cname))
        self.update(v, ll=cname.split('_'))

    def fits_numhead(self, v):
        """ Extension number of header (one for each headnum, starting with 01)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v <= 0:
            msgs.error("The number of fits headers must be >= 1")
        self.update(v)

    def fits_numlamps(self, v):
        """ How many lamps are listed in the header

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The number of lamps must be >= 0")
        self.update(v)

    def fits_timeunit(self, v):
        """ The unit of keyword time

        (s=seconds, m=minutes, h=hours, or any of the astropy Time formats)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['s', 'm', 'h']
        astropy_allowed = Time.FORMATS.keys()
        if v not in allowed and v not in astropy_allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed) + "or one of the astropy Time formats:" + msgs.newline() +
                       ", ".join(astropy_allowed))
        self.update(v)

    def keyword_ra(self, v):
        """ Right Ascension of the telescope pointing

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_dec(self, v):
        """ Declination of the telescope pointing

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_target(self, v):
        """ Header keyword for the name given by the observer to a given frame

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_airmass(self, v):
        """ Airmass at start of observation

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_binning(self, v):
        """ The binning of the data

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_binningspatial(self, v):
        """ Spatial binning

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_binningspectral(self, v):
        """ Spectral binning

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_date(self, v):
        """ The date of the observation (in the format YYYY-MM-DD  or  YYYY-MM-DDTHH:MM:SS.SS)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_decker(self, v):
        """ Which decker is being used

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_detrot(self, v):
        """ Detector Rotation angle

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_dichroic(self, v):
        """ Dichroic used for the observation

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_dispname(self, v):
        """ Disperser name

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_dispangle(self, v):
        """ Disperser angle

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_equinox(self, v):
        """ The equinox to use

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_exptime(self, v):
        """ Exposure time

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_filter1(self, v):
        """ Filter 1

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_filter2(self, v):
        """ Filter 2

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_hatch(self, v):
        """ Hatch open/close

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_idname(self, v):
        """ The keyword that identifies the frame type (i.e. bias, flat, etc.)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_lamps(self, v):
        """ Lamps being used

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_lampname(self, v, bnmbr=1):
        """ Name of a lamp. Multiple lamp nams can be specified by appending a
        two digit number (starting with 01) after lampname. There must be a
        corresponding keyword set for 'keyword lampstat'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(bnmbr=bnmbr)
        if "." not in v:
            # User must have passed the name of the lamp
            pass
        else:
            # Get the lamp status from the header
            v = key_keyword(v)
        self.update(v, ll=cname.split('_'))

    def keyword_lampstat(self, v, bnmbr=1):
        """ Status of a lamp. Multiple lamp statuses  can be specified by appending a
        two digit number (starting with 01) after lampstat. There must be a corresponding
        keyword set for 'keyword lampname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(bnmbr=bnmbr)
        if "." not in v:
            # User must have specified the status
            pass
        else:
            # Get the lamp status from the header
            v = key_keyword(v)
        self.update(v, ll=cname.split('_'))

    def keyword_naxis0(self, v):
        """ Number of pixels along the zeroth axis

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_naxis1(self, v):
        """ Number of pixels along the first axis

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_slitwid(self, v):
        """ Slit Width

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_slitlen(self, v):
        """ Slit Length

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def keyword_time(self, v):
        """ The time stamp of the observation (i.e. decimal MJD)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
        self.update(v)

    def mosaic_camera(self, v):
        """ Set the name of the instrument used (this will be used in the QA).

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def mosaic_elevation(self, v):
        """ Elevation of the telescope (in m)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def mosaic_latitude(self, v):
        """ Latitude of the telescope

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def mosaic_longitude(self, v):
        """ Longitude of the telescope

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def mosaic_ndet(self, v):
        """ Number of detectors in the mosaic

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        self.update(v)

    def mosaic_minexp(self, v):
        """ Minimum exposure time of the instrument (s)

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_float(v)
        self.update(v)

    def mosaic_reduction(self, v):
        """ Which reduction pipeline should be used to reduce data taken with this instrument

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['ARMLSD', 'ARMED']
        v = key_allowed(v, allowed, upper=True)
        self.update(v.upper())

    def pixelflat_canbe(self, v):
        """ If there are frames that will be a pixel flat in addition to other frame types,
        include the other frame types here.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none_list(v)
        self.update(v)

    def pixelflat_check_condition(self, v, cnmbr=1):
        """ Check that a frame satisfies a series of conditions before it is
        labelled as a pixel flat frame. Multiple conditions can be specified,
        where each new condition has a different integer suffix appended to
        the condition variable.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        v = key_check(v)
        self.update(v, ll=cname.split('_'))

    def pixelflat_idname(self, v):
        """ Header key value of pixel flat frames for header keyword: 'keyword idname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def pixelflat_lscomb(self, v):
        """ Combine frames with a different exposure time?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def pixelflat_number(self, v):
        """ Number of pixel flat frames to use

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        #assert key_min_val(v,-1)
        #if v < -1:
        #    msgs.error("The argument of {0:s} must be >= -1".format(get_current_name()))
        self.update(v)

    def science_canbe(self, v):
        """ If there are frames that will be a science frame in addition to other frame types,
        include the other frame types here.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none_list(v)
        self.update(v)

    def science_check_condition(self, v, cnmbr=1):
        """ Check that a frame satisfies a series of conditions before it is
        labelled as a science frame. Multiple conditions can be specified,
        where each new condition has a different integer suffix appended to
        the condition variable.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        v = key_check(v)
        self.update(v, ll=cname.split('_'))

    def science_idname(self, v):
        """ Header key value of science frames for header keyword: 'keyword idname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def set_arc(self, v):
        """ Manually force a given frame to be an arc frame. For example,
         'set arc filename1.fits,filename2.fits' will force filename1.fits
         and filename2.fits to be arc frames.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        v = self._spect['set']['arc'] + v
        self.update(v)

    def set_bias(self, v):
        """ Manually force a given frame to be a bias frame. For example,
         'set bias filename1.fits,filename2.fits' will force filename1.fits
         and filename2.fits to be bias frames.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        v = self._spect['set']['bias'] + v
        self.update(v)

    def set_pixelflat(self, v):
        """ Manually force a given frame to be a pixel flat frame. For example,
         'set pixelflat filename1.fits,filename2.fits' will force filename1.fits
         and filename2.fits to be pixel flat frames.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        v = self._spect['set']['pixelflat'] + v
        self.update(v)

    def set_science(self, v):
        """ Manually force a given frame to be a science frame. For example,
         'set science filename1.fits,filename2.fits' will force filename1.fits
         and filename2.fits to be science frames.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        v = self._spect['set']['science'] + v
        self.update(v)

    def set_pinhole(self, v):
        """ Manually force a given frame to be a pinhole frame. For example,
         'set pinhole filename1.fits,filename2.fits' will force filename1.fits
         and filename2.fits to be pinhole frames.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        v = self._spect['set']['pinhole'] + v
        self.update(v)

    def set_standard(self, v):
        """ Manually force a given frame to be a standard star frame. For example,
         'set standard filename1.fits,filename2.fits' will force filename1.fits
         and filename2.fits to be standard star frames.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        v = self._spect['set']['standard'] + v
        self.update(v)

    def set_trace(self, v):
        """ Manually force a given frame to be a trace frame. For example,
         'set trace filename1.fits,filename2.fits' will force filename1.fits
         and filename2.fits to be trace frames.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_list(v)
        v = self._spect['set']['trace'] + v
        self.update(v)

    def pinhole_canbe(self, v):
        """ If there are frames that will be a pinhole in addition to other frame types,
        include the other frame types here.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none_list(v)
        self.update(v)

    def pinhole_check_condition(self, v, cnmbr=1):
        """ Check that a frame satisfies a series of conditions before it is
        labelled as a pinhole frame. Multiple conditions can be specified,
        where each new condition has a different integer suffix appended to
        the condition variable.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        v = key_check(v)
        self.update(v, ll=cname.split('_'))

    def pinhole_idname(self, v):
        """ Header key value of pinhole frames for header keyword: 'keyword idname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def pinhole_lscomb(self, v):
        """ Combine frames with a different exposure time?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def pinhole_number(self, v):
        """ Number of pinhole frames to use

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        #key_min_val(v, -1)
        self.update(v)

    def standard_canbe(self, v):
        """ If there are frames that will be a standard star frame in addition
        to other frame types, include the other frame types here.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none_list(v)
        self.update(v)

    def standard_check_condition(self, v, cnmbr=1):
        """ Check that a frame satisfies a series of conditions before it is
        labelled as a standard frame. Multiple conditions can be specified,
        where each new condition has a different integer suffix appended to
        the condition variable.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        v = key_check(v)
        self.update(v, ll=cname.split('_'))

    def standard_idname(self, v):
        """ Header key value of standard star frames for header keyword: 'keyword idname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def standard_number(self, v):
        """ Number of standard star frames to use

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        #key_min_val(v,-1)
        self.update(v)

    def trace_canbe(self, v):
        """ If there are frames that will be a trace flat in addition to other frame types,
        include the other frame types here.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_none_list(v)
        self.update(v)

    def trace_check_condition(self, v, cnmbr=1):
        """ Check that a frame satisfies a series of conditions before it is
        labelled as a trace frame. Multiple conditions can be specified,
        where each new condition has a different integer suffix appended to
        the condition variable.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        cname = get_nmbr_name(cnmbr=cnmbr)
        v = key_check(v)
        self.update(v, ll=cname.split('_'))

    def trace_idname(self, v):
        """ Header key value of a trace frame for header keyword: 'keyword idname'

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        self.update(v)

    def trace_lscomb(self, v):
        """ Combine frames with a different exposure time?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)

    def trace_number(self, v):
        """ Number of trace frames to use

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        #if v < 0:
        #    msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)


class ARMLSD(BaseArgFlag):

    def reduce_calibrate_flux(self, v):
        """ Should a flux calibration be performed?

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_bool(v)
        self.update(v)


    def reduce_flexure_maxshift(self, v):
        """ Maximum allowed flexure shift in pixels

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        self.update(v)

    def reduce_flexure_method(self, v):
        """ Perform flexure correction on objects using boxcar extraction.
        If 'slitcen' is used, the flexure correction is performed before
        the extraction of objects

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        allowed = ['none', 'boxcar', 'slitcen']
        v = key_none_allowed(v, allowed)
        self.update(v)

    def reduce_flexure_spectrum(self, v):
        """ Specify the archive sky spectrum to be used for the flexure correction

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
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


class ARMED(BaseArgFlag):

    def reduce_flatfield_2dpca(self, v):
        """ Perform a simple 2D PCA on the echelle blaze fits
         if the value of this argument is >1. The argument value
         is equal to the number of PCA components. 0 means that
         no PCA will be performed.

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_tilts_disporder(self, v):
        """ What is the order of the polynomial function to be used to fit the tilts along the dispersion direction

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        self.update(v)

    def trace_slits_tilts_order(self, v):
        """ What is the order of the polynomial function to be used for the tilt of an individual arc line

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_int(v)
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        if v != 1:
            msgs.error("The argument of {0:s} must be equal to 1 for echelle data".format(get_current_name()))
        self.update(v)


class ARMLSD_spect(BaseSpect):
    pass


class ARMED_spect(BaseSpect):

    def keyword_echangle(self, v):
        """ The angle of the echelle grating

        Parameters
        ----------
        v : str
          value of the keyword argument given by the name of this function
        """
        v = key_keyword(v)
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
        defname = glob(dirname(__file__))[0] + "/data/settings/settings." + init[0].lower()
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
        defname = glob(dirname(__file__))[0] + "/data/settings/settings." + init[1].lower()
        return eval(init[0]+"_spect(defname='{0:s}', savname='{1:s}.spect')".format(defname, init[2]))
    except RuntimeError:
        msgs.error("{0:s} is not implemented yet".format(init[1]))


def get_current_name():
    """ Return the name of the function that called this function
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    return "'" + " ".join(ll) + "'"


def get_nmbr_name(anmbr=None, bnmbr=None, cnmbr=None):
    """ Return the name of the function that called this function,
    and append a two digit number to the first (when anmbr!=None),
    the second (when bnmbr!=None), or third (when cnmbr!=None) element
    of the function name.
    """
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
    """ Convert a detector index into a string used by the settings dictionary

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


def key_allowed(v, allowed, upper=False):
    """ Check that a keyword argument is in an allowed list of parameters.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take
    upper : bool (optional)
      If True, the allowed list is expected to contain only
      uppercase strings. If False, the allowed list is expected
      to contain only lowercase strings.

    Returns
    -------
    v : str
      A string used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if upper:
        v = v.upper()
    else:
        v = v.lower()
    if v not in allowed:
        msgs.error("The argument of {0:s} must be one of".format(func_name) + msgs.newline() +
                   ", ".join(allowed))
    if v.lower() == "none":
        v = None
    return v


def key_allowed_filename(v, allowed):
    """ Check that a keyword argument is in an allowed list of parameters.
    If not, assume that it is a filename.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : str
      A string used by the settings dictionary
    """
    vt = v.lower()
    if vt not in allowed:
        msgs.warn("Assuming the following is the name of a file:" + msgs.newline() + v)
    else:
        v = vt
    return v


def key_bool(v):
    """ Check that a keyword argument is a boolean variable.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : str
      A string used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if v.lower() == "true":
        v = True
    elif v.lower() == "false":
        v = False
    else:
        msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(func_name))
    return v


def key_check(v):
    """ Check that a keyword argument satisfies the form required of a
    keyword argument that checks frame types.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : str, list
      A value used by the settings dictionary
    """
    text = v.strip().replace('_', ' ')
    if ',' in text and text[0:2] != '%,':
        # There are multiple possibilities - split the text
        v = text.split(',')
    else:
        v = text
    return v


def key_float(v):
    """ Check that a keyword argument is a float.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : float
      A value used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    try:
        v = float(v)
    except ValueError:
        msgs.error("The argument of {0:s} must be of type float".format(func_name))
    return v


def key_int(v):
    """ Check that a keyword argument is an int.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : int
      A value used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    try:
        v = int(v)
    except ValueError:
        msgs.error("The argument of {0:s} must be of type int".format(func_name))
    return v


def key_keyword(v):
    """ Check that a keyword argument satisfies the form required
    for specifying a header keyword.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : str
      A value used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if v.lower() == "none":
        v = None
    else:
        try:
            vspl = v.split(".")
            int(vspl[0])
        except ValueError:
            msgs.error("The argument of {0:s} must be of the form:".format(func_name) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
    return v


def key_list(strlist):
    """ Check that a keyword argument is a list. Set the
    appropriate type of the list based on the supplied values.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : list
      A value used by the settings dictionary
    """
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
    """ Check that a keyword argument is a list. Set the
    appropriate type of the list based on the supplied values.
    Then, check that each value in the list is also in the
    supplied 'allowed' list.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : list
      A value used by the settings dictionary
    """
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
    """ Check if a keyword argument is set to None. If not,
    check that its value is in the 'allowed' list.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : None, str
      A value used by the settings dictionary
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
    """ Check if a keyword argument is set to None. If not,
    check that its value is in the 'allowed' list. Finally,
    assume that the supplied value is the name of a filename.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : None, str
      A value used by the settings dictionary
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


def key_none_list(v):
    """ Check if a keyword argument is set to None. If not,
    assume the supplied value is a list.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : list
      A value used by the settings dictionary
    """
    if v.lower() == "none":
        v = None
    else:
        if "," in v:
            v = v.strip("()[]").split(",")
        else:
            v = [v]
    return v


def key_none(v):
    """ Check if a keyword argument is set to None.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : None, str
      A value used by the settings dictionary
    """
    if v.lower() == "none":
        v = None
    return v


def key_min_val(v, vmin):
    """ Check that the value exceeds a minimum
    Returns
    -------
    bool

    """
    if v < vmin:
        msgs.error("The argument of {0:s} must be >= -1".format(get_current_name()))
    else:
        return True

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
    """ The options that can be used to replace saturated pixels when combining a set of frames
    """
    methods = ['reject', 'force', 'nothing']
    return methods
