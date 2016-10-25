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
        self._defname = defname
        self._afout = open(savname, 'w')

    def load_file(self, filename=None):
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
        arr = self.load_lines(lines)
        return arr

    def load_lines(self, lines):
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
        super(BaseArgFlag, self).__init__(defname, savname)
        self._argflag = NestedDict()

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
        for ll in range(len(lstall)):
            lst = lstall[ll]
            cnt = 1
            succeed = False
            members = [x for x, y in inspect.getmembers(self, predicate=inspect.ismethod)]
            while cnt < len(lst):
                func = "_".join(lst[:-cnt])
                if func in members:
                    func = "self." + func + "('{0:s}')".format(" ".join(lst[-cnt:]))
                    eval(func)
                    succeed = True
                    break
                else:
                    cnt += 1
            if not succeed:
                msgs.error("There appears to be an error on the following input line:" + msgs.newline() +
                           " ".join(lst))
        return

    def update(self, v, ll=None):
        """
        Update an element in the nested dictionary
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
        return

    def arc_combine_match(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_combine_method(self, v):
        # Check that v is allowed
        allowed = ['mean', 'median', 'weightmean']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def arc_combine_reject_cosmics(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_combine_reject_lowhigh(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def arc_combine_reject_level(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def arc_combine_reject_replace(self, v):
        # Check that v is allowed
        allowed = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def arc_combine_satpix(self, v):
        # Check that v is allowed
        allowed = ['reject', 'force', 'nothing']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def arc_extract_binby(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        if v < 1.0:
            msgs.error("The argument of {0:s} must be >= 1.0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_load_extracted(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_load_calibrated(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_calibrate_detection(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_calibrate_IDpixels(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def arc_calibrate_IDwaves(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def arc_calibrate_lamps(self, v):
        # Check that v is allowed
        allowed = ['ArI', 'CdI', 'HgI', 'HeI', 'KrI', 'NeI', 'XeI', 'ZnI', 'ThAr']
        v = load_list(v)
        for ll in v:
            if ll not in allowed:
                msgs.error("The current list of arc lines does not include: {0:s}".format(ll) + msgs.newline() +
                           "Please choose one of the following:" + msgs.newline() +
                           ", ".join(allowed) + msgs.newline() +
                           "for the argument of {0:s}".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_calibrate_method(self, v):
        # Check that v is allowed
        allowed = ['fit', 'simple']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def arc_calibrate_nfitpix(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v % 2 == 0:
            msgs.warn("An odd integer is recommended for the argument of {0:s}".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_calibrate_numsearch(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_useframe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            msgs.info("Assuming the following is the name of an arc frame:" + msgs.newline() + v)
        # Update argument
        self.update(v)
        return

    def bias_combine_method(self, v):
        # Check that v is allowed
        allowed = ['mean', 'median', 'weightmean']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def bias_combine_reject_cosmics(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def bias_combine_reject_lowhigh(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def bias_combine_reject_level(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def bias_combine_reject_replace(self, v):
        # Check that v is allowed
        allowed = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def bias_combine_satpix(self, v):
        # Check that v is allowed
        allowed = ['reject', 'force', 'nothing']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def bias_useframe(self, v):
        # Check that v is allowed
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
        elif v.lower() == "none":
            v = None
        elif v.lower() not in allowed:
            msgs.info("Assuming the following is the name of a bias frame:" + msgs.newline() +
                      v)
        # Update argument
        self.update(v)
        return

    def output_overwrite(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def output_sorted(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            if v.find('.') > 0:
                msgs.error("The argument for keyword 'output sorted' should have no extension")
        # Update argument
        self.update(v)
        return

    def output_verbosity(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if (v < 0) or (v > 2):
            msgs.error("The verbosity can only take values between 0 (minimum) and 2 (maximum)" + msgs.newline() +
                       "Please change the argument of {0:s}".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def pixelflat_combine_match(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def pixelflat_combine_method(self, v):
        # Check that v is allowed
        allowed = ['mean', 'median', 'weightmean']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def pixelflat_combine_reject_cosmics(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def pixelflat_combine_reject_lowhigh(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def pixelflat_combine_reject_level(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def pixelflat_combine_reject_replace(self, v):
        # Check that v is allowed
        allowed = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def pixelflat_combine_satpix(self, v):
        # Check that v is allowed
        allowed = ['reject', 'force', 'nothing']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def pixelflat_norm_recnorm(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def pixelflat_useframe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        elif v.lower() in ['pixelflat']:
            v = v.lower()
        else:
            msgs.info("Assuming the following is the name of a pixelflat frame:" + msgs.newline() + v)
        # Update argument
        self.update(v)
        return

    def reduce_badpix(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_calibrate_nonlinear(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_calibrate_refframe(self, v):
        # Check that v is allowed
        allowed = ['geocentric', 'heliocentric', 'barycentric']
        if v.lower() not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def reduce_calibrate_wavelength(self, v):
        # Check that v is allowed
        allowed = ['air', 'vacuum', 'none']
        if v.lower() not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def reduce_flatfield_method(self, v):
        # Check that v is allowed
        allowed = ['polyscan']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def reduce_flatfield_params(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def reduce_flatfield_perform(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_flatfield_useframe(self, v):
        # Check that v is allowed
        allowed = ['pixelflat', 'slitflat']
        vt = v.lower()
        if vt not in allowed:
            msgs.warn("Assuming the following is the name of a master flatfield frame:" + msgs.newline() + v)
        else:
            v = vt
        # Update argument
        self.update(v)
        return

    def reduce_trace_useframe(self, v):
        # Check that v is allowed
        allowed = ['trace', 'slitflat', 'science']
        vt = v.lower()
        if vt not in allowed:
            msgs.warn("Assuming the following is the name of a master trace frame:" + msgs.newline() + v)
        else:
            v = vt
        # Update argument
        self.update(v)
        return

    def reduce_flexure_maxshift(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_flexure_method(self, v):
        # Check that v is allowed
        allowed = ['none', 'boxcar', 'slitcen']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        if v == 'none':
            v = None
        # Update argument
        self.update(v)
        return

    def reduce_flexure_spectrum(self, v):
        # Check that v is allowed
        if not pathexists(self._argflag['run']['pypitdir'] + 'data/sky_spec/' + v):
            files_sky = glob(self._argflag['run']['pypitdir'] + 'data/sky_spec/*.fits')
            skyfiles = ""
            for i in files_sky:
                skyfiles += msgs.newline() + "  - " + str(i.split("/")[-1])
            msgs.error("The following archive sky spectrum file does not exist:" + msgs.newline() +
                       "  " + v + msgs.newline() + msgs.newline() + "Please use one of the files listed below:" +
                       skyfiles)
        # Update argument
        self.update(v)
        return

    def reduce_masters_file(self, v):
        # Check that v is allowed
        if v.lower() == 'none':
            v = ''
        # Update argument
        self.update(v)
        return

    def reduce_masters_loaded(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def reduce_masters_reuse(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_masters_setup(self, v):
        # Check that v is allowed
        if v.lower() == 'none':
            v = ''
        # Update argument
        self.update(v)
        return

    def reduce_overscan_method(self, v):
        # Check that v is allowed
        allowed = ['polynomial', 'savgol']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def reduce_overscan_params(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def reduce_pixel_locations(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        elif v.split(".")[-1] == "fits":
            pass
        elif v.split(".")[-2] == "fits" and v.split(".")[-1] == "gz":
            pass
        else:
            msgs.error("The argument of {0:s} must be 'None' or a fits file".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_pixel_size(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_skysub_bspline_everyn(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_skysub_method(self, v):
        # Check that v is allowed
        allowed = ['bspline', 'polyscan']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def reduce_skysub_perform(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_trim(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def run_calcheck(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'").format(get_current_name())
        # Update argument
        self.update(v)
        return

    def run_directory_master(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)
        return

    def run_directory_qa(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)
        return

    def run_directory_science(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)
        return

    def run_load_settings(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        elif not isfile(v):
                msgs.error("The argument of {0:s} must be a PYPIT settings file".format(get_current_name()) +
                           msgs.newline() + "or 'None'. The following file does not exist:" + msgs.newline() + v)
        # Update argument
        self.update(v)
        return

    def run_load_spect(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        elif isfile(v):
            pass
        elif isfile(v+".spect"):
            v += ".spect"
        else:
             msgs.error("The argument of {0:s} must be a PYPIT spectrograph settings".format(get_current_name()) +
                        msgs.newline() + "file or 'None'. The following file does not exist:" + msgs.newline() + v)
        # Update argument
        self.update(v)
        return

    def run_ncpus(self, v):
        # Check that v is allowed
        if 'ncpus' in self._argflag['run'].keys():
            curcpu = self._argflag['run']['ncpus']
        else:
            curcpu = 0
        cpucnt = cpu_count()
        if v == 'all':
            v = cpucnt  # Use all available cpus
            if v != curcpu:
                msgs.info("Setting {0:d} CPUs".format(v))
        elif v is None:
            v = cpucnt-1  # Use all but 1 available cpus
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
        # Update argument
        self.update(v)
        return

    def run_preponly(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def run_progname(self, v):
        # Update argument
        self.update(v)

    def run_pypitdir(self, v):
        # Update argument
        self.update(v)

    def run_qa(self, v):
        """
        run qcontrol
        """
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def run_redname(self, v):
        # Update argument
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
        # Update argument
        self.update(v)
        return

    def run_stopcheck(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def run_useIDname(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def science_extraction_maxnumber(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def science_extraction_profile(self, v):
        # Check that v is allowed
        allowed = ['gaussian', 'gaussfunc', 'moffat', 'moffatfunc']
        if v.lower() not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def science_extraction_reuse(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def slitflat_combine_match(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def slitflat_combine_method(self, v):
        # Check that v is allowed
        allowed = ['mean', 'median', 'weightmean']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def slitflat_combine_reject_cosmics(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def slitflat_combine_reject_lowhigh(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def slitflat_combine_reject_level(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def slitflat_combine_reject_replace(self, v):
        # Check that v is allowed
        allowed = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def slitflat_combine_satpix(self, v):
        # Check that v is allowed
        allowed = ['reject', 'force', 'nothing']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def slitflat_norm_recnorm(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def slitflat_useframe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        elif v.lower() in ['slitflat']:
            v = v.lower()
        else:
            msgs.info("Assuming the following is the name of a pixelflat frame:" + msgs.newline() + v)
        # Update argument
        self.update(v)
        return

    def trace_combine_match(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_combine_method(self, v):
        # Check that v is allowed
        allowed = ['mean', 'median', 'weightmean']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_combine_reject_cosmics(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_combine_reject_lowhigh(self, v):
        # Check that v is allowed
        v = load_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] < 0 or v[1] < 0:
            msgs.error("The list values of argument {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_combine_reject_level(self, v):
        # Check that v is allowed
        v = load_list(v)
        if len(v) != 2:
            msgs.error("The argument of {0:s} must be a two element list".format(get_current_name()))
        if v[0] <= 0.0 or v[1] <= 0.0:
            msgs.error("The list values of argument {0:s} must be > 0.0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_combine_reject_replace(self, v):
        # Check that v is allowed
        allowed = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_combine_satpix(self, v):
        # Check that v is allowed
        allowed = ['reject', 'force', 'nothing']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_slits_diffpolyorder(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_dispersion_window(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                v = load_sections(v)
            except ValueError:
                msgs.error("The argument of {0:s} must be a 2D region of the form:".format(get_current_name()) +
                           msgs.newline() + "[x1:x2,y1:y2]")
        # Update argument
        self.update(v)
        return

    def trace_dispersion_direction(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                v = int(v)
            except ValueError:
                msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
            if v != 0 and v != 1:
                msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                           "0 or 1 (if the dispersion axis is along a row or column respectively)")
        # Update argument
        self.update(v)
        return

    def trace_slits_fracignore(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        if v < 0.0 or v > 1.0:
            msgs.error("The argument of {0:s} must be between 0 and 1".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_function(self, v):
        # Check that v is allowed
        allowed = ['polynomial', 'legendre', 'chebyshev']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_slits_idsonly(self, v):
        # Check that v is allowed
        if v.lower() == "true":
            v = True
        elif v.lower() == "false":
            v = False
        else:
            msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_maxgap(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        else:
            try:
                v = int(v)
            except ValueError:
                msgs.error("The argument of {0:s} must be of type int, or set to 'none'".format(get_current_name()))
            if v <= 1:
                msgs.error("The argument of {0:s} must be > 1 to set the maximum slit gap".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_number(self, v):
        # Check that v is allowed
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
        # Update argument
        self.update(v)
        return

    def trace_slits_polyorder(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_pca_type(self, v):
        # Check that v is allowed
        allowed = ['pixel', 'order']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_slits_pca_params(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def trace_slits_pca_extrapolate_pos(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_pca_extrapolate_neg(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_sigdetect(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_single(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def trace_slits_tilts_method(self, v):
        # Check that v is allowed
        allowed = ['PCA', 'spline', 'spca', 'interp', 'perp', 'zero']
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_slits_tilts_params(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def trace_slits_tilts_disporder(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_slits_tilts_order(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        if v < 0:
            msgs.error("The argument of {0:s} must be >= 0".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_useframe(self, v):
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        elif v.lower() in ['trace', 'science']:
            v = v.lower()
        else:
            msgs.info("Assuming the following is the name of a trace frame:" + msgs.newline() + v)
        # Update argument
        self.update(v)
        return


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
                    self.update(lst[-1], ll=lst[:-1])
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
        # Update argument
        self.update(v)
        return

    def set_bias(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def set_pixelflat(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def set_science(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def set_slitflat(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def set_standard(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def set_trace(self, v):
        # Check that v is allowed
        v = load_list(v)
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


def get_argflag(init=None):
    """
    Get the Arguments and Flags
    ----------
    init : tuple
      For instantiation

    Returns
    -------
    argflag : Arguments and Flags
    """
    global pypit_argflag

    # Instantiate??
    if init is not None:
        try:
            defname = glob(dirname(__file__))[0] + "/settings/settings." + init[0].lower()
            pypit_argflag = eval(init[0]+"(defname='{0:s}', savname='{1:s}.settings')".format(defname, init[1]))
        except RuntimeError:
            msgs.error("Reduction type '{0:s}' is not allowed".format(init))

    return pypit_argflag


def get_spect(init=None):
    """
    Get the Arguments and Flags
    ----------
    init : tuple
      For instantiation

    Returns
    -------
    argflag : Arguments and Flags
    """
    global pypit_spect

    # Instantiate??
    if init is not None:
        try:
            defname = glob(dirname(__file__))[0] + "/settings/settings." + init[1].lower()
            pypit_spect = eval(init[0]+"_spect(defname='{0:s}', savname='{1:s}.spect')".format(defname, init[2]))
        except RuntimeError:
            msgs.error("{0:s} is not implemented yet".format(init[1]))
    return pypit_spect


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


def load_list(strlist):
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
