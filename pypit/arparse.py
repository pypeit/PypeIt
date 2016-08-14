import collections
import inspect
from multiprocessing import cpu_count
from os.path import dirname, basename, isfile
from textwrap import wrap as wraptext
from glob import glob

# Logging
import ardebug
debug = ardebug.init()
import armsgs
msgs = armsgs.get_logger()


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


class BaseArgFlag:
    def __init__(self, defname, savname):
        self._argflag = NestedDict()
        self._defname = defname
        self._afout = open(savname, 'w')
        # Load the default settings
        self.load()

    def load(self):
        msgs.info("Loading the default settings")
        lines = open(self._defname, 'r').readlines()
        self.load_lines(lines)
        return

    def load_lines(self, lines):
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
            self.set_flag(ll.split())
        return

    def save(self):
        """
        Save the settings used for this reduction
        """
        def savedict(dct):
            for key, value in dct.iteritems():
                self._afout.write(str(key))
                if isinstance(value, dict):
                    savedict(value)
                else:
                    self._afout.write(" " + str(value) + "\n")
        savedict(self._argflag.copy())
        self._afout.close()
        return

    def set_flag(self, lst):
        cnt = 1
        succeed = False
        members = [x for x, y in inspect.getmembers(self, predicate=inspect.ismethod)]
        while cnt < len(lst):
            func = "_".join(lst[:-cnt])
            if func in members:
                func = "self." + func + "('{0:s}')".format(" ".join(lst[-cnt:]))
                print func
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
        Update an element in argflag
        """
        def ingest(dct, upd):
            """
            Ingest the upd dictionary into dct
            """
            for kk, vv in upd.iteritems():
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
        for ii in xrange(1, len(ll)):
            udct = dict({ll[-ii-1]: udct.copy()})
        # Update the master dictionary
        self._argflag = ingest(dstr, udct).copy()
        return

    def arc_combine_match(self, v):
        """
        reduce arcmatch
        """
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_combine_method(self, v):
        """
        arc comb method
        """
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
        """
        arc comb rej_cosmicray
        """
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def arc_combine_reject_lowhigh(self, v):
        """
        arc comb rej_lowhigh
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def arc_combine_reject_level(self, v):
        """
        arc comb rej_level
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def arc_combine_reject_replace(self, v):
        """
        arc comb set_allrej
        """
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
        """
        arc comb sat_pix
        """
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

    def arc_calibrate_lamps(self, v):
        """
        arc calibrate linelist
        """
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

    def arc_calibrate_detection(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
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

    def arc_useframe(self, v):
        """
        reduce usearc
        """
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        # Update argument
        self.update(v)
        return

    def bias_combine_method(self, v):
        """
        bias_comb_method
        """
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
        """
        bias comb rej_cosmicray
        """
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def bias_combine_reject_lowhigh(self, v):
        """
        bias comb rej_lowhigh
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def bias_combine_reject_level(self, v):
        """
        bias comb rej_level
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def bias_combine_reject_replace(self, v):
        """
        bias comb set_allrej
        """
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
        """
        bias comb sat_pix
        """
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
        """
        reduce usebias
        """
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

    def output_verbosity(self, v):
        """
        out verbose
        """
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

    def output_sorted(self, v):
        """
        out sorted
        """
        # Check that v is allowed
        if v.lower() == "none":
            v = None
        # Update argument
        self.update(v)
        return

    def output_overwrite(self, v):
        """
        out overwrite
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

    def pixflat_combine_match(self, v):
        """
        reduce flatmatch
        """
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def pixflat_combine_method(self, v):
        """
        pixflat comb method
        """
        # Check that v is allowed
        allowed = ['mean', 'median', 'weightmean']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def pixflat_combine_reject_cosmics(self, v):
        """
        pixflat comb rej_cosmicray
        """
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def pixflat_combine_reject_lowhigh(self, v):
        """
        pixflat comb rej_lowhigh
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def pixflat_combine_reject_level(self, v):
        """
        pixflat comb rej_level
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def pixflat_combine_reject_replace(self, v):
        """
        pixflat comb set_allrej
        """
        # Check that v is allowed
        allowed = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def pixflat_combine_satpix(self, v):
        """
        pixflat comb sat_pix
        """
        # Check that v is allowed
        allowed = ['reject', 'force', 'nothing']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def pixflat_norm_recnorm(self, v):
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

    def reduce_calibrate(self, v):
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

    def reduce_flatfield_method(self, v):
        """
        reduce FlatMethod
        """
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
        """
        reduce FlatParams
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def reduce_flatfield_perform(self, v):
        """
        reduce flatfield
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

    def reduce_flatfield_useframe(self, v):
        """
        reduce useflat
        """
        # Check that v is allowed
        v = v.lower()
        # Update argument
        self.update(v)
        return

    def reduce_nonlinear(self, v):
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

    def reduce_overscan_method(self, v):
        """
        reduce oscanMethod
        """
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
        """
        reduce oscanParams
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def reduce_pixellocations(self, v):
        """
        reduce locations
        """
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

    def reduce_pixelsize(self, v):
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def reduce_refframe(self, v):
        """
        reduce heliocorr
        """
        # Check that v is allowed
        allowed = ['geocentric', 'heliocentric', 'barycentric']
        if v.lower() not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def reduce_skysub_perform(self, v):
        """
        reduce bgsubtraction perform
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
        """
        run  masterdir
        """
        # Check that v is allowed

        # Update argument
        self.update(v)
        return

    def run_directory_qa(self, v):
        """
        run  plotsdir
        """
        # Check that v is allowed

        # Update argument
        self.update(v)
        return

    def run_directory_science(self, v):
        """
        run  scidir
        """
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
        elif not isfile(v):
                msgs.error("The argument of {0:s} must be a PYPIT spectrograph settings".format(get_current_name()) +
                           msgs.newline() + "file or 'None'. The following file does not exist:" + msgs.newline() + v)
        # Update argument
        self.update(v)
        return

    def run_ncpus(self, v):
        # Check that v is allowed
        curcpu = self._argflag['run']['ncpus']
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

    def run_spectrograph(self, v):
        # Check that v is allowed
        stgs_arm = glob(dirname(__file__)+"/settings.arm*")
        stgs_all = glob(dirname(__file__)+"/settings.*")
        stgs_spc = list(set(stgs_arm) ^ set(stgs_all))
        spclist = [basename(stgs_spc[0]).split(".")[-1].lower()]
        for i in xrange(1, len(stgs_spc)):
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
        """
        run use_idname
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

    def science_load_extracted(self, v):
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

    def science_extraction_method(self, v):
        # Check that v is allowed
        allowed = ['2D', 'mean']
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
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

    def science_extraction_centorder(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def science_extraction_widthorder(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def science_extraction_function(self, v):
        # Check that v is allowed
        allowed = ['polynomial', 'legendre', 'chebyshev']
        if v.lower() not in allowed:
            msgs.error("The argument of {0:s} must be one of:".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def science_extraction_pcacent(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def science_extraction_pcawidth(self, v):
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def science_extraction_bintrace(self, v):
        # Check that v is allowed
        try:
            v = int(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type int".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_combine_match(self, v):
        """
        reduce flatmatch
        """
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_combine_method(self, v):
        """
        trace comb method
        """
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
        """
        trace comb rej_cosmicray
        """
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(get_current_name()))
        # Update argument
        self.update(v)
        return

    def trace_combine_reject_lowhigh(self, v):
        """
        trace comb rej_lowhigh
        """
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
        """
        trace comb rej_level
        """
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
        """
        trace comb set_allrej
        """
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
        """
        trace comb sat_pix
        """
        # Check that v is allowed
        allowed = ['reject', 'force', 'nothing']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_dispersion_window(self, v):
        """
        trace disp window
        """
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
        """
        trace disp direction
        """
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

    def trace_slits_function(self, v):
        """
        trace orders function
        """
        # Check that v is allowed
        allowed = ['polynomial', 'legendre', 'chebyshev']
        v = v.lower()
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_slits_nslits(self, v):
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
        """
        trace orders polyorder  3
        """
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

    def trace_slits_diffpolyorder(self, v):
        """
        trace orders diffpolyorder  3
        """
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
        """
        trace orders sigdetect
        """
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

    def trace_slits_fracignore(self, v):
        """
        trace orders fracignore
        """
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

    def trace_slits_pca_form(self, v):
        """
        trace orders pca
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return


    def trace_slits_pca_extrapolate_pos(self, v):
        """
        trace orders pcxpos
        """
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
        """
        trace orders pcxneg
        """
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

    def trace_slits_tilts_method(self, v):
        """
        trace orders tilts
        """
        # Check that v is allowed
        allowed = ['PCA', 'spline', 'interp', 'perp', 'zero']
        if v not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v)
        return

    def trace_slits_tilts_pcaform(self, v):
        """
        trace orders pcatilt
        """
        # Check that v is allowed
        v = load_list(v)
        # Update argument
        self.update(v)
        return

    def trace_slits_tilts_order(self, v):
        """
        trace orders tiltorder
        """
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


class BaseSpect:
    def __init__(self):
        self._spect = NestedDict()

    def det_xgap(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        # Update argument
        self.update(v)

    def det_ygap(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        # Update argument
        self.update(v)

    def det_ysize(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0.0".format(cname))
        # Update argument
        self.update(v)

    def det_darkcurr(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        # Update argument
        self.update(v)


    def det_ronoise(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v < 0.0:
            msgs.error("The argument of {0:s} must be >= 0.0".format(cname))
        # Update argument
        self.update(v)

    def det_gain(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0.0".format(cname))
        # Update argument
        self.update(v)

    def det_saturation(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v <= 0.0:
            msgs.error("The argument of {0:s} must be > 0.0".format(cname))
        # Update argument
        self.update(v)

    def det_nonlinear(self, v, dnmbr=1):
        cname = get_det_name(dnmbr)
        # Check that v is allowed
        try:
            v = float(v)
        except ValueError:
            msgs.error("The argument of {0:s} must be of type float".format(cname))
        if v <= 0.0 or v > 1.0:
            msgs.error("The argument of {0:s} must be > 0.0 and <= 1.0".format(cname))
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

    def mosaic_reduction(self, v):
        # Check that v is allowed
        allowed = ['ARMLSD', 'ARMED']
        if v.upper() not in allowed:
            msgs.error("The argument of {0:s} must be one of".format(get_current_name()) + msgs.newline() +
                       ", ".join(allowed))
        # Update argument
        self.update(v.upper())

    def set_spect(self, lst):
        cnt = 1
        succeed = False
        members = [x for x, y in inspect.getmembers(self, predicate=inspect.ismethod)]
        while cnt < len(lst):
            func = "_".join(lst[:-cnt])
            # Determine if there are options that need to be passed to this function
            options = ""
            if func[:3] == "det":
                dettmp = func.split("_")[0]
                try:
                    detnum = int(dettmp.lstrip("det"))
                except ValueError:
                    msgs.error("There must be an integer suffix on the det keyword argument" +
                               msgs.newline() + " ".join(lst))
                options += ", dnmbr={0:d}".format(detnum)
                func = func.replace(dettmp, "det")
            lup = ["ampsec", "datasec", "oscansec"]
            for ll in lup:
                if ll in func:
                    lltmp = func.split("_")[1]
                    try:
                        llnum = int(lltmp.lstrip("det"))
                    except ValueError:
                        msgs.error("There must be an integer suffix on the {0:s} keyword argument:".format(ll) +
                                   msgs.newline() + " ".join(lst))
                    options += ", anmbr={0:d}".format(llnum)
                    func = func.replace(lltmp, ll)
            # Now test if this is a function
            if func in members:
                func = "self." + func + "('{0:s}'".format(" ".join(lst[-cnt:]))
                func += options
                func += ")"
                print func
                eval(func)
                succeed = True
                break
            else:
                cnt += 1
        if not succeed:
            msgs.error("There appears to be an error on the following input line:" + msgs.newline() +
                       " ".join(lst))
        return


class ARMLSD(BaseArgFlag):

    def reduce_flexure_maxshift(self, v):
        """
        reduce flexure max_shift
        """
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

    def reduce_fluxcal_perform(self, v):
        """
        reduce fluxcalibrate
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
            defname = glob(dirname(__file__))[0] + "/settings." + init[0].lower()
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
        if init == "ARMLSD":
            pypit_spect = ARMLSD_spect()

    return pypit_spect


def get_current_name():
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    return "'" + " ".join(ll) + "'"


def get_det_name(dnmbr, anmbr=None):
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    tmp = "'" + " ".join(ll) + "'"
    cspl = tmp.split("_")
    cspl[0] += "{0:02d}".format(dnmbr)
    if anmbr is not None:
        cspl[1] += "{0:02d}".format(anmbr)
    return "_".join(cspl)


def load_list(strlist):
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
        xyarrx = [0,0]
    else:
        xyarrx = xyrng[0].split(':')
    if xyrng[1] == ":":
        xyarry = [0,0]
    else:
        xyarry = xyrng[1].split(':')
    return [[int(xyarrx[0]), int(xyarrx[1])], [int(xyarry[0]), int(xyarry[1])]]


"""
af = get_argflag("ARMLSD")
af.set_flag(["run", "ncpus", "5"])
af.set_flag(["run", "qa", "5"])
import pdb
pdb.set_trace()
"""