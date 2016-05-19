import collections
import inspect
from multiprocessing import cpu_count
import pdb

# Logging
import ardebug
debug = ardebug.init()
import armsgs
msgs = armsgs.get_logger((None, debug, "now", "0.0", 1))


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
        # self._argflag["run"]["ncpus"] = 1
        # self._argflag["a"]["b"]["c"]["d"]["e"] = 5
        # self._argflag["a"]["b"]["c"]["d"]["f"] = 6
        # self._argflag["a"]["b"]["c"]["next"] = 3
        # self._argflag["a"]["b"]["c"]["s"] = "zing"
        # self._argflag["a"]["b"]["g"]["d"]["e"] = "foo"
        # self._argflag["a"]["h"]["c"]["d"]["e"] = "bar"
        # Load the default settings
        self.load()

    def load(self):
        msgs.info("Loading the default settings")
        lines = open(self._defname, 'r').readlines()
        self.load_lines(lines)
        return

    def load_lines(self, lines):
        for ll in lines:
            if len(ll.strip()) == 0:
                # Nothing on a line
                continue
            elif ll.strip()[0] == '#':
                # A comment line
                continue
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
            except:
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
        func = None
        while cnt < len(lst):
            try:
                func = "self." + "_".join(lst[:-cnt]) + "({0:s})".format(" ".join(lst[-cnt:]))
                break
            except:
                cnt += 1
                continue
        if func is not None:
            eval(func)
        else:
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


class BaseSpect:
    def __init__(self):
        self._spect = NestedDict()

    def mosaic_ndet(self, v):
        # Check that v is allowed

        # Update argument
        self.update(v)

    def set_spect(self, lst):
        func = "self." + "_".join(lst[:-1]) + "({0:s})".format(lst[-1])
        eval(func)


class ARMLSD(BaseArgFlag):
    def does_nothing(self):
        return

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
            pypit_argflag = eval(init+"()")
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

"""
af = get_argflag("ARMLSD")
af.set_flag(["run", "ncpus", "5"])
af.set_flag(["a", "b", "c", "d", "e", "47"])
pdb.set_trace()
"""
