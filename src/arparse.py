import collections
import inspect
from multiprocessing import cpu_count

# Logging
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
    def __init__(self):
        self._argflag = NestedDict()

    def run_ncpus(self, v):
        # Check that v is allowed
        curcpu = self._argflag['run']['ncpus']
        cpucnt = cpu_count()
        if v == 'all':
            v = cpucnt  # Use all available cpus
            if v != curcpu: msgs.info("Setting {0:d} CPUs".format(v))
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

    def update(self, u):
        inspect.currentframe().f_back.f_code.co_name.split('_')
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                r = self.update(self._argflag.get(k, {}), v)
                self._argflag[k] = r
            else:
                self._argflag[k] = u[k]
        return self._argflag

    def set_flag(self, lst):
        func = "self." + "_".join(lst[:-1]) + "({0:s})".format(lst[-1])
        eval(func)


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
