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
