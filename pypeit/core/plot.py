"""
Convenience functions for plotting.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""


from matplotlib import pyplot as plt
import matplotlib


def pyplot_rcparams():
    """
    params for pretty matplotlib plots
    """
    # set some plotting parameters
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["ytick.direction"] = 'in'
    plt.rcParams["xtick.direction"] = 'in'
    plt.rcParams["xtick.major.size"] = 6
    plt.rcParams["ytick.major.size"] = 6
    plt.rcParams["xtick.minor.size"] = 3
    plt.rcParams["ytick.minor.size"] = 3
    plt.rcParams["xtick.major.width"] = 1
    plt.rcParams["ytick.major.width"] = 1
    plt.rcParams["xtick.minor.width"] = 1
    plt.rcParams["ytick.minor.width"] = 1
    plt.rcParams["axes.linewidth"] = 1
    plt.rcParams["lines.linewidth"] = 3
    plt.rcParams["lines.markeredgewidth"] = 2
    plt.rcParams["patch.linewidth"] = 3
    plt.rcParams["hatch.linewidth"] = 3
    plt.rcParams["font.size"] = 13
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["legend.handletextpad"] = 1


def pyplot_rcparams_default():
    """
    restore default rcparams
    """
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

