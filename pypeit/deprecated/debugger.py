"""
Module to setup the PypeIt debugger
"""
import matplotlib.pyplot as plt
import numpy as np

# These need to be outside of the def's
try:
    from pypeit.ginga.ginga import show_image
except ImportError:  # Ginga is not yet required
    pass
else:
    from pypeit.ginga.ginga import clear_canvas

# ADD-ONs from xastropy
def plot1d(*args, **kwargs):
    """ Plot 1d arrays

    Parameters
    ----------
    outfil= : string
      Outfil
    xlbl,ylbl= : string
      Labels for x,y axes
    xrng= list
      Range of x limits
    yrng= list
      Range of y limits
    xtwo= : ndarray
      x-values for a second array
    ytwo= : ndarray
      y-values for a second array
    mtwo= : str
      marker for xtwo
    scatter= : Bool
      True for a scatter plot
    NOTE: Any extra parameters are fed as kwargs to plt.plot()
    """
    # Error checking
    if len(args) == 0:
        print('x_guis.simple_splot: No arguments!')
        return

    if not isinstance(args[0], np.ndarray):
        print('x_guis: Input array is not a numpy.ndarray!')
        return

    plt_dict = {}

    # Outfil
    if ('outfil' in kwargs):
        plt_dict['outfil'] = kwargs['outfil']
        kwargs.pop('outfil')
    else:
        plt_dict['outfil'] = None

    # Scatter plot?
    if ('scatter' in kwargs):
        kwargs.pop('scatter')
        plt_dict['flg_scatt'] = 1
    else:
        plt_dict['flg_scatt'] = 0

    # Second array?
    if ('xtwo' in kwargs) & ('ytwo' in kwargs):
        plt_dict['xtwo'] = kwargs['xtwo']
        kwargs.pop('xtwo')
        plt_dict['ytwo'] = kwargs['ytwo']
        kwargs.pop('ytwo')
        plt_dict['flg_two'] = 1
        # mtwo
        if 'mtwo' in kwargs:
            plt_dict['mtwo']=kwargs['mtwo']
            kwargs.pop('mtwo')
        else:
            plt_dict['mtwo']=''
    else:
        plt_dict['flg_two'] = 0

    # Limits
    for irng in ['xrng','yrng']:
        try:
            plt_dict[irng] = kwargs[irng]
        except KeyError:
            plt_dict[irng] = None
        else:
            kwargs.pop(irng)

    # Labels
    for ilbl in ['xlbl','ylbl']:
        try:
            plt_dict[ilbl] = kwargs[ilbl]
        except KeyError:
            plt_dict[ilbl] = None
        else:
            kwargs.pop(ilbl)

    # Clear
    plt.clf()
    # Plot it right up
    if len(args) == 1:
        plt.plot(args[0].flatten(), **kwargs)
    else:
        for kk in range(1,len(args)):
            if plt_dict['flg_scatt'] == 0:
                plt.plot(args[0].flatten(),args[kk].flatten(), **kwargs)
            else:
                plt.scatter(args[0].flatten(),args[kk].flatten(), marker='o', **kwargs)

    if plt_dict['flg_two'] == 1:
        plt.plot(plt_dict['xtwo'], plt_dict['ytwo'], plt_dict['mtwo'], color='red', **kwargs)

    # Limits
    if plt_dict['xrng'] is not None:
        plt.xlim(plt_dict['xrng'])
    if plt_dict['yrng'] is not None:
        plt.ylim(plt_dict['yrng'])

    # Label
    if plt_dict['xlbl'] is not None:
        plt.xlabel(plt_dict['xlbl'])
    if plt_dict['ylbl'] is not None:
        plt.ylabel(plt_dict['ylbl'])

    # Output?
    if plt_dict['outfil'] is not None:
        plt.savefig(plt_dict['outfil'])
        print('Wrote figure to {:s}'.format(plt_dict['outfil']))
    else: # Show
        plt.show()
    return

