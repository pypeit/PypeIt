import numpy as np

from pypeit.utils import all_subclasses
from pypeit.scripts import *

# Build the list of script classes
def script_classes():

    # Recursively collect all subclasses
    scr_c = np.array(list(all_subclasses(scriptbase.ScriptBase)))
    scr_n = np.array([c.name() for c in scr_c])
    # Construct a dictionary with the script name and class
    srt = np.argsort(scr_n)
    return dict([ (n,c) for n,c in zip(scr_n[srt],scr_c[srt])])

pypeit_scripts = list(script_classes().keys())

