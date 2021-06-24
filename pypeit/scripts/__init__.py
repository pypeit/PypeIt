
from pypeit.utils import all_subclasses
from pypeit.scripts import scriptbase

# The import of all the script modules here is what enables the dynamic
# compiling of all the available scripts below
from pypeit.scripts import chk_2dslits
from pypeit.scripts import chk_alignments
from pypeit.scripts import chk_edges
from pypeit.scripts import chk_flats
from pypeit.scripts import chk_for_calibs
from pypeit.scripts import chk_wavecalib
from pypeit.scripts import coadd_1dspec
from pypeit.scripts import coadd_2dspec
from pypeit.scripts import coadd_datacube
from pypeit.scripts import collate_1d
from pypeit.scripts import compare_sky
from pypeit.scripts import flux_calib
from pypeit.scripts import flux_setup
from pypeit.scripts import identify
from pypeit.scripts import lowrdx_skyspec
from pypeit.scripts import multislit_flexure
from pypeit.scripts import obslog
from pypeit.scripts import qa_html
from pypeit.scripts import ql_keck_mosfire
from pypeit.scripts import ql_keck_nires
from pypeit.scripts import ql_mos
from pypeit.scripts import run_pypeit
from pypeit.scripts import sensfunc
from pypeit.scripts import setup
from pypeit.scripts import show_1dspec
from pypeit.scripts import show_2dspec
from pypeit.scripts import show_arxiv
from pypeit.scripts import show_wvcalib
from pypeit.scripts import skysub_regions
from pypeit.scripts import tellfit
from pypeit.scripts import trace_edges
from pypeit.scripts import view_fits


# Build the list of script classes
def script_classes():
    import numpy as np

    # Recursively collect all subclasses
    scr_c = np.array(list(all_subclasses(scriptbase.ScriptBase)))
    scr_n = np.array([c.name() for c in scr_c])
    # Construct a dictionary with the script name and class
    srt = np.argsort(scr_n)
    return dict([ (n,c) for n,c in zip(scr_n[srt],scr_c[srt])])

pypeit_scripts = list(script_classes().keys())
