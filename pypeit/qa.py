""" Module for QA in PypeIt

.. include:: ../include/links.rst

"""
import io
import pathlib
import inspect

import numpy as np
import yaml

from matplotlib import pyplot as plt
from matplotlib import lines, colormaps
from matplotlib import gridspec
from matplotlib.lines import Line2D
from astropy.stats import sigma_clipped_stats

from IPython import embed

# CANNOT INCLUDE log IN THIS MODULE AS
#  THE HTML GENERATION OCCURS FROM log
#from pypeit import log

# TODO: Move these names to the appropriate class.  This always writes
# to QA directory, even if the user sets something else...
def set_qa_filename(
    root:str, method:str, det:str|None=None, slit:int|None=None, 
    prefix:str|None=None, mode:str|None=None, out_dir:str|None=None
) -> str:
    """
    Generate the filename for the QA file from the input parameters.

    Parameters
    ----------
    root
        Root name for the output file
    method
        Describes the QA routine
    det
        The name of the detector or mosaic (e.g., DET01)
    slit
        Name of the slit / order being plotted
    prefix
        Start the name of the QA file (used for multiple-PNG PCA plots)
    mode
        Additional differentiating information (*e.g.*, ``gloabl`` vs ``local``
        flexure correction)
    out_dir
        Path to the QA/ directory

    Returns
    -------
        Output filename
    """
    if out_dir is None:
        out_dir = pathlib.Path.cwd()

    match method:
        case 'slit_trace_qa':
            # outfile = 'QA/PNGs/Slit_Trace_{:s}.png'.format(root)
            outfile = 'PNGs/Slit_Trace_{:s}.png'.format(root)

        case 'slit_profile_qa':
            outfile = 'QA/PNGs/Slit_Profile_{:s}_'.format(root)
            # outfile = 'PNGs/Slit_Profile_{:s}_'.format(root)

        case 'arc_fit_qa':
            # outfile = 'QA/PNGs/Arc_1dfit_{:s}_S{:04d}.png'.format(root, slit)
            outfile = 'PNGs/Arc_1dfit_{:s}_S{:04d}.png'.format(root, slit)

        case 'arc_fwhm_qa':
            outfile = 'PNGs/Arc_FWHMfit_{:s}_S{:04d}.png'.format(root, slit)

        case 'plot_orderfits_Arc':  # This is root for multiple PNGs
            outfile = 'QA/PNGs/Arc_lines_{:s}_S{:04d}_'.format(root, slit)
            # outfile = 'PNGs/Arc_lines_{:s}_S{:04d}_'.format(root, slit)

        case 'arc_fit2d_global_qa':
            # outfile = 'QA/PNGs/Arc_2dfit_global_{:s}'.format(root)
            outfile = 'PNGs/Arc_2dfit_global_{:s}'.format(root)

        case 'arc_fit2d_orders_qa':
            # outfile = 'QA/PNGs/Arc_2dfit_orders_{:s}'.format(root)
            outfile = 'PNGs/Arc_2dfit_orders_{:s}'.format(root)

        case 'arc_tilts_spec_qa':
            # outfile = 'QA/PNGs/Arc_tilts_spec_{:s}_S{:04d}.png'.format(root, slit)
            outfile = 'PNGs/Arc_tilts_spec_{:s}_S{:04d}.png'.format(root, slit)

        case 'arc_tilts_spat_qa':
            # outfile = 'QA/PNGs/Arc_tilts_spat_{:s}_S{:04d}.png'.format(root, slit)
            outfile = 'PNGs/Arc_tilts_spat_{:s}_S{:04d}.png'.format(root, slit)

        case 'arc_tilts_2d_qa':
            # outfile = 'QA/PNGs/Arc_tilts_2d_{:s}_S{:04d}.png'.format(root, slit)
            outfile = 'PNGs/Arc_tilts_2d_{:s}_S{:04d}.png'.format(root, slit)

        case 'pca_plot':  # This is root for multiple PNGs
            outfile = 'QA/PNGs/{:s}_pca_{:s}_'.format(prefix, root)
            # outfile = 'PNGs/{:s}_pca_{:s}_'.format(prefix, root)

        case 'pca_arctilt':  # This is root for multiple PNGs
            outfile = 'QA/PNGs/Arc_pca_{:s}_'.format(root)
            # outfile = 'PNGs/Arc_pca_{:s}_'.format(root)

        case 'plot_orderfits_Blaze':  # This is root for multiple PNGs
            outfile = 'QA/PNGs/Blaze_{:s}_'.format(root)
            # outfile = 'PNGs/Blaze_{:s}_'.format(root)

        case 'obj_trace_qa':
            outfile = 'PNGs/{:s}_{:s}_S{:04d}_obj_trace.png'.format(root, det, slit)

        case 'obj_profile_qa':
            outfile = 'PNGs/{:s}_{:s}_S{:04d}_obj_prof.png'.format(root, det, slit)

        case 'spat_flexure_qa_corr':
            outfile = 'QA/PNGs/{:s}_spat_flex_corr.png'.format(root)
            # outfile = 'PNGs/{:s}_spat_flex_corr.png'.format(root)

        case 'spec_flexure_qa_corr':
            # outfile = 'QA/PNGs/{:s}_D{:02d}_S{:04d}_spec_flex_corr.png'.format(root, det, slit)
            outfile = 'PNGs/{:s}_{:s}_{:s}_S{:04d}_spec_flex_corr.png'.format(root, mode, det, slit)

        case 'spec_flexure_qa_sky':
            # outfile = 'QA/PNGs/{:s}_D{:02d}_S{:04d}_spec_flex_sky.png'.format(root, det, slit)
            outfile = 'PNGs/{:s}_{:s}_{:s}_S{:04d}_spec_flex_sky.png'.format(root, mode, det, slit)

        case 'spatillum_finecorr':
            outfile = 'PNGs/{:s}_S{:04d}_spatillum_finecorr.png'.format(root, slit)

        case 'detector_structure':
            outfile = 'PNGs/{:s}_{:s}_detector_structure.png'.format(root, det)

        case _:
            raise IOError("NOT READY FOR THIS QA: {:s}".format(method))

    # Return
    return str(pathlib.Path(out_dir) / outfile)


def get_dimen(x:int, maxp:int=25) -> tuple[list,list]:
    """ Assign the plotting dimensions to be the "most square"

    Parameters
    ----------
    x : int
      An integer that equals the number of panels to be plot
    maxp : int (optional)
      The maximum number of panels to plot on a single page

    Returns
    -------
    pages : list
      The number of panels in the x and y direction on each page
    npp : list
      The number of panels on each page
    """
    pages, npp = [], []
    xr = x
    while xr > 0:
        if xr > maxp:
            xt = maxp
        else:
            xt = xr
        ypg = int(np.sqrt(float(xt)))
        if int(xt) % ypg == 0:
            xpg = int(xt)/ypg
        else:
            xpg = 1 + int(xt)/ypg
        pages.append([int(xpg), int(ypg)])
        npp.append(int(xt))
        xr -= xt
    return pages, npp


def html_header(title:str) -> str:
    """
    Generate a simple HTML header
    
    Parameters
    ----------
    title : str
        Title for the header

    Returns
    -------
    head : str
        An HTML header as a long string

    """
    head = '<?xml version="1.0" encoding="UTF-8"?>\n'
    head += '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n'


    head += '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n'
    head += '\n'
    head += '<head>\n'
    head += '\n'
    head += '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />\n'
    head += '<title>{:s}</title>\n'.format(title)
    head += '<meta name="keywords" content="" />\n'
    head += '<meta name="description" content="" />\n'
    head += '<script type="text/javascript" src="jquery/jquery-1.4.2.min.js"></script>\n'
    head += '<script type="text/javascript" src="jquery/jquery.slidertron-0.1.js"></script>\n'
    head += '<link href="style.css" rel="stylesheet" type="text/css" media="screen" />\n'
    head += '\n'
    head += '</head>\n'

    # Begin the Body
    head += '<body>\n'
    head += '<h1>{:s}</h1>\n'.format(title)
    head += '<hr>\n'

    return head

def html_end(f:io.TextIOWrapper, body:str, links:str|None=None) -> str:
    """
    Fill in the HTML file with a proper ending

    Parameters
    ----------
    f : `io.TextIOWrapper`_
    body : str
    links : str, optional

    Returns
    -------
    end : str
        The text written to the end of the HTML file
    
    """
    # Write links
    if links is not None:
        f.write(links)
        f.write('</ul>\n')
        f.write('<hr>\n')
    # Write body
    f.write(body)
    # Finish
    end = '</body>\n'
    end += '</html>\n'
    f.write(end)

    return end


def html_init(f:io.TextIOWrapper, title:str) -> str:
    """
    Initialize the HTML file

    Args:
        f (`io.TextIOWrapper`_):
            file object to write to
        title (str):
            title

    Returns:
        str: Initial HTML text incluing the header and links
    """
    head = html_header(title)
    f.write(head)
    # Init links
    links = '<h2>Quick Links</h2>\n'
    links += '<ul>\n'
    return links


def html_mf_pngs(idval:str) -> tuple[str,str]:
    """ Generate HTML for QA PNGs

    Args:
        idval: str
            Key identifier of the calibration set

    Returns:
        tuple: 

          - links -- HTML links to the PNGs
          - body -- HTML edits for the main body

    """
    links = ''
    body = ''

    # Organize the outputs
    html_dict = {
        'strace': dict(
            fname='slit_trace_qa', ext='', href='strace', label='Slit Trace', slit=False
        ),
        'sprof': dict(
            fname='slit_profile_qa', ext='*.png', href='sprof', label='Slit Profile', slit=False
        ),
        'blaze': dict(
            fname='plot_orderfits_Blaze', ext='*.png', href='blaze', label='Blaze', slit=False
        ),
        'arc_fit': dict(
            fname='arc_fit_qa', ext='', href='arc_fit', label='Arc 1D Fit', slit=True
        ),
        'arc_tilts_spec': dict(
            fname='arc_tilts_spec_qa', ext='', href='arc_tilts_spec', label='Arc Tilts Spec', slit=True
        ),
        'arc_tilts_spat': dict(
            fname='arc_tilts_spat_qa', ext='', href='arc_tilts_spat', label='Arc Tilts Spat', slit=True
        ),
        'arc_tilts_2d': dict(
            fname='arc_tilts_2d_qa', ext='', href='arc_tilts_2d', label='Arc Tilts 2D', slit=True
        ),
        'arc_pca': dict(
            fname='pca_arctilt', ext='*.png', href='arc_pca', label='Arc Tilt PCA', slit=False
        ),
        'arc_fit2d_global': dict(
            fname='arc_fit2d_global_qa', ext='*.png', href='arc_fit2d_global', label='2D Arc Global', slit=False
        ),
        'arc_fit2d_orders': dict(
            fname='arc_fit2d_orders_qa', ext='*.png', href='arc_fit2d_orders', label='2D Arc Orders', slit=False
        ),
    }

    # Generate HTML
    for key in ['strace', 'sprof', 'blaze', 'arc_fit', 'arc_pca', 'arc_fit2d_global', 'arc_fit2d_orders',
                'arc_tilts_spec', 'arc_tilts_spat', 'arc_tilts_2d']:
        # PNG Root
        png_fileroot = set_qa_filename(idval, html_dict[key]['fname'], slit=9999, out_dir='QA')
        if html_dict[key]['slit']:  # Kludge to handle multiple slits
            png_fileroot = png_fileroot.replace('S9999', 'S*')
        # Find the PNGs
        png_path, png_stem = pathlib.Path(png_fileroot).parent, pathlib.Path(png_fileroot).name
        pngs = sorted(pathlib.Path(png_path).glob(f"{png_stem}{html_dict[key]['ext']}"))
        if len(pngs) > 0:
            href="{:s}_{:s}".format(html_dict[key]['href'], idval)
            # Link
            links += '<li><a class="reference internal" href="#{:s}">{:s} {:s}</a></li>\n'.format(
                href, html_dict[key]['label'], idval)
            # Body
            body += '<hr>\n'
            body += '<div class="section" id="{:s}">\n'.format(href)
            body += '<h2> {:s} {:s} </h2>\n'.format(html_dict[key]['label'], idval)
            for png in pngs:
                # Remove QA
                if 'QA' not in [p.name for p in png.parents]:
                    raise ValueError("QA is expected to be in the path!")
                parent_dir = pathlib.Path(png.parent.name)
                if html_dict[key]['slit']:  # Kludge to handle multiple slits
                    slit_name = png.name[png.name.find(f"{idval}_S"):]
                    href = f"{html_dict[key]['href']}_{slit_name}"

                    body += '<img class ="research" src="{:s}" width="100%" id={:s} height="auto"/>\n'.format(
                        str(parent_dir / png.name), href)
                    links += '<li><a class="reference internal" href="#{:s}">{:s} {:s}</a></li>\n'.format(
                        href, html_dict[key]['label'], pathlib.Path(slit_name).stem)
                else:
                    body += '<img class ="research" src="{:s}" width="100%" height="auto"/>\n'.format(str(parent_dir / png.name))
            body += '</div>\n'

    # Return
    return links, body


def html_exp_pngs(exp_name:str, det:int) -> tuple[str,str]:
    """
    Generate HTML for Exposure PNGs

    Parameters
    ----------
    exp_name
        PypeIt-standard exposure name
    det
        Detector number

    Returns
    -------
    links
        Links to the individual images
    body
        Body HTML for the page showing the images

    """
    det_str = f"DET{det:02d}"
    links = ''
    body = ''

    # Organize the outputs
    html_dict = {
        'trace': dict(
            fname='obj_trace_qa', ext='', slit=True, href='otrace', label='Object Traces'
        ),
        'prof': dict(
            fname='obj_profile_qa', ext='', slit=True, href='oprofile', label='Object Profiles'
        ),
        'flex_corr': dict(
            fname='spec_flexure_qa_corr', ext='', slit=True, href='flex_corr', label='Flexure Cross Correlation'
        ),
        'flex_sky': dict(
            fname='spec_flexure_qa_sky', ext='', slit=True, href='flex_sky', label='Flexure Sky'
        ),
    }

    # Generate HTML
    for key in ['trace', 'prof', 'flex_corr', 'flex_sky']:
        # PNG Root
        png_fileroot = set_qa_filename(exp_name, html_dict[key]['fname'], det=det_str, slit=9999, mode="*", out_dir='QA')
        if html_dict[key]['slit']:  # Kludge to handle multiple slits
            png_fileroot = png_fileroot.replace('S9999', 'S*')
        # Find the PNGs
        png_path, png_stem = pathlib.Path(png_fileroot).parent, pathlib.Path(png_fileroot).name
        pngs = sorted(pathlib.Path(png_path).glob(f"{png_stem}{html_dict[key]['ext']}"))
        if len(pngs) > 0:
            href="{:s}_{:02d}".format(html_dict[key]['href'], det)
            # Link
            links += '<li><a class="reference internal" href="#{:s}">{:s} {:02d}</a></li>\n'.format(href, html_dict[key]['label'], det)
            # Body
            body += '<hr>\n'
            body += '<div class="section" id="{:s}">\n'.format(href)
            body += '<h2> {:s} {:02d} </h2>\n'.format(html_dict[key]['label'], det)
            for png in pngs:
                # Remove QA
                if 'QA' not in [p.name for p in png.parents]:
                    raise ValueError("QA is expected to be in the path!")
                parent_dir = pathlib.Path(png.parent.name)
                body += '<img class ="research" src="{:s}" width="100%" height="auto"/>\n'.format(str(parent_dir / png.name))
            body += '</div>\n'

    # Return
    return links, body


def gen_qa_dir(qa_path:str):
    """ Make the QA directory if it doesn't already exist

    Args:
        qa_path (str):
            Path to the QA folder
    """
    if not (the_path := pathlib.Path(qa_path)).exists():
        the_path.mkdir(parents=True, exist_ok=True)

# TODO: Need to revisit this...
def gen_mf_html(pypeit_file:str, qa_path:str):
    """ Generate the HTML for QA

    Args:
        pypeit_file (:obj:`str`):
            Name of the PypeIt file, no path
        qa_path (:obj:`str`):
            Path to the QA folder
    """
    # TODO: Can this instead just use the pypeit file?
    # Read calib file
    calib_file = pypeit_file.replace('.pypeit', '.calib')
    with open(calib_file, 'r') as infile:
        calib_dict = yaml.load(infile, Loader=yaml.FullLoader)
    # Parse
    setup = list(calib_dict.keys())[0]
    cbsets = []
    for key in calib_dict[setup].keys():
        if key == '--':
            continue
        #if isinstance(key,str):
        #    dets.append(int(key))
        else:
            cbsets.append(key)
    # TODO -- Read in spectograph from .pypeit file and then use spectrograph.ndet
    dets = (1+np.arange(20)).tolist()
    mscs = (1+np.arange(5)).tolist()
    # Generate MF file
    MF_filename = pathlib.Path(qa_path) / f"MF_{setup}.html"
    body = ''
    with open(MF_filename,'w') as f:
        # Start
        links = html_init(f, 'QA Setup {:s}: Calibration files'.format(setup))
        # Loop on calib_sets
        for cbset in cbsets:
            for det in dets:
                # Run
                idval = '{:s}_{:d}_DET{:02d}'.format(setup, cbset, det)
                new_links, new_body = html_mf_pngs(idval)
                # Save
                links += new_links
                body += new_body
            for msc in mscs:
                # Run
                idval = '{:s}_{:d}_MSC{:02d}'.format(setup, cbset, msc)
                new_links, new_body = html_mf_pngs(idval)
                # Save
                links += new_links
                body += new_body
        # End
        html_end(f, body, links)
    #
    print(f"Wrote: {MF_filename}")


def gen_exp_html():
    """ Generate the HTML for an Exposure set
    """
    # Find all obj_trace files -- Not fool proof but ok
    # NOTE: At some point, the obj_trace QA was removed from the repo.  Adding
    #       it back reactivates this code.  (TEB, 21-Oct-2025)
    obj_files = sorted((pathlib.Path("QA") / "PNGs").glob("*obj_trace.png"))
    # Parse for names
    uni_names = np.unique([obj_file.name.split("_DET")[0] for obj_file in obj_files])
    # Loop
    for uni_name in uni_names:
        # Generate MF file
        exp_filename = f"QA/{uni_name}.html"
        body = ""
        with open(exp_filename, "w", encoding="utf-8") as f_obj:
            # Start
            links = html_init(f_obj, f"QA for {uni_name}")
            # Loop on detector
            for det in range(1,99):
                # Run
                new_links, new_body = html_exp_pngs(uni_name, det)
                # Save
                links += new_links
                body += new_body
            # End
            html_end(f_obj, body, links)
        print(f"Wrote: {exp_filename}")


def close_qa(pypeit_file:str, qa_path:str):
    """
    Tie off QA under a crash

    Args:
        pypeit_file (str):
            PypeIt file name
        qa_path (str):
            Path to QA directory
    """
    if pypeit_file is None:
        return
    try:
        gen_mf_html(pypeit_file, qa_path)
    except:  # Likely crashed real early
        pass
    else:
        gen_exp_html()


# This method needs to match the name in set_qa_filename()
def arc_tilts_2d_qa(tilts_dspat, tilts, tilts_model, tot_mask, rej_mask, spat_order, spec_order, rms, fwhm,
                 slitord_id=0, setup='A', outfile=None, show_QA=False, out_dir=None):
    """

    ..todo.. this method needs docs

    Args:
        tilts_dspat:
        tilts:
        tilts_model:
        tot_mask:
        rej_mask:
        spat_order:
        spec_order:
        rms:
        fwhm:
        slitord_id:
        setup:
        outfile:
        show_QA:
        out_dir:

    Returns:

    """
    plt.rcdefaults()
    plt.rcParams['font.family'] = 'sans-serif'

    # Outfile
    method = inspect.stack()[0][3]
    if (outfile is None):
        outfile = set_qa_filename(setup, method, slit=slitord_id, out_dir=out_dir)

    # Show the fit
    fig, ax = plt.subplots(figsize=(12, 18))
    ax.cla()
    ax.plot(tilts_dspat[tot_mask], tilts[tot_mask], color='black', linestyle=' ', mfc='None', marker='o',
            markersize=9.0, markeredgewidth=1.0, zorder=4, label='Good Tilt')
    ax.plot(tilts_dspat[rej_mask], tilts[rej_mask], color='red', linestyle=' ', mfc='None', marker='o',
            markersize=9.0, markeredgewidth=2.0, zorder=5, label='Rejected')
    ax.plot(tilts_dspat[tot_mask], tilts_model[tot_mask], color='black', linestyle=' ', marker='o',
            markersize=2.0, markeredgewidth=1.0, zorder=1, label='2D Model')

    xmin = 1.1 * tilts_dspat[tot_mask].min()
    xmax = 1.1 * tilts_dspat[tot_mask].max()
    ax.set_xlim((xmin, xmax))
    ax.set_xlabel('Spatial Offset from Central Trace (pixels)', fontsize=15)
    ax.set_ylabel('Spectral Pixel', fontsize=15)
    ax.legend()
    ax.set_title('Tilts vs Fit (spat_order, spec_order)=({:d},{:d}) for slit={:d}: RMS = {:5.3f}, '
                 'RMS/FWHM={:5.3f}'.format(spat_order, spec_order, slitord_id, rms, rms / fwhm), fontsize=15)

    # Finish
    # plt.tight_layout(pad=1.0, h_pad=1.0, w_pad=1.0)

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    if show_QA:
        plt.show()

    plt.close()
    plt.rcdefaults()


# This method needs to match the name in set_qa_filename()
def arc_tilts_spec_qa(tilts_spec_fit, tilts, tilts_model, tot_mask, rej_mask, rms, fwhm,
                   slitord_id=0, setup='A', outfile=None, show_QA=False, out_dir=None):
    """ Generate a QA plot of the residuals for the fit to the tilts in the spectral direction one slit at a time

    Parameters
    ----------
    """

    plt.rcdefaults()
    plt.rcParams['font.family'] = 'sans-serif'

    # Outfil
    method = inspect.stack()[0][3]
    if (outfile is None):
        outfile = set_qa_filename(setup, method, slit=slitord_id, out_dir=out_dir)

    # Setup
    plt.figure(figsize=(14, 6))
    plt.clf()
    ax = plt.gca()

    # Scatter plot
    res = (tilts - tilts_model)

    nspat, nuse = tilts.shape
    # Show the fit residuals as a function of spatial position
    line_indx = np.outer(np.ones(nspat), np.arange(nuse))

    xmin = 0.90 * (tilts_spec_fit.min())
    xmax = 1.10 * (tilts_spec_fit.max())

    ax.hlines(0.0, xmin, xmax, linestyle='--', color='green')

    for iline in range(nuse):
        iall = (line_indx == iline) & tot_mask
        igd = (line_indx == iline) & tot_mask & (rej_mask == False)
        irej = (line_indx == iline) & tot_mask & rej_mask

        ax.plot(tilts_spec_fit[igd], (res[igd]), 'ko', mfc='k', markersize=4.0)
        ax.plot(tilts_spec_fit[irej], (res[irej]), 'ro', mfc='r', markersize=4.0)
        # Compute the RMS for this line
        all_rms = np.std(res[iall])
        good_rms = np.std(res[igd])
        # ToDo show the mean here as well
        if np.any(igd):
            ax.plot(tilts_spec_fit[igd][0], all_rms, marker='s', linestyle=' ', color='g', mfc='g', markersize=7.0)
            ax.plot(tilts_spec_fit[igd][0], good_rms, marker='^', linestyle=' ', color='orange', mfc='orange',
                    markersize=7.0)

    ax.text(0.90, 0.90, 'Slit {:d}:  Residual (pixels) = {:0.5f}'.format(slitord_id, rms),
            transform=ax.transAxes, ha='right', color='black', fontsize=16)
    ax.text(0.90, 0.80, ' Slit {:d}:  RMS/FWHM = {:0.5f}'.format(slitord_id, rms / fwhm),
            transform=ax.transAxes, ha='right', color='black', fontsize=16)
    # Label
    ax.set_xlabel('Spectral Pixel')
    ax.set_ylabel('RMS (pixels)')
    ax.set_title('RMS of Each Arc Line Traced')
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((-5.0 * rms, 5.0 * rms))
    # Legend
    legend_elements = [lines.Line2D([0], [0], linestyle=' ', color='k', marker='o', mfc='k',
                                    markersize=4.0, label='good'),
                       lines.Line2D([0], [0], linestyle=' ', color='r', marker='o', mfc='r',
                                    markersize=4.0, label='rejected'),
                       lines.Line2D([0], [0], linestyle=' ', color='g', marker='s', mfc='g',
                                    markersize=7.0, label='all RMS'),
                       lines.Line2D([0], [0], linestyle=' ', color='orange', marker='^',
                                    mfc='orange', markersize=7.0, label='good RMS')]
    ax.legend(handles=legend_elements)

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    if show_QA:
        plt.show()

    plt.close()
    plt.rcdefaults()


def arc_tilts_spat_qa(tilts_dspat, tilts, tilts_model, tilts_spec_fit, tot_mask, rej_mask,
                      spat_order, spec_order, rms, fwhm, setup='A', slitord_id=0, outfile=None,
                      show_QA=False, out_dir=None):
    """
    NEEDS A DOC STRING!
    """
    plt.rcdefaults()
    plt.rcParams['font.family'] = 'sans-serif'

    # Output file
    method = inspect.stack()[0][3]
    if outfile is None:
        outfile = set_qa_filename(setup, method, slit=slitord_id, out_dir=out_dir)

    nspat, nuse = tilts_dspat.shape
    # Show the fit residuals as a function of spatial position
    line_indx = np.outer(np.ones(nspat), np.arange(nuse))
    lines_spec = tilts_spec_fit[0, :]
    cmap = colormaps['coolwarm'].resampled(nuse)

    fig, ax = plt.subplots(figsize=(14, 12))
    # dummy mappable shows the spectral pixel
    dummie_cax = ax.scatter(lines_spec, lines_spec, c=lines_spec, cmap=cmap)
    ax.cla()

    for iline in range(nuse):
        iall = (line_indx == iline) & tot_mask
        irej = (line_indx == iline) & tot_mask & rej_mask
        this_color = cmap(iline)
        # plot the residuals
        ax.plot(tilts_dspat[iall], tilts[iall] - tilts_model[iall], color=this_color,
                linestyle='-', linewidth=3.0, marker='None', alpha=0.5)
        ax.plot(tilts_dspat[irej], tilts[irej] - tilts_model[irej], linestyle=' ',
                marker='o', color='limegreen', mfc='limegreen', markersize=5.0)

    xmin = 1.1 * tilts_dspat[tot_mask].min()
    xmax = 1.1 * tilts_dspat[tot_mask].max()
    ax.hlines(0.0, xmin, xmax, linestyle='--', linewidth=2.0, color='k', zorder=10)

    ax.set_xlim((xmin, xmax))
    ax.set_xlabel('Spatial Offset from Central Trace (pixels)')
    ax.set_ylabel('Arc Line Tilt Residual (pixels)')

    legend_elements = [lines.Line2D([0], [0], color='cornflowerblue', linestyle='-', linewidth=3.0,
                                    label='residual'),
                       lines.Line2D([0], [0], color='limegreen', linestyle=' ', marker='o',
                                    mfc='limegreen', markersize=7.0, label='rejected')]
    ax.legend(handles=legend_elements)
    ax.set_title('Tilts vs Fit (spat_order, spec_order)=({:d},{:d}) for slit={:d}: RMS = {:5.3f}, '
                 'RMS/FWHM={:5.3f}'.format(spat_order, spec_order, slitord_id, rms, rms / fwhm), fontsize=15)
    cb = fig.colorbar(dummie_cax, ax=ax, ticks=lines_spec)
    cb.set_label('Spectral Pixel')

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    if show_QA:
        plt.show()

    plt.close()
    plt.rcdefaults()


# TODO: With Python 3.14's deferred evaluation of annotations, may be able
#       to annotate `specobjs`; however, should really remove PypeIt-specific
#       objects from `core`.
def spec_flexure_qa(slitords:np.ndarray, bpm:np.ndarray, basename:str,
                    flex_list:list[dict], specobjs=None,
                    out_dir:str|None=None):
    """
    Generate QA for the spectral flexure calculation

    Args:
        slitords (`numpy.ndarray`_):
            Array of slit/order numbers
        bpm (`numpy.ndarray`_):
            Boolean mask; True = masked slit
        basename (str):
            Used to generate the output file name
        flex_list (list):
            list of :obj:`dict` objects containing the flexure information
        specobjs (:class:`~pypeit.specobjs.SpecObjs`, optional):
            Spectrally extracted objects
        out_dir (str, optional):
            Path to the output directory for the QA plots.  If None, the current
            is used.
    """
    # Extract the mode and detector from the ``basename``
    *_, mode, det = basename.split("_")

    plt.rcdefaults()
    plt.rcParams['font.family'] = 'serif'

    # What type of QA are we doing
    slit_cen = specobjs is None

    # Grab the named of the method
    method = inspect.stack()[0][3]

    # Mask
    gdslits = np.where(np.logical_not(bpm))[0]

    # Loop over slits, and then over objects here
    for islit in gdslits:
        # Slit/order number
        slitord = slitords[islit]

        this_flex_dict = flex_list[islit]
        # Check that the default was overwritten
        if len(this_flex_dict['spec_flexure']) == 0 or \
                (len(this_flex_dict['spec_flexure']) > 0 and np.all([ss is None for ss in this_flex_dict['spec_flexure']])):
            continue

        # Parse and Setup
        if slit_cen:
            nobj = 1
            ncol = 1
        else:
            indx = specobjs.slitorder_indices(slitord)
            this_specobjs = specobjs[indx]
            nobj = np.sum(indx)
            if nobj == 0:
                continue
            ncol = min(3, nobj)

        nrow = nobj // ncol + ((nobj % ncol) > 0)
        # Outfile, one QA file per slit
        outfile = set_qa_filename(
            basename, method + '_corr', slit=slitord, det=det, mode=mode, out_dir=out_dir
        )
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)
        # Correlation QA
        if slit_cen:
            ax = plt.subplot(gs[0, 0])
            spec_flexure_corrQA(ax, this_flex_dict, 0, 'Slit Center')
        else:
            iplt = 0
            for ss, specobj in enumerate(this_specobjs):
                if specobj is None or (specobj.BOX_WAVE is None and specobj.OPT_WAVE is None):
                    continue
                ax = plt.subplot(gs[iplt//ncol, iplt % ncol])
                spec_flexure_corrQA(ax, this_flex_dict, ss, '{:s}'.format(specobj.NAME))
                iplt += 1
        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile)#, dpi=400)
        plt.close()

        # Sky line QA (just one object)
        if slit_cen:
            iobj = 0
        else:
            # only show the first object in this slit that does not have None shift
            iobj = np.where([ss is not None for ss in this_flex_dict['spec_flexure']])[0][0]
            specobj = this_specobjs[iobj]

        # Repackage
        sky_spec = this_flex_dict['sky_spec'][iobj]
        arx_spec = this_flex_dict['arx_spec'][iobj]
        min_wave = max(np.amin(arx_spec.wave), np.amin(sky_spec.wave))
        max_wave = min(np.amax(arx_spec.wave), np.amax(sky_spec.wave))

        # Sky lines
        # TODO: Should these be defined / identified somewhere else?  Then they
        #       could more easily be included in the documentation.
        sky_lines = np.array([3370.0, 3914.0, 4046.56, 4358.34, 5577.338, 6300.304,
                              7340.885, 7993.332, 8430.174, 8919.610, 9439.660,
                              10013.99, 10372.88])
        dwv = 20.
        gdsky = np.where((sky_lines > min_wave) & (sky_lines < max_wave))[0]
        if len(gdsky) == 0:
            #log.warning("No sky lines for Flexure QA")
            continue
        if len(gdsky) > 6:
            idx = np.array([0, 1, len(gdsky)//2, len(gdsky)//2+1, -2, -1])
            gdsky = gdsky[idx]

        # Outfile
        outfile = set_qa_filename(
            basename, method+'_sky', slit=slitord, det=det, mode=mode, out_dir=out_dir
        )
        # Figure
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        nrow, ncol = 2, 3
        gs = gridspec.GridSpec(nrow, ncol)
        if slit_cen:
            plt.suptitle('Sky Comparison for Slit Center', y=0.99)
        else:
            plt.suptitle('Sky Comparison for {:s}'.format(specobj.NAME), y=0.99)

        for ii, igdsky in enumerate(gdsky):
            skyline = sky_lines[igdsky]
            ax = plt.subplot(gs[ii//ncol, ii % ncol])
            # Norm
            pix1 = np.where(np.abs(sky_spec.wave-skyline) < dwv)[0]
            pix2 = np.where(np.abs(arx_spec.wave-skyline) < dwv)[0]
            f1 = np.sum(sky_spec.flux[pix1])
            f2 = np.sum(arx_spec.flux[pix2])
            norm = f1/f2
            # Plot
            ax.plot(sky_spec.wave[pix1], sky_spec.flux[pix1], 'k-', label='Obj',
                    drawstyle='steps-mid')
            ax.plot(arx_spec.wave[pix2], arx_spec.flux[pix2]*norm, 'r-', label='Arx',
                    drawstyle='steps-mid')
            # Axes
            ax.xaxis.set_major_locator(plt.MultipleLocator(dwv))
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Counts')

        # Legend
        plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                   handletextpad=0.3, fontsize='small', numpoints=1)

        # Finish
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile)#, dpi=400)
        plt.close()
        #log.info("Wrote spectral flexure QA: {}".format(outfile))

    plt.rcdefaults()


def spec_flexure_corrQA(ax:plt.Axes, this_flex_dict:dict, cntr:int, name:str):
    """Spectral Flexure QA Plot

    Creates one panel of the spectral felxure QA plot, with the overall figure
    container being handled by the calling function.

    Parameters
    ----------
    ax
        Axes onto which to draw the plot
    this_flex_dict
        Dictionary of flexure-related information needed for the plot
    cntr
        The index into ``this_flex_dict``'s arrays corresponding to the
        particular object, trace, or location of interest.
    name
        Object, trace, or location name to be printed in the plot
    """
    # Fit
    fit = this_flex_dict['polyfit'][cntr]
    if fit is not None:
        xval = np.linspace(-10., 10, 100) + this_flex_dict['corr_cen'][cntr] + this_flex_dict['spec_flexure'][cntr]
        # model = (fit[2]*(xval**2.))+(fit[1]*xval)+fit[0]
        model = fit.eval(xval)
        # model = utils.func_val(fit, xval, 'polynomial')
        mxmod = np.max(model)
        ylim_min = np.min(model / mxmod) if np.isfinite(np.min(model / mxmod)) else 0.0
        ylim = [ylim_min, 1.3]
        ax.plot(xval - this_flex_dict['corr_cen'][cntr], model / mxmod, 'k-')
        # Measurements
        ax.scatter(this_flex_dict['subpix'][cntr] - this_flex_dict['corr_cen'][cntr],
                   this_flex_dict['corr'][cntr] / mxmod, marker='o')
        # Final spectral flexure shift
        ax.plot([this_flex_dict['spec_flexure'][cntr]] * 2, ylim, 'g:')
        # Label
        ax.text(0.5, 0.25, name, transform=ax.transAxes, size='large', ha='center')
        ax.text(0.5, 0.15, 'spec_flexure = {:g}'.format(this_flex_dict['spec_flexure'][cntr]),
                transform=ax.transAxes, size='large', ha='center')  # , bbox={'facecolor':'white'})
        # Axes
        ax.set_ylim(ylim)
        ax.set_xlabel('Lag')
    else:
        ax.text(0.5, 0.25, name, transform=ax.transAxes, size='large', ha='center')
        ax.text(0.5, 0.15, 'spectral flex_shift calculation failed', transform=ax.transAxes, size='large', ha='center')
        # Axes
        ax.set_xlabel('Lag')


def spat_flexure_qa(img, slits, spat_flexure, gpm=None, vrange=None, outfile=None):
    """
    Generate QA for the spatial flexure

    Args:
        img (`numpy.ndarray`_):
            Image of the detector
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
            Slits object
        spat_flexure (`numpy.ndarray`_):
            Shift in pixels for each slit edge. Shape is (nslits, 2), where the 2 refers to
            the left [0] and right[1] slits.
        gpm (`numpy.ndarray`_, optional):
            Good pixel mask (True = Bad)
        vrange (:obj:`tuple`, optional):
            Tuple with the min and max values for the imshow plot
        outfile (:obj:`str`, optional):
            Path to the output file where the QA is saved.  If None, the QA is
            shown on screen and not saved.
    """
    debug = True if outfile is None else False

    # check that vrange is a tuple
    if vrange is not None and not isinstance(vrange, tuple):
        #log.warning('vrange must be a tuple with the min and max values for the imshow plot. Ignoring vrange.')
        vrange = None

    # Spatial flexure is *always* relative to the initial slits, otherwise it pushes the tilts
    # to be out of sync with the slit edges.
    left_slits, right_slits, mask_slits = slits.select_edges(initial=True, spat_flexure=None)
    left_flex, right_flex, mask = slits.select_edges(initial=True, spat_flexure=spat_flexure)

    if debug:
        # where to start and end the plot in the spatial&spectral direction
        nxsnip = 1
        spat_starts = [0]
        spat_ends = [img.shape[1]]
        upper_ystart = 0
        upper_yend = img.shape[0]

    else:
        # where to start and end the plot in the spatial direction
        xstart = int(np.floor(np.min([left_slits, left_flex]) - 20))
        xend = int(np.ceil(np.max([right_slits, right_flex]) + 20))

        # how many snippets to plot in the spatial direction
        if slits.nslits == 1:
            # if longslit plot 2 snippets, one for the left edge and one for the right edge
            nxsnip = 2
            snippet = int((xend - xstart) // nxsnip)
            spat_starts = [xstart, xstart + snippet]
            spat_ends = [xend - snippet, xend]
        elif slits.nslits <= 12:
            # if 12 or less slits plot 3-4 snippets equally spaced
            nxsnip = 3 if slits.nslits <= 6 else 4
            snippet = int((xend - xstart) // nxsnip)
            spat_starts = [xstart, xstart + snippet, xstart + 2*snippet]
            spat_ends = [xend - 2*snippet, xend - snippet, xend]
            if slits.nslits > 6:
                # add the 4th snippet
                spat_starts.append(xstart + 3*snippet)
                spat_ends.insert(0, xend - 3*snippet)
        else:
            # if more than 12 slits plot 4 snippets
            nxsnip = 4
            # approximately, we want 3 slits in each snippet
            snippet = int(3 * (xend - xstart)/slits.nslits)
            # this would give nx many snippets
            nx = int((xend - xstart) // snippet)
            # but we want to plot only nxsnip of those snippets
            spat_starts = [xstart + i * snippet for i in np.linspace(0, nx - 1, nxsnip, dtype=int)]
            spat_ends = [xstart + i * snippet for i in np.linspace(1, nx, nxsnip, dtype=int)]

        # where to start and end the plot in the spectral direction for both the upper and lower sections
        lower_ystart = 0
        lower_yend = int(snippet)
        upper_ystart = int(img.shape[0] - snippet)
        upper_yend = img.shape[0]

    # plot the spatial flexure
    rows = 1 if debug else 2
    fig = plt.figure(figsize=(9, 8) if debug else (nxsnip*4, 8))
    gs = gridspec.GridSpec(rows, nxsnip, figure=fig)
    # spectral vector for plotting the slits
    spec = np.tile(np.arange(slits.nspec), (slits.nslits, 1)).T
    thin = 10
    # legend elements
    legend_elements = [Line2D([0], [0], color='C3', lw=1, ls='dashed', dashes=(5,5), label='initial left edges'),
                       Line2D([0], [0], color='C1', lw=1, ls='dashed', dashes=(5,5), label='initial right edges'),
                       Line2D([0], [0], color='C3', lw=1, ls='solid', label='shifted left edges'),
                       Line2D([0], [0], color='C1', lw=1, ls='solid', label='shifted right edges')]
    # loop over the 2 rows if we save the plot in the output directory, otherwise plot the whole detector
    for r in range(rows):
        _ystar, _yend = (upper_ystart, upper_yend) if r == 0 else (lower_ystart, lower_yend)
        # loop over the snippets
        for s in range(nxsnip):
            ax = fig.add_subplot(gs[r, s])
            if vrange is None:
                # get vmin and vmax for imshow
                _xstart = spat_starts[s] if spat_starts[s] >= 0 else 0
                _xend = spat_ends[s] if spat_ends[s] <= img.shape[1] else img.shape[1]
                _img = img[_ystar:_yend, _xstart:_xend]
                _gpm = gpm[_ystar:_yend, _xstart:_xend] if gpm is not None else np.ones_like(_img, dtype=bool)
                m, med, sig = sigma_clipped_stats(_img[_gpm], sigma_lower=5.0, sigma_upper=5.0)
                vmin = m - 1.0 * sig
                vmax = m + 4.0 * sig
            else:
                vmin, vmax = vrange
            # imshow img instead of _img to show the actual pixel values in each snippet
            ax.imshow(img, origin='lower', vmin=vmin, vmax=vmax)
            ax.set_ylim(_ystar, _yend)
            ax.set_xlim(spat_starts[s], spat_ends[s])

            # plot the slits
            for i in range(slits.nslits):
                plt.plot(left_slits[::thin, i], spec[::thin, i], color='C3', lw=1, ls='dashed', dashes=(5,5), zorder=5)
                plt.plot(right_slits[::thin, i], spec[::thin, i], color='C1', lw=1, ls='dashed', dashes=(5,5), zorder=5)
                plt.plot(left_flex[::thin, i], spec[::thin, i], color='C3', lw=1, ls='solid', zorder=6)
                plt.plot(right_flex[::thin, i], spec[::thin, i], color='C1', lw=1, ls='solid', zorder=6)
            ax.tick_params(axis='both', labelsize=6)
            if r == 0 and s == 0:
                plt.suptitle(f'L/R spatial flexure = {spat_flexure[i,0]:.1f}/{spat_flexure[i,1]:.1f} pixels', fontsize=18)
                ax.legend(handles=legend_elements, fontsize=7)
                if not debug:
                    ax.set_ylabel('Upper snippets', fontsize=18)
            elif r == 1 and s == 0:
                ax.set_ylabel('Lower snippets', fontsize=18)
    plt.tight_layout()
    if debug:
        plt.show()
    else:
        fig.savefig(outfile, dpi=200)
        plt.close(fig)
