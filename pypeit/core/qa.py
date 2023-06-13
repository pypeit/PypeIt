""" Module for QA in PypeIt
"""
import os
import datetime
import getpass
import glob
import numpy as np
import yaml

from IPython import embed

# CANNOT INCLUDE msgs IN THIS MODULE AS
#  THE HTML GENERATION OCCURS FROM msgs
#from pypeit import msgs

# TODO: Move these names to the appropriate class.  This always writes
# to QA directory, even if the user sets something else...
def set_qa_filename(root, method, det=None, slit=None, prefix=None, out_dir=None):
    """
    Generate the filename for the QA file from the input parameters.
    
    Parameters
    ----------
    root : str
      Root name
    method : str
      Describes the QA routine
    det : str, optional
        The name of the detector or mosaic (e.g., DET01)
    slit : int, optional
    prefix : str, optional
      start the name of the QA file
    out_dir : str, optional
      Path to QA/

    Returns
    -------
    outfile : str
      Filename
    """
    if out_dir is None:
        out_dir = os.getcwd()
    #
    if method == 'slit_trace_qa':
#        outfile = 'QA/PNGs/Slit_Trace_{:s}.png'.format(root)
        outfile = 'PNGs/Slit_Trace_{:s}.png'.format(root)
    elif method == 'slit_profile_qa':
        outfile = 'QA/PNGs/Slit_Profile_{:s}_'.format(root)
    elif method == 'arc_fit_qa':
#        outfile = 'QA/PNGs/Arc_1dfit_{:s}_S{:04d}.png'.format(root, slit)
        outfile = 'PNGs/Arc_1dfit_{:s}_S{:04d}.png'.format(root, slit)
    elif method == 'plot_orderfits_Arc':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/Arc_lines_{:s}_S{:04d}_'.format(root, slit)
    elif method == 'arc_fit2d_global_qa':
#        outfile = 'QA/PNGs/Arc_2dfit_global_{:s}'.format(root)
        outfile = 'PNGs/Arc_2dfit_global_{:s}'.format(root)
    elif method == 'arc_fit2d_orders_qa':
#        outfile = 'QA/PNGs/Arc_2dfit_orders_{:s}'.format(root)
        outfile = 'PNGs/Arc_2dfit_orders_{:s}'.format(root)
    elif method == 'arc_tilts_spec_qa':
#        outfile = 'QA/PNGs/Arc_tilts_spec_{:s}_S{:04d}.png'.format(root, slit)
        outfile = 'PNGs/Arc_tilts_spec_{:s}_S{:04d}.png'.format(root, slit)
    elif method == 'arc_tilts_spat_qa':
#        outfile = 'QA/PNGs/Arc_tilts_spat_{:s}_S{:04d}.png'.format(root, slit)
        outfile = 'PNGs/Arc_tilts_spat_{:s}_S{:04d}.png'.format(root, slit)
    elif method == 'arc_tilts_2d_qa':
#        outfile = 'QA/PNGs/Arc_tilts_2d_{:s}_S{:04d}.png'.format(root, slit)
        outfile = 'PNGs/Arc_tilts_2d_{:s}_S{:04d}.png'.format(root, slit)
    elif method == 'pca_plot':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/{:s}_pca_{:s}_'.format(prefix, root)
    elif method == 'pca_arctilt':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/Arc_pca_{:s}_'.format(root)
    elif method == 'plot_orderfits_Blaze':  # This is root for multiple PNGs
        outfile = 'QA/PNGs/Blaze_{:s}_'.format(root)
    elif method == 'obj_trace_qa':
        outfile = 'QA/PNGs/{:s}_{:s}_obj_trace.png'.format(root, det)
    elif method == 'obj_profile_qa':
        outfile = 'PNGs/{:s}_{:s}_S{:04d}_obj_prof.png'.format(root, det, slit)
    elif method == 'spec_flexure_qa_corr':
#        outfile = 'QA/PNGs/{:s}_D{:02d}_S{:04d}_spec_flex_corr.png'.format(root, det, slit)
        outfile = 'PNGs/{:s}_S{:04d}_spec_flex_corr.png'.format(root, slit)
    elif method == 'spec_flexure_qa_sky':
#        outfile = 'QA/PNGs/{:s}_D{:02d}_S{:04d}_spec_flex_sky.png'.format(root, det, slit)
        outfile = 'PNGs/{:s}_S{:04d}_spec_flex_sky.png'.format(root, slit)
    elif method == 'spatillum_finecorr':
        outfile = 'PNGs/{:s}_S{:04d}_spatillum_finecorr.png'.format(root, slit)
    elif method == 'detector_structure':
        outfile = 'PNGs/{:s}_{:s}_detector_structure.png'.format(root, det)
    else:
        raise IOError("NOT READY FOR THIS QA: {:s}".format(method))
    # Return
    return os.path.join(out_dir, outfile)


def get_dimen(x, maxp=25):
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


def gen_timestamp():
    """ Generate a simple time stamp including the current user

    Returns
    -------
    timestamp : str
      user_datetime
    """
    tstamp = datetime.datetime.today().strftime('%Y-%m-%d-T%Hh%Mm%Ss')
    try:
        user = getpass.getuser()
    except ModuleNotFoundError:
        # there appears to be a bug in getpass in windows systems where the pwd module doesn't load
        user = os.getlogin()
    # Return
    return '{:s}_{:s}'.format(user, tstamp)


def html_header(title):
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

def html_end(f, body, links=None):
    """ Fill in the HTML file with a proper ending

    Parameters
    ----------
    f : file
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


def html_init(f, title):
    """
    Initialize the HTML file

    Args:
        f (fileobj): file object to write to
        title (str): title

    Returns:
        str: Initial HTML text incluing the header and links
    """
    head = html_header(title)
    f.write(head)
    # Init links
    links = '<h2>Quick Links</h2>\n'
    links += '<ul>\n'
    return links


def html_mf_pngs(idval):
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
    # QA root
    # Search for PNGs

    # Organize the outputs
    html_dict = {}
    html_dict['strace'] = dict(fname='slit_trace_qa', ext='',
                               href='strace', label='Slit Trace', slit=False)
    html_dict['sprof'] = dict(fname='slit_profile_qa', ext='*.png',
                              href='sprof', label='Slit Profile', slit=False)
    html_dict['blaze'] = dict(fname='plot_orderfits_Blaze', ext='*.png',
                              href='blaze', label='Blaze', slit=False)
    html_dict['arc_fit'] = dict(fname='arc_fit_qa', ext='',
                                href='arc_fit', label='Arc 1D Fit', slit=True)
    html_dict['arc_tilts_spec'] = dict(fname='arc_tilts_spec_qa', ext='',
                                       href='arc_tilts_spec', label='Arc Tilts Spec', slit=True)
    html_dict['arc_tilts_spat'] = dict(fname='arc_tilts_spat_qa', ext='',
                                 href='arc_tilts_spat', label='Arc Tilts Spat', slit=True)
    html_dict['arc_tilts_2d'] = dict(fname='arc_tilts_2d_qa', ext='',
                                       href='arc_tilts_2d', label='Arc Tilts 2D', slit=True)
    html_dict['arc_pca'] = dict(fname='pca_arctilt', ext='*.png',
                                href='arc_pca', label='Arc Tilt PCA', slit=False)
    html_dict['arc_fit2d_global'] = dict(fname='arc_fit2d_global_qa', ext='*.png',
                                         href='arc_fit2d_global', label='2D Arc Global', slit=False)
    html_dict['arc_fit2d_orders'] = dict(fname='arc_fit2d_orders_qa', ext='*.png',
                                         href='arc_fit2d_orders', label='2D Arc Orders', slit=False)

    # Generate HTML
    for key in ['strace', 'sprof', 'blaze', 'arc_fit', 'arc_pca', 'arc_fit2d_global', 'arc_fit2d_orders',
                'arc_tilts_spec', 'arc_tilts_spat', 'arc_tilts_2d']:
        # PNG Root
        png_root = set_qa_filename(idval, html_dict[key]['fname'], slit=9999, out_dir='QA')
        if html_dict[key]['slit']:  # Kludge to handle multiple slits
            png_root = png_root.replace('S9999', 'S*')
        # Find the PNGs
        pngs = glob.glob(png_root+html_dict[key]['ext'])
        pngs.sort()
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
                ifnd = png.find('QA' + os.path.sep)
                if ifnd < 0:
                    raise ValueError("QA is expected to be in the path!")
                if html_dict[key]['slit']:  # Kludge to handle multiple slits
                    i0 = png.find('{:s}_S'.format(idval))
                    href="{:s}_{:s}".format(html_dict[key]['href'], png[i0:])
                    body += '<img class ="research" src="{:s}" width="100%" id={:s} height="auto"/>\n'.format(
                        png[ifnd+3:], href)
                    links += '<li><a class="reference internal" href="#{:s}">{:s} {:s}</a></li>\n'.format(
                        href, html_dict[key]['label'], png[i0:-4])
                else:
                    body += '<img class ="research" src="{:s}" width="100%" height="auto"/>\n'.format(png[ifnd+3:])
            body += '</div>\n'

    # Return
    return links, body


def html_exp_pngs(exp_name, det):
    """
    Generate HTML for Exposure PNGs

    Parameters
    ----------
    exp_name : str
    det : int

    Returns
    -------
    links : str
    body : str

    """
    links = ''
    body = ''
    # QA root

    # Organize the outputs
    html_dict = {}
    html_dict['trace'] = dict(fname='obj_trace_qa', ext='', slit=False,
                               href='otrace', label='Object Traces')
    html_dict['prof'] = dict(fname='obj_profile_qa', ext='', slit=True,
                              href='oprofile', label='Object Profiles')
    html_dict['flex_corr'] = dict(fname='flexure_qa_corr', ext='', slit=True,
                             href='flex_corr', label='Flexure Cross Correlation')
    html_dict['flex_sky'] = dict(fname='flexure_qa_sky', ext='', slit=True,
                                  href='flex_sky', label='Flexure Sky')

    # Generate HTML
    for key in ['trace', 'prof', 'flex_corr', 'flex_sky']:
        png_root = set_qa_filename(exp_name, html_dict[key]['fname'], det=det, slit=9999)
        if html_dict[key]['slit']:  # Kludge to handle multiple slits
            png_root = png_root.replace('S9999', 'S*')
        pngs = glob.glob(png_root+html_dict[key]['ext'])
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
                ifnd = png.find('QA/')
                if ifnd < 0:
                    raise ValueError("QA is expected to be in the path!")
                body += '<img class ="research" src="{:s}" width="100%" height="auto"/>\n'.format(png[ifnd+3:])
            body += '</div>\n'

    # Return
    return links, body


def gen_qa_dir(qa_path):
    """ Make the QA directory if it doesn't already exist

    Args:
        qa_path (str):
            Path to the QA folder
    """
    if not os.path.exists(qa_path):
        os.makedirs(qa_path)

# TODO: Need to revisit this...
def gen_mf_html(pypeit_file, qa_path):
    """ Generate the HTML for QA

    Args:
        pypeit_file (str):
            Name of the PypeIt file, no path
        qa_path (str):
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
    MF_filename = os.path.join('{:s}'.format(qa_path), 'MF_{:s}.html'.format(setup))
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
    print("Wrote: {:s}".format(MF_filename))


def gen_exp_html():
    """ Generate the HTML for an Exposure set
    """
    # Find all obj_trace files -- Not fool proof but ok
    obj_files = glob.glob('QA/PNGs/*obj_trace.png')
    # Parse for names
    names = []
    for obj_file in obj_files:
        i0 = obj_file.rfind('/')+1
        i1 = obj_file.rfind('_D')
        name = obj_file[i0:i1]
        names.append(name)
    uni_names = np.unique(names)
    # Loop
    for uni_name in uni_names:
        # Generate MF file
        exp_filename = 'QA/{:s}.html'.format(uni_name)
        body = ''
        with open(exp_filename,'w') as f:
            # Start
            links = html_init(f, 'QA for {:s}'.format(uni_name))
            # Loop on detector
            for det in range(1,99):
                # Run
                new_links, new_body = html_exp_pngs(uni_name, det)
                # Save
                links += new_links
                body += new_body
            # End
            html_end(f, body, links)
        print("Wrote: {:s}".format(exp_filename))


def close_qa(pypeit_file, qa_path):
    """Tie off QA under a crash

    Args:
        pypeit_file (_type_): _description_
        qa_path (_type_): _description_
    """
    if pypeit_file is None:
        return
    try:
        gen_mf_html(pypeit_file, qa_path)
    except:  # Likely crashed real early
        pass
    else:
        gen_exp_html()


