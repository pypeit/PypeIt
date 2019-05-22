# Module to compare PYPIT vs. LowRedux
# v0.1 -- First look [Kast_red only]
# v0.2 -- With improved trace
# v0.3 -- With improved skysub and airtovac 9 Nov 2015
import numpy as np
import sys, os
import yaml
import pdb

from scipy import stats

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
mpl.rc('xtick', labelsize=18) 
mpl.rc('ytick', labelsize=18) 
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

import astropy.units as u
from astropy.io import fits
from astropy.table import Table

from linetools.spectra import io as lsio
from linetools.spectra.xspectrum1d import XSpectrum1D

sys.path.append('../src/')
import pypit
import arwave as arwv

try:
    from xastropy.xutils import xdebug as xdb
except:
    pass

def compare_s2n(pp,lrdx_sciobj,pypit_boxfile, iso):
    '''Compare boxcar S/N
    '''
    # Read/Load
    pypit_boxspec = lsio.readspec(pypit_boxfile)
    # Read LowRedux
    sig = np.sqrt(lrdx_sciobj['MASK_BOX']/(lrdx_sciobj['SIVAR_BOX'] + (lrdx_sciobj['MASK_BOX']==0)))
    lwrdx_boxspec = XSpectrum1D.from_tuple( (lrdx_sciobj['WAVE_BOX'], lrdx_sciobj['FLUX_BOX'], sig) )

    # Plot
    plt.clf()
    fig = plt.figure(figsize=(16,7))
    fig.suptitle("Instr={:s}, Setup={:s} :: Boxcar S/N for {:s} :: PYPIT ({:s})".format(iso[0], iso[1], iso[2], pypit.version), fontsize=18.)
    ax = plt.gca()
    ymax = np.median(pypit_boxspec.flux)*2.
    # PYPIT
    gdpy = pypit_boxspec.sig > 0.
    pys2n =  pypit_boxspec.flux[gdpy]/pypit_boxspec.sig[gdpy]
    ax.plot(pypit_boxspec.dispersion[gdpy],pys2n, 'k-', drawstyle='steps', label='PYPIT')
    # LowRedux
    gdlx = lwrdx_boxspec.sig > 0.
    ax.plot(lwrdx_boxspec.dispersion[gdlx], lwrdx_boxspec.flux[gdlx]/lwrdx_boxspec.sig[gdlx], 
            '-', color='blue', label='LowRedux')
    # Axes
    ax.set_xlim(np.min(pypit_boxspec.dispersion.value), np.max(pypit_boxspec.dispersion.value))
    ax.set_ylim(0.,np.median(pys2n)*2.)
    ax.set_xlabel('Wavelength',fontsize=17.)
    ax.set_ylabel('S/N per pixel',fontsize=17.)
    # Legend
    legend = plt.legend(loc='upper right', borderpad=0.3,
                handletextpad=0.3, fontsize='x-large')
    # Finish
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1,rect=[0, 0.03, 1, 0.95])
    pp.savefig(bbox_inches='tight')
    plt.close()
#

def compare_boxcar(pp,lrdx_sciobj,pypit_boxfile, iso):
    '''Compare boxcar extractions
    '''
    # Read/Load
    pypit_boxspec = lsio.readspec(pypit_boxfile)
    # Read LowRedux
    sig = np.sqrt(lrdx_sciobj['MASK_BOX']/(lrdx_sciobj['SIVAR_BOX'] + (lrdx_sciobj['MASK_BOX']==0)))
    lwrdx_boxspec = XSpectrum1D.from_tuple( (lrdx_sciobj['WAVE_BOX'], lrdx_sciobj['FLUX_BOX'], sig) )

    # Plot
    plt.clf()
    fig = plt.figure(figsize=(16,11))
    gs = gridspec.GridSpec(2, 1)
    fig.suptitle("Instr={:s}, Setup={:s} :: Boxcar Extractions for {:s} :: PYPIT ({:s})".format(iso[0], iso[1], iso[2], pypit.version), fontsize=18.)

    for qq in range(2):
        ax = plt.subplot(gs[qq])
        if qq == 0:
            xlim = None
        else:
            xlim = (6700,7000)
        ymax = np.median(pypit_boxspec.flux)*2.
        # PYPIT
        ax.plot(pypit_boxspec.dispersion, pypit_boxspec.flux, 'k-', drawstyle='steps',label='PYPIT')
        ax.plot(pypit_boxspec.dispersion, pypit_boxspec.sig, 'g-', drawstyle='steps')
        # LowRedux
        ax.plot(lwrdx_boxspec.dispersion, lwrdx_boxspec.flux, '-', color='blue',label='LowRedux')
        ax.plot(lwrdx_boxspec.dispersion, lwrdx_boxspec.sig, '-', color='gray')
        # Axes
        if xlim is None:
            ax.set_xlim(np.min(pypit_boxspec.dispersion.value), np.max(pypit_boxspec.dispersion.value))
        else:
            ax.set_xlim(xlim)
        ax.set_ylim(0.,ymax)
        ax.set_xlabel('Wavelength',fontsize=19.)
        ax.set_ylabel('electrons',fontsize=19.)
        # Legend
        legend = plt.legend(loc='upper right', borderpad=0.3,
                    handletextpad=0.3, fontsize='x-large')

    # Finish
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1,rect=[0, 0.03, 1, 0.95])
    pp.savefig(bbox_inches='tight')
    plt.close()
#

def compare_chi2hist(pp,lrdx_scihdu,pypit_skyfile, pypit_skysubfile, pypit_varfile, pypit_objtrcfile, iso):
    '''Compare chi^2 histograms for skysub regions
    '''
    # Load PYPIT
    pypit_sky = fits.open(pypit_skyfile)[0].data
    pypit_skysub = fits.open(pypit_skysubfile)[0].data
    pypit_var = fits.open(pypit_varfile)[0].data
    pypit_resid = pypit_skysub / np.sqrt(pypit_var)
    pypit_objtrc = fits.open(pypit_objtrcfile)[0].data
    pypit_objmsk = fits.open(pypit_objtrcfile)[1].data.astype(int)
    # LowRedux
    lwrdx_proc = lrdx_scihdu[0].data
    lwrdx_ivar = lrdx_scihdu[1].data
    lwrdx_sky = lrdx_scihdu[2].data
    lwrdx_resid = (lwrdx_proc-lwrdx_sky) * np.sqrt(lwrdx_ivar)

    # Select regions
    # Median size of box car
    dx = np.median(np.sum(pypit_objmsk[:,:,0],axis=1))
    #
    pypit_skymask = np.zeros_like(pypit_sky)
    # Generate sky regions near the object trace
    for ii in xrange(pypit_sky.shape[0]):
        # Left
        idx = np.arange(int(pypit_objtrc[ii,0]-dx*3), int(pypit_objtrc[ii,0]-dx))
        pypit_skymask[ii,idx] = 1
        # Right
        idx = np.arange(int(pypit_objtrc[ii,0]+dx), int(pypit_objtrc[ii,0]+dx*3))
        pypit_skymask[ii,idx] = 1
    # 
    skypix = np.where(pypit_skymask==1)

    # Residuals
    # PYPIT
    pypit_chi = pypit_resid[skypix]
    # LowRedux
    lwrdx_chi = lwrdx_resid[skypix]

    # Histograms
    # Boundaries
    minv, maxv = -7., 7.
    binsz = 0.1

    # Set the boundaries sensibly given binsz
    i0 = int( minv / binsz) - 1
    i1 = int( maxv / binsz) + 1
    rng = tuple( binsz*np.array([i0,i1]) )
    nbin = i1-i0

    # Plot
    plt.clf()
    fig = plt.figure(figsize=(16,7))
    fig.suptitle("Instr={:s}, Setup={:s} :: Chi_Resid Histograms for {:s} :: PYPIT ({:s})".format(iso[0], iso[1], iso[2], pypit.version), fontsize=18.)
    gs = gridspec.GridSpec(1, 2)

    # PYPIT
    axp = plt.subplot(gs[0])
    # Histogram
    hist, edges = np.histogram(pypit_chi, range=rng, bins=nbin)
    axp.bar(edges[:-1], hist, width=binsz, color='black')#, alpha=kwargs['alpha'])
    # PDF for Gaussian
    area = pypit_chi.size * binsz
    xppf = np.linspace(stats.norm.ppf(0.0001), stats.norm.ppf(0.9999), 100)
    yppf = area*stats.norm.pdf(xppf)
    axp.plot(xppf, yppf, 'r-', alpha=1.0)
    # Median
    axp.plot([np.median(pypit_chi)]*2, [-9e9,9e9], 'g--')
    #
    axp.set_xlim(minv,maxv)  
    axp.set_ylim(0., np.max(yppf)*1.1)  
    axp.set_xlabel(r'$\chi$ (PYPIT)', fontsize=17)

    # LowRedux
    axl = plt.subplot(gs[1])
    # Histogram
    hist, edges = np.histogram(lwrdx_chi, range=rng, bins=nbin)
    axl.bar(edges[:-1], hist, width=binsz)#, alpha=kwargs['alpha'])
    # PDF for Gaussian
    area = lwrdx_chi.size * binsz
    yppf = area*stats.norm.pdf(xppf)
    axl.plot(xppf, yppf, 'r-', alpha=1.0)
    # Median
    axl.plot([np.median(lwrdx_chi)]*2, [-9e9,9e9], 'g--')
    #
    axl.set_xlim(minv,maxv)  
    axl.set_ylim(0., np.max(yppf)*1.1)  
    axl.set_xlabel(r'$\chi$ (LowRedux)', fontsize=17)

    # Finish
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1,rect=[0, 0.03, 1, 0.95])
    pp.savefig(bbox_inches='tight')
    plt.close()

def compare_chi2img(pp,lrdx_scihdu,pypit_skysubfile, pypit_varfile, iso):
    '''Compare chi^2 images for skysub
    '''
    # Load PYPIT
    pypit_skysub = fits.open(pypit_skysubfile)[0].data
    pypit_var = fits.open(pypit_varfile)[0].data
    pypit_resid = pypit_skysub / np.sqrt(pypit_var)
    # LowRedux
    lwrdx_proc = lrdx_scihdu[0].data
    lwrdx_ivar = lrdx_scihdu[1].data
    lwrdx_sky = lrdx_scihdu[2].data
    lwrdx_resid = (lwrdx_proc-lwrdx_sky) * np.sqrt(lwrdx_ivar)

    # Plot
    vmnx = (-3.,3.)
    plt.clf()
    fig = plt.figure(figsize=(16,8))
    fig.suptitle("Instr={:s}, Setup={:s} :: Chi_Resid Images for {:s} :: PYPIT ({:s})".format(iso[0], iso[1], iso[2], pypit.version), fontsize=18.)
    cm = plt.get_cmap('Greys') 
    wbox= {'facecolor':'white', 'edgecolor':'white'}
    gs = gridspec.GridSpec(2, 20)

    # PYPIT
    ax = plt.subplot(gs[0,:-1])
    mplt = ax.imshow(pypit_resid.T, origin='lower', cmap=cm)
    mplt.set_clim(vmin=vmnx[0], vmax=vmnx[1])
    ax.text(0.10, 0.80, 'PYPIT', transform=ax.transAxes, fontsize=19, ha='left',bbox=wbox, color='black')

    # LowRedux
    ax = plt.subplot(gs[1,:-1])
    mplt2 = ax.imshow(lwrdx_resid.T, origin='lower', cmap=cm)
    mplt2.set_clim(vmin=vmnx[0], vmax=vmnx[1])
    ax.text(0.10, 0.80, 'LowRedux', transform=ax.transAxes, fontsize=19, ha='left',bbox=wbox, color='blue')

    # Colorbar
    cbar_ax = plt.subplot(gs[:,-1])
    cb = plt.colorbar(mplt2, cax=cbar_ax)#,orientation='horizontal',)
    cb.set_label(r'$\chi$',fontsize=18.)
    # Finish
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig(bbox_inches='tight')
    plt.close()

def compare_skyspec(pp,sciobj,pypit_skyboxfile, iso):
    '''Compare traces
    Parameters:
    ----------
    pp: Pdf obj
    sciobj: Table
      LowRedux sciobj table 
    iso: tuple
      instr, setup, obj strings
    '''
    # Load
    pypit_skyspec = lsio.readspec(pypit_skyboxfile) 
    #
    lwrdx_skywv = sciobj['WAVE_BOX']
    lwrdx_skyfx = sciobj['SKY_BOX']
    lwrdx_skywv = arwv.vactoair(lwrdx_skywv*u.AA)
    #
    plt.clf()
    plt.figure(figsize=(16,8))
    gs = gridspec.GridSpec(2, 1)
    # Full
    ax = plt.subplot(gs[0])
    # PYPIT
    ax.plot(pypit_skyspec.dispersion, pypit_skyspec.flux, 'k-', label='PYPIT', drawstyle='steps')
    # LowRedux
    ax.plot(lwrdx_skywv, lwrdx_skyfx/(2.*sciobj['BOX_RAD']), '-', color='blue', label='LowRedux')
    # Axes
    ax.set_xlim(np.min(pypit_skyspec.dispersion.value),np.max(pypit_skyspec.dispersion.value))
    ax.set_ylim(0.,np.max(pypit_skyspec.flux))
    ax.set_xlabel('Wavelength',fontsize=17.)
    ax.set_ylabel('Sky (Counts/pix)',fontsize=17.)
    # Legend
    legend = plt.legend(loc='upper left', borderpad=0.3,
                handletextpad=0.3, fontsize='x-large')

    # ZOOM
    axz = plt.subplot(gs[1])
    # PYPIT
    axz.plot(pypit_skyspec.dispersion, pypit_skyspec.flux, 'k-', label='PYPIT', drawstyle='steps-mid')
    # LowRedux
    axz.plot(lwrdx_skywv, lwrdx_skyfx/(2.*sciobj['BOX_RAD']), '-', color='blue', label='LowRedux')
    # Axes
    zlim = np.array([7200., 7700.])*u.AA
    axz.set_xlim(zlim.value)
    ymx = np.max( pypit_skyspec.flux[np.where((pypit_skyspec.dispersion>zlim[0])&(pypit_skyspec.dispersion<zlim[1]))])
    axz.set_ylim(0.,ymx)
    axz.set_xlabel('Wavelength',fontsize=17.)
    axz.set_ylabel('Sky (Counts/pix)',fontsize=17.)
    # Legend
    legend = plt.legend(loc='upper right', borderpad=0.3,
                handletextpad=0.3, fontsize='x-large')
    ax.set_title("Instr={:s}, Setup={:s} :: Sky Spectra for {:s} :: PYPIT ({:s})".format(iso[0], iso[1], iso[2], pypit.version), fontsize=18.)

    # Finish
    plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
    pp.savefig(bbox_inches='tight')
    plt.close()
#

def compare_traces(pp,sciobj,pypit_objtrcfile, iso):
    '''Compare traces
    Parameters:
    ----------
    pp: Pdf obj
    sciobj: Table
      LowRedux sciobj table 
    iso: tuple
      instr, setup, obj strings
    '''
    # Read PYPIT
    pypit_objtrc = fits.open(pypit_objtrcfile)[0].data
    if pypit_objtrc.shape[1] > 1:
        raise ValueError('Not ready for multiple objects')
    # Read LowRedux
    lwrdx_objtrc = sciobj['XPOS']

    plt.clf()
    plt.figure(figsize=(16,8))
    ax = plt.gca()
    # PYPIT
    ax.plot(pypit_objtrc[:,0], 'k-', drawstyle='steps', label='PYPIT')
    # LowRedux
    ax.plot(lwrdx_objtrc, '-', color='blue', label='LowRedux')
    # Axes
    #ax.set_ylim(0.,np.median(pys2n)*2.)
    ax.set_xlabel('Row',fontsize=17.)
    ax.set_ylabel('Column',fontsize=17.)
    ax.set_title("Instr={:s}, Setup={:s} :: Object Traces for {:s} :: PYPIT ({:s})".format(iso[0], iso[1], iso[2], pypit.version), fontsize=18.)
    # Legend
    legend = plt.legend(loc='lower right', borderpad=0.3,
                handletextpad=0.3, fontsize='x-large')
    # Finish
    pp.savefig(bbox_inches='tight')
    plt.close()

def main():
    '''Loop through TEST_SUITES and perform comparison where applicable
    '''
    pypit_roots = dict(box='boxcar',sky='sky',skysub='skysub',var='var',skybox='skybox',objtrc='objtrc')
    # Point to all sub-folders
    raw_path = os.getenv('DROPBOX_DIR')+'PYPIT/TEST_SUITES/'
    walk = os.walk(raw_path) 
    instruments = next(walk)[1]
    # Loop on instruments
    for instr in instruments:
        # Setups
        setups = next(os.walk(raw_path+instr))[1]
        for setup in setups:
            wdir = os.path.join(os.getenv('TST_PYPIT'),instr,setup)
            # Look for LowRedux and MasterFrame folders
            low_rdx_path = os.path.join(raw_path,instr,setup,'LowRedux')
            yml_fil = low_rdx_path+'/objects.yaml'
            pypit_path = os.path.join(os.getenv('TST_PYPIT'),instr,setup,'MasterFrames') 
            if (os.path.exists(low_rdx_path) & os.path.exists(pypit_path) & os.path.isfile(yml_fil)):
                with open(yml_fil, 'r') as infile:
                    lrdx_dict = yaml.load(infile)
            else:
                print('No LowRedux for instr={:s}, setup={:s}'.format(instr,setup))
                continue
            # Loop on sources
            for obj in lrdx_dict.keys():
                iso = (instr, setup, obj)
                # Set LowRedux files
                lrdx_scifile = os.path.join(low_rdx_path,lrdx_dict[obj]['sci_file'])
                lrdx_scihdu = fits.open(lrdx_scifile)
                lrdx_sciobj =  Table(lrdx_scihdu[5].data)[0]
                lrdx_wvfil = lrdx_dict[obj]['wave_file']
                # Setup for PYPIT files
                pypit_prefix = pypit_path+'/'+obj+'_'
                pypit_wvimg_fil = pypit_path+'/mswvimg_red_000.fits'
                # Outfil
                outfil = 'PYPIT_vs_LRDX_'+instr+'_'+setup+'.pdf'
                pp = PdfPages(outfil)
                # Trace test
                compare_traces(pp,lrdx_sciobj,pypit_prefix+pypit_roots['objtrc']+'.fits', iso)
                # Sky
                compare_skyspec(pp,lrdx_sciobj,pypit_prefix+pypit_roots['skybox']+'.fits', iso)
                compare_chi2img(pp,lrdx_scihdu,pypit_prefix+pypit_roots['skysub']+'.fits', pypit_prefix+pypit_roots['var']+'.fits', iso)
                compare_chi2hist(pp,lrdx_scihdu,pypit_prefix+pypit_roots['sky']+'.fits', pypit_prefix+pypit_roots['skysub']+'.fits', pypit_prefix+pypit_roots['var']+'.fits', pypit_prefix+pypit_roots['objtrc']+'.fits', iso)
                # Spectra
                compare_boxcar(pp,lrdx_sciobj,pypit_prefix+pypit_roots['box']+'.fits', iso)
                compare_s2n(pp,lrdx_sciobj,pypit_prefix+pypit_roots['box']+'.fits', iso)
                # Finish
                pp.close()


# ################
if __name__ == "__main__":
    main()
