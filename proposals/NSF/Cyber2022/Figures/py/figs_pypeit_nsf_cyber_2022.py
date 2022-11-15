""" Figures for NSF Cyber Proposal """
from datetime import datetime
import os, sys
import numpy as np
import scipy
from scipy import stats
from urllib.parse import urlparse
import datetime

import argparse

import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import cartopy
from sympy import im

mpl.rcParams['font.family'] = 'stixgeneral'

import pandas
import seaborn as sns

import h5py

from pypeit.spectrographs import spectrograph_classes


from IPython import embed

def set_fontsize(ax, fsz):
    """
    Set the fontsize throughout an Axis

    Args:
        ax (Matplotlib Axis):
        fsz (float): Font size

    Returns:

    """
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)


def fig_spectrographs(outfile:str='fig_geo_spectrographs.png', 
                 debug=False): 
    """ Global geographic plot of PypeIt spectrographs

    Args:
        outfile (str): 
        debug (bool, optional): _description_. Defaults to False.
    """

    # Load up spectrographs 
    spectrographs = spectrograph_classes()

    # Grab the locations
    geo_dict = {}
    for key in spectrographs.keys():
        spectrograph = spectrographs[key]
        geo_dict[key] = spectrograph.telescope['longitude'], spectrograph.telescope['latitude']
        # Reset JWST
        if 'jwst' in key:
            geo_dict[key] = -76.615278, 39.289444
    all_lats = np.array([geo_dict[key][1] for key in geo_dict.keys()])
    specs = np.array(list(geo_dict.keys()))

    if debug:
        embed(header='70 of fig_geo_spectrographs')

    # Figure
    fig = plt.figure(figsize=(12,8))
    plt.clf()

    tformP = ccrs.PlateCarree()

    ax = plt.axes(projection=tformP)


    #cm = plt.get_cmap(color)
    uni_lats = np.unique(all_lats)

    # Plot em
    for uni_lat in uni_lats:
        idx = all_lats == uni_lat
        if np.abs(uni_lat-19.8283333333333) < 1e-5:
            # Gemini
            idx = idx | (np.abs(all_lats - 19.82380144722) < 1e-5)
        elif np.abs(uni_lat-19.82380144722) < 1e-5:
            continue
        # Magellan/CTIO
        if (np.abs(uni_lat+29.00333333333) < 1e-5) or (
            np.abs(uni_lat+29.256666666666) < 1e-5):
            # 
            idx = idx | (np.abs(all_lats + 30.240741666666672) < 1e-5)
            idx = idx | (np.abs(all_lats + 29.256666666666) < 1e-5)
        elif np.abs(uni_lat+30.240741666666672) < 1e-5:
            continue
        elif np.abs(uni_lat+29.256666666666) < 1e-5:
            continue
        lon, lat = geo_dict[specs[idx][0]][0], geo_dict[specs[idx][0]][1], 
        # Plot
        plt.plot(lon, lat, 'o', transform=tformP)
        # Label
        lbl = ''
        for spec in specs[idx]:
            lbl += spec + '\n'
        lbl = lbl[:-3]
        if 'vlt' in lbl or 'shane' in lbl:
            va = 'bottom'
        else:
            va = 'top'
        #lbl = specs[idx][0]
        ax.text(lon, lat, lbl, transform=tformP,
              fontsize=15, ha='left', va=va)

    
    # Zoom in
    ax.set_extent([-170, 10, -60, 60], 
                  crs=ccrs.PlateCarree())
    
    # Coast lines
    ax.coastlines(zorder=10)
    ax.add_feature(cartopy.feature.LAND, 
        facecolor='lightgray', edgecolor='black')

    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, 
        color='black', alpha=0.5, linestyle=':', draw_labels=True)
    #gl.xlabels_top = False
    #gl.ylabels_left = True
    #gl.ylabels_right=False
    #gl.xlines = True
    #gl.xformatter = LONGITUDE_FORMATTER
    #gl.yformatter = LATITUDE_FORMATTER
    #gl.xlabel_style = {'color': 'black'}# 'weight': 'bold'}
    #gl.ylabel_style = {'color': 'black'}# 'weight': 'bold'}

    set_fontsize(ax, 19.)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


def fig_slopevsDT(outfile='fig_slopevsDT.png', table=None,
                  local=False, vmax=None, 
                    cmap=None, cuts=None, scl = 1, debug=False):
    """ Bivariate of slope_min vs. DT

    Args:
        outfile (str, optional): [description]. Defaults to 'fig_slopevsDT.png'.
        local (bool, optional): [description]. Defaults to False.
        vmax ([type], optional): [description]. Defaults to None.
        cmap ([type], optional): [description]. Defaults to None.
        cuts ([type], optional): [description]. Defaults to None.
        scl (int, optional): [description]. Defaults to 1.
        debug (bool, optional): [description]. Defaults to False.
    """

    # Load table
    modis_tbl = ssl_paper_analy.load_modis_tbl(
        local=local, cuts=cuts, table=table)
    outfile = update_outfile(outfile, table)

    # Debug?
    if debug:
        modis_tbl = modis_tbl.loc[np.arange(1000000)].copy()

    # Plot
    fig = plt.figure(figsize=(12, 12))
    plt.clf()

    jg = sns.jointplot(data=modis_tbl, x='DT', y='min_slope', kind='hex',
                       bins='log', gridsize=250, xscale='log',
                       cmap=plt.get_cmap('winter'), mincnt=1,
                       marginal_kws=dict(fill=False, color='black', bins=100)) 
    jg.ax_joint.set_xlabel(r'$\Delta T$')
    jg.ax_joint.set_ylabel(metric_lbls['min_slope'])
    plt.colorbar()

    plotting.set_fontsize(jg.ax_joint, 15.)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


def fig_slopes(outfile='fig_slopes.png', 
               local=False, vmax=None, table=None,
               cmap=None, cuts=None, scl = 1, debug=False):

    # Load table
    modis_tbl = ssl_paper_analy.load_modis_tbl(
        local=local, cuts=cuts, table=table)
    outfile = update_outfile(outfile, table)

    # Debug?
    if debug:
        modis_tbl = modis_tbl.loc[np.arange(100000)].copy()

    # Check on isotropy
    diff = np.abs(modis_tbl.zonal_slope - modis_tbl.merid_slope)
    sig = np.sqrt(modis_tbl.zonal_slope_err**2 + modis_tbl.merid_slope**2)

    one_sig = diff < 1*sig
    frac = np.sum(one_sig) / len(diff)
    print(f"Fraction within 1 sigma = {frac}")

    # Plot
    fig = plt.figure(figsize=(12, 12))
    plt.clf()

    #ymnx = [-5000., 1000.]

    jg = sns.jointplot(data=modis_tbl, x='zonal_slope', y='merid_slope', 
                       kind='hex', #bins='log', xscale='log',
                       gridsize=100,
                       mincnt=1,
                       marginal_kws=dict(fill=False, 
                                         color='black', bins=100),
                       cmap=plt.get_cmap('YlGnBu')) 
                       #mincnt=1,
    plt.colorbar()
    
    jg.ax_joint.set_xlabel(metric_lbls['zonal_slope'])
    jg.ax_joint.set_ylabel(metric_lbls['merid_slope'])
    jg.ax_joint.plot([-5, 1.], [-5, 1.], 'k--')
    #jg.ax_joint.set_ylim(ymnx)

    plotting.set_fontsize(jg.ax_joint, 15.)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


def fig_2d_stats(outroot='fig_2dstats_', stat=None, table=None,
                local=False, vmax=None, nbins=40,
                cmap=None, cuts=None, scl = 1, debug=False):
    """ 2D histograms in the UMAP space

    Args:
        outroot (str, optional): [description]. Defaults to 'fig_2dstats_'.
        stat ([type], optional): [description]. Defaults to None.
        local (bool, optional): [description]. Defaults to False.
        vmax ([type], optional): [description]. Defaults to None.
        cmap ([type], optional): [description]. Defaults to None.
        cuts ([type], optional): [description]. Defaults to None.
        scl (int, optional): [description]. Defaults to 1.
        debug (bool, optional): [description]. Defaults to False.
    """

    # Load table
    modis_tbl = ssl_paper_analy.load_modis_tbl(local=local, cuts=cuts, table=table)

    # Debug?
    if debug:
        modis_tbl = modis_tbl.loc[np.arange(1000000)].copy()

    # Stat
    if stat is None:
        stat = 'min_slope'
    if cmap is None:
        cmap = 'hot'
    outfile = outroot+stat+'.png'

    # Decorate
    outfile = update_outfile(outfile, table)

    # Do it
    median_slope, x_edge, y_edge, ibins = scipy.stats.binned_statistic_2d(
        modis_tbl.U0, modis_tbl.U1, modis_tbl[stat],
        statistic='median', expand_binnumbers=True, bins=[nbins,nbins])

    # Plot
    fig = plt.figure(figsize=(12, 12))
    plt.clf()
    ax = plt.gca()


    cm = plt.get_cmap(cmap)
    mplt = ax.pcolormesh(x_edge, y_edge, 
                     median_slope.transpose(),
                     cmap=cm, 
                     vmax=None) 

    # Color bar
    cbaxes = plt.colorbar(mplt, pad=0., fraction=0.030)
    cbaxes.set_label(f'median({metric_lbls[stat]})', fontsize=17.)
    cbaxes.ax.tick_params(labelsize=15)

    ax.set_xlabel(r'$U_0$')
    ax.set_ylabel(r'$U_1$')

    plotting.set_fontsize(ax, 17.)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))


def fig_fit_metric(outroot='fig_fit_', metric=None, 
                   local=False, vmax=None, table=None,
                   distr='normal',
                   cmap=None, cuts=None, debug=False):

    # Load table
    modis_tbl = ssl_paper_analy.load_modis_tbl(local=local, cuts=cuts,
                                               table=table)

    # Debug?
    if debug:
        modis_tbl = modis_tbl.loc[np.arange(1000000)].copy()

    # Stat
    if metric is None:
        metric = 'DT'
    outfile = outroot+metric+'.png'
    # Decorate
    outfile = update_outfile(outfile, table)

    # Fit
    xmnx = modis_tbl[metric].min(), modis_tbl[metric].max()
    xval = np.linspace(xmnx[0], xmnx[1], 1000)
    dx = xval[1]-xval[0]
    if distr == 'normal':
        mean, sigma = stats.norm.fit(modis_tbl[metric])
        vals = stats.norm.pdf(xval, mean, sigma)
        print(f"Gaussian fit: mean={mean}, sigma={sigma}")
    elif distr == 'lognorm':
        shape,loc,scale = stats.lognorm.fit(modis_tbl[metric])
        vals = stats.lognorm.pdf(xval, shape, loc, scale)
        print(f"Log-norm fit: shape={shape}, loc={loc}, scale={scale}")
    else: 
        raise IOError(f"Bad distribution {distr}")

    # Normalize
    sum = dx * np.sum(vals)
    vals /= sum

    # Cumulative
    cumsum = np.cumsum(vals)
    cumsum /= cumsum[-1]
    
    # Plot
    fig = plt.figure(figsize=(10, 5))
    plt.clf()
    gs = gridspec.GridSpec(1,2)

    # Histogram
    ax_hist = plt.subplot(gs[0])

    _ = sns.histplot(modis_tbl, x=metric, ax=ax_hist,
                     stat='density')
    ax_hist.plot(xval, vals, 'k-')
    


    # CDF
    ax_cdf = plt.subplot(gs[1])
    _ = sns.ecdfplot(modis_tbl, x=metric, ax=ax_cdf)
    ax_cdf.plot(xval, cumsum, 'k--')

    for ax in [ax_hist, ax_cdf]:
        ax.set_xlabel(metric_lbls[metric])
        plotting.set_fontsize(ax, 17.)
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))

    # KS test
    #embed(header='778 of figs')
    #print(stats.kstest(modis_tbl[metric], distr))


def fig_learn_curve(outfile='fig_learn_curve.png'):
    # Grab the data
    #valid_losses_file = 's3://modis-l2/SSL/models/MODIS_R2019_96/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm/learning_curve/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm_losses_valid.h5'
    #valid_losses_file = 's3://modis-l2/SSL/models/MODIS_R2019_96/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm/learning_curve/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm_losses_valid.h5'
    valid_losses_file = 's3://modis-l2/SSL/models/MODIS_R2019_v4/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_256_temp_0.07_trial_5_cosine_warm/learning_curve/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_256_temp_0.07_trial_5_cosine_warm_losses_valid.h5'
    with ulmo_io.open(valid_losses_file, 'rb') as f:
        valid_hf = h5py.File(f, 'r')
    loss_avg_valid = valid_hf['loss_avg_valid'][:]
    loss_step_valid = valid_hf['loss_step_valid'][:]
    loss_valid = valid_hf['loss_valid'][:]
    valid_hf.close()

    #train_losses_file = 's3://modis-l2/SSL/models/MODIS_R2019_96/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm/learning_curve/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm_losses_train.h5'
    #train_losses_file = 's3://modis-l2/SSL/models/MODIS_R2019_96/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm/learning_curve/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_128_temp_0.07_trial_5_cosine_warm_losses_train.h5'
    train_losses_file = 's3://modis-l2/SSL/models/MODIS_R2019_v4/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_256_temp_0.07_trial_5_cosine_warm/learning_curve/SimCLR_resnet50_lr_0.05_decay_0.0001_bsz_256_temp_0.07_trial_5_cosine_warm_losses_train.h5'
    with ulmo_io.open(train_losses_file, 'rb') as f:
        train_hf = h5py.File(f, 'r')
    loss_train = train_hf['loss_train'][:]
    train_hf.close()

    # Plot
    fig = plt.figure(figsize=(10, 10))
    plt.clf()
    gs = gridspec.GridSpec(1,1)

    ax = plt.subplot(gs[0])

    ax.plot(loss_valid, label='valid', lw=3)
    ax.plot(loss_train, c='red', label='train', lw=3)

    ax.legend(fontsize=19.)

    # Label
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")

    plotting.set_fontsize(ax, 21.)
    
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))

def fig_DT_vs_U0(outfile='fig_DT_vs_U0.png',
                 local=None, table=None, nbins=40):
    # Grab the data
    modis_tbl = ssl_paper_analy.load_modis_tbl(local=local, table=table)

    median, x_edge, y_edge, ibins = scipy.stats.binned_statistic_2d(
        modis_tbl.U0, modis_tbl.U1, modis_tbl['DT'],
        statistic='median', expand_binnumbers=True, bins=[nbins,1])

    xvals = []
    for kk in range(len(x_edge)-1):
        xvals.append(np.mean(x_edge[kk:kk+2]))
        
    # Plot
    fig = plt.figure(figsize=(10, 10))
    plt.clf()
    gs = gridspec.GridSpec(1,1)

    ax = plt.subplot(gs[0])

    ax.plot(xvals, median.flatten(), 'o')

    #ax.legend(fontsize=15.)

    # Label
    ax.set_xlabel("U0")
    ax.set_ylabel("Median DT")

    plotting.set_fontsize(ax, 17.)
    
    plt.savefig(outfile, dpi=300)
    plt.close()
    print('Wrote {:s}'.format(outfile))

        
#### ########################## #########################
def main(pargs):

    # UMAP gallery
    if pargs.figure == 'geo_spec':
        fig_spectrographs(debug=pargs.debug)


def parse_option():
    """
    This is a function used to parse the arguments in the training.
    
    Returns:
        args: (dict) dictionary of the arguments.
    """
    parser = argparse.ArgumentParser("SSL Figures")
    parser.add_argument("figure", type=str, 
                        help="function to execute: 'slopes, 2d_stats, slopevsDT, umap_LL, learning_curve'")
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Debug?')
    args = parser.parse_args()
    
    return args

# Command line execution
if __name__ == '__main__':

    pargs = parse_option()
    main(pargs)

# Figures

# Geographic location of Spectrographs
# python py/figs_pypeit_nsf_cyber_2022.py geo_spec
