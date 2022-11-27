""" Figures for NSF Cyber Proposal """
from datetime import datetime
import os, sys
import numpy as np
import scipy
from scipy import stats
from urllib.parse import urlparse
import datetime
import geopy

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

# User Institutions

user_inst = [
    'Harvard University',
    'MPIA',
    'UC Berkeley',
    'Princeton University',
    'Cidade University, Brazil',
    'Ohio State University',
    'Vanderbilt University',
    'AIP, Germany',
    'Univeristy of Copenhagen',
    'Washington State University',
    'UC Santa Cruz',
    'University of Southhampton',
    'UC Los Angeles',
    'University of Illinois',
    'University of Melbourne',
    'University of Hong Kong',
    'University of Toronto',
    'Carnegie Observatories',
    'Northwestern University',
    'University of Dublin',
    'Johns Hopkins University',
    'STScI',
    'Las Cumbres Observatory',
    'UC Santa Barbara',
    'University of Hawaii',
    'Caltech',
    'Univeristy of Southern Queensland',
    'Univeristy of Texas',
    'Flatiron Institute',
    'Weizmann Institute',
    'Tel Aviv University',
    'Observatoire de Paris',
    'MIT',
    'Stockholm University',
    'University of Cambridge',
    'University of Maryland',
    'NASA Goddard',
    'University of Washington',
    'University of Portsmouth, UK',
    'Humboldt University, Berlin',
    'Universite Lyon',
    'Kunkoly Observatory, Hungary',
    'University of Arizona',
    'ESO Garching',
    'Gemini Observatory, AZ',
    'Leiden University',
    'University of Birmingham',
    'Technical Univeristy of Denmark',
    'Warsaw Univeristy',
    'Univeristy of Turku, Finland',
    'Radboud University',
    'SRON, Netherlands',
    'Stockholm University',
    'University of Edinburgh',
    'INAF Rome',
    'Queens University',
    'UC Davis',
    'UC Riverside',
    'York University',
    'Tufts University',
    'UC Irvine',
    'CSIC Madrid',
    'INAF Tieste',
    'INAF Napoli',
    'Universidad Andres Bello, Santiago',
    'Univeristy of Wisconsin',
    'INAF Padova',
    'San Jose State',
    'Univeristy of Waterloo',
    'Univeristy of Oulu',
    'Michigan State Univeristy',
    'Swinburne Univeristy',
    'RIT',
    'IAS, Princeton',
    'Queens Univeristy',
    'IAC, Canary Islands',
    'University of North Carolina',
    'Yale University',
    'CSIC Granada',
    'University of Manchester',
    'NAOJ',
    'Monash University',
    'Universidad de Zaragoza, Spain',
    'European University Cyprus, Spain',
    'NASA JPL',
    'Macquarie University',
    'UNIST, South Korea',
    'University of Firenze',
    'INFN Fiorentino',
    'University of Oslo',
    'INAF Bologna',
    'GEPI Meudon',
    'University of Pisa',
    'University of Minnesota',
    'LBNL',
    'Royal Observatory, Edinburgh',
    'UCL, London',
    'University of Tokyo',
    'University of Hertfordshire',
    'ASTRON',
    'MPIR Bonn',
    'CSIRO, Australia',
    'Curtin University',
    'SKA Observatory, UK',
    'University of Sydney',
    'Pontificia Universidad Catolica de Valparaiso',
    'University of Oxford',
    'University of Chicago',
    'INAF Naples',
    'CNRS Marseille',
    'Peking University',
    'Kyungpook National University, South Korea',
    'Pusan National University, South Korea',
    'Kyung Hee University, South Korea',
    'Korea Institute for Advanced Study',
    'SRI, Moscow',
    'MPA, Garching',
    'University of Michigan',
    'Karl Remeis-Observatory and Erlangen Centre for Astroparticle Physics',
    'South African Radio Astronomy Observatory',
    'IUCAA, Pune',
    'NRAO Socorro',
    'University of Geneva',
    'IAP, Paris',
    'Universidad de Chile',
    'University of the Western Cape',
    'University of New South Wales',
    'Arkansas Tech',
    'Japan Aerospace Exploration Agency',
    'Australian National University',
    'Maria Mitchell Observatory',
    'Universita degli Studi di Milano Bicocca',
    'Universita di Firenze',
    'INAF Florence',
    'San Diego State University',
    'W.M. Keck Observatory',
    ]

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
        # Magellan/CTIO
        if (np.abs(uni_lat+29.00333333333) < 1e-5): 
            # 
            idx = idx | (np.abs(all_lats + 30.240741666666672) < 1e-5)
            idx = idx | (np.abs(all_lats + 29.256666666666) < 1e-5)
        # Arziona
        if np.abs(uni_lat-31.6809444444) < 1e-5: 
            idx = idx | (np.abs(all_lats - 32.7015999999) < 1e-5) # LBT
            idx = idx | (np.abs(all_lats - 34.744305000) < 1e-5) # LDT
            idx = idx | (np.abs(all_lats - 31.963333333) < 1e-5)
        lon, lat = geo_dict[specs[idx][0]][0], geo_dict[specs[idx][0]][1], 
        # Plot
        plt.plot(lon, lat, 'o', transform=tformP)
        if (np.abs(uni_lat-32.) < 3.) & (
            np.abs(uni_lat-31.68094444) > 1e-5) & (
                np.abs(uni_lat-33.356000000) > 1e-5): # Arizona
            continue
        if np.abs(uni_lat+30.240741666666672) < 1e-5:
            continue
        if np.abs(uni_lat+29.256666666666) < 1e-5:
            continue
        if np.abs(uni_lat-19.82380144722) < 1e-5:
            continue
        # Label
        lbl = ''
        for spec in specs[idx]:
            lbl += spec + '\n'
        lbl = lbl[:-3]
        if 'vlt' in lbl or 'shane' in lbl:
            va = 'bottom'
        else:
            va = 'top'
        if 'p200' in lbl:
            ha = 'right'
        else:
            ha = 'left'
        #lbl = specs[idx][0]
        if np.abs(uni_lat-31.6809444) < 1e-5: 
            lon += 2.
        ax.text(lon, lat, lbl, transform=tformP,
              fontsize=15, ha=ha, va=va)

    
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


def fig_users(outfile:str='fig_geo_users.png', 
                 debug=False): 
    """ Global geographic plot of PypeIt users

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
        # Magellan/CTIO
        if (np.abs(uni_lat+29.00333333333) < 1e-5): 
            # 
            idx = idx | (np.abs(all_lats + 30.240741666666672) < 1e-5)
            idx = idx | (np.abs(all_lats + 29.256666666666) < 1e-5)
        # Arziona
        if np.abs(uni_lat-31.6809444444) < 1e-5: 
            idx = idx | (np.abs(all_lats - 32.7015999999) < 1e-5) # LBT
            idx = idx | (np.abs(all_lats - 34.744305000) < 1e-5) # LDT
            idx = idx | (np.abs(all_lats - 31.963333333) < 1e-5)
        lon, lat = geo_dict[specs[idx][0]][0], geo_dict[specs[idx][0]][1], 
        # Plot
        plt.plot(lon, lat, 'o', transform=tformP)
        if (np.abs(uni_lat-32.) < 3.) & (
            np.abs(uni_lat-31.68094444) > 1e-5) & (
                np.abs(uni_lat-33.356000000) > 1e-5): # Arizona
            continue
        if np.abs(uni_lat+30.240741666666672) < 1e-5:
            continue
        if np.abs(uni_lat+29.256666666666) < 1e-5:
            continue
        if np.abs(uni_lat-19.82380144722) < 1e-5:
            continue
        # Label
        lbl = ''
        for spec in specs[idx]:
            lbl += spec + '\n'
        lbl = lbl[:-3]
        if 'vlt' in lbl or 'shane' in lbl:
            va = 'bottom'
        else:
            va = 'top'
        if 'p200' in lbl:
            ha = 'right'
        else:
            ha = 'left'
        #lbl = specs[idx][0]
        if np.abs(uni_lat-31.6809444) < 1e-5: 
            lon += 2.
        ax.text(lon, lat, lbl, transform=tformP,
              fontsize=15, ha=ha, va=va)

    
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
