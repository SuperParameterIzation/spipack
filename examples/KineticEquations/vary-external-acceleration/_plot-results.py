import sys
import os, os.path

import numpy as np
from math import *

import random

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py
warnings.resetwarnings()

import matplotlib as mpl
#mpl.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams
from matplotlib import animation
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.backends.backend_pdf
import pylab

def MakeFigure(totalWidthPts, fraction, presentationVersion = False):
    fig_width_pt  = totalWidthPts * fraction

    inches_per_pt = 1.0 / 72.27
    golden_ratio  = (np.sqrt(5) - 1.0) / 2.0  # because it looks good

    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in * golden_ratio   # figure height in inches
    fig_dims      = [fig_width_in, fig_height_in] # fig dims as a list
    fig = plt.figure(figsize = fig_dims)

    greyColor = '#525252'
    whiteColor = '#ffffff'
    if not presentationVersion:
        rcParams['axes.labelsize'] = 9
        rcParams['xtick.labelsize'] = 9
        rcParams['ytick.labelsize'] = 9
        rcParams['legend.fontsize'] = 9
    else:
        rcParams['axes.labelsize'] = 12
        rcParams['xtick.labelsize'] = 12
        rcParams['ytick.labelsize'] = 12
        rcParams['legend.fontsize'] = 12
    rcParams['axes.edgecolor'] = greyColor
    rcParams['axes.facecolor'] = whiteColor
    rcParams['figure.facecolor'] = whiteColor
    rcParams['axes.labelcolor'] = greyColor
    rcParams['text.color'] = greyColor
    rcParams['xtick.color'] = greyColor
    rcParams['ytick.color'] = greyColor

    rcParams['lines.antialiased'] = True
    if not presentationVersion:
        rcParams['font.family'] = 'serif'
        #rcParams['font.serif'] = ['Computer Modern Roman']
        rcParams['text.usetex'] = True
#        rcParams['axes.formatter.use_mathtext'] = True
    else:
        rcParams['text.usetex'] = False
        rcParams['font.family'] = 'sans-serif'
#        rcParams['font.sans-serif'] = ['Helvetica']
#        rcParams['mathtext.fontset'] = 'stixsans'
#        rcParams['mathtext.fallback_to_cm'] = False
        rcParams['lines.linewidth'] = 1.5
#        rcParams['axes.formatter.use_mathtext'] = True

    #rcParams['backend'] = 'pdf'

    return fig

baseFilename = 'output'

outputDir = './output/'
filenames = [name for name in os.listdir(outputDir) if os.path.isfile(os.path.join(outputDir, name))]

time = list()
energy = list()
for filename in filenames:
    figname = 'figures/'+filename[:-3]+'.png'

    # load the data file
    hdf5file = h5py.File(outputDir+filename, 'r')

    time = time + hdf5file['/time'] [()] [0].tolist()
    samples = hdf5file['/samples'] [()].T
    acceleration = hdf5file['/acceleration'] [()]
    energy = energy + hdf5file['/expected energy'] [()] [0].tolist()
    sampleMean = hdf5file['/sample mean'] [()].T [0]
    vel = hdf5file['/expected velocity'] [()].T [0]
    acc = hdf5file['/expected acceleration'] [()].T [0]
    cov = hdf5file['/covariance'] [()]

    eigvals, eigvecs = np.linalg.eig(cov)

    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    ax.arrow(vel[0], vel[1], eigvals[0]*eigvecs[0, 0], eigvals[0]*eigvecs[1, 0])
    ax.arrow(vel[0], vel[1], eigvals[1]*eigvecs[0, 1], eigvals[1]*eigvecs[1, 1])
    ax.arrow(vel[0], vel[1], acc[0], acc[1], color='#08519c')
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, color='#525252')
    scale = 0.1 # scale the acceleration for plotting purposes
    #for samp, acc in zip(samples, acceleration):
    #    ax.arrow(samp[0], samp[1], scale*acc[0], scale*acc[1], color='#252525', alpha=0.8)
    ax.plot(sampleMean[0], sampleMean[1], 'o', color='#a50f15', markersize=3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    plt.savefig(figname, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)


sortIndex = np.argsort(time)
time = [time[i] for i in sortIndex]
energy = [energy[i] for i in sortIndex]

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
ax.plot(time, energy, color='#737373')
ax.set_xlabel(r'$\epsilon$')
ax.set_ylabel(r'$\Sigma_l^{\prime}$')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('figures/ExpectedEnergy.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)