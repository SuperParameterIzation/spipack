import sys

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

# load the data file
hdf5file = h5py.File('outputData.h5', 'r')

samples = hdf5file['/samples'] [()].T
rhs = hdf5file['/right hand side'] [()].T [0]
inverse = hdf5file['/weighted Poisson solution'] [()].T [0]
gradient = hdf5file['/weighted Poisson solution gradient'] [()]

pdf = matplotlib.backends.backend_pdf.PdfPages("figures/Kolmogorov_inverse.pdf")

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=rhs, vmin=min(rhs), vmax=max(rhs))
cbar = plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.ax.set_ylabel(r'Right hand side magnitude')
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=4, c=inverse, vmin=min(inverse), vmax=max(inverse))
cbar = plt.colorbar(scatter)
scale = 0.1 # scale the gradient for plotting purposes
for samp, grad in zip(samples, gradient):
    ax.arrow(samp[0], samp[1], scale*grad[0], scale*grad[1], color='#252525', alpha=0.8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.ax.set_ylabel(r'Weighted Poisson solution magnitude')
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)
pdf.close()
