import sys

import numpy as np

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
hdf5file = h5py.File('samples.h5', 'r')

samples = hdf5file['/samples'] [()].T
logBandwidth = hdf5file['/log bandwidth'] [()].T [0]
bandwidth = hdf5file['/bandwidth parameter candidate'] [()].T [0]
sigmaprime = hdf5file['/sigma prime'] [()].T [0]
#eigenvalues = hdf5file['/heat matrix eigenvalues'] [()].T [0]

sigmaprimeMaxInd = np.argmax(sigmaprime)
print('key bandwidth', bandwidth[sigmaprimeMaxInd])

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=logBandwidth, vmin=min(logBandwidth), vmax=max(logBandwidth))
plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('figures/Samples.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
ax.plot([min(bandwidth), max(bandwidth)], [1, 1], '--', color='#737373')
ax.semilogx(bandwidth, sigmaprime)
ax.plot(bandwidth[sigmaprimeMaxInd], sigmaprime[sigmaprimeMaxInd], 'o', markersize=3, markeredgecolor='#a50f15', markerfacecolor='#a50f15')
ax.set_xlim([min(bandwidth), max(bandwidth)])
ax.set_ylim([0, 1])
ax.set_xlabel(r'$\epsilon_l$')
ax.set_ylabel(r'$\Sigma_l^{\prime}$')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('figures/SigmaPrime.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)

# fig = MakeFigure(425, 0.9, False)
# ax = plt.gca()
# ax.plot(eigenvalues, color='#525252')
# ax.set_xlim([0, len(eigenvalues)])
# ax.set_ylim([min(eigenvalues), max(eigenvalues)])
# ax.set_xlabel(r'Index')
# ax.set_ylabel(r'Eigenvalue')
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# ax.yaxis.set_ticks_position('left')
# ax.xaxis.set_ticks_position('bottom')
# plt.savefig('figures/HeatMatrixEigenvalues.pdf', format='pdf', bbox_inches='tight')
# plt.close(fig)
