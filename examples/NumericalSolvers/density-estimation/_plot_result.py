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
squaredBandwidth = hdf5file['/squared bandwidth'] [()].T [0]
trueDensity = hdf5file['/true density'] [()].T [0]
densityEstimate = hdf5file['/density estimate'] [()].T [0]

#densityEstimate = np.sqrt(2)*densityEstimate

print(densityEstimate-trueDensity)

bandwidthPara = hdf5file['/tune/candidate bandwidth parameters'] [()].T [0]
logKernelAvgDerivative = hdf5file['/tune/log kernel average derivative'] [()].T [0]

optInd = np.argmax(logKernelAvgDerivative)
print('key bandwidth', bandwidthPara[optInd])

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=squaredBandwidth, norm=mcolors.LogNorm(), vmin=min(squaredBandwidth), vmax=max(squaredBandwidth))
cbar = plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.ax.set_ylabel(r'Squared bandwidth $r_i^2$')
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
plt.savefig('figures/SquaredBandwidth.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=densityEstimate, vmin=0.0, vmax=0.16)
plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('figures/DensityEstimation.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)

#(\boldsymbol{x}^{(i)}; \boldsymbol{0}, \boldsymbol{I})

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=trueDensity, vmin=0.0, vmax=0.16)
cbar = plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.ax.set_ylabel(r'True density evaluation $\mathcal{N}(\mathbf{x}^{(i)}; \mathbf{0}, \mathbf{I})$')
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
plt.savefig('figures/TrueDensity.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=trueDensity-densityEstimate)
plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('figures/DensityDifference.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
#ax.plot([min(bandwidth), max(bandwidth)], [1, 1], '--', color='#737373')
ax.semilogx(bandwidthPara[:len(bandwidthPara)-1], logKernelAvgDerivative, color='#d73027')
#ax.set_xlim([min(bandwidth), max(bandwidth)])
ax.set_ylim([0, 1])
ax.set_xlabel(r'$\epsilon_l$')
ax.set_ylabel(r'$\Sigma_l^{\prime}$')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('figures/LogKernelAvgDerivative.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)
