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
        rcParams['axes.labelsize'] = 18
        rcParams['xtick.labelsize'] = 14
        rcParams['ytick.labelsize'] = 14
        rcParams['legend.fontsize'] = 14
    else:
        rcParams['axes.labelsize'] = 18
        rcParams['xtick.labelsize'] = 14
        rcParams['ytick.labelsize'] = 14
        rcParams['legend.fontsize'] = 14
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
rhs = hdf5file['/true right hand side'] [()].T [0]
rhsExpanded = hdf5file['/expanded right hand side'] [()].T [0]
inverse = hdf5file['/weighted Poisson solution'] [()].T [0]
gradient = hdf5file['/weighted Poisson solution gradient'] [()]


eigvecsLhat = hdf5file['/eigenvectors Qhat'] [()].T

pdf = matplotlib.backends.backend_pdf.PdfPages("figures/Kolmogorov_inverse.pdf")

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=rhs, vmin=min(rhs), vmax=max(rhs), cmap='jet')
cbar = plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.set_label(label=r'Function magnitude', size=16)
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
plt.savefig('figures/fig_WeightedLaplace_TrueRHS.pdf', format='pdf', bbox_inches='tight')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=rhsExpanded, vmin=min(rhs), vmax=max(rhs), cmap='jet')
cbar = plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.set_label(label=r'Function magnitude', size=16)
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
plt.savefig('figures/fig_WeightedLaplace_ExpandedRHS.pdf', format='pdf', bbox_inches='tight')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=rhs-rhsExpanded, vmin=-max(abs(min(rhs-rhsExpanded)), max(rhs-rhsExpanded)), vmax=max(abs(min(rhs-rhsExpanded)), max(rhs-rhsExpanded)), cmap='jet')
cbar = plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.set_label(label=r'Difference', size=16)
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
plt.savefig('figures/fig_WeightedLaplace_RHSDifference.pdf', format='pdf', bbox_inches='tight')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)

trueSolution_c0 = [0.0]*len(samples)
for i in range(len(samples)):
    trueSolution_c0[i] = samples[i][0]*samples[i][0]*samples[i][0]/6.0

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=4, c=inverse, vmin=min(inverse), vmax=max(inverse), cmap='jet')
cbar = plt.colorbar(scatter)
scale = 0.1 # scale the gradient for plotting purposes
#for i in range(0, min(1000, len(samples))):
for i in range(1000, 2000):
    ax.arrow(samples[i][0], samples[i][1], scale*gradient[i][0], scale*gradient[i][1], color='#252525', alpha=0.8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.ax.set_ylabel(r'Solution magnitude')
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
plt.savefig('figures/fig_WeightedLaplace_Solution.pdf', format='pdf', bbox_inches='tight')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)
pdf.close()

pdf = matplotlib.backends.backend_pdf.PdfPages("figures/Eigenfunctions_Lhat.pdf")

for i in range(min([25, len(eigvecsLhat)])):
    mn = min(eigvecsLhat[i])
    mx = max(eigvecsLhat[i])
    mntkx = floor(log10(abs(mn)))
    mxtkx = floor(log10(abs(mx)))
    if mx<0.0:
        mxtk = -10**mxtkx
    else:
        if mn>0.0:
            mxtkx = max([mntkx+1, mxtkx])
            mxtk = 10**mxtkx
        else:
            mxtk = 10**mntkx
    if mn<0.0:
        if mx<0.0:
            mntkx = max([mxtkx+1, mntkx])
            mntk = -10**mntkx
        else:
            mntk = -10**mntkx
    else:
        mntk = 10**mntkx

    # normalize the eigenvector to be sqrt(N)
    eigvec = np.sqrt(len(eigvecsLhat[i]))*eigvecsLhat[i]/np.linalg.norm(eigvecsLhat[i])

    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scale = 0.001*max(abs(mn), abs(mx))
    #scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvecsLhat[i], norm=mcolors.SymLogNorm(linthresh=scale, linscale=scale, vmin=min(mn, mntk), vmax=max(mx, mxtk)))
    if mxtk<0.0:
        mxtkm = 1.4*mxtk
    else:
        mxtkm = 0.6*mxtk
    if mntk<0.0:
        mntkm = 0.6*mntk
    else:
        mntkm = 1.4*mntk
    #scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvecsLhat[i], vmin=mntkm, vmax=mxtkm)
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvec, vmin=max([-5.0, min(eigvec)]), vmax=min([5.0, max(eigvec)]), cmap='jet')
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(r'Eigenfunction '+str(i+1)+' ($\lambda_{'+str(i+1)+'}$)')
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
pdf.close()
