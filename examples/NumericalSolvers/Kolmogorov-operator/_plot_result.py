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
densityEstimate = hdf5file['/density estimate'] [()].T [0]

bandwidthPara = hdf5file['/tune/candidate bandwidth parameters'] [()].T [0]
logKernelAvgDerivative = hdf5file['/tune/log kernel average change'] [()].T [0]

optInd = np.argmax(logKernelAvgDerivative)
print('Optimal bandwidth parameter:', bandwidthPara[optInd])

eigvecs = hdf5file['/eigenvectors'] [()].T
eigvals = hdf5file['/eigenvalues'] [()].T [0]
#eigvals = np.sort(eigvals)

f = hdf5file['/apply function'] [()].T
Lf = hdf5file['/applied Kolmogorov operator'] [()].T
Linvf = hdf5file['/applied inverse Kolmogorov operator'] [()].T

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=densityEstimate, vmin=0.0, vmax=0.16)
cbar = plt.colorbar(scatter)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
cbar.ax.set_ylabel(r'Density estimtaion $\psi^{(i)}$')
ax.set_xlabel(r'$x_0$')
ax.set_ylabel(r'$x_1$')
plt.savefig('figures/DensityEstimation.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)

pdf = matplotlib.backends.backend_pdf.PdfPages("figures/AppliedKolmogorovOperator.pdf")

for i in range(len(f)):
    ftitle = ''
    if i==0:
        ftitle = r'$f(\mathbf{x})=1$'
    if i==1:
        ftitle = r'$f(\mathbf{x})=x_0$'
    if i==2:
        ftitle = r'$f(\mathbf{x})=x_1$'

    Lftitle = r'$\mathcal{L}_{\psi, 1} f$ (' + ftitle + ')'
    Linvftitle = r'$\mathcal{L}_{\psi, 1}^{-1} f$ (' + ftitle + ')'

    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=f[i])
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(ftitle)
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=Lf[i])
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(Lftitle)
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=Linvf[i])
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(Linvftitle)
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
pdf.close()

pdf = matplotlib.backends.backend_pdf.PdfPages("figures/Eigenfunctions.pdf")
fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
ax.plot(range(1, len(eigvals)+1), eigvals, color='#525252')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_xlim([1, len(eigvals)])
ax.set_xlabel(r'Index $i$')
ax.set_ylabel(r'Eigenvalues $\lambda_i$')
pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)

for i in range(max([25, len(eigvecs)])):
    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvecs[i])
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(r'Eigenfunction '+str(i+1)+' ($\lambda_{'+str(i+1)+'}='+str(round(eigvals[i], 2))+'$)')
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
pdf.close()

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
ax.semilogx(bandwidthPara, logKernelAvgDerivative, color='#525252')
ax.plot(bandwidthPara[optInd], logKernelAvgDerivative[optInd], 'o', markersize=4, markerfacecolor='#a50f15', markeredgecolor='#a50f15')
ax.set_xlim([min(bandwidthPara), max(bandwidthPara)])
ax.set_xlabel(r'$\epsilon$')
ax.set_ylabel(r'$\Sigma_l^{\prime}$')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig('figures/LogKernelAvgDerivative.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)
