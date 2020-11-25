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
hdf5file = h5py.File('outputData_eigendecomposition.h5', 'r')

samples = hdf5file['/samples'] [()].T
densityEstimate = hdf5file['/density estimate'] [()].T [0]

bandwidthPara = hdf5file['/tune/candidate bandwidth parameters'] [()].T [0]
logKernelAvgDerivative = hdf5file['/tune/log kernel average change'] [()].T [0]

optInd = np.argmax(logKernelAvgDerivative)
print('Optimal bandwidth parameter:', bandwidthPara[optInd])

eigvecsL = hdf5file['/eigenvectors (L)'] [()].T
eigvecsLhat = hdf5file['/eigenvectors (Lhat)'] [()].T
eigvals = hdf5file['/eigenvalues'] [()].T [0]

f = hdf5file['/apply function'] [()].T
Lf = hdf5file['/applied Kolmogorov operator (L)'] [()].T
SinvLhatSf = hdf5file['/applied Kolmogorov operator (Sinv Lhat S)'] [()].T
Lhatinvf = hdf5file['/applied inverse transfromed Kolmogorov operator'] [()].T
Linvf = hdf5file['/applied inverse Kolmogorov operator'] [()].T


# find the point that is farthest from 0 (to identify the tails)
dist = 0.0
tailind = -1
for i in range(len(samples)):
    d = np.linalg.norm(samples[i])
    if d>dist:
        dist = d
        tailind = i

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

    Lftitle = r'$\mathbf{L}_{\psi, 1} \mathbf{f}$ (' + ftitle + ')'
    SinvLhatSftitle = r'$\mathbf{S}^{-1} \mathbf{\hat{L}}_{\psi, 1} \mathbf{S} \mathbf{f}$ (' + ftitle + ')'
    Lhatinvftitle = r'$\mathbf{\hat{L}}_{\psi, 1}^{-\dagger} \mathbf{f} \approx \mathbf{\hat{Q}} \mathbf{\Lambda}^{-\dagger} \mathbf{\hat{Q}}^T \mathbf{f}$ (' + ftitle + ')'
    Linvftitle = r'$\mathbf{L}_{\psi, 1}^{-\dagger} \mathbf{f} \approx \mathbf{S}^{-1} \mathbf{\hat{Q}} \mathbf{\Lambda}^{-\dagger} \mathbf{\hat{Q}}^T \mathbf{S} \mathbf{f}$ (' + ftitle + ')'

    # plot f
    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=f[i], vmin=min(f[i]), vmax=max(f[i]))
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

    # plot L f
    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=Lf[i], vmin=min(Lf[i]), vmax=max(Lf[i]))
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

    # plot Sinv Lhat S f
    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=SinvLhatSf[i], vmin=min(SinvLhatSf[i]), vmax=max(SinvLhatSf[i]))
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(SinvLhatSftitle)
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # plot Lhatinv f = Uhat lambdainv Uhat^T f
    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=Lhatinvf[i], vmin=min(Lhatinvf[i]), vmax=max(Lhatinvf[i]))
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(Lhatinvftitle)
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # plot Linv f = S^{-1} Uhat lambdainv Uhat^T S f
    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=Linvf[i], vmin=min(Linvf[i]), vmax=max(Linvf[i]))
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

pdf = matplotlib.backends.backend_pdf.PdfPages("figures/Eigenfunctions_L.pdf")
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

for i in range(min([25, len(eigvecsL)])):
    mn = min(eigvecsL[i])
    mx = max(eigvecsL[i])
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
    eigvec = np.sqrt(len(eigvecsL[i]))*eigvecsL[i]/np.linalg.norm(eigvecsL[i])

    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    scale = 0.001*max(abs(mn), abs(mx))
    #scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvecsL[i], norm=mcolors.SymLogNorm(linthresh=scale, linscale=scale, vmin=min(mn, mntk), vmax=max(mx, mxtk)))
    if mxtk<0.0:
        mxtkm = 1.5*mxtk
    else:
        mxtkm = 0.5*mxtk
    if mntk<0.0:
        mntkm = 0.5*mntk
    else:
        mntkm = 1.5*mntk
    #scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvecsL[i], vmin=mntkm, vmax=mxtkm)
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvec, vmin=max([-5.0, min(eigvec)]), vmax=min([5.0, max(eigvec)]))
    #cbar = plt.colorbar(scatter, ticks=np.sort([mntk, 0.0, mxtk]))
    cbar = plt.colorbar(scatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    #cbar.ax.set_ylabel(r'Magnitude')
    ax.set_title(r'Eigenfunction '+str(i+1)+' ($\lambda_{'+str(i+1)+'}='+str(round(eigvals[i], 2))+'$)')
    ax.set_xlabel(r'$x_0$')
    ax.set_ylabel(r'$x_1$')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
pdf.close()

pdf = matplotlib.backends.backend_pdf.PdfPages("figures/Eigenfunctions_Lhat.pdf")

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
    scatter = ax.scatter(samples.T[0], samples.T[1], s=3, c=eigvec, vmin=max([-5.0, min(eigvec)]), vmax=min([5.0, max(eigvec)]))
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
