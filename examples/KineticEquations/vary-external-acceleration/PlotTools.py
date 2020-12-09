import os, os.path

import numpy as np

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

def PlotFrames(figureDir, outputDir):
    filenames = [name for name in os.listdir(outputDir) if os.path.isfile(os.path.join(outputDir, name))]

    time = list()
    energy = list()
    for filename in filenames:
        if filename[-3:]!='.h5':
            continue

        figname = figureDir+'/frames/'+filename[:-3]+'.png'

        # load the data file
        hdf5file = h5py.File(outputDir+filename, 'r')

        time = time + hdf5file['/time'] [()] [0].tolist()
        samples = hdf5file['/samples'] [()].T
        acceleration = hdf5file['/acceleration'] [()]
        energy = energy + hdf5file['/expected energy'] [()] [0].tolist()
        sampleMean = hdf5file['/sample mean'] [()].T [0]

        #fig = MakeFigure(425, 0.9, False)
        #ax = plt.gca()
        #scatter = ax.scatter(samples.T[0], samples.T[1], s=3, color='#525252')
        #ax.spines['right'].set_visible(False)
        #ax.spines['top'].set_visible(False)
        #ax.yaxis.set_ticks_position('left')
        #ax.xaxis.set_ticks_position('bottom')
        #ax.set_ylim([-3.5, 3.5])
        #ax.set_xlim([-3.5, 3.5])
        #ax.set_xlabel(r'$V_0$')
        #ax.set_ylabel(r'$V_1$')
        #plt.savefig(figname, format='png', bbox_inches='tight', dpi=300)
        #plt.close(fig)

    sortIndex = np.argsort(time)
    time = [time[i] for i in sortIndex]
    energy = [energy[i] for i in sortIndex]

    return time, energy

def PlotExpectedEnergy(figureDir, time, energy):
    fig = MakeFigure(425, 0.9, False)
    ax = plt.gca()
    ax.plot(time, energy, color='#737373')
    ax.set_xlim([time[0], time[-1]])
    ax.set_ylim([0.0, 1.01*max(energy)])
    ax.set_xlabel(r'Time $t$')
    ax.set_ylabel(r'Expected energy $e_{\psi}$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.savefig(figureDir+'ExpectedEnergy.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)
