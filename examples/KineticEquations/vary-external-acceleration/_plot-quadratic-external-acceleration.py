import sys
import os, os.path

import numpy as np
from math import *

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

from PlotTools import *

figureDir = 'figures/quadratic-external-acceleration/'
outputDir = './output/quadratic-external-acceleration/'

time, energy = PlotFrames(figureDir, outputDir)

fig = MakeFigure(425, 0.9, False)
ax = plt.gca()
ax.plot(time, energy, color='#737373')
ax.set_xlabel(r'$\epsilon$')
ax.set_ylabel(r'$\Sigma_l^{\prime}$')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig(figureDir+'ExpectedEnergy.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)
