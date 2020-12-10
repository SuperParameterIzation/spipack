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

figureDir = 'figures/inelastic-collisions/'
outputDir = './output/inelastic-collisions/'

time, energy = PlotFrames(figureDir, outputDir)

PlotExpectedEnergy(figureDir, time, energy)
