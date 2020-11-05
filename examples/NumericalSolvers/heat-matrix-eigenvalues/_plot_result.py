import sys

import numpy as np
import scipy.stats

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

# load the data file
file = h5py.File('samples.h5', 'r')
